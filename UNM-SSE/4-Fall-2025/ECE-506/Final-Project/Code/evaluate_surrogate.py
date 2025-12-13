import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.time import Time, TimeDelta
from astropy import units as u

import torch
import torch.nn as nn

from porkchop_run import delta_v_objective

try:
    from scipy.optimize import differential_evolution, minimize
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


CLEAN_CSV_PATH = "porkchop_earth_mars_2030_2032_daily_clean.csv"
MODEL_CHECKPOINT = "dv_surrogate_model_2030_2032_supervised.pt"
NORM_STATS_PATH = "normalization_stats.json"


class DVSurrogate(nn.Module):
    def __init__(self, in_dim=2, hidden_dim=128, dropout_p=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_p),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_p),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x):
        return self.net(x)


def _window_tag(dep_start, dep_end, min_tof_days, max_tof_days):
    t_dep0 = Time(dep_start, scale="tdb").tdb
    t_dep1 = Time(dep_end, scale="tdb").tdb
    t_arr0 = (t_dep0 + TimeDelta(float(min_tof_days) * u.day)).tdb
    t_arr1 = (t_dep1 + TimeDelta(float(max_tof_days) * u.day)).tdb

    dep0 = t_dep0.to_datetime().strftime("%Y%m%d")
    dep1 = t_dep1.to_datetime().strftime("%Y%m%d")
    arr0 = t_arr0.to_datetime().strftime("%Y%m%d")
    arr1 = t_arr1.to_datetime().strftime("%Y%m%d")

    return f"dep{dep0}-{dep1}_arr{arr0}-{arr1}"


def rebuild_grid_from_df(df, dep_decimals=9, tof_decimals=6):
    dep_julians = np.sort(df["dep_julian"].unique()).astype(np.float64)
    tofs_days = np.sort(df["tof_days"].unique()).astype(np.float64)
    dep_times = Time(dep_julians, format="jd", scale="tdb").tdb

    dep_keys = np.round(dep_julians, dep_decimals)
    tof_keys = np.round(tofs_days, tof_decimals)

    dep_index = {float(k): i for i, k in enumerate(dep_keys)}
    tof_index = {float(k): j for j, k in enumerate(tof_keys)}

    grid = (
        df.assign(
            dep_key=np.round(df["dep_julian"].astype(np.float64), dep_decimals),
            tof_key=np.round(df["tof_days"].astype(np.float64), tof_decimals),
        )
        .pivot(index="dep_key", columns="tof_key", values="dv_total_kms")
        .reindex(index=dep_keys, columns=tof_keys)
        .to_numpy()
        .astype(np.float64)
    )

    return dep_times, tofs_days.astype(np.float64), grid, dep_julians, dep_index, tof_index, dep_decimals, tof_decimals


def plot_porkchop_triptych(dep_times, tofs_days, dv_true, dv_pred, dv_err, tag):
    os.makedirs("plots", exist_ok=True)

    dep_mpl = mdates.date2num(dep_times.to_datetime())
    D_dep, D_tof = np.meshgrid(dep_mpl, tofs_days, indexing="ij")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    vmin = np.nanmin(dv_true)
    vmax = np.nanmax(dv_true)
    err_abs_max = np.nanmax(np.abs(dv_err))

    configs = [
        ("True Lambert ΔV", dv_true, (vmin, vmax), r"Total $\Delta V$ [km/s]"),
        ("Surrogate Predicted ΔV", dv_pred, (vmin, vmax), r"Total $\Delta V$ [km/s]"),
        ("Error (Pred − True)", dv_err, (-err_abs_max, err_abs_max), r"Error [km/s]"),
    ]

    for ax, (title, grid, (vmin_c, vmax_c), cbar_label) in zip(axes, configs):
        Z = np.ma.masked_invalid(grid)
        cs = ax.contourf(D_dep, D_tof, Z, levels=30, vmin=vmin_c, vmax=vmax_c)
        ax.set_title(title)

        try:
            cs_lines = ax.contour(D_dep, D_tof, Z, levels=10, colors="k", linewidths=0.5)
            ax.clabel(cs_lines, inline=True, fontsize=7, fmt="%.1f")
        except Exception:
            pass

        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.tick_params(axis="x", rotation=30)
        ax.set_xlabel("Departure Date")
        ax.set_ylabel("Time of Flight [days]")

        cbar = fig.colorbar(cs, ax=ax)
        cbar.set_label(cbar_label)

    save_path = os.path.join("plots", f"porkchop_triptych_{tag}.png")
    plt.savefig(save_path, dpi=300)
    print(f"Saved 1x3 triptych: {save_path}")

def plot_error_histogram(errors, tag):
    os.makedirs("plots", exist_ok=True)

    plt.figure(figsize=(6, 4))
    plt.hist(errors, bins=50)
    plt.xlabel(r"Error $\Delta V_{\mathrm{pred}} - \Delta V_{\mathrm{true}}$ [km/s]")
    plt.ylabel("Count")
    plt.title("Surrogate ΔV Error Histogram")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    hist_path = os.path.join("plots", f"porkchop_error_histogram_{tag}.png")
    plt.savefig(hist_path, dpi=300)
    print(f"Saved error histogram: {hist_path}")

def plot_porkchop_four_mins(
    dep_times,
    tofs_days,
    dv_true_grid,
    grid_min_true,
    true_opt_min,
    sur_grid_min,
    sur_opt_min,
    tag,
):
    os.makedirs("plots", exist_ok=True)

    dep_mpl = mdates.date2num(dep_times.to_datetime())
    D_dep, D_tof = np.meshgrid(dep_mpl, tofs_days, indexing="ij")
    Z = np.ma.masked_invalid(dv_true_grid)

    plt.figure(figsize=(10, 5))
    cs = plt.contourf(D_dep, D_tof, Z, levels=30)
    cbar = plt.colorbar(cs)
    cbar.set_label(r"Total $\Delta V$ [km/s]")

    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.xticks(rotation=30, ha="right")
    plt.xlabel("Departure Date")
    plt.ylabel("Time of Flight [days]")

    def pt(dep_time_tdb, tof_days):
        return mdates.date2num(dep_time_tdb.to_datetime()), float(tof_days)

    xg, yg = pt(grid_min_true["dep_time"], grid_min_true["tof_days"])
    plt.scatter([xg], [yg], marker="*", s=180, label=f"Grid true min: {grid_min_true['dv_kms']:.3f}")

    if true_opt_min is not None:
        xo, yo = pt(true_opt_min["dep_time"], true_opt_min["tof_days"])
        plt.scatter([xo], [yo], marker="X", s=140, label=f"True opt min: {true_opt_min['dv_kms']:.3f}")

    xsg, ysg = pt(sur_grid_min["dep_time"], sur_grid_min["tof_days"])
    plt.scatter([xsg], [ysg], marker="D", s=90, label=f"Sur grid min (pred): {sur_grid_min['dv_kms']:.3f}")

    xso, yso = pt(sur_opt_min["dep_time"], sur_opt_min["tof_days"])
    plt.scatter([xso], [yso], marker="o", s=90, label=f"Sur opt min (pred): {sur_opt_min['dv_pred_kms']:.3f}")

    plt.title("Porkchop minima overlay")
    plt.legend(loc="best")
    plt.tight_layout()

    out_path = os.path.join("plots", f"porkchop_four_mins_{tag}.png")
    plt.savefig(out_path, dpi=300)
    print(f"Saved four-mins plot: {out_path}")

def load_norm_stats(path):
    with open(path, "r") as f:
        norm_info = json.load(f)
    X_mean = np.array(norm_info["X_mean"], dtype=np.float32)
    X_std = np.array(norm_info["X_std"], dtype=np.float32)
    return X_mean, X_std


def load_model(path, device):
    ckpt = torch.load(path, map_location=device)
    in_dim = int(ckpt.get("in_dim", 2))
    hidden_dim = int(ckpt.get("hidden_dim", 128))
    dropout_p = float(ckpt.get("dropout_p", 0.1))
    y_mean = float(ckpt.get("y_mean", 0.0))
    y_std = float(ckpt.get("y_std", 1.0))

    model = DVSurrogate(in_dim=in_dim, hidden_dim=hidden_dim, dropout_p=dropout_p)
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(device)
    model.eval()

    return model, in_dim, hidden_dim, dropout_p, y_mean, y_std


def predict_surrogate_kms(model, dep_jd, tof_days, X_mean, X_std, y_mean, y_std, device):
    x = np.array([[dep_jd, tof_days]], dtype=np.float32)
    x_norm = (x - X_mean) / X_std
    x_t = torch.from_numpy(x_norm).to(device)
    with torch.no_grad():
        y_norm = float(model(x_t).detach().cpu().numpy()[0, 0])
    return float(y_norm * y_std + y_mean)


def build_surrogate_grid(
    df_clean,
    dv_true_grid,
    dep_index,
    tof_index,
    model,
    X_mean,
    X_std,
    y_mean,
    y_std,
    device,
    dep_decimals=9,
    tof_decimals=6,
):
    dv_pred_grid = np.full_like(dv_true_grid, np.nan, dtype=np.float64)

    dep_vals64 = df_clean["dep_julian"].values.astype(np.float64)
    tof_vals64 = df_clean["tof_days"].values.astype(np.float64)

    X = np.stack([dep_vals64, tof_vals64], axis=1).astype(np.float32)
    X_norm = (X - X_mean.reshape(1, -1)) / X_std.reshape(1, -1)

    X_t = torch.from_numpy(X_norm).to(device)
    with torch.no_grad():
        y_norm = model(X_t).detach().cpu().numpy().reshape(-1)

    y_pred = y_norm * y_std + y_mean

    dep_keys = np.round(dep_vals64, dep_decimals)
    tof_keys = np.round(tof_vals64, tof_decimals)

    for k in range(len(df_clean)):
        i = dep_index[float(dep_keys[k])]
        j = tof_index[float(tof_keys[k])]
        dv_pred_grid[i, j] = float(y_pred[k])

    return dv_pred_grid


def find_grid_min_from_grid(dep_times, tofs_days, dv_grid):
    Z = np.array(dv_grid, dtype=np.float64)
    if not np.isfinite(Z).any():
        raise RuntimeError("No finite values in grid.")
    idx = np.nanargmin(Z)
    i, j = np.unravel_index(idx, Z.shape)
    dep_time = dep_times[i].tdb
    tof = float(tofs_days[j])
    arr_time = (dep_time + TimeDelta(tof * u.day)).tdb
    dv = float(Z[i, j])
    return {"i": int(i), "j": int(j), "dep_time": dep_time, "tof_days": tof, "arr_time": arr_time, "dv_kms": dv}


def find_surrogate_grid_min(dep_times, tofs_days, dv_pred_grid):
    Z = np.array(dv_pred_grid, dtype=np.float64)
    if not np.isfinite(Z).any():
        raise RuntimeError("No finite values in surrogate grid.")
    idx = np.nanargmin(Z)
    i, j = np.unravel_index(idx, Z.shape)
    dep_time = dep_times[i].tdb
    tof = float(tofs_days[j])
    arr_time = (dep_time + TimeDelta(tof * u.day)).tdb
    dv = float(Z[i, j])
    return {"i": int(i), "j": int(j), "dep_time": dep_time, "tof_days": tof, "arr_time": arr_time, "dv_kms": dv}


def run_surrogate_optimizer(
    dep_start,
    dep_end,
    min_tof_days,
    max_tof_days,
    model,
    X_mean,
    X_std,
    y_mean,
    y_std,
    device,
    method="de",
    seed=2,
):
    if not _HAS_SCIPY:
        raise RuntimeError("SciPy is required for surrogate optimization")

    t0 = Time(dep_start, scale="tdb").tdb
    t1 = Time(dep_end, scale="tdb").tdb
    dep_bounds = (float(t0.jd), float(t1.jd))
    tof_bounds = (float(min_tof_days), float(max_tof_days))

    BIG = 1e6
    history = []

    def f(x):
        dep_jd = float(x[0])
        tof = float(x[1])
        if dep_jd < dep_bounds[0] or dep_jd > dep_bounds[1] or tof < tof_bounds[0] or tof > tof_bounds[1]:
            return BIG
        return predict_surrogate_kms(model, dep_jd, tof, X_mean, X_std, y_mean, y_std, device)

    x0 = np.array([(dep_bounds[0] + dep_bounds[1]) / 2.0, (tof_bounds[0] + tof_bounds[1]) / 2.0], dtype=np.float64)

    if method.lower() == "de":
        bounds = [dep_bounds, tof_bounds]

        def cb(xk, convergence=None):
            history.append((float(xk[0]), float(xk[1]), float(f(xk))))
            return False

        res = differential_evolution(
            f,
            bounds=bounds,
            seed=seed,
            callback=cb,
            polish=True,
            updating="deferred",
            workers=1,
        )
        xbest = res.x
    else:
        bounds = [dep_bounds, tof_bounds]

        def cb(xk):
            history.append((float(xk[0]), float(xk[1]), float(f(xk))))

        res = minimize(
            f,
            x0=x0,
            method="L-BFGS-B",
            bounds=bounds,
            callback=cb,
        )
        xbest = res.x

    dep_jd_best = float(xbest[0])
    tof_best = float(xbest[1])
    dep_time_best = Time(dep_jd_best, format="jd", scale="tdb").tdb
    arr_time_best = (dep_time_best + TimeDelta(tof_best * u.day)).tdb
    dv_pred_best = float(f([dep_jd_best, tof_best]))

    return {
        "method": method,
        "success": bool(getattr(res, "success", True)),
        "message": str(getattr(res, "message", "")),
        "nfev": int(getattr(res, "nfev", -1)) if getattr(res, "nfev", None) is not None else -1,
        "dep_time": dep_time_best,
        "tof_days": tof_best,
        "arr_time": arr_time_best,
        "dv_pred_kms": dv_pred_best,
        "history": history,
    }


def evaluate_surrogate(
    dep_start,
    dep_end,
    min_tof_days,
    max_tof_days,
    clean_csv_path=CLEAN_CSV_PATH,
    model_checkpoint=MODEL_CHECKPOINT,
    norm_stats_path=NORM_STATS_PATH,
    true_opt_result=None,
    dv_cutoff=20.0,
    surrogate_opt_method="de",
    surrogate_opt_seed=2,
    dep_decimals=9,
    tof_decimals=6,
):
    if not os.path.exists(clean_csv_path):
        raise FileNotFoundError(clean_csv_path)
    if not os.path.exists(model_checkpoint):
        raise FileNotFoundError(model_checkpoint)
    if not os.path.exists(norm_stats_path):
        raise FileNotFoundError(norm_stats_path)

    tag = _window_tag(dep_start, dep_end, min_tof_days, max_tof_days)

    df_clean = pd.read_csv(clean_csv_path)

    (
        dep_times,
        tofs_days,
        dv_true_grid,
        _,
        dep_index,
        tof_index,
        dep_decimals_used,
        tof_decimals_used,
    ) = rebuild_grid_from_df(df_clean, dep_decimals=dep_decimals, tof_decimals=tof_decimals)

    X_mean, X_std = load_norm_stats(norm_stats_path)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model, in_dim, hidden_dim, dropout_p, y_mean, y_std = load_model(model_checkpoint, device)

    dv_pred_grid = build_surrogate_grid(
        df_clean,
        dv_true_grid,
        dep_index,
        tof_index,
        model,
        X_mean,
        X_std,
        y_mean,
        y_std,
        device,
        dep_decimals=dep_decimals_used,
        tof_decimals=tof_decimals_used,
    )

    dv_err_grid = dv_pred_grid - dv_true_grid

    plot_porkchop_triptych(dep_times, tofs_days, dv_true_grid, dv_pred_grid, dv_err_grid, tag)

    valid_mask = np.isfinite(dv_true_grid) & np.isfinite(dv_pred_grid)
    errors = dv_err_grid[valid_mask]
    plot_error_histogram(errors, tag)

    mae = float(np.mean(np.abs(errors)))
    rmse = float(np.sqrt(np.mean(errors ** 2)))

    grid_min_true = find_grid_min_from_grid(dep_times, tofs_days, dv_true_grid)
    sur_grid_min = find_surrogate_grid_min(dep_times, tofs_days, dv_pred_grid)

    sur_opt = run_surrogate_optimizer(
        dep_start,
        dep_end,
        min_tof_days,
        max_tof_days,
        model,
        X_mean,
        X_std,
        y_mean,
        y_std,
        device,
        method=surrogate_opt_method,
        seed=surrogate_opt_seed,
    )

    dv_true_at_sur = float(
        delta_v_objective(
            sur_opt["dep_time"],
            sur_opt["tof_days"],
            dv_cutoff=dv_cutoff,
            return_components=False,
        )
    )

    dv_true_at_sur_grid = float(
        delta_v_objective(
            sur_grid_min["dep_time"],
            sur_grid_min["tof_days"],
            dv_cutoff=dv_cutoff,
            return_components=False,
        )
    )

    plot_porkchop_four_mins(
        dep_times,
        tofs_days,
        dv_true_grid,
        grid_min_true=grid_min_true,
        true_opt_min=true_opt_result,
        sur_grid_min=sur_grid_min,
        sur_opt_min=sur_opt,
        tag=tag,
    )

    print("SURROGATE EVAL")
    print(f"mae_kms {mae:.6f}")
    print(f"rmse_kms {rmse:.6f}")

    print("GRID TRUE MIN")
    print(grid_min_true["dep_time"].iso)
    print(f"{grid_min_true['tof_days']:.6f}")
    print(grid_min_true["arr_time"].iso)
    print(f"{grid_min_true['dv_kms']:.6f}")

    if true_opt_result is not None:
        print("TRUE OBJECTIVE OPT MIN")
        print(true_opt_result["dep_time"].iso)
        print(f"{true_opt_result['tof_days']:.6f}")
        print(true_opt_result["arr_time"].iso)
        print(f"{true_opt_result['dv_kms']:.6f}")
        print(len(true_opt_result.get("history", [])))

    print("SURROGATE GRID MIN")
    print(sur_grid_min["dep_time"].iso)
    print(f"{sur_grid_min['tof_days']:.6f}")
    print(sur_grid_min["arr_time"].iso)
    print(f"{sur_grid_min['dv_kms']:.6f}")
    print(f"{dv_true_at_sur_grid:.6f}")

    print("SURROGATE OPT MIN")
    print(sur_opt["dep_time"].iso)
    print(f"{sur_opt['tof_days']:.6f}")
    print(sur_opt["arr_time"].iso)
    print(f"{sur_opt['dv_pred_kms']:.6f}")
    print(f"{dv_true_at_sur:.6f}")
    print(len(sur_opt.get("history", [])))

    print("SUMMARY")
    print(f"grid_true_dv_kms {grid_min_true['dv_kms']:.6f}")
    if true_opt_result is not None:
        print(f"trueopt_true_dv_kms {true_opt_result['dv_kms']:.6f}")
    print(f"surggrid_pred_dv_kms {sur_grid_min['dv_kms']:.6f}")
    print(f"surggrid_true_dv_kms {dv_true_at_sur_grid:.6f}")
    print(f"suropt_pred_dv_kms {sur_opt['dv_pred_kms']:.6f}")
    print(f"suropt_true_dv_kms {dv_true_at_sur:.6f}")

    return {
        "mae_kms": mae,
        "rmse_kms": rmse,
        "grid_min_true": grid_min_true,
        "true_opt_result": true_opt_result,
        "sur_grid_min": sur_grid_min,
        "sur_grid_min_true_dv": dv_true_at_sur_grid,
        "sur_opt": sur_opt,
        "sur_opt_true_dv": dv_true_at_sur,
        "device": device,
        "model_meta": {
            "in_dim": in_dim,
            "hidden_dim": hidden_dim,
            "dropout_p": dropout_p,
            "y_mean": y_mean,
            "y_std": y_std,
        },
        "rounding": {"dep_decimals": dep_decimals_used, "tof_decimals": tof_decimals_used},
        "tag": tag,
    }


if __name__ == "__main__":
    evaluate_surrogate(
        dep_start="2030-01-01",
        dep_end="2032-01-01",
        min_tof_days=120.0,
        max_tof_days=360.0,
        clean_csv_path=CLEAN_CSV_PATH,
        model_checkpoint=MODEL_CHECKPOINT,
        norm_stats_path=NORM_STATS_PATH,
        true_opt_result=None,
        dv_cutoff=20.0,
        surrogate_opt_method="de",
        surrogate_opt_seed=2,
        dep_decimals=9,
        tof_decimals=6,
    )
