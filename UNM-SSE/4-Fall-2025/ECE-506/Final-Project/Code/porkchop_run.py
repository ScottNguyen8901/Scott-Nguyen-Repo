#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import argparse
import numpy as np
import pandas as pd

from astropy import units as u
from astropy.time import Time, TimeDelta

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from poliastro.bodies import Sun, Earth, Mars
from poliastro.twobody import Orbit
from poliastro.iod import izzo

from tqdm.auto import tqdm

try:
    from joblib import Parallel, delayed
    _HAS_JOBLIB = True
except Exception:
    _HAS_JOBLIB = False

try:
    from poliastro.ephem import Ephem
    _HAS_EPHEM = True
except Exception:
    _HAS_EPHEM = False

try:
    from scipy.optimize import differential_evolution, minimize
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

DV_CUTOFF = 20.0
DATA_CSV_PATH = "porkchop_earth_mars_2030_2032_daily.csv"
N_JOBS = -1


def get_heliocentric_state(body, epoch):
    orb = Orbit.from_body_ephem(body, epoch)
    return orb.r, orb.v


def lambert_delta_v(r1, v1_planet, r2, v2_planet, tof_seconds):
    (v1_trans, v2_trans), = izzo.lambert(Sun.k, r1, r2, tof_seconds)
    dv1_vec = v1_trans - v1_planet
    dv2_vec = v2_planet - v2_trans
    dv1 = np.linalg.norm(dv1_vec.to(u.km / u.s).value)
    dv2 = np.linalg.norm(dv2_vec.to(u.km / u.s).value)
    return dv1 + dv2, dv1, dv2


def delta_v_objective(dep_time, tof_days, dv_cutoff=DV_CUTOFF, return_components=False):
    try:
        tof_days = float(tof_days)
        if tof_days <= 0.0:
            return (np.nan, np.nan, np.nan) if return_components else np.nan

        tof_seconds = (tof_days * u.day).to(u.s)
        r_E, v_E = get_heliocentric_state(Earth, dep_time)
        arr_time = dep_time + TimeDelta(tof_days * u.day)
        r_M, v_M = get_heliocentric_state(Mars, arr_time)
        dv_total, dv1, dv2 = lambert_delta_v(r_E, v_E, r_M, v_M, tof_seconds)

        if (not np.isfinite(dv_total)) or (dv_total > dv_cutoff):
            return (np.nan, np.nan, np.nan) if return_components else np.nan

        return (dv_total, dv1, dv2) if return_components else dv_total

    except Exception:
        return (np.nan, np.nan, np.nan) if return_components else np.nan


def _mars_states_vectorized(arr_times):
    if _HAS_EPHEM:
        try:
            eph = Ephem.from_body(Mars, arr_times)
            r_M, v_M = eph.rv(arr_times)
            r_list = [r_M[i] for i in range(len(arr_times))]
            v_list = [v_M[i] for i in range(len(arr_times))]
            return r_list, v_list
        except Exception:
            pass

    r_list, v_list = [], []
    for t in arr_times:
        r, v = get_heliocentric_state(Mars, t)
        r_list.append(r)
        v_list.append(v)
    return r_list, v_list


def _compute_one_departure_fast(t_dep, tofs_days, dv_cutoff):
    n_tof = len(tofs_days)

    dv_row = np.full(n_tof, np.nan, dtype=np.float64)
    dv1_row = np.full(n_tof, np.nan, dtype=np.float64)
    dv2_row = np.full(n_tof, np.nan, dtype=np.float64)

    r_E, v_E = get_heliocentric_state(Earth, t_dep)

    arr_times = t_dep + TimeDelta(tofs_days * u.day)
    arr_tdb = arr_times.tdb

    arr_jd_row = arr_tdb.jd.astype(np.float64)
    arr_iso_row = arr_tdb.iso.astype(object)

    r_M_list, v_M_list = _mars_states_vectorized(arr_times)

    for j, tof_days in enumerate(tofs_days):
        try:
            tof_seconds = (float(tof_days) * u.day).to(u.s)
            dv_total, dv1, dv2 = lambert_delta_v(r_E, v_E, r_M_list[j], v_M_list[j], tof_seconds)

            if (not np.isfinite(dv_total)) or (dv_total > dv_cutoff):
                continue

            dv_row[j] = dv_total
            dv1_row[j] = dv1
            dv2_row[j] = dv2
        except Exception:
            continue

    return dv_row, dv1_row, dv2_row, arr_jd_row, arr_iso_row


def generate_porkchop_data_daily(
    dep_start="2030-01-01",
    dep_end="2032-01-01",
    arr_start=None,
    arr_end=None,
    departure_step_days=1,
    min_tof_days=100.0,
    max_tof_days=400.0,
    n_tof=60,
    csv_path=DATA_CSV_PATH,
    dv_cutoff=DV_CUTOFF,
    n_jobs=N_JOBS,
):
    t_dep_start = Time(dep_start, scale="tdb").tdb
    t_dep_end = Time(dep_end, scale="tdb").tdb

    if arr_start is not None and arr_end is not None:
        t_arr_start = Time(arr_start, scale="tdb").tdb
        t_arr_end = Time(arr_end, scale="tdb").tdb
        min_tof = (t_arr_start - t_dep_end).to(u.day).value
        max_tof = (t_arr_end - t_dep_start).to(u.day).value
        min_tof_days = float(max(min_tof, 1.0))
        max_tof_days = float(max(max_tof, min_tof_days + 1.0))
        print(f"Derived TOF range: [{min_tof_days:.3f}, {max_tof_days:.3f}] days")

    total_days = int((t_dep_end - t_dep_start).to(u.day).value)
    day_offsets = np.arange(0, total_days + 1, int(departure_step_days), dtype=int)
    dep_times = (t_dep_start + TimeDelta(day_offsets * u.day)).tdb
    n_departures = len(dep_times)

    tofs_days = np.linspace(float(min_tof_days), float(max_tof_days), int(n_tof), dtype=np.float64)
    n_tof = len(tofs_days)

    print(f"Grid: {n_departures} departures × {n_tof} TOFs")
    print(f"Departure window: {t_dep_start.iso} to {t_dep_end.iso}")
    if arr_start is not None and arr_end is not None:
        print(f"Arrival window: {Time(arr_start, scale='tdb').tdb.iso} to {Time(arr_end, scale='tdb').tdb.iso}")
    print(f"DV cutoff: {dv_cutoff} km/s")
    if not _HAS_EPHEM:
        print("Ephem unavailable; Mars ephemeris will be slower.")

    dep_jd = dep_times.jd.astype(np.float64)
    dep_iso = dep_times.iso.astype(object)

    dv_grid = np.full((n_departures, n_tof), np.nan, dtype=np.float64)
    dv1_grid = np.full((n_departures, n_tof), np.nan, dtype=np.float64)
    dv2_grid = np.full((n_departures, n_tof), np.nan, dtype=np.float64)

    arr_jd_grid = np.empty((n_departures, n_tof), dtype=np.float64)
    arr_iso_grid = np.empty((n_departures, n_tof), dtype=object)

    do_parallel = _HAS_JOBLIB and isinstance(n_jobs, int) and (n_jobs != 1) and (n_departures > 1)

    if do_parallel:
        results = Parallel(n_jobs=n_jobs, prefer="processes")(
            delayed(_compute_one_departure_fast)(dep_times[i], tofs_days, dv_cutoff)
            for i in range(n_departures)
        )
        for i, (dv_row, dv1_row, dv2_row, arr_jd_row, arr_iso_row) in enumerate(results):
            dv_grid[i, :] = dv_row
            dv1_grid[i, :] = dv1_row
            dv2_grid[i, :] = dv2_row
            arr_jd_grid[i, :] = arr_jd_row
            arr_iso_grid[i, :] = arr_iso_row
    else:
        for i, t_dep in enumerate(tqdm(dep_times, desc="Departures", unit="dep")):
            dv_row, dv1_row, dv2_row, arr_jd_row, arr_iso_row = _compute_one_departure_fast(
                t_dep, tofs_days, dv_cutoff
            )
            dv_grid[i, :] = dv_row
            dv1_grid[i, :] = dv1_row
            dv2_grid[i, :] = dv2_row
            arr_jd_grid[i, :] = arr_jd_row
            arr_iso_grid[i, :] = arr_iso_row

    n_total = n_departures * n_tof
    dep_jd_flat = np.repeat(dep_jd, n_tof)
    dep_iso_flat = np.repeat(dep_iso, n_tof)
    tof_flat = np.tile(tofs_days, n_departures)

    df = pd.DataFrame(
        {
            "dep_julian": dep_jd_flat,
            "arr_julian": arr_jd_grid.reshape(n_total),
            "dep_iso": dep_iso_flat,
            "arr_iso": arr_iso_grid.reshape(n_total),
            "tof_days": tof_flat.astype(np.float64),
            "dv_total_kms": dv_grid.reshape(n_total),
            "dv_dep_kms": dv1_grid.reshape(n_total),
            "dv_arr_kms": dv2_grid.reshape(n_total),
        }
    )

    if arr_start is not None and arr_end is not None:
        t_arr_start = Time(arr_start, scale="tdb").tdb
        t_arr_end = Time(arr_end, scale="tdb").tdb
        arr_t = Time(df["arr_julian"].values.astype(np.float64), format="jd", scale="tdb").tdb
        in_window = (arr_t.jd >= t_arr_start.jd) & (arr_t.jd <= t_arr_end.jd)
        df.loc[~in_window, ["dv_total_kms", "dv_dep_kms", "dv_arr_kms"]] = np.nan

    df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path} ({len(df)} rows)")
    return df, dep_times, tofs_days, dv_grid


def plot_porkchop(dep_times, tofs_days, dv_grid, title_suffix="run"):
    os.makedirs("plots", exist_ok=True)

    dep_datetimes = dep_times.to_datetime()
    dep_mpl = mdates.date2num(dep_datetimes)

    D_dep, D_tof = np.meshgrid(dep_mpl, tofs_days, indexing="ij")
    Z = np.ma.masked_invalid(dv_grid)

    plt.figure(figsize=(10, 5))
    cs = plt.contourf(D_dep, D_tof, Z, levels=30)
    cbar = plt.colorbar(cs)
    cbar.set_label(r"Total $\Delta V$ [km/s]")

    try:
        cs_lines = plt.contour(D_dep, D_tof, Z, levels=10, colors="k", linewidths=0.5)
        plt.clabel(cs_lines, inline=True, fontsize=8, fmt="%.1f")
    except Exception:
        pass

    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    plt.xticks(rotation=30, ha="right")
    plt.xlabel("Departure Date")
    plt.ylabel("Time of Flight [days]")
    plt.title(rf"Earth–Mars Pork-Chop (Total $\Delta V$) {title_suffix}")

    plt.tight_layout()
    save_path = f"plots/porkchop_{title_suffix}.png"
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot: {save_path}")

def rebuild_grid_from_csv(df):
    dep_julians = np.sort(df["dep_julian"].unique())
    tofs_days = np.sort(df["tof_days"].unique())
    dep_times = Time(dep_julians, format="jd", scale="tdb").tdb

    grid = (
        df.pivot(index="dep_julian", columns="tof_days", values="dv_total_kms")
        .reindex(index=dep_julians, columns=tofs_days)
        .to_numpy()
    )
    return dep_times, tofs_days, grid


def find_global_min_from_grid(dep_times, tofs_days, dv_grid):
    Z = np.array(dv_grid, dtype=np.float64)
    if not np.isfinite(Z).any():
        raise RuntimeError("dv_grid has no finite values")
    idx = np.nanargmin(Z)
    i, j = np.unravel_index(idx, Z.shape)
    dep_time = dep_times[i].tdb
    tof = float(tofs_days[j])
    arr_time = (dep_time + TimeDelta(tof * u.day)).tdb
    dv = float(Z[i, j])
    return {"i": int(i), "j": int(j), "dep_time": dep_time, "tof_days": tof, "arr_time": arr_time, "dv_kms": dv}


def _continuous_objective_factory(dv_cutoff, dep_jd_bounds, tof_bounds):
    dep_lo, dep_hi = dep_jd_bounds
    tof_lo, tof_hi = tof_bounds
    BIG = 1e6

    def f(x):
        dep_jd = float(x[0])
        tof_days = float(x[1])

        if dep_jd < dep_lo or dep_jd > dep_hi or tof_days < tof_lo or tof_days > tof_hi:
            return BIG

        dep_time = Time(dep_jd, format="jd", scale="tdb").tdb
        dv = delta_v_objective(dep_time, tof_days, dv_cutoff=dv_cutoff, return_components=False)
        if (dv is None) or (not np.isfinite(dv)):
            return BIG
        return float(dv)

    return f


def run_optimizer_with_history(
    dep_start, dep_end,
    min_tof_days, max_tof_days,
    dv_cutoff=DV_CUTOFF,
    method="de",
    seed=1,
):
    if not _HAS_SCIPY:
        raise RuntimeError("SciPy not available")

    t0 = Time(dep_start, scale="tdb").tdb
    t1 = Time(dep_end, scale="tdb").tdb

    dep_bounds = (float(t0.jd), float(t1.jd))
    tof_bounds = (float(min_tof_days), float(max_tof_days))

    f = _continuous_objective_factory(dv_cutoff, dep_bounds, tof_bounds)
    history = []

    x0 = np.array([(dep_bounds[0] + dep_bounds[1]) / 2.0, (tof_bounds[0] + tof_bounds[1]) / 2.0], dtype=float)

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
    dv_best = float(delta_v_objective(dep_time_best, tof_best, dv_cutoff=dv_cutoff))

    return {
        "method": method,
        "success": bool(getattr(res, "success", True)),
        "message": str(getattr(res, "message", "")),
        "nfev": int(getattr(res, "nfev", -1)) if getattr(res, "nfev", None) is not None else -1,
        "dep_time": dep_time_best,
        "tof_days": tof_best,
        "arr_time": arr_time_best,
        "dv_kms": dv_best,
        "history": history,
    }


def plot_porkchop_with_trajectory(dep_times, tofs_days, dv_grid, global_min, opt_result, title_suffix="opt"):
    os.makedirs("plots", exist_ok=True)

    dep_datetimes = dep_times.to_datetime()
    dep_mpl = mdates.date2num(dep_datetimes)

    D_dep, D_tof = np.meshgrid(dep_mpl, tofs_days, indexing="ij")
    Z = np.ma.masked_invalid(dv_grid)

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

    g_dep = mdates.date2num(global_min["dep_time"].to_datetime())
    g_tof = float(global_min["tof_days"])
    plt.scatter([g_dep], [g_tof], marker="*", s=180, label=f"Grid min: {global_min['dv_kms']:.3f} km/s")

    hist = opt_result.get("history", [])
    if len(hist) > 0:
        hist_dep_dt = [Time(h[0], format="jd", scale="tdb").to_datetime() for h in hist]
        hist_dep_mpl = mdates.date2num(hist_dep_dt)
        hist_tof = [h[1] for h in hist]
        plt.plot(hist_dep_mpl, hist_tof, linewidth=2, label=f"Opt path ({opt_result['method']})")
        plt.scatter([hist_dep_mpl[-1]], [hist_tof[-1]], marker="o", s=60)

    o_dep = mdates.date2num(opt_result["dep_time"].to_datetime())
    o_tof = float(opt_result["tof_days"])
    plt.scatter([o_dep], [o_tof], marker="X", s=120, label=f"Opt min: {opt_result['dv_kms']:.3f} km/s")

    plt.title(rf"Pork-Chop + Opt Path ({title_suffix})")
    plt.legend(loc="best")
    plt.tight_layout()

    save_path = f"plots/porkchop_{title_suffix}_optpath.png"
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot: {save_path}")

def prepare_ml_npz_from_csv(csv_path, out_npz="ml_data_all.npz", out_norm_json="normalization_stats.json"):
    df = pd.read_csv(csv_path)
    df_clean = df.dropna(subset=["dv_total_kms"]).copy()

    X = df_clean[["dep_julian", "tof_days"]].values.astype(np.float32)
    y = df_clean["dv_total_kms"].values.astype(np.float32).reshape(-1, 1)

    dep_julian = df_clean["dep_julian"].values.astype(np.float64)
    dep_year = df_clean["dep_iso"].str.slice(0, 4).astype(np.int16).values

    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std == 0.0] = 1.0
    X_norm = (X - X_mean) / X_std

    with open(out_norm_json, "w") as f:
        json.dump({"X_mean": X_mean.tolist(), "X_std": X_std.tolist()}, f)

    np.savez(
        out_npz,
        X_norm=X_norm,
        y=y,
        dep_julian=dep_julian,
        dep_year=dep_year,
    )

    print(f"Saved: {out_norm_json}")
    print(f"Saved: {out_npz}")
    print(f"Rows: {len(df_clean)} / {len(df)}")
    return out_npz


def parse_args():
    p = argparse.ArgumentParser(description="Earth->Mars porkchop generator")

    p.add_argument("--dep-start", type=str, default="2030-01-01")
    p.add_argument("--dep-end", type=str, default="2032-01-01")
    p.add_argument("--arr-start", type=str, default=None)
    p.add_argument("--arr-end", type=str, default=None)

    p.add_argument("--dep-step", type=int, default=1)
    p.add_argument("--min-tof", type=float, default=120.0)
    p.add_argument("--max-tof", type=float, default=360.0)
    p.add_argument("--n-tof", type=int, default=60)

    p.add_argument("--dv-cutoff", type=float, default=DV_CUTOFF)
    p.add_argument("--csv", type=str, default=DATA_CSV_PATH)
    p.add_argument("--n-jobs", type=int, default=N_JOBS)

    p.add_argument("--no-generate", action="store_true")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--opt-plot", action="store_true")
    p.add_argument("--opt-method", type=str, default="de", choices=["de", "lbfgsb"])
    p.add_argument("--opt-seed", type=int, default=1)
    p.add_argument("--prep-ml", action="store_true")
    p.add_argument("--title", type=str, default="run")

    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    dep_times = None
    tofs_days = None
    dv_grid = None

    if not args.no_generate:
        _, dep_times, tofs_days, dv_grid = generate_porkchop_data_daily(
            dep_start=args.dep_start,
            dep_end=args.dep_end,
            arr_start=args.arr_start,
            arr_end=args.arr_end,
            departure_step_days=args.dep_step,
            min_tof_days=args.min_tof,
            max_tof_days=args.max_tof,
            n_tof=args.n_tof,
            csv_path=args.csv,
            dv_cutoff=args.dv_cutoff,
            n_jobs=args.n_jobs,
        )

    if args.plot:
        if dep_times is None or tofs_days is None or dv_grid is None:
            df_plot = pd.read_csv(args.csv)
            dep_times, tofs_days, dv_grid = rebuild_grid_from_csv(df_plot)
        plot_porkchop(dep_times, tofs_days, dv_grid, title_suffix=args.title)

    if args.opt_plot:
        if dep_times is None or tofs_days is None or dv_grid is None:
            df_plot = pd.read_csv(args.csv)
            dep_times, tofs_days, dv_grid = rebuild_grid_from_csv(df_plot)

        global_min = find_global_min_from_grid(dep_times, tofs_days, dv_grid)
        print(f"Grid min dep: {global_min['dep_time'].iso}")
        print(f"Grid min tof: {global_min['tof_days']:.6f}")
        print(f"Grid min arr: {global_min['arr_time'].iso}")
        print(f"Grid min dv:  {global_min['dv_kms']:.6f}")

        opt_result = run_optimizer_with_history(
            dep_start=args.dep_start,
            dep_end=args.dep_end,
            min_tof_days=args.min_tof,
            max_tof_days=args.max_tof,
            dv_cutoff=args.dv_cutoff,
            method=args.opt_method,
            seed=args.opt_seed,
        )

        print(f"Opt min dep:  {opt_result['dep_time'].iso}")
        print(f"Opt min tof:  {opt_result['tof_days']:.6f}")
        print(f"Opt min arr:  {opt_result['arr_time'].iso}")
        print(f"Opt min dv:   {opt_result['dv_kms']:.6f}")
        print(f"History pts:  {len(opt_result['history'])}")

        plot_porkchop_with_trajectory(
            dep_times, tofs_days, dv_grid,
            global_min=global_min,
            opt_result=opt_result,
            title_suffix=args.title,
        )

    if args.prep_ml:
        prepare_ml_npz_from_csv(args.csv, out_npz="ml_data_all.npz", out_norm_json="normalization_stats.json")
