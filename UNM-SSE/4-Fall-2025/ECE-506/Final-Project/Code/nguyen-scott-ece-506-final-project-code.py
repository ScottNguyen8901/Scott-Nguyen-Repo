import os
import sys
import atexit
from datetime import datetime
import pandas as pd
from tqdm import tqdm

from porkchop_run import (
    generate_porkchop_data_daily,
    plot_porkchop,
    find_global_min_from_grid,
    run_optimizer_with_history,
    plot_porkchop_with_trajectory,
    prepare_ml_npz_from_csv,
)

from train_surrogate import train_from_npz
from evaluate_surrogate import evaluate_surrogate


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            try:
                s.write(data)
                s.flush()
            except Exception:
                pass

    def flush(self):
        for s in self.streams:
            try:
                s.flush()
            except Exception:
                pass


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def combo_tag(dep_start, dep_end, arr_start, arr_end):
    return f"dep{dep_start}_{dep_end}_arr{arr_start}_{arr_end}"


run_tag = datetime.now().strftime("%Y%m%d_%H%M%S")
base_out_dir = os.path.abspath(f"runs_{run_tag}")
ensure_dir(base_out_dir)

log_name = os.path.abspath(f"run_log_{run_tag}.txt")
log_file = open(log_name, "w", encoding="utf-8")

sys.stdout = Tee(sys.__stdout__, log_file)
sys.stderr = Tee(sys.__stderr__, log_file)

print(f"Logging console output to: {log_name}")


def _close_log():
    try:
        log_file.flush()
        log_file.close()
    except Exception:
        pass


atexit.register(_close_log)


min_tof = 120.0
max_tof = 360.0
n_tof = 60
dv_cutoff = 20.0

years = [2030, 2032, 2034, 2036, 2038]

for y0 in tqdm(years, desc="Running porkchop + surrogate experiments", unit="window", ncols=100):
    dep_start = f"{y0}-01-01"
    dep_end = f"{y0 + 2}-01-01"
    arr_start = f"{y0}-05-01"
    arr_end = f"{y0 + 2}-12-26"

    tag = combo_tag(dep_start, dep_end, arr_start, arr_end)
    out_dir = os.path.abspath(os.path.join(base_out_dir, tag))
    ensure_dir(out_dir)

    print("\n" + "=" * 80)
    print(f"RUN COMBO: {tag}")
    print("=" * 80)

    csv_name = f"porkchop_{tag}.csv"
    clean_csv_name = f"porkchop_{tag}_clean.csv"
    npz_name = f"ml_data_{tag}.npz"
    norm_name = f"normalization_stats_{tag}.json"
    model_name = f"dv_surrogate_model_{tag}.pt"

    csv_path = os.path.join(out_dir, csv_name)
    clean_csv_path = os.path.join(out_dir, clean_csv_name)
    npz_path = os.path.join(out_dir, npz_name)
    norm_path = os.path.join(out_dir, norm_name)
    model_path = os.path.join(out_dir, model_name)

    df, dep_times, tofs_days, dv_grid = generate_porkchop_data_daily(
        dep_start=dep_start,
        dep_end=dep_end,
        arr_start=arr_start,
        arr_end=arr_end,
        min_tof_days=min_tof,
        max_tof_days=max_tof,
        n_tof=n_tof,
        csv_path=csv_path,
        dv_cutoff=dv_cutoff,
        n_jobs=-1,
    )

    old_cwd = os.getcwd()
    os.chdir(out_dir)
    try:
        plot_porkchop(dep_times, tofs_days, dv_grid, title_suffix=f"{dep_start}_{dep_end}")
    finally:
        os.chdir(old_cwd)

    global_min = find_global_min_from_grid(dep_times, tofs_days, dv_grid)
    print("GRID MIN")
    print(global_min["dep_time"].iso)
    print(f"{global_min['tof_days']:.6f}")
    print(global_min["arr_time"].iso)
    print(f"{global_min['dv_kms']:.6f}")

    opt_result = run_optimizer_with_history(
        dep_start=dep_start,
        dep_end=dep_end,
        min_tof_days=min_tof,
        max_tof_days=max_tof,
        dv_cutoff=dv_cutoff,
        method="de",
        seed=1,
    )

    print("OPT MIN")
    print(opt_result["dep_time"].iso)
    print(f"{opt_result['tof_days']:.6f}")
    print(opt_result["arr_time"].iso)
    print(f"{opt_result['dv_kms']:.6f}")
    print(len(opt_result["history"]))

    os.chdir(out_dir)
    try:
        plot_porkchop_with_trajectory(
            dep_times,
            tofs_days,
            dv_grid,
            global_min=global_min,
            opt_result=opt_result,
            title_suffix=f"{dep_start}_{dep_end}",
        )
    finally:
        os.chdir(old_cwd)

    if not os.path.exists(clean_csv_path):
        df_raw = pd.read_csv(csv_path)
        df_clean = df_raw.dropna(subset=["dv_total_kms"]).copy()
        df_clean.to_csv(clean_csv_path, index=False)
        print(f"Saved clean CSV: {clean_csv_path}")
    else:
        print(f"Clean CSV exists: {clean_csv_path}")

    prepare_ml_npz_from_csv(csv_path, out_npz=npz_path, out_norm_json=norm_path)
    train_from_npz(npz_path, save_path=model_path)

    os.chdir(out_dir)
    try:
        evaluate_surrogate(
            dep_start=dep_start,
            dep_end=dep_end,
            min_tof_days=min_tof,
            max_tof_days=max_tof,
            clean_csv_path=clean_csv_name,
            model_checkpoint=model_name,
            norm_stats_path=norm_name,
            true_opt_result=opt_result,
            dv_cutoff=dv_cutoff,
            surrogate_opt_method="de",
            surrogate_opt_seed=2,
        )
    finally:
        os.chdir(old_cwd)

print("\nALL COMBINATIONS COMPLETE")
print(f"Results saved under: {base_out_dir}")
