import time, warnings, random, os
import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy import units as u
from astropy.constants import G, M_earth, R_earth
from poliastro.bodies import Sun, Earth, Mars
from poliastro.ephem import get_body_ephem
from poliastro import iod
from scipy.optimize import minimize
from numpy.linalg import norm
from tqdm import tqdm

warnings.filterwarnings("ignore")

# =========================
# Constants
# =========================
mu_sun = Sun.k.to(u.km**3 / u.s**2).value
mu_earth = (G * M_earth).to(u.km**3 / u.s**2).value
R_E = R_earth.to(u.km).value
alt_earth = 300
mu_mars = Mars.k.to(u.km**3 / u.s**2).value
R_mars = Mars.R.to(u.km).value
alt_mars = 300

JD_DEP_MIN = Time("2025-01-01", scale="tdb").jd
JD_DEP_MAX = Time("2030-01-01", scale="tdb").jd
JD_ARR_MAX = Time("2035-01-01", scale="tdb").jd
TOF_MIN_DAYS = 150
TOF_MAX_DAYS = 305


# =========================
# Helper functions
# =========================
def get_state(body: str, jd: float):
    t = Time(jd, format="jd", scale="tdb")
    pos, vel = get_body_ephem(body, t)
    return pos.to(u.km).value, vel.to(u.km / u.s).value


def earth_mars_total_delta_v(jd_start: float, jd_end: float) -> float:
    """Compute total delta-V for Earth–Mars transfer."""
    if jd_end <= jd_start:
        return np.inf
    r_earth, v_earth = get_state("earth", jd_start)
    r_mars, v_mars = get_state("mars", jd_end)
    tof = (jd_end - jd_start) * 86400.0
    try:
        (v_dep, v_arr), = iod.lambert(Sun.k, r_earth * u.km, r_mars * u.km, tof * u.s)
    except Exception:
        return np.inf
    v_dep = v_dep.to(u.km / u.s).value
    v_arr = v_arr.to(u.km / u.s).value
    v_inf_dep = norm(v_dep - v_earth)
    v_inf_arr = norm(v_arr - v_mars)
    r_pE, r_pM = R_E + alt_earth, R_mars + alt_mars
    dv_dep = abs(np.sqrt(v_inf_dep**2 + 2*mu_earth/r_pE) - np.sqrt(mu_earth/r_pE))
    dv_arr = abs(np.sqrt(v_inf_arr**2 + 2*mu_mars/r_pM) - np.sqrt(mu_mars/r_pM))
    return dv_dep + dv_arr


def total_delta_v_opt(x):
    jd_start, jd_end = x
    return earth_mars_total_delta_v(jd_start, jd_end)


def constraints(x):
    jd_start, jd_end = x
    tof = jd_end - jd_start
    return np.array([
        tof - TOF_MIN_DAYS,
        jd_start - JD_DEP_MIN,
        JD_DEP_MAX - jd_start,
        TOF_MAX_DAYS - tof
    ])


# =========================
# Optimizer
# =========================
def run_optimizer(num_guesses=6):
    """Randomized initial guesses."""
    jd_end_lower = JD_DEP_MIN + TOF_MIN_DAYS
    bounds = [(JD_DEP_MIN, JD_DEP_MAX), (jd_end_lower, JD_ARR_MAX)]
    nonlin = {"type": "ineq", "fun": lambda x: constraints(x)}
    random.seed(42)
    initial_guesses = []
    for _ in range(num_guesses):
        jd_start_guess = random.uniform(JD_DEP_MIN, JD_DEP_MAX)
        tof_guess = random.uniform(TOF_MIN_DAYS, 600)
        jd_end_guess = min(jd_start_guess + tof_guess, JD_ARR_MAX)
        initial_guesses.append([jd_start_guess, jd_end_guess])

    results = []
    for x0 in initial_guesses:
        res = minimize(total_delta_v_opt, x0, method="SLSQP", bounds=bounds, constraints=[nonlin])
        if res.success:
            jd_start_opt, jd_end_opt = res.x
            results.append({
                "initial_guess": (x0[0], x0[1]),
                "dep_jd": jd_start_opt, "arr_jd": jd_end_opt,
                "total_delta_v": res.fun
            })
        else:
            results.append({
                "initial_guess": (x0[0], x0[1]),
                "dep_jd": None, "arr_jd": None, "total_delta_v": None
            })
    return results, initial_guesses


def print_results_table(results):
    print("\n=== Optimizer Results Table ===")
    print(f"{'Initial Guess (Dep,Arr)':>40} | {'Optimized Dep':>15} | {'Optimized Arr':>15} | {'ΔV (km/s)':>10}")
    print("-" * 95)
    for r in results:
        ig = r["initial_guess"]
        ig_dep = Time(ig[0], format="jd").to_datetime().strftime("%Y-%m-%d")
        ig_arr = Time(ig[1], format="jd").to_datetime().strftime("%Y-%m-%d")
        if r["dep_jd"]:
            dep = Time(r["dep_jd"], format="jd").to_datetime().strftime("%Y-%m-%d")
            arr = Time(r["arr_jd"], format="jd").to_datetime().strftime("%Y-%m-%d")
            dv = f"{r['total_delta_v']:.3f}"
        else:
            dep, arr, dv = "Failed", "Failed", "—"
        print(f"{f'({ig_dep},{ig_arr})':>40} | {dep:>15} | {arr:>15} | {dv:>10}")


# =========================
# Porkchop Plot
# =========================
def build_porkchop(n=60):
    dep_dates = np.linspace(JD_DEP_MIN, JD_DEP_MAX, n)
    arr_dates = np.linspace(JD_DEP_MIN + TOF_MIN_DAYS, JD_ARR_MAX, n)
    dv_matrix = np.full((n, n), np.nan)
    for i, t_dep in enumerate(tqdm(dep_dates, desc="Porkchop grid")):
        for j, t_arr in enumerate(arr_dates):
            if t_arr <= t_dep:
                continue
            dv_matrix[i, j] = earth_mars_total_delta_v(t_dep, t_arr)
    i_min, j_min = np.unravel_index(np.nanargmin(dv_matrix), dv_matrix.shape)
    return dep_dates, arr_dates, dv_matrix, dv_matrix[i_min, j_min], dep_dates[i_min], arr_dates[j_min]


def plot_porkchop(dep, arr, dv, min_dv, min_dep, min_arr,
                  results=None, guesses=None, clean=False, savepath=None):
    dep_grid, arr_grid = np.meshgrid(dep, arr, indexing="ij")
    fig, ax = plt.subplots(figsize=(12, 9))
    valid = dv[~np.isnan(dv)]
    levels = np.linspace(np.nanmin(valid), np.nanpercentile(valid, 99), 40)
    cp = ax.contourf(dep_grid, arr_grid, dv, levels=levels, cmap="viridis")
    cbar = fig.colorbar(cp)
    cbar.set_label("Total ΔV (km/s)")
    ax.set_xlabel("Departure (JD)")
    ax.set_ylabel("Arrival (JD)")
    title = "Earth → Mars Porkchop Plot" + ("" if clean else " (Annotated)")
    ax.set_title(title)

    if not clean:
        ax.plot(min_dep, min_arr, "r*", markersize=14, label=f"Min {min_dv:.2f} km/s")
        if guesses:
            g_dep = [g[0] for g in guesses]
            g_arr = [g[1] for g in guesses]
            ax.scatter(g_dep, g_arr, c="cyan", edgecolors="k", s=60, label="Initial guesses")
        if results:
            opt_dep = [r["dep_jd"] for r in results if r["dep_jd"]]
            opt_arr = [r["arr_jd"] for r in results if r["arr_jd"]]
            ax.scatter(opt_dep, opt_arr, c="#e600ff", s=100, marker="*", label="Optimizer solutions")
            for r, (ig_dep, ig_arr) in zip(results, guesses):
                if r["dep_jd"]:
                    ax.plot([ig_dep, r["dep_jd"]], [ig_arr, r["arr_jd"]], "w--", lw=1.2, alpha=0.8)
        ax.legend()

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=300)
        print(f"Saved: {savepath}")
    plt.show()


# =========================
# Main
# =========================
if __name__ == "__main__":
    t0 = time.time()
    print("=== Running Randomized SLSQP Optimizer ===")
    results, guesses = run_optimizer()
    print_results_table(results)

    print("\n=== Building Porkchop Grid ===")
    dep, arr, mat, min_dv, min_dep, min_arr = build_porkchop(60)
    print(f"\nPorkchop Minimum ΔV = {min_dv:.3f} km/s at "
          f"Dep={Time(min_dep,format='jd').to_datetime().strftime('%Y-%m-%d')} "
          f"Arr={Time(min_arr,format='jd').to_datetime().strftime('%Y-%m-%d')}")

    out_dir = os.getcwd()
    annotated_path = os.path.join(out_dir, "porkchop_annotated.png")
    clean_path = os.path.join(out_dir, "porkchop.png")

    plot_porkchop(dep, arr, mat, min_dv, min_dep, min_arr,
                  results, guesses, clean=False, savepath=annotated_path)
    plot_porkchop(dep, arr, mat, min_dv, min_dep, min_arr,
                  clean=True, savepath=clean_path)

    print(f"\nTotal runtime: {time.time()-t0:.2f} s")
