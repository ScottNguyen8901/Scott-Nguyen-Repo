from __future__ import annotations

import csv
import math
import random
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import matplotlib.pyplot as plt

from filterpy.kalman import UnscentedKalmanFilter as UKF
from filterpy.kalman import MerweScaledSigmaPoints

try:
    from tqdm.auto import tqdm
except Exception:
    tqdm = None


# =========================
# Files
# =========================
FULL_TRAJ_CSV = "iss_trajectory_full_1hz.csv"
PASS_MEAS_CSV = "iss_pass_measurements_1hz.csv"

# =========================
# Sensor / noise
# =========================
LAT_DEG = 35.0844
LON_DEG = -106.6504
ELEV_M = 0.0

SIGMA_AZ_DEG = 0.01
SIGMA_EL_DEG = 0.01
RNG_SEED = 7

# =========================
# UKF / dynamics
# =========================
MU_E = 398600.4418  # km^3/s^2

ALPHA = 1e-3
BETA = 2.0
KAPPA = 0.0

# process noise per second (scaled by dt)
Q_POS_KM2_PER_S = 1e-5
Q_VEL_KM2_S2_PER_S = 1e-5

# init from truth at first measurement time (angles-only ambiguity)
INIT_POS_SIGMA_KM = 10.0
INIT_VEL_SIGMA_KM_S = 0.05

# speed-up when not visible (used only in AFTER window, if any)
COAST_STEP_S = 10.0

# pass selection
PASS_CHOICE = 0
PASS_AFTER_S = 0.0  # set to 30.0, etc.


# =========================
# Time parsing
# =========================
def parse_time_utc(s: str) -> datetime:
    s = s.strip()
    if s.endswith(" UTC"):
        s = s[:-4]
    dt = datetime.strptime(s, "%Y-%m-%d %H:%M:%S")
    return dt.replace(tzinfo=timezone.utc)


# =========================
# Earth/site geometry (minimal)
# =========================
def julian_date(dt_utc: datetime) -> float:
    dt = dt_utc.astimezone(timezone.utc)
    y, m = dt.year, dt.month
    D = dt.day + (dt.hour + (dt.minute + dt.second / 60.0) / 60.0) / 24.0
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + (A // 4)
    JD = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + D + B - 1524.5
    return JD


def gmst_rad(dt_utc: datetime) -> float:
    JD = julian_date(dt_utc)
    T = (JD - 2451545.0) / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (JD - 2451545.0)
        + 0.000387933 * T * T
        - (T * T * T) / 38710000.0
    ) % 360.0
    return math.radians(gmst_deg)


def rot_z(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[ c,  s, 0.0],
                     [-s,  c, 0.0],
                     [0.0, 0.0, 1.0]], dtype=float)


def wgs84_site_ecef_km(lat_deg: float, lon_deg: float, elev_m: float) -> np.ndarray:
    a = 6378.137  # km
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)

    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    h = elev_m / 1000.0

    s = math.sin(lat)
    c = math.cos(lat)
    N = a / math.sqrt(1.0 - e2 * s * s)

    x = (N + h) * c * math.cos(lon)
    y = (N + h) * c * math.sin(lon)
    z = (N * (1.0 - e2) + h) * s
    return np.array([x, y, z], dtype=float)


def ecef_to_enu_matrix(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    return np.array([
        [-slon,       clon,      0.0],
        [-slat*clon, -slat*slon, clat],
        [ clat*clon,  clat*slon, slat],
    ], dtype=float)


def wrap_to_pi(rad: np.ndarray) -> np.ndarray:
    return (rad + np.pi) % (2.0 * np.pi) - np.pi


# =========================
# Two-body dynamics (RK4)
# =========================
def f_2body(x: np.ndarray) -> np.ndarray:
    r = x[0:3]
    v = x[3:6]
    rn = np.linalg.norm(r)
    a = -(MU_E / (rn**3)) * r
    return np.hstack((v, a))


def rk4_step(x: np.ndarray, dt: float) -> np.ndarray:
    k1 = f_2body(x)
    k2 = f_2body(x + 0.5 * dt * k1)
    k3 = f_2body(x + 0.5 * dt * k2)
    k4 = f_2body(x + dt * k3)
    return x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def propagate(x: np.ndarray, dt: float) -> np.ndarray:
    if abs(dt) < 1e-12:
        return x.copy()

    n = int(round(dt))
    if abs(dt - n) < 1e-9 and n > 0:
        xx = x.copy()
        for _ in range(n):
            xx = rk4_step(xx, 1.0)
        return xx

    return rk4_step(x, dt)


# =========================
# PD-safe covariance
# =========================
def make_pd(P: np.ndarray, min_jitter: float = 1e-14) -> np.ndarray:
    P = 0.5 * (P + P.T)
    try:
        np.linalg.cholesky(P)
        return P
    except np.linalg.LinAlgError:
        pass

    w = np.linalg.eigvalsh(P)
    min_eig = float(np.min(w))
    jitter = max(min_jitter, -min_eig + min_jitter)
    P2 = P + jitter * np.eye(P.shape[0])

    for _ in range(6):
        try:
            np.linalg.cholesky(P2)
            return P2
        except np.linalg.LinAlgError:
            jitter *= 10.0
            P2 = P + jitter * np.eye(P.shape[0])

    return P + (1e-6) * np.eye(P.shape[0])


# =========================
# Measurement prediction (NO GLOBALS)
# =========================
def z_from_sigmas(sigmas: np.ndarray, R_step: np.ndarray, site_ecef: np.ndarray, ecef2enu: np.ndarray) -> np.ndarray:
    r_eci = sigmas[:, 0:3]
    r_ecef = (r_eci @ R_step.T)
    rho = r_ecef - site_ecef[None, :]
    enu = rho @ ecef2enu.T
    e = enu[:, 0]
    n = enu[:, 1]
    u = enu[:, 2]
    az = np.arctan2(e, n)
    az = np.where(az < 0.0, az + 2.0*np.pi, az)
    el = np.arctan2(u, np.hypot(e, n))
    return np.column_stack((az, el))


def ukf_update_vectorized(
    ukf: UKF,
    z: np.ndarray,
    R_step: np.ndarray,
    site_ecef: np.ndarray,
    ecef2enu: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ukf.P = make_pd(ukf.P)

    sigmas = ukf.points_fn.sigma_points(ukf.x, ukf.P)
    Zsig = z_from_sigmas(sigmas, R_step, site_ecef, ecef2enu)

    Wm = ukf.Wm
    Wc = ukf.Wc

    sin_az = np.sin(Zsig[:, 0])
    cos_az = np.cos(Zsig[:, 0])
    az_mean = math.atan2(float(np.dot(Wm, sin_az)), float(np.dot(Wm, cos_az)))
    if az_mean < 0.0:
        az_mean += 2.0*np.pi
    el_mean = float(np.dot(Wm, Zsig[:, 1]))
    z_mean = np.array([az_mean, el_mean], dtype=float)

    dz = Zsig - z_mean[None, :]
    dz[:, 0] = wrap_to_pi(dz[:, 0])

    dx = sigmas - ukf.x[None, :]

    S = (dz.T * Wc) @ dz + ukf.R
    Pxz = (dx.T * Wc) @ dz

    K = Pxz @ np.linalg.solve(S, np.eye(2))

    y = z - z_mean
    y[0] = float(wrap_to_pi(np.array([y[0]]))[0])

    ukf.x = ukf.x + K @ y
    ukf.P = make_pd(ukf.P - K @ S @ K.T)

    return y, z_mean


# =========================
# CSV helpers
# =========================
def load_pass_measurements(pass_csv: Path) -> dict[int, tuple[float, float]]:
    d: dict[int, tuple[float, float]] = {}
    with pass_csv.open("r", newline="", encoding="utf-8") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        idx = {name: i for i, name in enumerate(header)}

        t_i = idx["time_since_epoch_s"]
        az_i = idx["az_deg"]
        el_i = idx["el_deg"]
        vis_i = idx["visible"]

        for row in rdr:
            if int(row[vis_i]) != 1:
                continue
            t = int(round(float(row[t_i])))
            d[t] = (float(row[az_i]), float(row[el_i]))
    return d


def detect_passes_from_times(times_sorted: np.ndarray) -> list[tuple[int, int]]:
    if len(times_sorted) == 0:
        return []
    passes: list[tuple[int, int]] = []
    start = int(times_sorted[0])
    prev = int(times_sorted[0])
    for t in times_sorted[1:]:
        t = int(t)
        if t == prev + 1:
            prev = t
            continue
        passes.append((start, prev))
        start = t
        prev = t
    passes.append((start, prev))
    return passes


# =========================
# Main
# =========================
def main():
    random.seed(RNG_SEED)
    np.random.seed(RNG_SEED)

    full_path = Path(FULL_TRAJ_CSV)
    pass_path = Path(PASS_MEAS_CSV)

    if not full_path.exists():
        raise FileNotFoundError(f"Missing {FULL_TRAJ_CSV}")
    if not pass_path.exists():
        raise FileNotFoundError(f"Missing {PASS_MEAS_CSV}")

    pass_meas = load_pass_measurements(pass_path)
    meas_times = np.array(sorted(pass_meas.keys()), dtype=int)
    passes = detect_passes_from_times(meas_times)

    if len(passes) == 0:
        raise RuntimeError("No passes found in PASS_MEAS_CSV.")

    if PASS_CHOICE < 0 or PASS_CHOICE >= len(passes):
        raise ValueError(f"PASS_CHOICE={PASS_CHOICE} out of range (found {len(passes)} passes).")

    t_start, t_end = passes[PASS_CHOICE]
    t_stop = int(round(t_end + PASS_AFTER_S))

    print(f"Using pass #{PASS_CHOICE}: t={t_start}..{t_end} plus {PASS_AFTER_S:.0f}s -> stop at t={t_stop}")

    # measurement covariance (rad^2)
    Rmeas = np.diag([
        math.radians(SIGMA_AZ_DEG) ** 2,
        math.radians(SIGMA_EL_DEG) ** 2,
    ])

    site_ecef = wgs84_site_ecef_km(LAT_DEG, LON_DEG, ELEV_M)
    ecef2enu = ecef_to_enu_matrix(LAT_DEG, LON_DEG)

    # Load truth states + time_utc strings for this window
    times: list[int] = []
    timeutc: list[str] = []
    X_true: list[np.ndarray] = []
    Z_meas: list[np.ndarray] = []
    visible: list[bool] = []

    with full_path.open("r", newline="", encoding="utf-8") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        idx = {name: i for i, name in enumerate(header)}

        t_i = idx["time_since_epoch_s"]
        timeutc_i = idx["time_utc"]
        rx_i = idx["r_x_km"]; ry_i = idx["r_y_km"]; rz_i = idx["r_z_km"]
        vx_i = idx["v_x_km_s"]; vy_i = idx["v_y_km_s"]; vz_i = idx["v_z_km_s"]

        for row in rdr:
            t_int = int(round(float(row[t_i])))
            if t_int < t_start:
                continue
            if t_int > t_stop:
                break

            r = np.array([float(row[rx_i]), float(row[ry_i]), float(row[rz_i])], dtype=float)
            v = np.array([float(row[vx_i]), float(row[vy_i]), float(row[vz_i])], dtype=float)
            xtruth = np.hstack((r, v))

            meas_avail = (t_int in pass_meas)
            if meas_avail:
                az_deg, el_deg = pass_meas[t_int]
                az_meas = az_deg + random.gauss(0.0, SIGMA_AZ_DEG)
                el_meas = el_deg + random.gauss(0.0, SIGMA_EL_DEG)
                z = np.array([math.radians(az_meas), math.radians(el_meas)], dtype=float)
            else:
                z = np.array([np.nan, np.nan], dtype=float)

            times.append(t_int)
            timeutc.append(row[timeutc_i])
            X_true.append(xtruth)
            Z_meas.append(z)
            visible.append(bool(meas_avail))

    times = np.array(times, dtype=int)
    X_true = np.vstack(X_true)
    Z_meas = np.vstack(Z_meas)
    visible = np.array(visible, dtype=bool)

    if len(times) < 2:
        raise RuntimeError("Window too small. Check PASS_CHOICE / PASS_AFTER_S.")

    # Find first measurement inside window (init like your generator)
    first_meas_idx = int(np.argmax(visible))
    if not visible[first_meas_idx]:
        raise RuntimeError("No measurements in selected window.")

    # Initialize at first measurement epoch (matched)
    x0 = X_true[first_meas_idx].copy()
    x0[0:3] += np.random.normal(0.0, INIT_POS_SIGMA_KM, size=3)
    x0[3:6] += np.random.normal(0.0, INIT_VEL_SIGMA_KM_S, size=3)

    P0 = np.zeros((6, 6), dtype=float)
    P0[0:3, 0:3] = (INIT_POS_SIGMA_KM**2) * np.eye(3)
    P0[3:6, 3:6] = (INIT_VEL_SIGMA_KM_S**2) * np.eye(3)
    P0 = make_pd(P0)

    points = MerweScaledSigmaPoints(n=6, alpha=ALPHA, beta=BETA, kappa=KAPPA)

    ukf = UKF(dim_x=6, dim_z=2, fx=lambda x, dt: propagate(x, dt), hx=lambda x: x[:2], dt=1.0, points=points)
    ukf.x = x0
    ukf.P = P0
    ukf.R = Rmeas
    ukf.Q = np.zeros((6, 6), dtype=float)

    # Prepare outputs (start from first measurement index)
    times_run = times[first_meas_idx:]
    timeutc_run = timeutc[first_meas_idx:]
    X_true_run = X_true[first_meas_idx:]
    Z_meas_run = Z_meas[first_meas_idx:]
    vis_run = visible[first_meas_idx:]

    N = len(times_run)
    X_hat = np.zeros((N, 6), dtype=float)
    P_diag = np.zeros((N, 6), dtype=float)

    X_hat[0] = ukf.x
    P_diag[0] = np.diag(ukf.P)

    t_plot = (times_run - times_run[0]).astype(float)

    pbar = None
    total_sim = float(t_plot[-1] - t_plot[0])
    if tqdm is not None and total_sim > 0:
        pbar = tqdm(total=total_sim, desc="Sim time (matched)", unit="s", dynamic_ncols=True, smoothing=0.0)

    accum_dt = 0.0

    # Run filter
    for k in range(1, N):
        dt = float(times_run[k] - times_run[k - 1])
        if dt <= 0.0:
            X_hat[k] = ukf.x
            P_diag[k] = np.diag(ukf.P)
            continue

        # rotation for this epoch
        dt_utc = parse_time_utc(timeutc_run[k])
        R_step = rot_z(gmst_rad(dt_utc))

        if vis_run[k]:
            # process noise scaled by dt
            Q = np.zeros((6, 6), dtype=float)
            Q[0:3, 0:3] = (Q_POS_KM2_PER_S * dt) * np.eye(3)
            Q[3:6, 3:6] = (Q_VEL_KM2_S2_PER_S * dt) * np.eye(3)
            ukf.Q = Q

            ukf.P = make_pd(ukf.P)
            ukf.predict(dt=dt)
            ukf.P = make_pd(ukf.P)

            ukf_update_vectorized(ukf, Z_meas_run[k], R_step, site_ecef, ecef2enu)
            accum_dt = 0.0
        else:
            accum_dt += dt
            if accum_dt >= COAST_STEP_S:
                Qc = np.zeros((6, 6), dtype=float)
                Qc[0:3, 0:3] = (Q_POS_KM2_PER_S * accum_dt) * np.eye(3)
                Qc[3:6, 3:6] = (Q_VEL_KM2_S2_PER_S * accum_dt) * np.eye(3)
                ukf.Q = Qc

                ukf.P = make_pd(ukf.P)
                ukf.predict(dt=accum_dt)
                ukf.P = make_pd(ukf.P)

                accum_dt = 0.0

        X_hat[k] = ukf.x
        P_diag[k] = np.diag(ukf.P)

        if pbar is not None:
            pbar.update(dt)

    if pbar is not None:
        pbar.close()

    # Plot error + shaded ±3σ
    err = X_true_run - X_hat
    sig3 = 3.0 * np.sqrt(np.maximum(P_diag, 0.0))

    fig, axs = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
    labels = ["x", "y", "z"]

    for i in range(3):
        e_r = err[:, i]
        b_r = sig3[:, i]
        axs[i, 0].plot(t_plot, e_r, label="error")
        axs[i, 0].fill_between(t_plot, -b_r, b_r, alpha=0.25, label="±3σ")
        axs[i, 0].set_ylabel(f"r_{labels[i]} err [km]")

        e_v = err[:, i + 3]
        b_v = sig3[:, i + 3]
        axs[i, 1].plot(t_plot, e_v, label="error")
        axs[i, 1].fill_between(t_plot, -b_v, b_v, alpha=0.25, label="±3σ")
        axs[i, 1].set_ylabel(f"v_{labels[i]} err [km/s]")

    for ax in axs.flat:
        ax.set_xlabel("Time [s]")
        ax.grid(True)

    handles, leglabels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, leglabels, loc="upper right")

    plt.suptitle(f"Matched UKF Errors with 3σ Bands — Pass #{PASS_CHOICE} + {PASS_AFTER_S:.0f}s After", fontsize=14)
    plt.show()


if __name__ == "__main__":
    main()
