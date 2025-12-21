from __future__ import annotations

import csv
import math
import random
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
from filterpy.kalman import UnscentedKalmanFilter as UKF
from filterpy.kalman import MerweScaledSigmaPoints

try:
    from tqdm import tqdm
except Exception:
    tqdm = None


# =========================
# Files
# =========================
FULL_TRAJ_CSV = "iss_trajectory_full_1hz.csv"
PASS_MEAS_CSV = "iss_pass_measurements_1hz.csv"
OUTPUT_CSV = "ukf_estimates_angles_only.csv"

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

# process noise per second
Q_POS_KM2_PER_S = 1e-10
Q_VEL_KM2_S2_PER_S = 1e-12

# init from truth at first measurement time (angles-only ambiguity)
INIT_POS_SIGMA_KM = 10.0
INIT_VEL_SIGMA_KM_S = 0.05

# progress bar update cadence (seconds)
PBAR_UPDATE_EVERY_S = 50


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


SITE_ECEF = wgs84_site_ecef_km(LAT_DEG, LON_DEG, ELEV_M)
ECEF2ENU = ecef_to_enu_matrix(LAT_DEG, LON_DEG)


# =========================
# Two-body dynamics (RK4, 1Hz optimized)
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
    # Your timeline is 1Hz; keep this fast.
    if abs(dt) < 1e-12:
        return x.copy()

    # If dt is an integer number of seconds, do dt times 1-second RK4 steps
    # (keeps accuracy and avoids expensive substepping logic).
    n = int(round(dt))
    if abs(dt - n) < 1e-9 and n > 0:
        xx = x.copy()
        for _ in range(n):
            xx = rk4_step(xx, 1.0)
        return xx

    # Otherwise, one RK4 step
    return rk4_step(x, dt)


# FilterPy expects fx(x, dt) and hx(x); we’ll still provide them, but hx will be cheap.
_R_STEP = None  # cached ECI->ECEF rotation at current epoch (3x3)

def fx(x: np.ndarray, dt: float) -> np.ndarray:
    return propagate(x, dt)

def hx(x: np.ndarray) -> np.ndarray:
    # Single-point hx (FilterPy may call this; we avoid using it in the fast update)
    global _R_STEP
    r_ecef = _R_STEP @ x[0:3]
    rho_ecef = r_ecef - SITE_ECEF
    e, n, u = (ECEF2ENU @ rho_ecef)
    az = math.atan2(e, n)
    if az < 0.0:
        az += 2.0 * math.pi
    el = math.atan2(u, math.hypot(e, n))
    return np.array([az, el], dtype=float)


# =========================
# PD-safe covariance
# =========================
def make_pd(P: np.ndarray, min_jitter: float = 1e-14) -> np.ndarray:
    """
    Ensure P is symmetric positive definite (SPD) by:
    - symmetrizing
    - adding diagonal jitter if needed
    """
    P = 0.5 * (P + P.T)

    # quick try
    try:
        np.linalg.cholesky(P)
        return P
    except np.linalg.LinAlgError:
        pass

    # add jitter based on smallest eigenvalue
    w = np.linalg.eigvalsh(P)
    min_eig = float(np.min(w))
    jitter = max(min_jitter, -min_eig + min_jitter)
    P2 = P + jitter * np.eye(P.shape[0])

    # if still not SPD, escalate a bit
    for _ in range(6):
        try:
            np.linalg.cholesky(P2)
            return P2
        except np.linalg.LinAlgError:
            jitter *= 10.0
            P2 = P + jitter * np.eye(P.shape[0])

    # last resort
    return P + (1e-6) * np.eye(P.shape[0])


# =========================
# Vectorized measurement prediction for sigma points
# =========================
def z_from_sigmas(sigmas: np.ndarray) -> np.ndarray:
    """
    sigmas: (Ns, 6)
    returns Z: (Ns, 2) in radians [az, el]
    """
    global _R_STEP
    r_eci = sigmas[:, 0:3]                 # (Ns,3)
    r_ecef = (r_eci @ _R_STEP.T)           # (Ns,3) faster than looping (note transpose)
    rho = r_ecef - SITE_ECEF[None, :]      # (Ns,3)
    enu = rho @ ECEF2ENU.T                 # (Ns,3)
    e = enu[:, 0]
    n = enu[:, 1]
    u = enu[:, 2]
    az = np.arctan2(e, n)
    az = np.where(az < 0.0, az + 2.0*np.pi, az)
    el = np.arctan2(u, np.hypot(e, n))
    return np.column_stack((az, el))


def ukf_update_vectorized(ukf: UKF, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorized UKF measurement update using FilterPy sigma points + weights, but
    computing Zsig in one NumPy batch.

    Returns (innovation y, z_pred_mean)
    """
    # Ensure SPD before sigma points
    ukf.P = make_pd(ukf.P)

    sigmas = ukf.points_fn.sigma_points(ukf.x, ukf.P)  # (Ns,6)
    Zsig = z_from_sigmas(sigmas)                       # (Ns,2)

    Wm = ukf.Wm
    Wc = ukf.Wc

    # circular mean for az, linear mean for el
    sin_az = np.sin(Zsig[:, 0])
    cos_az = np.cos(Zsig[:, 0])
    az_mean = math.atan2(float(np.dot(Wm, sin_az)), float(np.dot(Wm, cos_az)))
    if az_mean < 0.0:
        az_mean += 2.0*np.pi
    el_mean = float(np.dot(Wm, Zsig[:, 1]))
    z_mean = np.array([az_mean, el_mean], dtype=float)

    # innovation sigmas
    dz = Zsig - z_mean[None, :]
    dz[:, 0] = wrap_to_pi(dz[:, 0])

    dx = sigmas - ukf.x[None, :]

    # S = sum Wc * dz dz^T + R
    S = np.zeros((2, 2), dtype=float)
    for i in range(Zsig.shape[0]):
        S += Wc[i] * np.outer(dz[i], dz[i])
    S += ukf.R

    # Pxz = sum Wc * dx dz^T
    Pxz = np.zeros((6, 2), dtype=float)
    for i in range(Zsig.shape[0]):
        Pxz += Wc[i] * np.outer(dx[i], dz[i])

    # Kalman gain
    K = Pxz @ np.linalg.inv(S)

    # innovation
    y = z - z_mean
    y[0] = float(wrap_to_pi(np.array([y[0]]))[0])

    # state update
    ukf.x = ukf.x + K @ y

    # Joseph-ish stabilized covariance update: P = P - K S K^T, then sym + PD fix
    ukf.P = ukf.P - K @ S @ K.T
    ukf.P = make_pd(ukf.P)

    return y, z_mean


# =========================
# Fast CSV helpers
# =========================
def load_pass_measurements(pass_csv: Path) -> dict[int, tuple[float, float]]:
    """
    {time_since_epoch_s_int: (az_deg, el_deg)} for rows where visible==1
    """
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

def scan_duration_seconds(full_csv: Path) -> tuple[int, int]:
    # fast-ish: read first and last data rows
    with full_csv.open("r", newline="", encoding="utf-8") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        idx = {name: i for i, name in enumerate(header)}
        t_i = idx["time_since_epoch_s"]

        first = next(rdr, None)
        if first is None:
            raise RuntimeError("Full trajectory CSV is empty.")
        t0 = int(round(float(first[t_i])))

        last = first
        for row in rdr:
            last = row
        tf = int(round(float(last[t_i])))

    return t0, tf


# =========================
# Main
# =========================
def main():
    random.seed(RNG_SEED)
    np.random.seed(RNG_SEED)

    full_path = Path(FULL_TRAJ_CSV)
    pass_path = Path(PASS_MEAS_CSV)
    out_path = Path(OUTPUT_CSV)

    if not full_path.exists():
        raise FileNotFoundError(f"Missing {FULL_TRAJ_CSV}")
    if not pass_path.exists():
        raise FileNotFoundError(f"Missing {PASS_MEAS_CSV}")

    pass_meas = load_pass_measurements(pass_path)

    t0, tf = scan_duration_seconds(full_path)
    total_sim_s = max(0, tf - t0)

    pbar = None
    if tqdm is not None:
        pbar = tqdm(total=total_sim_s, desc="Simulation time", unit="s", smoothing=0.0)

    # measurement covariance (rad^2)
    R = np.diag([
        math.radians(SIGMA_AZ_DEG) ** 2,
        math.radians(SIGMA_EL_DEG) ** 2,
    ])

    ukf: UKF | None = None
    have_filter = False
    prev_t_int: int | None = None
    last_pbar_update_at: int | None = None

    # Open full CSV (fast reader)
    with full_path.open("r", newline="", encoding="utf-8") as f_in, out_path.open("w", newline="", encoding="utf-8") as f_out:
        rdr = csv.reader(f_in)
        header = next(rdr)
        idx = {name: i for i, name in enumerate(header)}

        # indices we need
        t_i = idx["time_since_epoch_s"]
        timeutc_i = idx["time_utc"]
        rx_i = idx["r_x_km"]
        ry_i = idx["r_y_km"]
        rz_i = idx["r_z_km"]
        vx_i = idx["v_x_km_s"]
        vy_i = idx["v_y_km_s"]
        vz_i = idx["v_z_km_s"]

        w = csv.writer(f_out)
        w.writerow([
            "time_utc", "time_since_epoch_s",
            "meas_available",
            "z_az_deg_noisy", "z_el_deg_noisy",
            "innov_az_deg", "innov_el_deg",
            "xhat_rx_km", "xhat_ry_km", "xhat_rz_km",
            "xhat_vx_km_s", "xhat_vy_km_s", "xhat_vz_km_s",
            "P_rr_x", "P_rr_y", "P_rr_z",
            "P_vv_x", "P_vv_y", "P_vv_z",
        ])

        for row in rdr:
            t_int = int(round(float(row[t_i])))
            dt = 0.0 if prev_t_int is None else float(t_int - prev_t_int)

            # progress bar: update only every PBAR_UPDATE_EVERY_S seconds
            if pbar is not None:
                if last_pbar_update_at is None:
                    last_pbar_update_at = t_int
                if (t_int - last_pbar_update_at) >= PBAR_UPDATE_EVERY_S:
                    pbar.n = max(0, t_int - t0)
                    pbar.refresh()
                    last_pbar_update_at = t_int

            # cache rotation ONCE per epoch
            dt_utc = parse_time_utc(row[timeutc_i])
            global _R_STEP
            _R_STEP = rot_z(gmst_rad(dt_utc))

            # truth state (for init)
            r_truth = np.array([float(row[rx_i]), float(row[ry_i]), float(row[rz_i])], dtype=float)
            v_truth = np.array([float(row[vx_i]), float(row[vy_i]), float(row[vz_i])], dtype=float)

            meas_available = (t_int in pass_meas)

            # Predict step if filter exists
            if have_filter:
                assert ukf is not None

                # process noise scaled by dt
                if dt > 0.0:
                    Q = np.zeros((6, 6), dtype=float)
                    Q[0:3, 0:3] = (Q_POS_KM2_PER_S * dt) * np.eye(3)
                    Q[3:6, 3:6] = (Q_VEL_KM2_S2_PER_S * dt) * np.eye(3)
                    ukf.Q = Q
                else:
                    ukf.Q = np.zeros((6, 6), dtype=float)

                # PD-fix before predict (prevents your crash)
                ukf.P = make_pd(ukf.P)
                ukf.predict(dt=dt)
                ukf.P = make_pd(ukf.P)

            z_az_deg_noisy = ""
            z_el_deg_noisy = ""
            innov_az_deg = ""
            innov_el_deg = ""

            if meas_available:
                az_deg, el_deg = pass_meas[t_int]

                # add noise to CSV measurements
                az_meas_deg = az_deg + random.gauss(0.0, SIGMA_AZ_DEG)
                el_meas_deg = el_deg + random.gauss(0.0, SIGMA_EL_DEG)

                z_az_deg_noisy = f"{az_meas_deg:.6f}"
                z_el_deg_noisy = f"{el_meas_deg:.6f}"
                z = np.array([math.radians(az_meas_deg), math.radians(el_meas_deg)], dtype=float)

                if not have_filter:
                    # initialize at first measurement epoch
                    x0 = np.hstack((r_truth, v_truth)).astype(float)
                    x0[0:3] += np.random.normal(0.0, INIT_POS_SIGMA_KM, size=3)
                    x0[3:6] += np.random.normal(0.0, INIT_VEL_SIGMA_KM_S, size=3)

                    P0 = np.zeros((6, 6), dtype=float)
                    P0[0:3, 0:3] = (INIT_POS_SIGMA_KM**2) * np.eye(3)
                    P0[3:6, 3:6] = (INIT_VEL_SIGMA_KM_S**2) * np.eye(3)
                    P0 = make_pd(P0)

                    points = MerweScaledSigmaPoints(n=6, alpha=ALPHA, beta=BETA, kappa=KAPPA)
                    ukf = UKF(dim_x=6, dim_z=2, fx=fx, hx=hx, dt=1.0, points=points)
                    ukf.x = x0
                    ukf.P = P0
                    ukf.R = R
                    ukf.Q = np.zeros((6, 6), dtype=float)

                    have_filter = True

                    # if dt>0 (starting mid-file), step to this epoch
                    if dt > 0.0:
                        Q = np.zeros((6, 6), dtype=float)
                        Q[0:3, 0:3] = (Q_POS_KM2_PER_S * dt) * np.eye(3)
                        Q[3:6, 3:6] = (Q_VEL_KM2_S2_PER_S * dt) * np.eye(3)
                        ukf.Q = Q
                        ukf.P = make_pd(ukf.P)
                        ukf.predict(dt=dt)
                        ukf.P = make_pd(ukf.P)

                assert ukf is not None

                # FAST vectorized update (also PD-stabilized)
                y, _zmean = ukf_update_vectorized(ukf, z)
                innov_az_deg = f"{math.degrees(float(y[0])):.6f}"
                innov_el_deg = f"{math.degrees(float(y[1])):.6f}"

            # output
            if have_filter:
                assert ukf is not None
                x = ukf.x
                P = ukf.P
                w.writerow([
                    row[timeutc_i], f"{t_int:.3f}",
                    int(meas_available),
                    z_az_deg_noisy, z_el_deg_noisy,
                    innov_az_deg, innov_el_deg,
                    f"{x[0]:.6f}", f"{x[1]:.6f}", f"{x[2]:.6f}",
                    f"{x[3]:.9f}", f"{x[4]:.9f}", f"{x[5]:.9f}",
                    f"{P[0,0]:.6e}", f"{P[1,1]:.6e}", f"{P[2,2]:.6e}",
                    f"{P[3,3]:.6e}", f"{P[4,4]:.6e}", f"{P[5,5]:.6e}",
                ])
            else:
                w.writerow([
                    row[timeutc_i], f"{t_int:.3f}",
                    int(meas_available),
                    z_az_deg_noisy, z_el_deg_noisy,
                    "", "",
                    "", "", "",
                    "", "", "",
                    "", "", "",
                    "", "", "",
                ])

            prev_t_int = t_int

    if pbar is not None:
        pbar.n = total_sim_s
        pbar.refresh()
        pbar.close()

    print(f"Done. Wrote: {out_path.resolve()}")
    print(f"Used {len(pass_meas)} visible measurement epochs from {PASS_MEAS_CSV}")


if __name__ == "__main__":
    main()
