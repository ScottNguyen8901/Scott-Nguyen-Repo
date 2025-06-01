import numpy as np
from astropy.time import Time
from astropy import units as u
from poliastro.bodies import Sun, Earth, Mars
from poliastro.ephem import get_body_ephem
from poliastro import iod
from scipy.optimize import minimize
from numpy.linalg import norm
from astropy.constants import G, M_earth, R_earth, g0

# Constants
mu_sun = Sun.k.to(u.km**3 / u.s**2).value  # Sun gravitational parameter (km^3/s^2)
mu_earth = (G * M_earth).to(u.km**3 / u.s**2).value  # Earth's gravitational parameter (km^3/s^2) from astropy
R_E = R_earth.to(u.km).value  # Earth radius (km) from astropy
alt_earth = 300  # Altitude of flyby (km)

# Define time window for ephemerides queries
start_time = Time("2025-01-01", scale="tdb")
end_time = Time("2035-12-31", scale="tdb")

def get_state(body_name, t_jd):
    """Get heliocentric position and velocity vectors of a body at Julian Date t_jd."""
    t_astropy = Time(t_jd, format='jd', scale='tdb')
    pos, vel = get_body_ephem(body_name, t_astropy)  # unpack tuple here
    r = pos.to(u.km).value
    v = vel.to(u.km / u.s).value
    return r, v

def delta_v_interplanetary(x):
    """Calculate total delta-V including Earth departure and Mars capture burns."""
    jd_start, jd_end = x

    if jd_end <= jd_start:
        return np.inf

    # Get state vectors
    r1, v_earth = get_state("earth", jd_start)
    r2, v_mars = get_state("mars", jd_end)

    tof = (jd_end - jd_start) * 86400  # seconds

    try:
        (v_dep, v_arr), = iod.lambert(Sun.k, r1 * u.km, r2 * u.km, tof * u.s)
    except Exception:
        return np.inf

    v_dep = v_dep.to(u.km / u.s).value
    v_arr = v_arr.to(u.km / u.s).value

    # --- Earth departure ΔV ---
    v_inf_dep = norm(v_dep - v_earth)
    r_p_earth = R_E + alt_earth
    v_circ_earth = np.sqrt(mu_earth / r_p_earth)
    v_hyp_earth = np.sqrt(v_inf_dep**2 + 2 * mu_earth / r_p_earth)
    delta_v_dep = abs(v_hyp_earth - v_circ_earth)

    # --- Mars arrival ΔV (MOI burn) ---
    mu_mars = Mars.k.to(u.km**3 / u.s**2).value
    R_mars = Mars.R.to(u.km).value
    alt_mars = 300  # target Mars orbit altitude (km)

    v_inf_arr = norm(v_arr - v_mars)
    r_p_mars = R_mars + alt_mars
    v_circ_mars = np.sqrt(mu_mars / r_p_mars)
    v_hyp_mars = np.sqrt(v_inf_arr**2 + 2 * mu_mars / r_p_mars)
    delta_v_arr = abs(v_hyp_mars - v_circ_mars)

    # --- Total ΔV ---
    total_delta_v = delta_v_dep + delta_v_arr

    return total_delta_v  # Minimize total ΔV

def constraints(x):
    jd_start, jd_end = x
    tof = jd_end - jd_start

    min_tof = 150  # days
    max_tof = 305  # days

    c1 = tof - min_tof  # tof >= min_tof
    c2 = jd_start - Time("2025-01-01", scale="tdb").jd  # jd_start >= lower bound
    c3 = Time("2030-01-01", scale="tdb").jd - jd_start  # jd_start <= upper bound
    c4 = max_tof - tof  # tof <= max_tof

    return [c1, c2, c3, c4]

jd_start_lower = Time("2025-01-01", scale="tdb").jd
jd_start_upper = Time("2030-01-01", scale="tdb").jd
jd_end_lower = jd_start_lower + 150
jd_end_upper = Time("2035-01-01", scale="tdb").jd

bounds = [(jd_start_lower, jd_start_upper), (jd_end_lower, jd_end_upper)]

x0 = [Time("2028-01-01", scale="tdb").jd, Time("2030-01-01", scale="tdb").jd]

nonlinear_constraint = {
    'type': 'ineq',
    'fun': lambda x: np.array(constraints(x))
}

result = minimize(delta_v_interplanetary, x0, method='SLSQP', bounds=bounds, constraints=[nonlinear_constraint])

if result.success:
    jd_start_opt, jd_end_opt = result.x
    start_date = Time(jd_start_opt, format='jd').to_datetime()
    end_date = Time(jd_end_opt, format='jd').to_datetime()
    print(f"Optimal departure date: {start_date.strftime('%Y-%m-%d')}")
    print(f"Optimal arrival date: {end_date.strftime('%Y-%m-%d')}")
    print(f"Total Delta-v (km/s): {result.fun:.2f}")
else:
    print("Optimization failed:", result.message)