import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy import units as u
from poliastro.bodies import Sun, Earth, Mars
from poliastro.ephem import get_body_ephem
from poliastro import iod
from numpy.linalg import norm

# Constants
mu_sun = Sun.k.to(u.km**3 / u.s**2).value  # Sun gravitational parameter (km^3/s^2)
mu_earth = Earth.k.to(u.km**3 / u.s**2).value  # Earth's gravitational parameter (km^3/s^2)
R_E = Earth.R.to(u.km).value  # Earth radius (km)
alt_earth = 300  # Altitude of flyby (km)

mu_mars = Mars.k.to(u.km**3 / u.s**2).value
R_mars = Mars.R.to(u.km).value
alt_mars = 300  # target Mars orbit altitude (km)

# Date ranges for departure and arrival (Julian Dates)
dep_start = Time("2025-01-01", scale="tdb").jd
dep_end = Time("2030-01-01", scale="tdb").jd
arr_start = dep_start + 150  # minimum 150 days transfer time
arr_end = Time("2035-01-01", scale="tdb").jd

# Create grids of departure and arrival dates (in JD)
num_points = 60
dep_dates = np.linspace(dep_start, dep_end, num_points)
arr_dates = np.linspace(arr_start, arr_end, num_points)

# Initialize delta-V matrix
delta_v_matrix = np.full((num_points, num_points), np.nan)

total_steps = num_points * num_points
step_count = 0

print("Computing porkchop plot delta-V values:")

for i, t_dep in enumerate(dep_dates):
    # Get Earth state once per departure date
    r_earth_astropy, v_earth_astropy = get_body_ephem("earth", Time(t_dep, format='jd', scale='tdb'))
    r_earth = r_earth_astropy.to(u.km).value
    v_earth = v_earth_astropy.to(u.km / u.s).value

    for j, t_arr in enumerate(arr_dates):
        step_count += 1
        progress = (step_count / total_steps) * 100
        print(f"\rProgress: {progress:.1f}%", end='', flush=True)

        if t_arr <= t_dep:
            continue  # invalid case

        tof = (t_arr - t_dep) * 86400  # seconds

        # Get Mars state once per arrival date (not inside Lambert)
        r_mars_astropy, v_mars_astropy = get_body_ephem("mars", Time(t_arr, format='jd', scale='tdb'))
        r_mars = r_mars_astropy.to(u.km).value
        v_mars = v_mars_astropy.to(u.km / u.s).value

        try:
            (v_dep, v_arr), = iod.lambert(Sun.k, r_earth * u.km, r_mars * u.km, tof * u.s)
        except Exception:
            continue  # skip if Lambert solver fails

        v_dep = v_dep.to(u.km / u.s).value
        v_arr = v_arr.to(u.km / u.s).value

        # Earth departure ΔV
        v_inf_dep = norm(v_dep - v_earth)
        r_p_earth = R_E + alt_earth
        v_circ_earth = np.sqrt(mu_earth / r_p_earth)
        v_hyp_earth = np.sqrt(v_inf_dep**2 + 2 * mu_earth / r_p_earth)
        delta_v_dep = abs(v_hyp_earth - v_circ_earth)

        # Mars arrival ΔV
        v_inf_arr = norm(v_arr - v_mars)
        r_p_mars = R_mars + alt_mars
        v_circ_mars = np.sqrt(mu_mars / r_p_mars)
        v_hyp_mars = np.sqrt(v_inf_arr**2 + 2 * mu_mars / r_p_mars)
        delta_v_arr = abs(v_hyp_mars - v_circ_mars)

        total_delta_v = delta_v_dep + delta_v_arr
        
        delta_v_matrix[i, j] = total_delta_v

print("\nPorkchop plot generation complete!")

# Find minimum ΔV and its indices
min_index = np.nanargmin(delta_v_matrix)
min_i, min_j = np.unravel_index(min_index, delta_v_matrix.shape)
min_dv = delta_v_matrix[min_i, min_j]

min_dep_jd = dep_dates[min_i]
min_arr_jd = arr_dates[min_j]

# Convert JD to calendar date strings
min_dep_date_str = Time(min_dep_jd, format='jd').to_datetime().strftime('%Y-%m-%d')
min_arr_date_str = Time(min_arr_jd, format='jd').to_datetime().strftime('%Y-%m-%d')

print(f"Minimum total ΔV: {min_dv:.4f} km/s")
print(f"Departure date: {min_dep_date_str}")
print(f"Arrival date: {min_arr_date_str}")

# Plotting
dep_grid, arr_grid = np.meshgrid(dep_dates, arr_dates, indexing='ij')

# Convert JD to datetime for axis ticks
dep_dates_dt = Time(dep_dates, format='jd').to_datetime()
arr_dates_dt = Time(arr_dates, format='jd').to_datetime()

plt.figure(figsize=(12, 9))

# Levels: 0 to 6 km/s by 0.01 for smoother contour
levels = np.arange(0, 100.00, 0.01)

cp = plt.contourf(dep_grid, arr_grid, delta_v_matrix, levels=levels, cmap='viridis', extend='both')

cbar = plt.colorbar(cp)
cbar.set_label('Total ΔV (km/s)')

plt.xlabel('Departure Date')
plt.ylabel('Arrival Date')
plt.title('Earth to Mars Porkchop Plot (Total ΔV)')

tick_step = 10
plt.xticks(dep_dates[::tick_step], [d.strftime('%Y-%m-%d') for d in dep_dates_dt[::tick_step]], rotation=45)
plt.yticks(arr_dates[::tick_step], [d.strftime('%Y-%m-%d') for d in arr_dates_dt[::tick_step]])

# Plot red star for minimum delta-V with calendar dates
plt.plot(min_dep_jd, min_arr_jd, 'r*', markersize=15, label=f"Min ΔV = {min_dv:.4f} km/s")
plt.text(min_dep_jd, min_arr_jd, f"{min_dv:.4f} km/s\n{min_dep_date_str} → {min_arr_date_str}",
         color='red', fontsize=10, verticalalignment='bottom', horizontalalignment='right')

plt.legend()
plt.tight_layout()
plt.show()