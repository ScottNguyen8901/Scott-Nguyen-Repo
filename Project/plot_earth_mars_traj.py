import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy import units as u
from poliastro.bodies import Sun
from poliastro.ephem import get_body_ephem
from poliastro import iod
from poliastro.twobody import Orbit
import os

# --- Parameters ---
dep_date_str = "2026-11-12"
arr_date_str = "2027-09-08"

dep_time = Time(dep_date_str, scale='tdb')
arr_time = Time(arr_date_str, scale='tdb')

tof = (arr_time - dep_time).to(u.s)

# --- Get planetary states ---
r_earth, v_earth = get_body_ephem("earth", dep_time)
r_mars, v_mars = get_body_ephem("mars", arr_time)

r_earth = r_earth.to(u.km)
v_earth = v_earth.to(u.km / u.s)

r_mars = r_mars.to(u.km)
v_mars = v_mars.to(u.km / u.s)

# --- Solve Lambert's problem for transfer ---
(v_dep, v_arr), = iod.lambert(Sun.k, r_earth, r_mars, tof)

# Create transfer orbit from departure state around the Sun
transfer_orbit = Orbit.from_vectors(Sun, r_earth, v_dep, dep_time)

# --- Times to sample orbits ---
num_points = 300
times = dep_time + np.linspace(0, tof.value, num_points) * u.s

# --- Precompute positions ---
earth_positions = np.array([get_body_ephem("earth", t)[0].to(u.km).value for t in times])
mars_positions = np.array([get_body_ephem("mars", t)[0].to(u.km).value for t in times])
spacecraft_positions = np.array([transfer_orbit.propagate(t - dep_time).r.to(u.km).value for t in times])

# --- Create directory for saving trajectories ---
save_dir = "trajectory"
os.makedirs(save_dir, exist_ok=True)

# --- Save to CSV ---
np.savetxt(os.path.join(save_dir, "earth_trajectory.csv"), earth_positions, delimiter=",", header="X_km,Y_km,Z_km", comments="")
np.savetxt(os.path.join(save_dir, "mars_trajectory.csv"), mars_positions, delimiter=",", header="X_km,Y_km,Z_km", comments="")
np.savetxt(os.path.join(save_dir, "spacecraft_trajectory.csv"), spacecraft_positions, delimiter=",", header="X_km,Y_km,Z_km", comments="")

# --- Plot ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f"Earth-Mars Transfer from {dep_date_str} to {arr_date_str}")
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_zlabel("Z (km)")

ax.set_xlim(-1.7e8, 1.7e8)
ax.set_ylim(-1.7e8, 1.7e8)
ax.set_zlim(-1.7e8, 1.7e8)

# Sun at origin
ax.scatter(0, 0, 0, color='yellow', s=300, label='Sun')

# Plot orbits as lines
ax.plot(earth_positions[:, 0], earth_positions[:, 1], earth_positions[:, 2], 'b-', label='Earth Orbit')
ax.plot(mars_positions[:, 0], mars_positions[:, 1], mars_positions[:, 2], 'r-', label='Mars Orbit')
ax.plot(spacecraft_positions[:, 0], spacecraft_positions[:, 1], spacecraft_positions[:, 2], 'g-', label='Transfer Trajectory')

# Mark start and end points with dots
ax.scatter(*earth_positions[0], color='blue', s=50, label='Earth at Departure')
ax.scatter(*mars_positions[-1], color='red', s=50, label='Mars at Arrival')
ax.scatter(*spacecraft_positions[0], color='green', s=50, label='Spacecraft Start')
ax.scatter(*spacecraft_positions[-1], color='green', s=50, marker='X', label='Spacecraft End')

ax.legend()
plt.show()