import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rpm_predictor import RPMPredictor

# Set global font to Times New Roman and bold
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# --- Load and process data ---
df = pd.read_csv("data.csv", skiprows=3)
df.columns = df.columns.str.strip()

cols = ['pitch_angle[deg]', 'vel_x[m/s]', 'vel_y[m/s]', 'prop_speed[rpm]', 'est_wind_x[m/s]', 'est_wind_y[m/s]']
df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')

pitch = df['pitch_angle[deg]'].to_numpy()
vx = df['vel_x[m/s]'].to_numpy()
vy = df['vel_y[m/s]'].to_numpy()
rpm_true = df['prop_speed[rpm]'].to_numpy()
wx = df['est_wind_x[m/s]'].to_numpy()
wy = df['est_wind_y[m/s]'].to_numpy()

# --- Prediction ---
model = RPMPredictor()
rpm_pred, rpm_flag = model.predict(pitch, vx, vy, wx, wy)
rpm_error = rpm_true - rpm_pred

# --- Percent error ---
percent_error = 100 * rpm_error / np.where(rpm_true != 0, rpm_true, np.nan)

# --- Recompute AoA and Airspeed for visualization ---
ax = vx - wx
ay = vy - wy
airspeed = np.hypot(ax, ay)
fpa = np.degrees(np.arctan2(ay, ax))
aoa_data = pitch - fpa

# --- Error reporting before plots ---
mean_error = np.nanmean(rpm_error)
std_error = np.nanstd(rpm_error)
mean_pct_error = np.nanmean(percent_error)
std_pct_error = np.nanstd(percent_error)
max_over_pct_error = np.nanmax(percent_error)
max_under_pct_error = np.nanmin(percent_error)

print(f"Mean RPM Error: {mean_error:.4f}")
print(f"Standard Deviation of RPM Error: {std_error:.4f}")
print(f"Mean Percent Error: {mean_pct_error:.2f}%")
print(f"Standard Deviation of Percent Error: {std_pct_error:.2f}%")
print(f"\nMaximum Over-Prediction Percent Error: {max_over_pct_error:.2f}%")
print(f"Maximum Under-Prediction Percent Error: {max_under_pct_error:.2f}%")

print(f"\nRPM Quality Flags: GOOD = {np.sum(rpm_flag == 'GOOD')}, BAD = {np.sum(rpm_flag == 'BAD')}")

bad_indices = np.where(rpm_flag == 'BAD')[0]
print("\nDetailed Errors for BAD RPM Predictions:")
if bad_indices.size == 0:
    print("All predictions within 3σ bounds.")
else:
    for i in bad_indices:
        print(f"Index {i:>4}: True = {rpm_true[i]:.2f}, Predicted = {rpm_pred[i]:.2f}, Error = {rpm_error[i]:+.2f}, Percent Error = {percent_error[i]:+.2f}%")

# --- Plot: Predicted, True, and Error RPM ---
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

good_mask = rpm_flag == 'GOOD'
bad_mask = rpm_flag == 'BAD'

sc0 = axs[0].scatter(aoa_data[good_mask], airspeed[good_mask], c=rpm_pred[good_mask], cmap='viridis', edgecolor='k', label='GOOD')
axs[0].scatter(aoa_data[bad_mask], airspeed[bad_mask], marker='x', color='red', label='BAD')
axs[0].set_title('Predicted RPM')
axs[0].set_xlabel('AoA [deg]')
axs[0].set_ylabel('Airspeed [m/s]')
axs[0].grid(True)
plt.colorbar(sc0, ax=axs[0], label='RPM')
axs[0].legend()

sc1 = axs[1].scatter(aoa_data, airspeed, c=rpm_true, cmap='viridis', edgecolor='k')
axs[1].set_title('True RPM')
axs[1].set_xlabel('AoA [deg]')
axs[1].set_ylabel('Airspeed [m/s]')
axs[1].grid(True)
plt.colorbar(sc1, ax=axs[1], label='RPM')

sc2 = axs[2].scatter(aoa_data[good_mask], airspeed[good_mask], c=rpm_error[good_mask], cmap='coolwarm', edgecolor='k', label='GOOD')
axs[2].scatter(aoa_data[bad_mask], airspeed[bad_mask], marker='x', color='red', label='BAD')
axs[2].set_title('RPM Error (True - Pred)')
axs[2].set_xlabel('AoA [deg]')
axs[2].set_ylabel('Airspeed [m/s]')
axs[2].grid(True)
plt.colorbar(sc2, ax=axs[2], label='RPM Error')
axs[2].legend()

plt.tight_layout()

# --- Histogram of RPM Error ---
plt.figure(figsize=(10, 6))
plt.hist(rpm_error, bins=50, color='lightblue', edgecolor='k')
plt.axvline(mean_error, color='red', linestyle='dashed', linewidth=2, label=f'Mean = {mean_error:.2f}')
plt.axvline(mean_error + std_error, color='green', linestyle='dashed', linewidth=2, label=f'+1σ = {mean_error + std_error:.2f}')
plt.axvline(mean_error - std_error, color='green', linestyle='dashed', linewidth=2, label=f'-1σ = {mean_error - std_error:.2f}')
plt.title('Histogram of RPM Prediction Error', fontweight='bold')
plt.xlabel('RPM Error (True - Pred)', fontweight='bold')
plt.ylabel('Frequency', fontweight='bold')
plt.grid(True)
plt.legend(fontsize='medium')
plt.tight_layout()

# --- Histogram of Percent Error ---
plt.figure(figsize=(10, 6))
plt.hist(percent_error[good_mask], bins=50, color='lightcoral', edgecolor='k')
plt.axvline(mean_pct_error, color='red', linestyle='dashed', linewidth=2, label=f'Mean = {mean_pct_error:.2f}%')
plt.axvline(mean_pct_error + std_pct_error, color='green', linestyle='dashed', linewidth=2, label=f'+1σ = {mean_pct_error + std_pct_error:.2f}%')
plt.axvline(mean_pct_error - std_pct_error, color='green', linestyle='dashed', linewidth=2, label=f'-1σ = {mean_pct_error - std_pct_error:.2f}%')
plt.title('Histogram of RPM Percent Error', fontweight='bold')
plt.xlabel('Percent Error (%)', fontweight='bold')
plt.ylabel('Frequency', fontweight='bold')
plt.grid(True)
plt.legend(fontsize='medium')
plt.tight_layout()

plt.show()