import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.ensemble import RandomForestRegressor

# --- Settings ---
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.weight'] = 'bold'

# --- Load Data ---
filename = 'data.csv'
df_all = pd.read_csv(filename, header=None)
aoa, clcd = df_all.iloc[0, 1:].astype(float), df_all.iloc[1, 1:].astype(float)

df_data = pd.read_csv(filename, skiprows=3)
df_data.columns = df_data.columns.str.strip()
cols = ['pitch_angle[deg]', 'vel_x[m/s]', 'vel_y[m/s]', 'prop_speed[rpm]', 'est_wind_x[m/s]', 'est_wind_y[m/s]']
df_data[cols] = df_data[cols].apply(pd.to_numeric, errors='coerce')
pitch, vx, vy, rpm, wx, wy = (df_data[col].to_numpy() for col in cols)

print("\n--- Mean and Standard Deviation of Raw Flight Data ---")
for col in ['vel_x[m/s]', 'vel_y[m/s]', 'est_wind_x[m/s]', 'est_wind_y[m/s]', 'pitch_angle[deg]']:
    mean_val = df_data[col].mean(skipna=True)
    std_val = df_data[col].std(skipna=True)
    print(f"{col}: Mean = {mean_val:.4f}, Std = {std_val:.4f}")

# --- Derived Quantities ---
ax, ay = vx - wx, vy - wy
airspeed = np.hypot(ax, ay)
fpa = np.degrees(np.arctan2(ay, ax))
aoa_data = pitch - fpa

A, cd, rho = 0.625, 0.13, 1.225
deg_clcd = 10
poly_clcd = np.poly1d(np.polyfit(aoa, clcd, deg_clcd))
cl_lookup = poly_clcd(aoa_data) * cd

# --- Print Cl/Cd Terms ---
print("\nCl/Cd Polynomial Terms:")
for i, coef in enumerate(poly_clcd.coefficients[::-1]):
    print(f"x^{i}: {coef:.6e}")

# --- Aerodynamic Forces ---
L = 0.5 * rho * airspeed**2 * A * cl_lookup
D = 0.5 * rho * airspeed**2 * A * cd
f = D / (rpm**2)

# --- Filter Valid ---
X = np.column_stack((aoa_data, airspeed))
valid = np.isfinite(f) & np.all(np.isfinite(X), axis=1)
X, f = X[valid], f[valid]
aoa_data, airspeed, rpm, L, D = [arr[valid] for arr in [aoa_data, airspeed, rpm, L, D]]

# --- Print AoA and Airspeed Distribution ---
print("\n--- AoA and Airspeed Distribution ---")
print(f"AoA [deg]: Mean = {np.mean(aoa_data):.4f}, Std = {np.std(aoa_data):.4f}")
print(f"Airspeed [m/s]: Mean = {np.mean(airspeed):.4f}, Std = {np.std(airspeed):.4f}")

# --- Standardize and Polynomial Expansion ---
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("\n--- Feature Scaling Parameters (StandardScaler) ---")
print(f"Mean (aoa, v): {scaler.mean_}")
print(f"Variance (aoa, v): {scaler.var_}")

poly = PolynomialFeatures(degree=4)
X_poly = poly.fit_transform(X_scaled)
feature_names = poly.get_feature_names_out(['aoa', 'v'])[1:]

# --- Fit Linear Regression Model ---
lin_model = LinearRegression().fit(X_poly, f)
pred_lin = lin_model.predict(X_poly)
error_f = f - pred_lin

print("\n--- Linear Regression Coefficients ---")
print(f"Intercept: {lin_model.intercept_:.6e}")
for name, coef in zip(feature_names, lin_model.coef_[1:]):
    print(f"{name}: {coef:.6e}")
print(f"R² (Linear): {lin_model.score(X_poly, f):.4f}")

# --- Fit Ridge and Random Forest Models ---
ridge_model = Ridge(alpha=1.0).fit(X_poly, f)
pred_ridge = ridge_model.predict(X_poly)
rf_model = RandomForestRegressor(n_estimators=100, random_state=42).fit(X_poly, f)
pred_rf = rf_model.predict(X_poly)

print("\n--- Ridge Regression Coefficients ---")
print(f"Intercept: {ridge_model.intercept_:.6e}")
for name, coef in zip(feature_names, ridge_model.coef_[1:]):
    print(f"{name}: {coef:.6e}")
print(f"R² (Ridge): {ridge_model.score(X_poly, f):.4f}")
print(f"R² (Random Forest): {rf_model.score(X_poly, f):.4f}")

# --- Plot: Distributions of State Variables ---
fig, axs = plt.subplots(2, 3, figsize=(18, 10))
axs = axs.ravel()

variables = [
    ('vel_x[m/s]', vx, 'Velocity X [m/s]'),
    ('vel_y[m/s]', vy, 'Velocity Y [m/s]'),
    ('est_wind_x[m/s]', wx, 'Wind X [m/s]'),
    ('est_wind_y[m/s]', wy, 'Wind Y [m/s]'),
    ('pitch_angle[deg]', pitch, 'Pitch Angle [deg]'),
    ('aoa [deg]', aoa_data, 'Angle of Attack [deg]'),
]

for i, (label, data, title) in enumerate(variables):
    axs[i].hist(data, bins=50, color='skyblue', edgecolor='k')
    axs[i].set_title(title, fontweight='bold')
    axs[i].set_xlabel(label, fontweight='bold')
    axs[i].set_ylabel('Frequency', fontweight='bold')
    axs[i].grid(True)

plt.suptitle('Distributions of Flight State Variables', fontsize=16, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])

# --- Surface Grid for Plotting (Linear Model) ---
aoa_vals = np.linspace(aoa_data.min(), aoa_data.max(), 50)
v_vals = np.linspace(airspeed.min(), airspeed.max(), 50)
AOA_grid, V_grid = np.meshgrid(aoa_vals, v_vals)
X_grid_poly = poly.transform(scaler.transform(np.column_stack((AOA_grid.ravel(), V_grid.ravel()))))
Z_grid = lin_model.predict(X_grid_poly).reshape(AOA_grid.shape)

# --- Plot: Cl/Cd and Cl vs AoA ---
x_fit = np.linspace(aoa.min(), aoa.max(), 200)
y_fit_clcd = poly_clcd(x_fit)
poly_cl = np.poly1d(np.polyfit(aoa, clcd * cd, deg_clcd))
y_fit_cl = poly_cl(x_fit)

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.scatter(aoa, clcd, label='Data')
plt.plot(x_fit, y_fit_clcd, label=f'Poly Fit {deg_clcd}')
plt.xlabel('Angle of Attack [deg]', fontweight='bold')
plt.ylabel('Cl/Cd', fontweight='bold')
plt.title('AOA vs Cl/Cd', fontweight='bold')
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
plt.scatter(aoa, clcd * cd, label='Data')
plt.plot(x_fit, y_fit_cl, label=f'Poly Fit {deg_clcd}')
plt.xlabel('Angle of Attack [deg]', fontweight='bold')
plt.ylabel('Cl', fontweight='bold')
plt.title('AOA vs Cl', fontweight='bold')
plt.grid(True)
plt.legend()
plt.tight_layout()

# --- Plot: Lift/Drag vs RPM ---
fig1, ax = plt.subplots(1, 3, figsize=(18, 5))
ax[0].plot(L, rpm, 'o')
ax[0].set_xlabel('Lift [N]', fontweight='bold')
ax[0].set_ylabel('Propeller Speed [rpm]', fontweight='bold')

ax[1].plot(D, rpm, 's')
ax[1].set_xlabel('Drag [N]', fontweight='bold')
ax[1].set_ylabel('Propeller Speed [rpm]', fontweight='bold')

lift_lin = np.linspace(min(L), max(L), 100)
drag_lin = np.linspace(min(D), max(D), 100)
L_grid, D_grid = np.meshgrid(lift_lin, drag_lin)
rpm_grid = griddata((L, D), rpm, (L_grid, D_grid), method='cubic')
contour = ax[2].contourf(L_grid, D_grid, rpm_grid, levels=20, cmap='viridis', alpha=0.7)
ax[2].scatter(L, D, c=rpm, cmap='viridis', edgecolor='k')
fig1.colorbar(contour, ax=ax[2], label='Propeller Speed [rpm]')
ax[2].set_xlabel('Lift [N]', fontweight='bold')
ax[2].set_ylabel('Drag [N]', fontweight='bold')

labels = ['Lift vs RPM', 'Drag vs RPM', 'Prop Speed (rpm) by Lift and Drag']
for i in range(3):
    ax[i].set_title(labels[i], fontweight='bold')
    ax[i].grid(True)
plt.tight_layout()

# --- Plot: AoA & Airspeed vs RPM ---
fig2, ax2 = plt.subplots(1, 3, figsize=(18, 5))
ax2[0].plot(aoa_data, rpm, 'o')
ax2[0].set_xlabel('AoA [deg]', fontweight='bold')
ax2[0].set_ylabel('RPM', fontweight='bold')

ax2[1].plot(airspeed, rpm, 's')
ax2[1].set_xlabel('Airspeed [m/s]', fontweight='bold')
ax2[1].set_ylabel('RPM', fontweight='bold')

xi = np.linspace(aoa_data.min(), aoa_data.max(), 100)
yi = np.linspace(airspeed.min(), airspeed.max(), 100)
XI, YI = np.meshgrid(xi, yi)
ZI = griddata((aoa_data, airspeed), rpm, (XI, YI), method='cubic')
ax2[2].contourf(XI, YI, ZI, levels=20, cmap='viridis')
ax2[2].scatter(aoa_data, airspeed, c=rpm, cmap='viridis', edgecolor='k')
ax2[2].set_xlabel('AoA [deg]', fontweight='bold')
ax2[2].set_ylabel('Airspeed [m/s]', fontweight='bold')

for i, label in enumerate(['AoA vs RPM', 'Airspeed vs RPM', 'RPM by AoA & V']):
    ax2[i].set_title(label, fontweight='bold')
    ax2[i].grid(True)
plt.tight_layout()

# --- Plot: Model Comparison ---
models = [
    ('Linear Regression', pred_lin),
    ('Ridge Regression', pred_ridge),
    ('Random Forest', pred_rf)
]

fig3, axs3 = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
for ax, (name, pred) in zip(axs3, models):
    ax.scatter(f, pred, alpha=0.6)
    ax.plot([min(f), max(f)], [min(f), max(f)], 'r--')
    ax.set_title(f'{name}\n(f vs f_pred)', fontweight='bold')
    ax.set_xlabel('Actual f (T / RPM²)', fontweight='bold')
    ax.set_ylabel('Predicted f', fontweight='bold')
    ax.grid(True)
plt.tight_layout()

# --- Final Error Plot ---
fig4 = plt.figure(figsize=(18, 5))

ax1 = fig4.add_subplot(131)
sc = ax1.scatter(aoa_data, airspeed, c=f, cmap='viridis', edgecolor='k')
plt.colorbar(sc, ax=ax1, label='T / RPM²')
ax1.set_xlabel('AoA [deg]', fontweight='bold')
ax1.set_ylabel('Airspeed [m/s]', fontweight='bold')
ax1.set_title('T/RPM² Data', fontweight='bold')
ax1.grid(True)

ax2 = fig4.add_subplot(132, projection='3d')
ax2.scatter(aoa_data, airspeed, f, color='k', s=10, label='Actual Data')
ax2.plot_surface(AOA_grid, V_grid, Z_grid, cmap='viridis', alpha=0.6, zorder=0)
ax2.set_xlabel('AoA [deg]')
ax2.set_ylabel('Airspeed [m/s]')
ax2.set_zlabel('T / RPM²')
ax2.set_title('Linear Fit Surface', fontweight='bold')

ax3 = fig4.add_subplot(133)
err_sc = ax3.scatter(aoa_data, airspeed, c=error_f, cmap='coolwarm', edgecolor='k')
plt.colorbar(err_sc, ax=ax3, label='Error in f')
ax3.set_xlabel('AoA [deg]', fontweight='bold')
ax3.set_ylabel('Airspeed [m/s]', fontweight='bold')
ax3.set_title('Error (Actual - Predicted)', fontweight='bold')
ax3.grid(True)

plt.tight_layout()
plt.show()