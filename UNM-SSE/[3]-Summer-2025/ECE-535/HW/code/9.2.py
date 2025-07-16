import math
import numpy as np
import satcom as sc

# -----------------------------
# Problem 12.15 – Effective Noise Temperature
# -----------------------------
T_amp = 200         # Amplifier noise temperature in K
atten_dB = 4        # Attenuator loss in dB

T_effective = sc.effective_noise_temp(T_amp, atten_dB)
print("Problem 12.15 – Effective Noise Temperature")
print(f"Effective noise temperature referred to the attenuator input: {T_effective:.2f} K\n")

# -----------------------------
# Problem 12.19 – Carrier-to-Noise and Carrier-to-Noise Density Ratios
# -----------------------------
Pc = 400e-12        # Carrier power in W (400 pW)
T_sys = 450         # System noise temperature in K
B = 36e6            # Bandwidth in Hz (36 MHz)

cn0, cn = sc.cn_calc(Pc, T_sys, B)
print("Problem 12.19 – C/N0 and C/N")
print(f"C/N0 = {cn0:.2f} dBHz")
print(f"C/N  = {cn:.2f} dB\n")

# -----------------------------
# Problem 12.21 – Received C/N Using Link Budget
# -----------------------------
EIRP = 45           # dBW
GT = 11             # dB/K
L_total = 203       # Total link losses in dB (200 dB + 3 dB margin)

cn_link = sc.cn_dB(EIRP, GT, L_total, B)
print("Problem 12.21 – Received C/N from Link Budget")
print(f"C/N = {cn_link:.2f} dB\n")

# -----------------------------
# Problem 12.23 – Power Flux Density at Antenna
# -----------------------------
P_r = 250e-12       # Received power in W (250 pW)
D = 1.8             # Antenna diameter in meters

F_W_m2, F_dBW_m2 = sc.calc_flux_density(P_r, D)
print("Problem 12.23 – Power Flux Density")
print(f"(a) Power flux density = {F_W_m2:.4e} W/m²")
print(f"(b) Power flux density = {F_dBW_m2:.2f} dBW/m²\n")

# -----------------------------
# Problem 12.25 – Required EIRP for Saturation Flux
# -----------------------------
F_sat = -110        # Saturation flux density in dBW/m²
L_total = 200       # Total link losses in dB

eirp_needed = sc.required_eirp(F_sat, L_total)
print("Problem 12.25 – Required EIRP for Saturation Flux")
print(f"Required EIRP: {eirp_needed:.2f} dBW")