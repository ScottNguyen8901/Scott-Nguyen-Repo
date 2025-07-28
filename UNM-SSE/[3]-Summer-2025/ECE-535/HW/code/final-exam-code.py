import satcom as sc

# Given values
G1_dB = 15
G2_dB = 15
Te1 = 190  # K
Te2 = 190  # K

# Compute
gain_dB, gain_linear, T_effective = sc.cascaded_amp_noise(G1_dB, G2_dB, Te1, Te2)

# Print results
print(f"Overall Gain: {gain_dB} dB or {gain_linear:.0f} linear")
print(f"Effective Noise Temperature referred to input: {T_effective:.2f} K")

EIRP_sat = 25       # dBW
G_T = 30            # dB/K
L = 196             # dB
B = 3e6             # Hz
backoff = 6         # dB
total_bw = 36e6     # Hz

cn_dB, n_max, n_backoff = sc.fdma_downlink_cn(EIRP_sat, G_T, L, B, backoff, total_bw)

print(f"Single-carrier C/N: {cn_dB:.2f} dB")
print(f"Max carriers without backoff: {n_max}")
print(f"Carriers with 6 dB backoff: {n_backoff}")

# Inputs
P_c = 400e-12         # Received carrier power in watts (400 pW)
T_sys = 450           # System noise temperature in K
B = 36e6              # Bandwidth in Hz (36 MHz)

# Call the function
cn0, cn = sc.calc_cn(P_c, T_sys, B)

# Print results
print(f"C/N₀ = {cn0:.2f} dBHz")
print(f"C/N  = {cn:.2f} dB")

# Inputs from the problem
eirp = 45               # dBW
fs_loss = 206           # dB
pointing_loss = 1       # dB
atm_loss = 2            # dB
feeder_loss = 1         # dB
gt_ratio = 19.5         # dB/K

# Calculate C/N0
cn0_result = sc.calculate_cn0(eirp, fs_loss, pointing_loss, atm_loss, feeder_loss, gt_ratio)
print(f"C/N₀ = {cn0_result:.2f} dBHz")