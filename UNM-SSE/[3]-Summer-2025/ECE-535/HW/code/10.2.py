import satcom as sc

# -----------------------------
# Problem 14.9 – Number of FDMA Carriers with Backoff
# -----------------------------
sat_eirp = 25         # Saturation EIRP in dBW
backoff = 7           # Output backoff in dB
carrier_eirp = 12     # Required EIRP per carrier in dBW
transponder_bw = 27   # Transponder bandwidth in MHz
carrier_bw = 3        # Carrier bandwidth in MHz

num_carriers = sc.calc_fdma_carriers(sat_eirp, backoff, carrier_eirp, transponder_bw, carrier_bw)
print("Problem 14.9 – Number of FDMA Carriers with Backoff")
print(f"Number of carriers that can be accommodated: {num_carriers}\n")

# -----------------------------
# Problem 14.14 – Carriers With and Without Output Backoff
# -----------------------------
transponder_bw = 80     # Transponder bandwidth in MHz
carrier_bw = 6          # Carrier bandwidth in MHz
backoff = 6.5           # Output backoff in dB

no_backoff, with_backoff = sc.fdma_carrier_capacity(transponder_bw, carrier_bw, backoff)
print("Problem 14.14 – FDMA Carriers With and Without Backoff")
print(f"Carriers without backoff: {no_backoff}")
print(f"Carriers with 6.5 dB backoff: {with_backoff}\n")