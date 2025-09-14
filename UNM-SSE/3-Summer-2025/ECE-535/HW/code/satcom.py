import math

def effective_noise_temp(T_amp_K, atten_dB, T0_K=290):
    """
    Calculates the effective noise temperature referred to the input of an attenuator.

    Parameters:
        T_amp_K (float): Noise temperature of the amplifier in Kelvin.
        atten_dB (float): Attenuator loss in dB.
        T0_K (float): Reference temperature, default is 290 K.

    Returns:
        float: Effective noise temperature referred to the attenuator input in Kelvin.
    """
    # Convert attenuator loss from dB to linear scale
    L = 10 ** (atten_dB / 10)

    # Noise temperature contribution from the attenuator
    T_att = T0_K * (L - 1)

    # Amplifier noise temperature referred to the attenuator input
    T_amp_eff = T_amp_K / L

    # Total effective noise temperature
    T_eq = T_att + T_amp_eff

    return T_eq

def cn_calc(Pc_watts, T_sys_K, B_Hz):
    """
    Compute C/N0 (dBHz) and C/N (dB) from carrier power, system noise temp, and bandwidth.
    
    Parameters:
        Pc_watts (float): Carrier power in watts.
        T_sys_K (float): System noise temperature in Kelvin.
        B_Hz (float): Bandwidth in Hz.
    
    Returns:
        tuple: (C/N0 in dBHz, C/N in dB)
    """
    k = 1.38e-23  # Boltzmann constant

    cn0_dBHz = 10 * math.log10(Pc_watts / (k * T_sys_K))
    cn_dB = cn0_dBHz - 10 * math.log10(B_Hz)

    return cn0_dBHz, cn_dB

def cn_dB(EIRP_dBW, G_T_dB, L_dB, B_Hz):
    """
    Computes carrier-to-noise ratio [C/N] in dB for a satellite link.

    Parameters:
        EIRP_dBW (float): Effective Isotropic Radiated Power in dBW.
        G_T_dB (float): Receiver G/T in dB/K.
        L_dB (float): Total link loss in dB.
        B_Hz (float): System bandwidth in Hz.

    Returns:
        float: C/N in dB.
    """
    k_dB = 10 * math.log10(1.38e-23)  # Boltzmann constant in dB
    B_dBHz = 10 * math.log10(B_Hz)

    CN_dB = EIRP_dBW + G_T_dB - L_dB - k_dB - B_dBHz
    return CN_dB

def calc_flux_density(P_r_watts, D_m, to_dBW=True):
    """
    Calculates power flux density (W/m² and dBW/m²) at a parabolic antenna.

    Parameters:
        P_r_watts (float): Received power in watts.
        D_m (float): Antenna diameter in meters.
        to_dBW (bool): If True, returns dBW/m² as well.

    Returns:
        tuple: (flux density in W/m², flux density in dBW/m² if to_dBW is True)
    """
    # Antenna aperture area (A = πD²/4)
    A_m2 = math.pi * (D_m ** 2) / 4

    # Power flux density in W/m²
    flux_W_m2 = P_r_watts / A_m2

    if to_dBW:
        flux_dBW_m2 = 10 * math.log10(flux_W_m2)
        return flux_W_m2, flux_dBW_m2
    else:
        return flux_W_m2
    
def required_eirp(flux_density_dBW_m2, total_losses_dB):
    """
    Calculates required EIRP to achieve a given flux density over a link with known losses.

    Parameters:
        flux_density_dBW_m2 (float): Desired flux density in dBW/m².
        total_losses_dB (float): Total path and system losses in dB.

    Returns:
        float: Required EIRP in dBW.
    """
    return flux_density_dBW_m2 + total_losses_dB

def calc_fdma_carriers(sat_eirp_dBW, backoff_dB, carrier_eirp_dBW, transponder_bw_MHz, carrier_bw_MHz):
    """
    Calculates the number of FDMA carriers a satellite transponder can support.

    Parameters:
        sat_eirp_dBW (float): Saturation EIRP of the transponder in dBW.
        backoff_dB (float): Output backoff in dB.
        carrier_eirp_dBW (float): EIRP required per carrier in dBW.
        transponder_bw_MHz (float): Total transponder bandwidth in MHz.
        carrier_bw_MHz (float): Bandwidth per carrier in MHz.

    Returns:
        int: Number of carriers the transponder can accommodate.
    """
    available_eirp_dBW = sat_eirp_dBW - backoff_dB
    num_carriers_power = 10 ** ((available_eirp_dBW - carrier_eirp_dBW) / 10)
    num_carriers_bw = transponder_bw_MHz / carrier_bw_MHz
    return int(min(num_carriers_power, num_carriers_bw))

def fdma_carrier_capacity(transponder_bw_MHz, carrier_bw_MHz, backoff_dB=0):
    """
    Calculates number of FDMA carriers with and without backoff.

    Parameters:
        transponder_bw_MHz (float): Total transponder bandwidth in MHz.
        carrier_bw_MHz (float): Bandwidth required per carrier in MHz.
        backoff_dB (float): Output backoff in dB (default is 0).

    Returns:
        tuple: (carriers_without_backoff, carriers_with_backoff)
    """
    carriers_no_backoff = int(transponder_bw_MHz // carrier_bw_MHz)

    if backoff_dB > 0:
        linear_loss = 10 ** (backoff_dB / 10)
        effective_bw_per_carrier = carrier_bw_MHz * linear_loss
        carriers_with_backoff = int(transponder_bw_MHz // effective_bw_per_carrier)
    else:
        carriers_with_backoff = carriers_no_backoff

    return carriers_no_backoff, carriers_with_backoff

def cascaded_amp_noise(G1_dB, G2_dB, Te1, Te2):
    # Convert gains from dB to linear
    G1 = 10**(G1_dB / 10)
    G2 = 10**(G2_dB / 10)
    
    # Overall gain in dB and linear
    G_total_dB = G1_dB + G2_dB
    G_total = G1 * G2

    # Friis formula for effective noise temperature
    T_total = Te1 + Te2 / G1

    return G_total_dB, G_total, T_total

def fdma_downlink_cn(EIRP_sat_dB, G_T_dB, L_dB, B_Hz, backoff_dB, total_bw_Hz):
    k_dB = -228.6
    EIRP_actual = EIRP_sat_dB - backoff_dB
    B_dBHz = 10 * math.log10(B_Hz)
    
    CN_dB = EIRP_actual - L_dB + G_T_dB - k_dB - B_dBHz

    # Number of carriers
    carrier_bw = B_Hz
    N_total = total_bw_Hz / carrier_bw
    power_ratio = 10**(-backoff_dB / 10)
    N_with_backoff = int(N_total * power_ratio)
    
    return CN_dB, int(N_total), N_with_backoff

def calc_cn(P_c_watts, T_sys_K, bandwidth_Hz):
    k = 1.38e-23  # Boltzmann constant (J/K)
    
    cn0_linear = P_c_watts / (k * T_sys_K)
    cn0_dBHz = 10 * math.log10(cn0_linear)
    
    B_dBHz = 10 * math.log10(bandwidth_Hz)
    cn_dB = cn0_dBHz - B_dBHz
    
    return cn0_dBHz, cn_dB

def calculate_cn0(eirp_dBW, fs_loss_dB, pointing_loss_dB, atm_loss_dB, feeder_loss_dB, gt_dB):
    k_dB = -228.6  # 10*log10(k), Boltzmann's constant in dB
    total_loss = fs_loss_dB + pointing_loss_dB + atm_loss_dB + feeder_loss_dB
    cn0 = eirp_dBW - total_loss + gt_dB - k_dB
    return cn0