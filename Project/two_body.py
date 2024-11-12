import numpy as np
import math

def kep_eqn(M: float, e: float):
    error = 1
    tol = 1e-8
    max_iter = 1000
    iteration = 0

    # Initial guess for E
    if -np.pi < M < 0 or M > np.pi:
        E_old = M - e
    else:
        E_old = M + e 

    # Iteratively solve for E using Newton's method
    while error > tol and iteration < max_iter:
        E_new = E_old + (M - E_old + e * np.sin(E_old)) / (1 - e * np.cos(E_old))
        error = np.abs(E_new - E_old)
        E_old = E_new
        iteration += 1

    E = E_new

    return E

def kep_eqn_par(dt: float, p: float, mu: float):
    n_p = 2 * math.sqrt(mu / p ** 3)
    cot_2s = 3 / 2 * n_p * dt
    s = np.arctan(1 / cot_2s) / 2
    tan_w = np.tan(s)
    w = np.arctan(tan_w ** (1 / 3))
    B = 2 / np.tan(2 * w)

    return B

def kep_eqn_hyp(M: float, e: float):
    error = 1
    tol = 1e-6
    max_iter = 1000
    iter = 0

    if e < 1.6:
        if -np.pi < M < 0 or M > np.pi:
            H_old = M - e
        else:
            H_old = M + e
    else:
        if e < 3.6 and np.abs(M) > np.pi:
            H_old = M - np.sign(M) * e
        else: 
            H_old = M / (e - 1)

    while error > tol and iter < max_iter:
        H_new = H_old + (M - e * np.sinh(H_old) + H_old) / (e * np.cosh(H_old) - 1)
        error = np.abs(H_new - H_old)
        H_old = H_new
        iter += 1

    H = H_new

    return H

def ta_to_anom(e: float, nu: float):
    tol = 1e-6

    if e < 1:
        y = (np.sin(nu) * math.sqrt(1 - e ** 2)) / (1 + e * np.cos(nu))
        x = (e + np.cos(nu)) / (1 + e * np.cos(nu))

        E = np.atan2(y, x)
        anom = E
    elif np.abs(e - 1) < tol:
        B = np.tan(nu / 2)
        anom = B
    elif e > 1:
        H = np.asinh((np.sin(nu) * math.sqrt(e ** 2 - 1)) / (1 + e * np.cos(nu)))
        anom = H
    
    return anom

def anom_to_ta(e, anom, p=None, r=None):
    tol = 1e-6

    if e < 1:
        E = anom
        y = (np.sin(E) * math.sqrt(1 - e ** 2)) / (1 - e * np.cos(E))
        x = (np.cos(E) - e) / (1 - e * np.cos(E))

    elif np.abs(e - 1) < tol:
        B = anom
        y = (p * B) / r
        x = (p - r) / r

    elif e > 1:
        H = anom
        y = (-np.sinh(H) * math.sqrt(e ** 2 - 1)) / (1 - e * np.cosh(H))
        x = (np.cosh(H) - e) / (1 - e * np.cosh(H))
    
    nu = np.atan2(y, x)

    return nu

def rv_to_koe(r_vec, v_vec, mu):
    tol = 1e-6

    # Calculate norms of position and velocity vectors
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    # Define the z-axis unit vector
    K_hat = np.array([0, 0, 1])

    # Calculate the specific angular momentum vector and its norm
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    
    # Calculate the node vector (pointing towards ascending node) and its norm
    n_vec = np.cross(K_hat, h_vec)
    n = np.linalg.norm(n_vec)

    # Calculate eccentricity vector and eccentricity
    e_vec = ((v ** 2 - mu / r) * r_vec - np.dot(r_vec, v_vec) * v_vec) / mu
    e = np.linalg.norm(e_vec)

    # Calculate specific orbital energy
    eps = v ** 2 / 2 - mu / r

    # Semi-major axis (a) or semi-latus rectum (p) for parabolic case
    if np.abs(e - 1) < tol:
        p = h ** 2 / mu  # Parabolic orbit
    else:
        a = -mu / (2 * eps)  # Non-parabolic orbit

    # Inclination (i)
    i = np.arccos(h_vec[2] / h)

    # Right Ascension of the Ascending Node (RAAN), W
    if n != 0:
        W = np.arccos(n_vec[0] / n)
        if n_vec[1] < 0:
            W = 2 * math.pi - W
    else:
        W = 0  # Circular equatorial orbit case, RAAN is undefined

    # Argument of Perigee (w)
    if e != 0:
        w = np.arccos(np.dot(n_vec, e_vec) / (n * e))
        if e_vec[2] < 0:
            w = 2 * math.pi - w
    else:
        w = 0  # Circular orbit case, Argument of Perigee is undefined

    # True anomaly (nu)
    nu = np.arccos(np.dot(e_vec, r_vec) / (e * r))
    if np.dot(r_vec, v_vec) < 0:
        nu = 2 * math.pi - nu

    # Special case handling: elliptical equatorial, circular inclined, and circular equatorial
    w_true = np.arccos(e_vec[0] / e)
    if e_vec[1] < 0:
        w_true = 2 * math.pi - w_true
    
    u = np.arccos(np.dot(n_vec, r_vec) / (n * r))
    if r_vec[2] < 0:
        u = 2 * math.pi - u

    lambda_true = np.arccos(r_vec[0] / r)  # True longitude
    if r_vec[1] < 0:
        lambda_true = 2 * math.pi - lambda_true

    return [
        p if np.abs(e - 1) < tol else a,  # Return p for parabolic, a for others
        e,            # Eccentricity (e)
        i,            # Inclination (i)
        W,            # Right Ascension of the Ascending Node (RAAN, W)
        w,            # Argument of Perigee (w)
        nu,           # True Anomaly (nu)
        lambda_true,  # True Longitude for circular equatorial orbits (not calculated)
        w_true        # True Argument of Perigee for elliptical equatorial cases
    ]

def koe_2_rv(a, e, i, W, w, nu, mu):
    tol = 1e-6  # Tolerance for checking near-zero values
    p = a * (1 - e**2)  # Semi-latus rectum

    # Define nu based on orbital configurations
    if np.abs(e - 0) < tol and np.abs(i - 0) < tol:
        # Circular equatorial case
        W = 0
        w = 0
    elif np.abs(e - 0) < tol and i > 0:
        # Circular inclined case
        w = 0
    elif 0 < e < 1 and np.abs(i - 0) < tol:
        # Elliptic equatorial case
        W = 0

    # Position vector in PQW frame
    r_PQW = np.array([
        [p * np.cos(nu) / (1 + e * np.cos(nu))],
        [p * np.sin(nu) / (1 + e * np.cos(nu))],
        [0]
    ])

    # Velocity vector in PQW frame
    v_PQW = np.array([
        [-np.sqrt(mu / p) * np.sin(nu)],
        [np.sqrt(mu / p) * (e + np.cos(nu))],
        [0]
    ])

    # Rotation matrix from IJK to PQW frame
    R_IJK_to_PQW = np.array([
        [np.cos(W) * np.cos(w) - np.sin(W) * np.sin(w) * np.cos(i),
         -np.cos(W) * np.sin(w) - np.sin(W) * np.cos(w) * np.cos(i),
         np.sin(W) * np.sin(i)],
        
        [np.sin(W) * np.cos(w) + np.cos(W) * np.sin(w) * np.cos(i),
         -np.sin(W) * np.sin(w) + np.cos(W) * np.cos(w) * np.cos(i),
         -np.cos(W) * np.sin(i)],
        
        [np.sin(w) * np.sin(i),
         np.cos(w) * np.sin(i),
         np.cos(i)]
    ])

    # Convert PQW vectors to IJK frame
    r_IJK = R_IJK_to_PQW @ r_PQW
    v_IJK = R_IJK_to_PQW @ v_PQW

    return [
        r_IJK.flatten(),
        v_IJK.flatten()
    ]


