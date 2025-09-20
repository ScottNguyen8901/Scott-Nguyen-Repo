import numpy as np

class RPMPredictor:
    def __init__(self):
        # Constants for drag calculation
        self.A = 0.625
        self.cd = 0.13
        self.rho = 1.225

        # Feature scaling parameters
        self.mean_aoa = -2.53270792
        self.var_aoa = 0.89707128
        self.mean_v = 25.83081732
        self.var_v = 3.87410052

        # Polynomial coefficients for f = D / RPM^2 (15 terms)
        self.coefs = np.array([
            1.277587e-06, -4.437110e-07, -4.640848e-07,
           -3.563343e-07, -1.260215e-06, -4.866570e-07,
            4.018992e-07,  1.742271e-06,  2.098580e-06,
            8.238411e-07,  3.148649e-07,  1.435648e-06,
            2.638828e-06,  2.153819e-06,  5.652796e-07
        ])

        # 3-sigma bounds for input sanity checks
        self.mean_std_bounds = {
            'pitch_angle': (-0.2808, 3.1342),
            'vel_x': (25.4471, 3.5163),
            'vel_y': (0.9641, 1.0884),
            'est_wind_x': (-0.3382, 2.8489),
            'est_wind_y': (0.0049, 0.2622)
        }

    def predict(self, pitch_angle, vel_x, vel_y, wind_x, wind_y):
        input_vars = {
            'pitch_angle': np.asarray(pitch_angle),
            'vel_x': np.asarray(vel_x),
            'vel_y': np.asarray(vel_y),
            'est_wind_x': np.asarray(wind_x),
            'est_wind_y': np.asarray(wind_y)
        }

        N = len(pitch_angle)
        rpm_flag = np.full(N, 'GOOD', dtype='<U4')

        try:
            for key, (mu, sigma) in self.mean_std_bounds.items():
                val = input_vars[key]
                lower, upper = mu - 3 * sigma, mu + 3 * sigma
                bad_mask = (val < lower) | (val > upper)
                if np.any(bad_mask):
                    print(f"\nWarning: {key} has {np.sum(bad_mask)} values out of bounds:")
                    print(f"  Expected range: [{lower:.2f}, {upper:.2f}]")
                    print(f"  Out-of-bounds values: {val[bad_mask]}")
                    rpm_flag[bad_mask] = 'BAD'
        except Exception as e:
            print(f"3σ check failed with error: {e}")

        ax = input_vars['vel_x'] - input_vars['est_wind_x']
        ay = input_vars['vel_y'] - input_vars['est_wind_y']
        airspeed = np.hypot(ax, ay)
        fpa = np.degrees(np.arctan2(ay, ax))
        aoa = input_vars['pitch_angle'] - fpa

        # AoA and airspeed bounds check
        lower_aoa, upper_aoa = self.mean_aoa - 3 * np.sqrt(self.var_aoa), self.mean_aoa + 3 * np.sqrt(self.var_aoa)
        lower_v, upper_v = self.mean_v - 3 * np.sqrt(self.var_v), self.mean_v + 3 * np.sqrt(self.var_v)

        aoa_bad = (aoa < lower_aoa) | (aoa > upper_aoa)
        v_bad = (airspeed < lower_v) | (airspeed > upper_v)

        if np.any(aoa_bad):
            print(f"\nWarning: aoa has {np.sum(aoa_bad)} values out of bounds:")
            print(f"  Expected range: [{lower_aoa:.2f}, {upper_aoa:.2f}]")
            print(f"  Out-of-bounds values: {aoa[aoa_bad]}")
            rpm_flag[aoa_bad] = 'BAD'

        if np.any(v_bad):
            print(f"\nWarning: airspeed has {np.sum(v_bad)} values out of bounds:")
            print(f"  Expected range: [{lower_v:.2f}, {upper_v:.2f}]")
            print(f"  Out-of-bounds values: {airspeed[v_bad]}")
            rpm_flag[v_bad] = 'BAD'

        aoa_s = (aoa - self.mean_aoa) / np.sqrt(self.var_aoa)
        v_s = (airspeed - self.mean_v) / np.sqrt(self.var_v)

        X = np.column_stack([
            np.ones_like(aoa_s), aoa_s, v_s,
            aoa_s**2, aoa_s * v_s, v_s**2,
            aoa_s**3, aoa_s**2 * v_s, aoa_s * v_s**2, v_s**3,
            aoa_s**4, aoa_s**3 * v_s, aoa_s**2 * v_s**2,
            aoa_s * v_s**3, v_s**4
        ])

        f_pred = np.dot(X, self.coefs)
        D = 0.5 * self.rho * airspeed**2 * self.A * self.cd
        rpm_pred = np.sqrt(D / f_pred)

        return rpm_pred, rpm_flag