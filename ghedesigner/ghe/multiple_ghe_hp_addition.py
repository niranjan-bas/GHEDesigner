# This module is developed to enable ghedesigner to handle multiple ground heat exchanger and heat pump layouts. Most
# part of the code is custom-built for 3ghe-6hp layout.

import pandas as pd
from pathlib import Path
from ghedesigner.ghe.gfunction import calc_g_func_for_multiple_lengths
from scipy.interpolate import interp1d
import numpy as np
from ghedesigner.enums import BHPipeType
from pygfunction.boreholes import Borehole
from ghedesigner.media import Grout, Pipe, Soil
from ghedesigner.ghe.simulation import SimulationParameters
from ghedesigner.ghe.coaxial_borehole import get_bhe_object
from ghedesigner.constants import VERSION


repo_root = Path(__file__).resolve().parents[2]
file_path_hp_loads = repo_root / "buildingloads" / "MultipleGHE-HP_heatpumploads_DULUTH_3ghe-6hp.csv"
file_path_cases = repo_root / "cases" / "3ghe-6hp_case.csv"

selected_configuration = "3ghe-6hp_1-pipe"
selected_case = "Case: 1-pipe"


def read_hp_loads(hp_loads_path):
    """Reads HP loads from CSVs with MultiIndex columns."""
    df1 = pd.read_csv(hp_loads_path, header=[0, 1, 2, 3])
    df1.columns = pd.MultiIndex.from_tuples([(int(l1), int(l2), int(l3),  l4) for l1, l2, l3, l4 in df1.columns])

    return df1


df1 = read_hp_loads(file_path_hp_loads)


def extract_time_array():
    """
    Extracts the time values from the specified location in the dataframe
    and converts them into a NumPy array.

    Returns:
    A NumPy array containing the extracted time values.
    """
    array = df1.iloc[:, 0].to_numpy()
    return array


def time_array_size():
    size = len(extract_time_array())
    return size


time_array = extract_time_array()
n = time_array_size()


def read_case_data(cases_path, selected_configuration: str, selected_case: str):
    if selected_configuration == "3ghe-6hp_1-pipe":

        """Reads selected case data and extracts parameters."""
        df_cases = pd.read_csv(cases_path)
        case_data = df_cases[df_cases["Case"] == selected_case].iloc[0]

        # Borehole layout
        n_rows1, n_cols1 = int(case_data["n_rows1"]), int(case_data["n_columns1"])
        n_rows2, n_cols2 = int(case_data["n_rows2"]), int(case_data["n_columns2"])
        n_rows3, n_cols3 = int(case_data["n_rows3"]), int(case_data["n_columns3"])
        b_spacing = case_data["b_spacing"]

        # Heat pump design flowrates
        (m_hp1_design, m_hp2_design, m_hp3_design,
        m_hp4_design, m_hp5_design, m_hp6_design) = (
        case_data["m_hp1"], case_data["m_hp2"], case_data["m_hp3"],
        case_data["m_hp4"], case_data["m_hp5"], case_data["m_hp6"]
        )

        # Boreholes per GHE and GHE flowrates
        nbh1, nbh2, nbh3 = case_data["nbh_1GHE"], case_data["nbh_2GHE"], case_data["nbh_3GHE"]
        m_ghe1, m_ghe2, m_ghe3 = case_data["m_ghe1"], case_data["m_ghe2"], case_data["m_ghe3"]

        # Single HP values
        (m_hp1_single_as, m_hp2_single_as, m_hp3_single_as, m_hp_single_wtw) = (
            case_data["m_hp1_single_WTA"], case_data["m_hp2_single_WTA"],
            case_data["m_hp3_single_WTA"], case_data["m_hp_single_WTW"]
        )

        return (n_rows1, n_cols1, n_rows2, n_cols2, n_rows3, n_cols3, b_spacing,
                m_hp1_design, m_hp2_design, m_hp3_design, m_hp4_design, m_hp5_design, m_hp6_design,
                nbh1, nbh2, nbh3, m_ghe1, m_ghe2, m_ghe3,
                m_hp1_single_as, m_hp2_single_as, m_hp3_single_as, m_hp_single_wtw)
    else:
        print("Configuration not recognized")


if selected_configuration == "3ghe-6hp_1-pipe":
    (n_rows1, n_cols1, n_rows2, n_cols2, n_rows3, n_cols3, b_spacing, m_hp1_design, m_hp2_design, m_hp3_design,
    m_hp4_design, m_hp5_design, m_hp6_design, nbh1, nbh2, nbh3, m_ghe1, m_ghe2, m_ghe3, m_hp1_single_as, m_hp2_single_as,
    m_hp3_single_as, m_hp_single_wtw) = (read_case_data(file_path_cases, selected_configuration, selected_case))

else:
    print("Configuration not recognized")


def generate_g_functions_for_all_ghe(
    bhe1, bhe2, bhe3,
    m_ghe1_borehole, m_ghe2_borehole, m_ghe3_borehole,
    bhe_type, log_time, b_spacing,
    n_rows1, n_cols1, n_rows2, n_cols2, n_rows3, n_cols3,
    calc_g_func_for_multiple_lengths
):
    h_values = [100.0]  # fixed h value

    r_b1 = bhe1.calc_effective_borehole_resistance()
    r_b2 = bhe2.calc_effective_borehole_resistance()
    r_b3 = bhe3.calc_effective_borehole_resistance()
    depth = bhe1.b.D

    coordinates_ghe1 = [(i * b_spacing, j * b_spacing) for i in range(n_rows1) for j in range(n_cols1)]
    coordinates_ghe2 = [(i * b_spacing, j * b_spacing) for i in range(n_rows2) for j in range(n_cols2)]
    coordinates_ghe3 = [(i * b_spacing, j * b_spacing) for i in range(n_rows3) for j in range(n_cols3)]

    gFunction1 = calc_g_func_for_multiple_lengths(
        b_spacing, h_values, r_b1, depth, m_ghe1_borehole, bhe_type, log_time,
        coordinates_ghe1, bhe1.fluid, bhe1.pipe, bhe1.grout, bhe1.soil
    )

    gFunction2 = calc_g_func_for_multiple_lengths(
        b_spacing, h_values, r_b2, depth, m_ghe2_borehole, bhe_type, log_time,
        coordinates_ghe2, bhe2.fluid, bhe2.pipe, bhe2.grout, bhe2.soil
    )

    gFunction3 = calc_g_func_for_multiple_lengths(
        b_spacing, h_values, r_b3, depth, m_ghe3_borehole, bhe_type, log_time,
        coordinates_ghe3, bhe3.fluid, bhe3.pipe, bhe3.grout, bhe3.soil
    )

    return gFunction1, gFunction2, gFunction3


def calculation_of_ghe_constant_c_n(g1, g2, g3, ts, k_soil, rb1, rb2, rb3, time_array, n):

    """
    Calculate C_n values for three GHEs based on their g-functions.

    Cn = 1 / (2 * pi * K_s) * g((tn - tn-1) / t_s) + R_b
    """

    two_pi_k = 2 * np.pi * k_soil
    c_n1 = np.zeros(n, dtype=float)
    c_n2 = np.zeros(n, dtype=float)
    c_n3 = np.zeros(n, dtype=float)

    for i in range(1, n):
        delta_log_time = np.log((time_array[i] - time_array[i - 1]) / (ts / 3600))
        g_val1 = g1(delta_log_time)
        g_val2 = g2(delta_log_time)
        g_val3 = g3(delta_log_time)

        c_n1[i] = (1 / two_pi_k * g_val1) + rb1
        c_n2[i] = (1 / two_pi_k * g_val2) + rb2
        c_n3[i] = (1 / two_pi_k * g_val3) + rb3

    return c_n1, c_n2, c_n3


# Read the file for Heat Pump Coefficients
def get_coeffs(param, df):
    """Returns a tuple of values for hp1 through hp6 for the given parameter."""
    return tuple(df.loc[param, f"hp{i}"] for i in range(1, 7))


def load_all_hp_coefficients():
    """
        Reads the heat pump coefficients CSV and returns a dictionary of coefficient tuples
        for all heating and cooling parameters.
        """
    file_path_hp_coefficients = repo_root / "HeatPumpData" / "HeatPumpCoefficients.csv"
    df = pd.read_csv(file_path_hp_coefficients).set_index("parameter")

    hp_coefficients = {
        "c1_htg": get_coeffs("c1_htg", df),
        "c2_htg": get_coeffs("c2_htg", df),
        "c3_htg": get_coeffs("c3_htg", df),
        "c1_clg": get_coeffs("c1_clg", df),
        "c2_clg": get_coeffs("c2_clg", df),
        "c3_clg": get_coeffs("c3_clg", df),
        "a_htg": get_coeffs("a_htg", df),
        "b_htg": get_coeffs("b_htg", df),
        "c_htg": get_coeffs("c_htg", df),
        "a_clg": get_coeffs("a_clg", df),
        "b_clg": get_coeffs("b_clg", df),
        "c_clg": get_coeffs("c_clg", df),
    }

    return hp_coefficients

def get_net_htg_loads(df1):
    """Returns a tuple of 6 arrays for net heating loads of all HPs."""
    return (
        df1.loc[:, (1, 1, 1, "HPHtgLd_W")] - df1.loc[:, (1, 1, 1, "HPClgLd_W")].to_numpy(),
        df1.loc[:, (1, 2, 2, "HPHtgLd_W")].to_numpy(),
        df1.loc[:, (2, 1, 3, "HPHtgLd_W")] - df1.loc[:, (2, 1, 3, "HPClgLd_W")].to_numpy(),
        df1.loc[:, (2, 2, 4, "HPHtgLd_W")].to_numpy(),
        df1.loc[:, (3, 1, 5, "HPHtgLd_W")] - df1.loc[:, (3, 1, 5, "HPClgLd_W")].to_numpy(),
        df1.loc[:, (3, 2, 6, "HPHtgLd_W")].to_numpy(),
    )

def get_hourly_h_c(i, df1):
    """Returns heating and cooling loads at time step i for all 6 HPs."""
    h = [
        df1.loc[i, (1, 1, 1, "HPHtgLd_W")],
        df1.loc[i, (1, 2, 2, "HPHtgLd_W")],
        df1.loc[i, (2, 1, 3, "HPHtgLd_W")],
        df1.loc[i, (2, 2, 4, "HPHtgLd_W")],
        df1.loc[i, (3, 1, 5, "HPHtgLd_W")],
        df1.loc[i, (3, 2, 6, "HPHtgLd_W")]
    ]
    c = [
        df1.loc[i, (1, 1, 1, "HPClgLd_W")],
        0,
        df1.loc[i, (2, 1, 3, "HPClgLd_W")],
        0,
        df1.loc[i, (3, 1, 5, "HPClgLd_W")],
        0
    ]
    return np.array(h), np.array(c)


def calculate_r1_r2(t_eft, h, c, a_htg, b_htg, c_htg, a_clg, b_clg, c_clg):
    """Calculate r1 and r2 arrays for all 6 heat pumps."""
    h_array = np.array(h)
    c_array = np.array(c)

    a_htg = np.array(a_htg)
    b_htg = np.array(b_htg)
    c_htg = np.array(c_htg)

    a_clg = np.array(a_clg)
    b_clg = np.array(b_clg)
    c_clg = np.array(c_clg)

    # Heating calculations
    slope_htg = 2 * a_htg * t_eft + b_htg
    ratio_htg = a_htg * t_eft**2 + b_htg * t_eft + c_htg
    u = ratio_htg - slope_htg * t_eft
    v = slope_htg

    # Cooling calculations
    slope_clg = 2 * a_clg * t_eft + b_clg
    ratio_clg = a_clg * t_eft**2 + b_clg * t_eft + c_clg
    a = ratio_clg - slope_clg * t_eft
    b = slope_clg

    # Final results
    r1 = v * h_array - b * c_array
    r2 = u * h_array - a * c_array

    return r1, r2


def calculate_hp_capacity(
    t_eft,
    c1_htg, c2_htg, c3_htg,
    c1_clg, c2_clg, c3_clg,
    m_design, m_single,
    q_net_htg
):
    """
    Returns the array of heat pump capacities for heating or cooling.
    """

    capacity_htg = c1_htg * t_eft**2 + c2_htg * t_eft + c3_htg
    capacity_clg = c1_clg * t_eft**2 + c2_clg * t_eft + c3_clg

    # Choose heating or cooling based on sign of load
    hp_capacity = np.where(q_net_htg > 0, capacity_htg, capacity_clg) * (m_design / m_single)

    return hp_capacity


def calculate_rtf_and_mass_flows(
    q_net_htg, t_eft,
    c1_htg, c2_htg, c3_htg,
    c1_clg, c2_clg, c3_clg,
    m_design, m_single,
    m_hp1_design, m_hp2_design, m_hp3_design,
    m_hp4_design, m_hp5_design, m_hp6_design,
    nbh_ghe1, nbh_ghe2, nbh_ghe3, nbh_total,
    beta=1.5
):
    """
    Computes runtime fractions, heat pump mass flow rates, loop flow, and GHE flows.

    Returns:
        rtf: array of runtime fractions [6]
        m_dot: array of HP mass flows [6]
        m_loop: total loop mass flow [scalar]
        m_ghe1, m_ghe2, m_ghe3: split GHE mass flows [scalars]
    """

    # Calculate heating and cooling capacity
    capacity_htg = c1_htg * t_eft**2 + c2_htg * t_eft + c3_htg
    capacity_clg = c1_clg * t_eft**2 + c2_clg * t_eft + c3_clg

    # Select appropriate capacity depending on heating or cooling
    hp_capacity = np.where(q_net_htg > 0, capacity_htg, capacity_clg) * (m_design / m_single)

    # Compute run-time fraction (clip to [0, 1])
    with np.errstate(divide='ignore', invalid='ignore'):
        rtf = np.abs(q_net_htg) / np.abs(hp_capacity)
        rtf = np.where(np.isnan(rtf), 0, rtf)
        rtf = np.clip(rtf, 0, 1)

    # Individual HP mass flow rates
    m_dot = np.array([
        m_hp1_design * rtf[0],
        m_hp2_design * rtf[1],
        m_hp3_design * rtf[2],
        m_hp4_design * rtf[3],
        m_hp5_design * rtf[4],
        m_hp6_design * rtf[5]
    ])

    # Loop mass flow rate
    m_loop = beta * np.sum(m_dot)

    # GHE mass flow rates
    m_ghe1 = m_loop * nbh_ghe1 / nbh_total
    m_ghe2 = m_loop * nbh_ghe2 / nbh_total
    m_ghe3 = m_loop * nbh_ghe3 / nbh_total

    return rtf, m_dot, m_loop, m_ghe1, m_ghe2, m_ghe3


def compute_ghe_history_terms(
    i, time_array, ts, two_pi_k,
    q_ghe1, q_ghe2, q_ghe3,
    g1, g2, g3,
    tg,
    total_values_ghe1, total_values_ghe2, total_values_ghe3,
    H_n_ghe1, H_n_ghe2, H_n_ghe3
):
    time_n = time_array[i]

    # Compute dimensionless time for all indices from 1 to i-1
    indices = np.arange(1, i)
    dim_less_time = np.log((time_n - time_array[indices - 1]) / (ts / 3600))

    # Compute contributions from all previous steps
    delta_q_ghe1 = (q_ghe1[indices] - q_ghe1[indices - 1]) / two_pi_k
    delta_q_ghe2 = (q_ghe2[indices] - q_ghe2[indices - 1]) / two_pi_k
    delta_q_ghe3 = (q_ghe3[indices] - q_ghe3[indices - 1]) / two_pi_k

    values_ghe1 = np.sum(delta_q_ghe1 * g1(dim_less_time))
    values_ghe2 = np.sum(delta_q_ghe2 * g2(dim_less_time))
    values_ghe3 = np.sum(delta_q_ghe3 * g3(dim_less_time))

    total_values_ghe1[i] = values_ghe1
    total_values_ghe2[i] = values_ghe2
    total_values_ghe3[i] = values_ghe3

    # Dimensionless time for the immediate previous step
    dim1_less_time = np.log((time_n - time_array[i - 1]) / (ts / 3600))
    H_n_ghe1[i] = tg - total_values_ghe1[i] + (q_ghe1[i - 1] / two_pi_k * g1(dim1_less_time))
    H_n_ghe2[i] = tg - total_values_ghe2[i] + (q_ghe2[i - 1] / two_pi_k * g2(dim1_less_time))
    H_n_ghe3[i] = tg - total_values_ghe3[i] + (q_ghe3[i - 1] / two_pi_k * g3(dim1_less_time))


def solve_temperature_matrix(
    i,
    r1_vals, r2_vals, m_loop, cp,
    m_ghe1, m_ghe2, m_ghe3,
    c_n_ghe1, c_n_ghe2, c_n_ghe3,
    H_n_ghe1, H_n_ghe2, H_n_ghe3,
    nbh_ghe1, nbh_ghe2, nbh_ghe3,
    bhe1, bhe2, bhe3,
    time_array, labels, units,
    t, q_ghe1, q_ghe2, q_ghe3, results
):
    # Unpack r1 and r2
    r1_hp1, r1_hp2, r1_hp3, r1_hp4, r1_hp5, r1_hp6 = r1_vals
    r2_hp1, r2_hp2, r2_hp3, r2_hp4, r2_hp5, r2_hp6 = r2_vals

    # Assemble matrix A
    A = np.array([
        [1 - (r1_hp1 / (m_loop * cp)), -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1 - (r1_hp2 / (m_loop * cp)), -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1 - (r1_hp3 / (m_loop * cp)), -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1 - (r1_hp4 / (m_loop * cp)), -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1 - (r1_hp5 / (m_loop * cp)), -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1 - (r1_hp6 / (m_loop * cp)), -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],

        [0, 0, 0, 0, 0, 0, (m_loop - m_ghe1) * cp, -m_loop * cp, 0, 0, 0, 0, m_ghe1 * cp, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, c_n_ghe1[i], 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0, 0, 2, 0, 0, -1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, m_ghe1 * cp, 0, 0, 0, 0, 0, -m_ghe1 * cp, 0, 0, bhe1.b.H * nbh_ghe1, 0, 0],

        [0, 0, 0, 0, 0, 0, 0, (m_loop - m_ghe2) * cp, -m_loop * cp, 0, 0, 0, 0, m_ghe2 * cp, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, c_n_ghe2[i], 0],
        [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 2, 0, 0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, m_ghe2 * cp, 0, 0, 0, 0, 0, -m_ghe2 * cp, 0, 0, bhe2.b.H * nbh_ghe2, 0],

        [-m_loop * cp, 0, 0, 0, 0, 0, 0, 0, (m_loop - m_ghe3) * cp, 0, 0, 0, 0, 0, m_ghe3 * cp, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, c_n_ghe3[i]],
        [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 2, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, m_ghe3 * cp, 0, 0, 0, 0, 0, -m_ghe3 * cp, 0, 0, bhe3.b.H * nbh_ghe3]
    ])

    # Vector B
    B = np.array([
        r2_hp1 / (m_loop * cp), r2_hp2 / (m_loop * cp), r2_hp3 / (m_loop * cp),
        r2_hp4 / (m_loop * cp), r2_hp5 / (m_loop * cp), r2_hp6 / (m_loop * cp),
        0, H_n_ghe1[i], 0, 0, 0, H_n_ghe2[i], 0, 0, 0, H_n_ghe3[i], 0, 0
    ])

    # Solve linear system
    X = np.linalg.solve(A, B)

    # Update solution arrays
    t[i, :] = X[:15]
    q_ghe1[i] = X[15]
    q_ghe2[i] = X[16]
    q_ghe3[i] = X[17]

    # Save current row in results
    row = {"Time Step": time_array[i]}
    row.update({label: value for label, value in zip(labels, X)})
    results.append(row)

    return X


class MultiGHEHP:
    def __init__(
            self,
            v_flow_system: float,
            b_spacing: float,
            bhe_type: BHPipeType,
            fluid,
            borehole: Borehole,
            pipe: Pipe,
            grout: Grout,
            soil: Soil,
            sim_params: SimulationParameters,
            hourly_extraction_ground_loads: list,
            field_type="N/A",
            field_specifier="N/A",
    ) -> None:
        self.fieldType = field_type
        self.fieldSpecifier = field_specifier
        self.V_flow_system = v_flow_system
        self.B_spacing = b_spacing

        # Additional simulation parameters
        self.sim_params = sim_params
        # Hourly ground extraction loads
        # Building cooling is negative, building heating is positive
        self.hourly_extraction_ground_loads = hourly_extraction_ground_loads
        self.times = np.empty((0,), dtype=np.float64)
        self.loading = None

        # This is layout specific (My addition)
        self.n_rows1, self.n_rows2, self.n_rows3 = n_rows1, n_rows2, n_rows3,
        self.n_cols1, self.n_cols2, self.n_cols3 = n_cols1, n_cols2, n_cols3
        self.nbh1, self.nbh2, self.nbh3 = nbh1, nbh2, nbh3
        self.m_ghe1, self.m_ghe2, self.m_ghe3 = m_ghe1, m_ghe2, m_ghe3

        self.m_ghe1_borehole = float(self.m_ghe1) / float(self.nbh1)
        self.m_ghe2_borehole = float(self.m_ghe2) / float(self.nbh2)
        self.m_ghe3_borehole = float(self.m_ghe3) / float(self.nbh3)
        self.log_time = np.linspace(-10, 4, 25).tolist()

        # Borehole Heat Exchanger
        self.bhe_type = bhe_type
        self.bhe1 = get_bhe_object(bhe_type, self.m_ghe1_borehole, fluid, borehole, pipe, grout, soil)
        self.bhe2 = get_bhe_object(bhe_type, self.m_ghe2_borehole, fluid, borehole, pipe, grout, soil)
        self.bhe3 = get_bhe_object(bhe_type, self.m_ghe3_borehole, fluid, borehole, pipe, grout, soil)

        # Equivalent borehole Heat Exchanger

        self.bhe1_eq = self.bhe1.to_single()
        self.bhe2_eq = self.bhe2.to_single()
        self.bhe3_eq = self.bhe3.to_single()

        # Radial numerical short time step
        self.bhe1_eq.calc_sts_g_functions()
        self.bhe2_eq.calc_sts_g_functions()
        self.bhe3_eq.calc_sts_g_functions()

        # Calculations of g-functions

        self.gFunction1, self.gFunction2, self.gFunction3 = generate_g_functions_for_all_ghe(
            self.bhe1, self.bhe2, self.bhe3,
            self.m_ghe1_borehole, self.m_ghe2_borehole, self.m_ghe3_borehole,
            self.bhe_type, self.log_time, self.B_spacing,
            self.n_rows1, self.n_cols1, self.n_rows2, self.n_cols2, self.n_rows3, self.n_cols3,
            calc_g_func_for_multiple_lengths)

        b_over_h1 = self.B_spacing / self.bhe1.b.H
        b_over_h2 = self.B_spacing / self.bhe2.b.H
        b_over_h3 = self.B_spacing / self.bhe3.b.H

        self.g1, _ = self.grab_g_function(self.gFunction1, b_over_h1, self.bhe1_eq)
        self.g2, _ = self.grab_g_function(self.gFunction2, b_over_h2, self.bhe2_eq)
        self.g3, _ = self.grab_g_function(self.gFunction3, b_over_h3, self.bhe3_eq)

        # Calculation of c_n constants

        ts = self.bhe1_eq.t_s  # (-)
        print(ts)
        two_pi_k = 2*np.pi * self.bhe1.soil.k
        rb1 = self.bhe1.calc_effective_borehole_resistance()  # (m.K/W)
        rb2 = self.bhe2.calc_effective_borehole_resistance()  # (m.K/W)
        rb3 = self.bhe3.calc_effective_borehole_resistance()  # (m.K/W)
        k_soil = self.bhe1.soil.k

        self.c_n1, self.c_n2, self.c_n3 = calculation_of_ghe_constant_c_n(
            self.g1, self.g2, self.g3,
            ts, k_soil,
            rb1, rb2, rb3,
            time_array, n)

    @staticmethod
    def combine_sts_lts(log_time_lts: list, g_lts: list, log_time_sts: list, g_sts: list) -> interp1d:
        # make sure the short time step doesn't overlap with the long time step
        max_log_time_sts = max(log_time_sts)
        min_log_time_lts = min(log_time_lts)

        if max_log_time_sts < min_log_time_lts:
            log_time = log_time_sts + log_time_lts
            g = g_sts + g_lts
        else:
            # find where to stop in sts
            i = 0
            value = log_time_sts[i]
            while value <= min_log_time_lts:
                i += 1
                value = log_time_sts[i]
            log_time = log_time_sts[0:i] + log_time_lts
            g = g_sts[0:i] + g_lts
        g = interp1d(log_time, g)

        return g

    # I have made changes in grab_g_function so that I can generate 3 different gfunctions for my 3 different GHEs
    def grab_g_function(self, g_function_obj, b_over_h, bhe_eq):
        # Interpolate for the Long Time Step (LTS) g-function
        g_function, rb_value, _, _ = g_function_obj.g_function_interpolation(b_over_h)

        # Correct the g-function for borehole radius
        g_function_corrected = g_function_obj.borehole_radius_correction(
            g_function, rb_value, bhe_eq.b.r_b
        )

        # Combine short and long time step g-functions
        g = self.combine_sts_lts(
            self.log_time,
            g_function_corrected,
            bhe_eq.lntts.tolist(),
            bhe_eq.g.tolist(),
        )

        g_bhw = self.combine_sts_lts(
            self.log_time,
            g_function_corrected,
            bhe_eq.lntts.tolist(),
            bhe_eq.g_bhw.tolist(),
        )

        return g, g_bhw

    def _simulate_detailed(self):
        # My part of code for multiple GHE systems begins here
        if selected_configuration == "3ghe-6hp_1-pipe":

            ts = self.bhe1_eq.t_s  # (-)
            tg = self.bhe1.soil.ugt  # (Celsius)
            two_pi_k = 2*np.pi * self.bhe1.soil.k  # (W/m.K)

            # Input of HP coefficients
            hp_coeffs = load_all_hp_coefficients()

            c1_hp1_htg, c1_hp2_htg, c1_hp3_htg, c1_hp4_htg, c1_hp5_htg, c1_hp6_htg = hp_coeffs["c1_htg"]
            c2_hp1_htg, c2_hp2_htg, c2_hp3_htg, c2_hp4_htg, c2_hp5_htg, c2_hp6_htg = hp_coeffs["c2_htg"]
            c3_hp1_htg, c3_hp2_htg, c3_hp3_htg, c3_hp4_htg, c3_hp5_htg, c3_hp6_htg = hp_coeffs["c3_htg"]

            c1_hp1_clg, c1_hp2_clg, c1_hp3_clg, c1_hp4_clg, c1_hp5_clg, c1_hp6_clg = hp_coeffs["c1_clg"]
            c2_hp1_clg, c2_hp2_clg, c2_hp3_clg, c2_hp4_clg, c2_hp5_clg, c2_hp6_clg = hp_coeffs["c2_clg"]
            c3_hp1_clg, c3_hp2_clg, c3_hp3_clg, c3_hp4_clg, c3_hp5_clg, c3_hp6_clg = hp_coeffs["c3_clg"]

            a_htg_hp1, a_htg_hp2, a_htg_hp3, a_htg_hp4, a_htg_hp5, a_htg_hp6 = hp_coeffs["a_htg"]
            b_htg_hp1, b_htg_hp2, b_htg_hp3, b_htg_hp4, b_htg_hp5, b_htg_hp6 = hp_coeffs["b_htg"]
            c_htg_hp1, c_htg_hp2, c_htg_hp3, c_htg_hp4, c_htg_hp5, c_htg_hp6 = hp_coeffs["c_htg"]

            a_clg_hp1, a_clg_hp2, a_clg_hp3, a_clg_hp4, a_clg_hp5, a_clg_hp6 = hp_coeffs["a_clg"]
            b_clg_hp1, b_clg_hp2, b_clg_hp3, b_clg_hp4, b_clg_hp5, b_clg_hp6 = hp_coeffs["b_clg"]
            c_clg_hp1, c_clg_hp2, c_clg_hp3, c_clg_hp4, c_clg_hp5, c_clg_hp6 = hp_coeffs["c_clg"]

            nbh_ghe1, nbh_ghe2, nbh_ghe3 = self.nbh1, self.nbh2, self.nbh3
            nbh_total = nbh_ghe1 + nbh_ghe2 + nbh_ghe3
            cp = self.bhe1.fluid.cp  # (J/kg.s)
            rho_fluid = self.bhe1.fluid.rho

            g1, g2, g3 = self.g1, self.g2, self.g3
            c_n_ghe1, c_n_ghe2, c_n_ghe3 = self.c_n1, self.c_n2, self.c_n3

            # Initializing the values
            q_ghe1, q_ghe2, q_ghe3 = np.zeros(n), np.zeros(n), np.zeros(n)
            H_n_ghe1, H_n_ghe2, H_n_ghe3 = np.zeros(n), np.zeros(n), np.zeros(n)
            total_values_ghe1, total_values_ghe2, total_values_ghe3 = np.zeros(n), np.zeros(n), np.zeros(n)

            # Initializing mass flow rates
            (m_hp1_array, m_hp2_array, m_hp3_array, m_hp4_array, m_hp5_array, m_hp6_array, m_loop_array, m_ghe1_array,
             m_ghe2_array, m_ghe3_array) = (np.zeros(n) for _ in range(10))

            # Initializing temperatures
            t = np.full((n, 15), tg)
            t[0, :] = tg  # T1 through T_ghe2_ext

            # Extracting net heating loads
            q_net_htg_hp1, q_net_htg_hp2, q_net_htg_hp3, q_net_htg_hp4, q_net_htg_hp5, q_net_htg_hp6 = (
                get_net_htg_loads(df1))

            results = []

            # Special case for the first iteration
            total_values_ghe1[0], total_values_ghe2[0], total_values_ghe3[0] = (0, 0, 0)
            H_n_ghe1[0], H_n_ghe2[0], H_n_ghe3[0] = (tg, tg, tg)
            q_ghe2[0], q_ghe2[0], q_ghe3[0] = (0, 0, 0)

            labels = [
                "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "Tf_ghe1", "Tf_ghe2", "Tf_ghe3",
                "T_ghe1_ext", "T_ghe2_ext", "T_ghe3_ext", "q_ghe1", "q_ghe2", "q_ghe3"
            ]
            units = ["°C"] * 16 + ["W/m", "W/m", "W/m"]

            for i in range(1, n):
                if (np.isclose(q_net_htg_hp1[i], 0) and np.isclose(q_net_htg_hp2[i], 0)
                        and np.isclose(q_net_htg_hp3[i], 0) and np.isclose(q_net_htg_hp4[i], 0)
                        and np.isclose(q_net_htg_hp5[i], 0) and np.isclose(q_net_htg_hp6[i], 0)):
                    t[i, :] = t[i - 1, :]
                    q_ghe1[i], q_ghe2[i], q_ghe3[i] = (0, 0, 0)

                # Import hourly heating and cooling loads
                h, c = get_hourly_h_c(i, df1)

                t_eft_hp1, t_eft_hp2, t_eft_hp3, t_eft_hp4, t_eft_hp5, t_eft_hp6 = t[i - 1, 1:7]

                # Calculating values for heat pump constants r1 and r2 for all heat pumps

                t_eft = np.array([t_eft_hp1, t_eft_hp2, t_eft_hp3, t_eft_hp4, t_eft_hp5, t_eft_hp6])

                # Calculating r1 and r2

                a_htg = np.array([a_htg_hp1, a_htg_hp2, a_htg_hp3, a_htg_hp4, a_htg_hp5, a_htg_hp6])
                b_htg = np.array([b_htg_hp1, b_htg_hp2, b_htg_hp3, b_htg_hp4, b_htg_hp5, b_htg_hp6])
                c_htg = np.array([c_htg_hp1, c_htg_hp2, c_htg_hp3, c_htg_hp4, c_htg_hp5, c_htg_hp6])

                a_clg = np.array([a_clg_hp1, a_clg_hp2, a_clg_hp3, a_clg_hp4, a_clg_hp5, a_clg_hp6])
                b_clg = np.array([b_clg_hp1, b_clg_hp2, b_clg_hp3, b_clg_hp4, b_clg_hp5, b_clg_hp6])
                c_clg = np.array([c_clg_hp1, c_clg_hp2, c_clg_hp3, c_clg_hp4, c_clg_hp5, c_clg_hp6])

                r1, r2 = calculate_r1_r2(t_eft, h, c, a_htg, b_htg, c_htg, a_clg, b_clg, c_clg)

                # Calculating heat pump capacity

                # Coefficients
                c1_htg = np.array([c1_hp1_htg, c1_hp2_htg, c1_hp3_htg, c1_hp4_htg, c1_hp5_htg, c1_hp6_htg])
                c2_htg = np.array([c2_hp1_htg, c2_hp2_htg, c2_hp3_htg, c2_hp4_htg, c2_hp5_htg, c2_hp6_htg])
                c3_htg = np.array([c3_hp1_htg, c3_hp2_htg, c3_hp3_htg, c3_hp4_htg, c3_hp5_htg, c3_hp6_htg])

                c1_clg = np.array([c1_hp1_clg, c1_hp2_clg, c1_hp3_clg, c1_hp4_clg, c1_hp5_clg, c1_hp6_clg])
                c2_clg = np.array([c2_hp1_clg, c2_hp2_clg, c2_hp3_clg, c2_hp4_clg, c2_hp5_clg, c2_hp6_clg])
                c3_clg = np.array([c3_hp1_clg, c3_hp2_clg, c3_hp3_clg, c3_hp4_clg, c3_hp5_clg, c3_hp6_clg])

                m_design = np.array(
                    [m_hp1_design, m_hp2_design, m_hp3_design, m_hp4_design, m_hp5_design, m_hp6_design])
                m_single = np.array(
                    [m_hp1_single_as, m_hp_single_wtw, m_hp2_single_as, m_hp_single_wtw, m_hp3_single_as,
                     m_hp_single_wtw])

                q_net_htg = np.array(
                    [q_net_htg_hp1[i], q_net_htg_hp2[i], q_net_htg_hp3[i], q_net_htg_hp4[i], q_net_htg_hp5[i],
                     q_net_htg_hp6[i]])

                # Get capacity
                hp_capacity = calculate_hp_capacity(t_eft, c1_htg, c2_htg, c3_htg, c1_clg, c2_clg, c3_clg, m_design,
                                                    m_single, q_net_htg)

                # Unpack
                # hp1_capacity, hp2_capacity, hp3_capacity, hp4_capacity, hp5_capacity, hp6_capacity = hp_capacity

                # I am introducing the scaling factor here (m_hp3_design/m_hp_single) because loads taken are for
                # whole building but heat pump used is only one. Obviously the heat pump is not able to handle all
                # the building loads. To tackle this I scaled up mass flow rate by ratio of htg/clg_loads to hp
                # capacity. To be consistent in this assumption, I need to scale up hp capacity also by the same factor
                # or else  I will always get rtf greater than 1.

                rtf, m_dot, m_loop, m_ghe1, m_ghe2, m_ghe3 = calculate_rtf_and_mass_flows(
                    q_net_htg, t_eft,
                    c1_htg, c2_htg, c3_htg,
                    c1_clg, c2_clg, c3_clg,
                    m_design, m_single,
                    m_hp1_design, m_hp2_design, m_hp3_design,
                    m_hp4_design, m_hp5_design, m_hp6_design,
                    nbh_ghe1, nbh_ghe2, nbh_ghe3, nbh_total
                )

                m_hp1, m_hp2, m_hp3, m_hp4, m_hp5, m_hp6 = m_dot

                # Storing these values
                values = [m_hp1, m_hp2, m_hp3, m_hp4, m_hp5, m_hp6, m_loop, m_ghe1, m_ghe2, m_ghe3]
                arrays = [m_hp1_array, m_hp2_array, m_hp3_array, m_hp4_array, m_hp5_array, m_hp6_array,
                          m_loop_array, m_ghe1_array, m_ghe2_array, m_ghe3_array]

                for arr, val in zip(arrays, values):
                    arr[i] = val
                time_n = time_array[i]

                compute_ghe_history_terms(
                    i, time_array, ts, two_pi_k,
                    q_ghe1, q_ghe2, q_ghe3,
                    g1, g2, g3,
                    tg,
                    total_values_ghe1, total_values_ghe2, total_values_ghe3,
                    H_n_ghe1, H_n_ghe2, H_n_ghe3
                )

                # from multiple_ghe_hp_addition import solve_temperature_matrix
                X = solve_temperature_matrix(
                    i, r1, r2, m_loop, cp, m_ghe1, m_ghe2, m_ghe3, c_n_ghe1, c_n_ghe2, c_n_ghe3, H_n_ghe1, H_n_ghe2,
                    H_n_ghe3, nbh_ghe1, nbh_ghe2, nbh_ghe3, self.bhe1, self.bhe2, self.bhe3, time_array, labels, units,
                    t, q_ghe1, q_ghe2, q_ghe3, results)

            # Create DataFrame and insert units as second row
            df = pd.DataFrame(results)
            unit_row = {"Time Step": "Units"}
            unit_row.update(dict(zip(labels, units)))

            # Insert units row at the top
            df_with_units = pd.concat([pd.DataFrame([unit_row]), df], ignore_index=True)

            # Save to CSV
            df_with_units.to_csv("detailed_simulation_results.csv", index=False)
        else:
            print("Configuration not recognized")
