import pandas as pd
import numpy as np
from ghedesigner.media import Grout, Soil, GHEFluid
from ghedesigner.media import Pipe as MediaPipe   # I am importing Pipe from media as MediaPipe to avoid name conflict with my Pipe class
from pygfunction.boreholes import Borehole
from ghedesigner.ghe.coaxial_borehole import get_bhe_object
from ghedesigner.ghe.gfunction import calc_g_func_for_multiple_lengths
from ghedesigner.ghe.simulation import SimulationParameters
from ghedesigner.enums import BHPipeType, TimestepType
from ghedesigner.ghe.gfunction import GFunction, calc_g_func_for_multiple_lengths
from ghedesigner.ghe.ground_heat_exchangers import BaseGHE

from OpenGL.GL import *
from OpenGL_2D_class_GLFW import gl2D, gl2DCircle, gl2DText,gl2DArrow, gl2DArc
from HersheyFont import HersheyFont
hf = HersheyFont()

import json


class GHX:
    def __init__(self):
        self.ID = None
        self.type = "GHX"
        self.nodeID = None
        self.input = None
        self.n_rows = None
        self.n_cols = None
        self.row_spacing = None
        self.col_spacing = None
        self.nbh = None
        self.height = None
        self.upstream_device = None
        self.downstream_device = None
        self.matrix_line = None
        self.height = None
        self.row_index = None

        # Parameters to be assigned later
        self.m_dot_total = None

        # Thermal object references (to be set during setup)
        self.pipe = Pipe
        self.soil = Soil
        self.grout = Grout
        self.borehole = Borehole
        self.fluid = None
        self.bhe_type = BHPipeType.SINGLEUTUBE
        self.split_ratio = None

        # Computed properties
        self.bhe = None
        self.r_b = None # borehole thermal resistance
        self.gFunction = None
        self.mass_flow_ghe = None
        self.mass_flow_ghe_borehole = None
        self.depth = None

        self.mass_flow_ghe_design = None
        self.mass_flow_ghe_borehole_design = None
        self.H_n_ghe = None
        self.total_values_ghe = None

        # for output
        self.t_eft = None
        self.t_mean = None
        self.q_ghe = None
        self.t_exit = None

    def generate_g_function_object(self, log_time, calc_g_func_for_multiple_lengths, h_values):
        self.r_b = self.bhe.calc_effective_borehole_resistance()
        self.depth = self.bhe.b.D
        self.mass_flow_ghe_borehole_design = self.mass_flow_ghe_design/self.nbh
        h_values = [self.height]
        coordinates_ghe = [(i * self.row_spacing, j * self.row_spacing) for i in range(int(self.n_rows)) for j in range(int(self.n_cols))]
        self.gFunction = calc_g_func_for_multiple_lengths(
            self.row_spacing, h_values, self.r_b, self.depth, self.mass_flow_ghe_borehole_design, self.bhe_type, log_time,
            coordinates_ghe, self.bhe.fluid, self.bhe.pipe, self.bhe.grout, self.bhe.soil
        )
        return self.gFunction

    def grab_g_function(self, log_time):
        """
        Interpolates g-function values using self.gFunction and self.bhe,
        and returns g and g_bhw arrays.
        """

        # Interpolate LTS g-function
        g_function, rb_value, _, _ = self.gFunction.g_function_interpolation(self.row_spacing / self.height)

        # Correct the g-function for borehole radius
        g_function_corrected = self.gFunction.borehole_radius_correction(
            g_function, rb_value, self.bhe.b.r_b
        )

        # Combine STS and LTS g-functions
        g = BaseGHE.combine_sts_lts(
            log_time,
            g_function_corrected,
            self.bhe.lntts.tolist(),
            self.bhe.g.tolist(),
        )

        g_bhw = BaseGHE.combine_sts_lts(
            log_time,
            g_function_corrected,
            self.bhe.lntts.tolist(),
            self.bhe.g_bhw.tolist(),
        )

        return g, g_bhw

    def calculation_of_ghe_constant_c_n(self, g, ts, time_array, n_timesteps):

        """
        Calculate C_n values for three GHEs based on their g-functions.

        Cn = 1 / (2 * pi * K_s) * g((tn - tn-1) / t_s) + R_b
        """

        two_pi_k = 2 * np.pi * self.soil.k
        c_n = np.zeros(n_timesteps, dtype=float)

        for i in range(1, n_timesteps):
            delta_log_time = np.log((time_array[i] - time_array[i - 1]) / (ts / 3600))
            g_val = g(delta_log_time)

            c_n[i] = (1 / two_pi_k * g_val) + self.r_b

        return c_n

    def compute_history_term(self, i, time_array, ts, two_pi_k, g, tg, H_n_ghe, total_values_ghe, q_ghe):
        """
        Computes the history term H_n for this GHX at time index `i`.
        Updates self.total_values_ghe and self.H_n_ghe in place.
        """
        if i == 0:
            H_n_ghe[i] = tg
            total_values_ghe[i] = 0
            return

        time_n = time_array[i]

        # Compute dimensionless time for all indices from 1 to i-1
        indices = np.arange(1, i)
        dim_less_time = np.log((time_n - time_array[indices - 1]) / (ts / 3600))

        # Compute contributions from all previous steps
        delta_q_ghe = (q_ghe[indices] - q_ghe[indices - 1]) / two_pi_k
        values = np.sum(delta_q_ghe * g(dim_less_time))

        total_values_ghe[i] = values

        # Contribution from the last time step only
        dim1_less_time = np.log((time_n - time_array[i - 1]) / (ts / 3600))
        H_n_ghe[i] = tg - total_values_ghe[i] + (
                q_ghe[i - 1] / two_pi_k * g(dim1_less_time)
        )
        return H_n_ghe[i]

    def generate_GHX_matrix_row(self, matrix_size, m_loop, mass_flow_ghe, cp, H_n_ghe, c_n):
        row1 = np.zeros(matrix_size)
        row2 = np.zeros(matrix_size)
        row3 = np.zeros(matrix_size)
        row4 = np.zeros(matrix_size)

        row_index = self.row_index
        neighbour_index = self.downstream_device.row_index

        row1[row_index] = (m_loop - mass_flow_ghe)*cp
        row1[row_index + 3] = mass_flow_ghe * cp
        row1[neighbour_index] = - m_loop * cp

        row2[row_index + 1] = 1
        row2[row_index + 2] = c_n

        row3[row_index] = -1
        row3[row_index + 1] = 2
        row3[row_index + 3] = -1

        row4[row_index] = mass_flow_ghe * cp
        row4[row_index + 2] = self.height * (self.n_rows * self.n_cols)
        row4[row_index + 3] = - mass_flow_ghe * cp

        rhs1, rhs2, rhs3, rhs4 = 0, H_n_ghe, 0, 0

        rows = [row1, row2, row3, row4]
        rhs = [rhs1, rhs2, rhs3, rhs4]
        return rows, rhs

class Building:
    def __init__(self):
        self.name = None
        self.ID = None
        self.zoneIDs = []  # list of zone ids
        self.zones = []  # list of zones


class Zone:
    def __init__(self):
        # values read from the file
        self.name = None
        self.ID = None
        self.type = "zone"
        self.nodeID = None
        self.node = None
        self.HPmodel = None
        self.HP = None
        self.loads_file = None
        self.matrix_line = None
        self.row_index = None
        self.index = None
        self.mass_flow_zone = None
        self.df_zone = None
        self.t_eft = None
        self.upstream_device = None
        self.downstream_device = None

    def q_net_htg(self):
        """
        Calculate net heat extracted/rejected each hour for the zone.
        If either column is missing, default to zeros.
        """
        if "HPHtgLd_W" in self.df_zone.columns:
            h = np.array(self.df_zone["HPHtgLd_W"])
        else:
            h = np.zeros(len(self.df_zone))

        if "HPClgLd_W" in self.df_zone.columns:
            c = np.array(self.df_zone["HPClgLd_W"])
        else:
            c = np.zeros(len(self.df_zone))

        return h - c

    def zone_mass_flow_rate(self, t_eft, q_net_htg, i):
        hp = self.HP
        cap_htg = hp.c1_htg * t_eft ** 2 + hp.c2_htg * t_eft + hp.c3_htg
        cap_clg = hp.c1_clg * t_eft ** 2 + hp.c2_clg * t_eft + hp.c3_clg
        m_single_hp = hp.m_single_hp

        q_i = q_net_htg[i]
        hp_capacity = cap_htg if q_i > 0 else cap_clg

        # compute mass flow rates
        self.mass_flow_zone = np.abs(q_i) / hp_capacity * m_single_hp

        return self.mass_flow_zone

    def calculate_r1_r2(self, t_eft, hour_index):
        """
        Calculate r1 and r2 for this zone based on entering fluid temperature and HP coefficients.
        """

        # Extract loads
        h = self.df_zone["HPHtgLd_W"].iloc[hour_index] if "HPHtgLd_W" in self.df_zone.columns else 0.0
        c = self.df_zone["HPClgLd_W"].iloc[hour_index] if "HPClgLd_W" in self.df_zone.columns else 0.0

        # Extract HP coefficients
        a_htg = self.HP.a_htg
        b_htg = self.HP.b_htg
        c_htg = self.HP.c_htg

        a_clg = self.HP.a_clg
        b_clg = self.HP.b_clg
        c_clg = self.HP.c_clg

        # Heating calculations
        slope_htg = 2 * a_htg * t_eft + b_htg
        ratio_htg = a_htg * t_eft ** 2 + b_htg * t_eft + c_htg
        u = ratio_htg - slope_htg * t_eft
        v = slope_htg

        # Cooling calculations
        slope_clg = 2 * a_clg * t_eft + b_clg
        ratio_clg = a_clg * t_eft ** 2 + b_clg * t_eft + c_clg
        a = ratio_clg - slope_clg * t_eft
        b = slope_clg

        # Final arrays
        r1 = v * h - b * c
        r2 = u * h - a * c

        return r1, r2

    def generate_zone_matrix_row(self, matrix_size, m_loop, cp, r1, r2):
        neighbour_index = self.downstream_device.row_index
        row = np.zeros(matrix_size)
        row[self.row_index] = 1 - r1/(m_loop * cp)
        row[neighbour_index] = -1
        rhs = r2/(m_loop * cp)
        return row, rhs


class Node:
    def __init__(self):
        self.ID = None
        self.type = None
        self.x = None
        self.y = None
        self.z = None
        self.input = None
        self.output = None
        self.diversion = None

class Pipe:
    def __init__(self):
        self.ID = None
        self.node_in_name = None
        self.node_out_name = None
        self.input = None
        self.output = None
        self.length = None
        self.type = None


class HPmodel:
    def __init__(self):
        self.name = None
        self.ID = None
        self.a_htg, self.b_htg, self.c_htg = None, None, None
        self.a_clg, self.b_clg, self.c_clg = None, None, None
        self.c1_htg, self.c2_htg, self.c3_htg = None, None, None
        self.c1_clg, self.c2_clg, self.c3_clg = None, None, None
        self.m_single_hp = None
        self.design_htg_cap = None
        self.design_clg_cap = None


class IsolationHX:
    def __init__(self):
        self.name = None
        self.type = "ISHX"
        self.ID = None
        self.node_network_inlet_ID = None
        self.node_HP_inlet_ID = None
        self.node_HP_outlet_ID = None
        self.beta_ISHX = None
        self.input = None
        self.HP_output = None
        self.HP_input = None
        self.upstream_device = None
        self.downstream_device = None
        self.upstream_device_HP = None
        self.downstream_device_HP = None

        self.zoneIDs = []  # list of zone ids
        self.zones = []  # list of zones


class GHEHPSystem:
    def __init__(self):
        self.title = None
        self.GHXs = []
        self.buildings = []
        self.zones = []
        self.nodes = []
        self.pipes = []
        self.HPmodels = []
        self.ISHXs = []
        self.current_row = 0
        self.m_loop = None
        self.bhe = None
        self.g_value = {}
        self.c_n = {}
        self.time_array = None
        self.time_array_size = None

        # Thermal object references (to be set during setup)
        self.pipe = None
        self.soil = None
        self.grout = None
        self.borehole = None
        self.fluid = None
        self.mass_flow_ghe_borehole = None
        self.nbh_total = None
        self.gFunction = None
        self.g = None
        self.log_time = None
        self.mass_flow_ghe = None
        self.bhe_eq = None
        self.c_n = None
        self.m_loop = None
        self.beta = None

        self.df = None
        self.current_frame = 0
        self.data = None


    def read_GHEHPSystem_data(self, data):
        next_matrix_line = 0
        for line in data:  # loop over all the lines
            cells = [c.strip() for c in line.strip().split(',')]
            keyword = cells[0].lower()

            if keyword == 'title':
                self.title = cells[1].replace("'", "")

            if keyword == 'ghx':
                thisghx = GHX()
                thisghx.ID = str(cells[1])
                thisghx.nodeID = str(cells[2])
                thisghx.n_rows = float(cells[3])
                thisghx.n_cols = float(cells[4])
                thisghx.row_spacing = float(cells[5])
                thisghx.col_spacing = float(cells[6])
                thisghx.ghe_height = float(cells[7])
                thisghx.mass_flow_ghe_design = float(cells[8])
                thisghx.matrix_line = next_matrix_line
                next_matrix_line += 4
                self.GHXs.append(thisghx)

            if keyword == 'building':
                thisbuilding = Building()
                thisbuilding.name = str(cells[1])
                thisbuilding.ID = str(cells[2])
                thisbuilding.zoneIDs = ([zones.strip() for zones in cells[3:]])
                self.buildings.append(thisbuilding)

            if keyword == 'zone':
                df = pd.read_csv(cells[5])
                self.time_array = df['Hours'].values
                self.time_array_size = len(self.time_array)

                thiszone = Zone()
                thiszone.name = str(cells[1])
                thiszone.ID = str(cells[2])
                thiszone.nodeID = str(cells[3])
                thiszone.HPmodel = str(cells[4])
                thiszone.loads_file = pd.read_csv(cells[5])
                thiszone.matrix_line = next_matrix_line
                next_matrix_line += 1
                self.zones.append(thiszone)

            if keyword == 'ishx':
                thisishx = IsolationHX()
                thisishx.name = str(cells[1])
                thisishx.ID = str(cells[2])
                thisishx.node_network_inlet_ID = str(cells[3])
                thisishx.node_HP_inlet_ID = str(cells[4])
                thisishx.node_HP_outlet_ID = str(cells[5])
                thisishx.beta_ISHX = float(cells[6])
                thisishx.zoneIDs = ([zones.strip() for zones in cells[7:]])
                self.ISHXs.append(thisishx)

            if keyword == 'node':
                thisnode = Node()
                thisnode.ID = str(cells[1])
                thisnode.type = str(cells[2])
                thisnode.x = float(cells[3])
                thisnode.y = float(cells[4])
                thisnode.z = float(cells[5])
                self.nodes.append(thisnode)

            if keyword == 'pipe':
                thispipe = Pipe()
                thispipe.ID = str(cells[1])
                thispipe.type = str(cells[2])
                thispipe.node_in_name = str(cells[3])
                thispipe.node_out_name = str(cells[4])
                thispipe.length = float(cells[5])
                self.pipes.append(thispipe)

            if keyword == 'hpmodel':
                thishpmodel = HPmodel()
                thishpmodel.name = str(cells[1])
                thishpmodel.ID = str(cells[2])
                thishpmodel.a_htg, thishpmodel.b_htg, thishpmodel.c_htg = (float(cells[3]), float(cells[4]),
                                                                           float(cells[5]))
                thishpmodel.a_clg, thishpmodel.b_clg, thishpmodel.c_clg = (float(cells[6]), float(cells[7]),
                                                                           float(cells[8]))
                thishpmodel.c1_htg, thishpmodel.c2_htg, thishpmodel.c3_htg = (float(cells[9]), float(cells[10]),
                                                                              float(cells[11]))
                thishpmodel.c1_clg, thishpmodel.c2_clg, thishpmodel.c3_clg = (float(cells[12]), float(cells[13]),
                                                                              float(cells[14]))
                thishpmodel.m_single_hp = float(cells[15])
                self.HPmodels.append(thishpmodel)

            if keyword == "beta":
                self.beta = float(cells[1])

        # end for line
        self.UpdateConnections()

    def read_data_from_json_file(self):
        with open("find_design_bi_rectangle_single_u_tube.json", 'r') as f:
            self.data = json.load(f)

        # Extract input values
        fluid_data = self.data["fluid"]
        soil_data = self.data["ground-heat-exchanger"]["ghe1"]["soil"]
        grout_data = self.data["ground-heat-exchanger"]["ghe1"]["grout"]
        pipe_data = self.data["ground-heat-exchanger"]["ghe1"]["pipe"]
        borehole_data = self.data["ground-heat-exchanger"]["ghe1"]["borehole"]
        geometric_data = self.data["ground-heat-exchanger"]["ghe1"]["geometric_constraints"]

        # Construct objects
        fluid = (
            GHEFluid(
                fluid_data["fluid_name"],
                fluid_data["concentration_percent"],
                fluid_data["temperature"]
            ))
        # Pipe object (Single U-tube)
        r_in = pipe_data["inner_diameter"] / 2.0
        r_out = pipe_data["outer_diameter"] / 2.0
        s = pipe_data["shank_spacing"]

        pipe_positions = MediaPipe.place_pipes(s, r_out, 1)

        pipe = MediaPipe(
            pipe_positions,
            r_in,
            r_out,
            s,
            pipe_data["roughness"],
            pipe_data["conductivity"],
            pipe_data["rho_cp"]
        )

        soil = Soil(soil_data["conductivity"], soil_data["rho_cp"], soil_data["undisturbed_temp"])
        grout = Grout(grout_data["conductivity"], grout_data["rho_cp"])
        borehole = Borehole(100.0, borehole_data["buried_depth"], borehole_data["diameter"] / 2.0, 0.0, 0.0)

        # Simulation parameters
        sim_params = SimulationParameters(num_months=12)
        sim_params.set_design_heights(geometric_data["max_height"], geometric_data["min_height"])

        return fluid, pipe, grout, soil, borehole, sim_params

    def solveSystem(self, fluid, pipe, grout, soil, borehole, sim_params):
        # precompute all time invariant constants

        time_array = self.time_array
        n_timesteps = self.time_array_size
        matrix_size = 4 * len(self.GHXs) + len(self.zones)
        self.log_time = np.linspace(-10, 4, 25).tolist()
        nbh_total = sum(GHX.n_rows * GHX.n_cols for GHX in self.GHXs)
        self.nbh_total = nbh_total

        for GHX in self.GHXs:
            GHX.fluid = fluid
            GHX.pipe = pipe
            GHX.grout = grout
            GHX.soil = soil
            GHX.borehole = borehole
            GHX.sim_params = sim_params

        # for getting g_functions and bhe object
        for GHX in self.GHXs:
            GHX.borehole = borehole
            GHX.height = GHX.borehole.H
            GHX.nbh = GHX.n_rows * GHX.n_cols
            GHX.mass_flow_ghe_borehole_design = GHX.mass_flow_ghe_design / GHX.nbh
            GHX.bhe = get_bhe_object(GHX.bhe_type, GHX.mass_flow_ghe_borehole_design, GHX.fluid, GHX.borehole,
                                     GHX.pipe, GHX.grout, GHX.soil)
            GHX.bhe_eq = GHX.bhe.to_single()
            GHX.bhe_eq.calc_sts_g_functions()
            ts = GHX.bhe_eq.t_s
            cp = GHX.bhe.fluid.cp
            tg = GHX.bhe.soil.ugt
            borehole.H = GHX.height
            h_values = [borehole.H]
            self.gFunction = GHX.generate_g_function_object(self.log_time, calc_g_func_for_multiple_lengths, h_values)
            self.g, _ = GHX.grab_g_function(self.log_time)
            self.c_n = GHX.calculation_of_ghe_constant_c_n(self.g, ts, time_array, n_timesteps)

            # Initializing the values
            for GHX in self.GHXs:
                GHX.H_n_ghe, GHX.total_values_ghe, GHX.q_ghe = np.full((n_timesteps), tg), np.zeros(
                    n_timesteps), np.zeros(n_timesteps)

            # Assigning indices to zones
            for idx, zone in enumerate(self.zones):
                zone.index = idx

            # Initializing t_eft, t__mena, q_ghe, t_exit
            for zone in self.zones:
                zone.t_eft = np.full(n_timesteps, tg)

            for GHX in self.GHXs:
                GHX.t_eft = np.full(n_timesteps, tg)
                GHX.t_mean = np.full(n_timesteps, tg)
                GHX.q_ghe = np.zeros(n_timesteps)
                GHX.t_exit = np.full(n_timesteps, tg)

            # Assigning row_indices
            for k, zone in enumerate(self.zones):
                zone.row_index = k
            for k, GHX in enumerate(self.GHXs):
                GHX.row_index = len(self.zones) + k * 4

        for i in range(1, n_timesteps):  # loop over all timestep
            matrix_rows = []
            matrix_rhs = []
            total_hp_flow = 0
            for zone in self.zones:
                t_eft = zone.t_eft[i - 1]
                zone.df_zone = zone.loads_file
                q_net_htg = zone.q_net_htg()
                m_zone = zone.zone_mass_flow_rate(t_eft, q_net_htg, i)
                total_hp_flow += m_zone

            m_loop = total_hp_flow * self.beta

            for zone in self.zones:
                t_eft = zone.t_eft[i - 1]
                r1, r2 = zone.calculate_r1_r2(t_eft, i)
                this_zone_row, rhs = zone.generate_zone_matrix_row(matrix_size, m_loop, cp, r1, r2)
                matrix_rows.append(this_zone_row)
                matrix_rhs.append(rhs)

            for j, GHX in enumerate(self.GHXs):
                q_ghe = GHX.q_ghe[:i]  # <--- FIXED: slice of all past values, it is an array
                two_pi_k = 2 * np.pi * GHX.soil.k
                nbh = GHX.n_rows * GHX.n_cols
                split_ratio = nbh / nbh_total
                mass_flow_ghe = m_loop * split_ratio
                g = self.g
                c_n = self.c_n[i]
                H_n_ghe = GHX.compute_history_term(i, time_array, ts, two_pi_k, g, tg, GHX.H_n_ghe,
                                                   GHX.total_values_ghe, q_ghe)
                rows, rhs_values = GHX.generate_GHX_matrix_row(matrix_size, m_loop, mass_flow_ghe, cp, H_n_ghe, c_n)
                for row, rhs in zip(rows, rhs_values):
                    matrix_rows.append(row)
                    matrix_rhs.append(rhs)

            # Solve the matrix
            A = np.array(matrix_rows, dtype=float)
            B = np.array(matrix_rhs, dtype=float)

            X = np.linalg.solve(A, B)

            for j, zone in enumerate(self.zones):
                zone.t_eft[i] = X[j]

            X_ghe = X[len(self.zones):]

            for j, GHX in enumerate(self.GHXs):
                base = 4 * j
                GHX.t_eft[i] = X_ghe[base]
                GHX.t_mean[i] = X_ghe[base + 1]
                GHX.q_ghe[i] = X_ghe[base + 2]
                GHX.t_exit[i] = X_ghe[base + 3]

    def createOutput(self):
        # create csv files
        n_timesteps = self.time_array_size
        data_rows = []

        for i in range(n_timesteps):
            row = []
            for zone in self.zones:
                row.append(zone.t_eft[i])

            for GHX in self.GHXs:
                row.append(GHX.t_eft[i])
                row.append(GHX.t_mean[i])
                row.append(GHX.q_ghe[i])
                row.append(GHX.t_exit[i])

            data_rows.append(row)

        # Step 2: Create column labels
        column_names = []

        for j, zone in enumerate(self.zones):
            column_names.append(f"Zone{j}_t_eft")

        for j, GHX in enumerate(self.GHXs):
            column_names += [
                f"GHX{j}_t_eft",
                f"GHX{j}_t_mean",
                f"GHX{j}_q_ghe",
                f"GHX{j}_t_exit"
            ]

        # Step 3: Create and save DataFrame
        self.df = pd.DataFrame(data_rows, columns=column_names)
        self.df.index.name = "Hour"

        # Drop timestep 0 and reindex starting from 1
        self.df = self.df.iloc[1:]
        self.df.index = range(1, len(self.df) + 1)

        # Save to CSV
        self.df.to_csv("output_results.csv")

    def UpdateConnections(self):

        for pipe in self.pipes:
            pipe.input = FindItemByID(pipe.node_in_name, self.nodes)
            pipe.output = FindItemByID(pipe.node_out_name, self.nodes)
            if pipe.type == "1way":
                pipe.input.output = pipe
                pipe.output.input = pipe
            else:
                pipe.input.diversion = pipe
                pipe.output.input = pipe

        for zone in self.zones:
            zone.HP = FindItemByID(zone.HPmodel, self.HPmodels)
            zone.input = FindItemByID(zone.nodeID, self.nodes)
            zone.input.output = zone

        for building in self.buildings:
            for zoneID in building.zoneIDs:
                zone = FindItemByID(zoneID, self.zones)
                building.zones.append(zone)

        for GHX in self.GHXs:
            GHX.input = FindItemByID(GHX.nodeID, self.nodes)
            GHX.input.output = GHX

        for ISHX in self.ISHXs:
            ISHX.input = FindItemByID(ISHX.node_network_inlet_ID, self.nodes)
            ISHX.HP_input = FindItemByID(ISHX.node_HP_inlet_ID, self.nodes)
            ISHX.HP_output = FindItemByID(ISHX.node_HP_outlet_ID, self.nodes)

            ISHX.input.output = ISHX
            ISHX.HP_output.input = ISHX
            ISHX.HP_input.output = ISHX

        for ISHX in self.ISHXs:
            for zoneID in ISHX.zoneIDs:
                zone = FindItemByID(zoneID, self.zones)
                ISHX.zones.append(zone)

        # finding upstream and downstream device for GHX

        for GHX in self.GHXs:

            # find the first upstream mixing node
            device = GHX.input
            while device.type != "mixing":
                device = device.input

            # find the second upstream mixing node
            device = device.input
            while device.type != "mixing":
                device = device.input

            # find the upstream device
            device = device.diversion
            while device.type != "GHX" and device.type != "zone":
                device = device.output

            GHX.upstream_device = device
            device.downstream_device = GHX

        # finding upstream and donwstream device for zones

        for zone in self.zones:
            # find the first upstream mixing node
            device = zone.input
            while device.type != "mixing":
                device = device.input

            # find the second upstream mixing node or if upstream device if it is connected to ISHX
            device = device.input
            while device.type != "mixing" and device.type != "device":
                device = device.input

            # find the upstream device
            if device.type == "mixing":
                device = device.diversion
                while device.type != "GHX" and device.type != "zone" and device.type != "ISHX":
                    device = device.output

            else:
                device = device.input

            zone.upstream_device = device
            device.downstream_device = zone

        # finding upstream and donwstream device for ISHX
        for ISHX in self.ISHXs:
            # find first upstream node
            device = ISHX.input
            while device.type != "mixing":
                device = device.input

            # find the second upstream mixing node
            device = device.input
            while device.type != "mixing":
                device = device.input

            # finding upstream device
            device = device.diversion
            while device.type != "GHX" and device.type != "zone":
                device = device.output

            ISHX.upstream_device = device
            device.downstream_device = ISHX

        #finding ISHX upstream device in HP side

        for ISHX in self.ISHXs:
            device = ISHX.HP_input
            # finding first upstream mixing node
            while device.type != "mixing":
                device = device.input

            # finding upstream device
            device = device.diversion
            while device.type != "zone":
                device = device.output

            ISHX.upstream_device_HP = device
            device.downstream_device = ISHX

        # finding ISHX downstream device in HP side

        for ISHX in self.ISHXs:
            device = ISHX.HP_output

            # finding first mixing node
            while device.type != "mixing":
                device = device.output

            # finding downstream device
            device = device.diversion
            while device.type != "zone":
                device = device.output

            ISHX.downstream_device_HP = device
            device.upstream_device = ISHX


    def drawnetwork(self):
        pipes = self.pipes
        nodes = self.nodes
        zones = self.zones
        GHXs = self.GHXs
        ISHXs = self.ISHXs


        # Drawing zones
        glLineWidth(5)
        glColor3f(0, 0, 0)

        for zone in zones:
            glBegin(GL_LINE_LOOP)  # begin drawing connected lines
            glVertex2f(zone.input.x, zone.input.y + 2)
            glVertex2f(zone.input.x + 10, zone.input.y + 2)
            glVertex2f(zone.input.x + 10, zone.input.y - 2)
            glVertex2f(zone.input.x, zone.input.y - 2)
            glEnd()

        # Drawing GHXs
        glLineWidth(5)
        glColor3f(1, 1, 1)

        for GHX in GHXs:
            glBegin(GL_LINE_LOOP)  # begin drawing connected lines
            glVertex2f(GHX.input.x, GHX.input.y + 2)
            glVertex2f(GHX.input.x - 10, GHX.input.y + 2)
            glVertex2f(GHX.input.x - 10, GHX.input.y - 2)
            glVertex2f(GHX.input.x, GHX.input.y - 2)
            glEnd()

        # Drawing ISHXs
        glLineWidth(5)
        glColor3f(0,0, 0)

        for ISHX in ISHXs:
            glBegin(GL_LINE_LOOP)  # begin drawing connected lines
            glVertex2f(ISHX.input.x, ISHX.input.y - 18)
            glVertex2f(ISHX.input.x, ISHX.input.y + 18)
            glVertex2f(ISHX.HP_output.x, ISHX.HP_output.y + 3)
            glVertex2f(ISHX.HP_input.x, ISHX.HP_input.y - 3)
            glEnd()

        # Drawing pipes
        glColor3f(0, 0, 1)
        glLineWidth(3)

        for pipe in pipes:
            if pipe.type == "1way":
                glColor3f(0,0,1)
            else:
                glColor3f(0,1,0)

            glBegin(GL_LINES)  # begin drawing connected lines
            glVertex2f(pipe.input.x, pipe.input.y)
            glVertex2f(pipe.output.x, pipe.output.y)
            glEnd()

        # Drawing arrows
        glLineWidth(3)
        for pipe in pipes:
            if pipe.type == "1way":
                glColor3f(0, 0, 1)
                xtip, ytip = pipe.output.x, pipe.output.y
                xstart, ystart = pipe.input.x, pipe.input.y
                angle = np.arctan2(ytip - ystart, xtip - xstart) * 180 / np.pi

            else:
                glColor3f(0, 1, 0)
                xtip, ytip = pipe.output.x, pipe.output.y
                xstart, ystart = pipe.input.x, pipe.input.y
                angle = np.arctan2(ytip - ystart, xtip - xstart) * 180 / np.pi
            gl2DArrow(xtip, ytip, size=1, angleDeg=angle, widthDeg=30, toCenter=False, fill=True)

        # Drawing nodes
        glLineWidth(3)
        radius = 0.5
        for node in nodes:
            if node.type == "mixing":
                glColor3f(1, 0, 0)
            elif node.type == "simple":
                glColor3f(0,1, 0)
            else:
                glColor3f(0, 0, 1)

            gl2DCircle(node.x, node.y, radius, fill=True)


def FindItemByID(ID, objectlist):
    # search a list of objects to find one with a particular name
    # of course, the objects must have a "name" member
    for item in objectlist:  # all objects in the list
        if item.ID == ID:  # does it have the ID I am seeking?
            return item  # then return this one
    # next item
    return None  # couldn't find it


System = GHEHPSystem()

def main():
    f1 = open("1-pipe_3ghe-6hp_system_w_pumping_station_input.txt", 'r')
    data = f1.readlines()  # read the entire file as a list of strings
    f1.close()  # close the file  ... very important

    System.read_GHEHPSystem_data(data)

    # Draw
    gl2d = gl2D(None, System.drawnetwork, width=2000, height=1500)
    gl2d.setViewSize(-10, 80, -10, 100, False)
    gl2d.glWait()  # wait for the user to close the window

    print("Finished drawing 1")


if __name__ == "__main__":
    main()






