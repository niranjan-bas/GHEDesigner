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


class Isolation_HX:
    def __init__(self):
        self.name = None
        self.ID = None
        self.node_network_inlet_ID = None
        self.node_network_outlet_ID = None
        self.node_HP_inlet_ID = None
        self.node_HP_outlet_ID = None
        self.zoneID = None


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
                thisishx = Isolation_HX()
                thisishx.name = str(cells[1])
                thisishx.ID = str(cells[2])
                thisishx.node_network_inlet_ID = str(cells[3])
                thisishx.node_network_outlet_ID = str(cells[4])
                thisishx.node_HP_inlet_ID = str(cells[5])
                thisishx.node_HP_outlet_ID = str(cells[6])
                thisishx.node_zoneIDs = ([zones.strip() for zones in cells[7:]])

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
                thishpmodel.m_design_htg_cap = float(cells[16])
                thishpmodel.m_design_clg_cap = float(cells[17])
                self.HPmodels.append(thishpmodel)

            if keyword == "beta":
                self.beta = float(cells[1])

        # end for line
        self.UpdateConnections()

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
            ISHX.network_input = FindItemByID(ISHX.node_network_inlet_ID, self.nodes)
            ISHX.network_output = FindItemByID(ISHX.node_network_outlet_ID, self.nodes)
            ISHX.HP_input = FindItemByID(ISHX.node_HP_inlet_ID, self.nodes)
            ISHX.HP_output = FindItemByID(ISHX.node_HP_outlet_ID, self.nodes)

        for GHX in self.GHXs:
            # find the upstream device

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

        #find the upstream device

        for zone in self.zones:
            # find the first upstream mixing node
            device = zone.input
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

            zone.upstream_device = device
            device.downstream_device = zone

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
        radius = 1
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






