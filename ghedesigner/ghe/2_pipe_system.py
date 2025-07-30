import numpy as np
import pandas as pd
from OpenGL.GL import *
from OpenGL_2D_class_GLFW import gl2D, gl2DCircle, gl2DText,gl2DArrow,gl2DArc
from HersheyFont import HersheyFont
hf = HersheyFont()
import os

# Define filename

file_path = "3ghe-6hp_system.txt"
file_name = os.path.splitext(os.path.basename(file_path))[0]  # Removes .txt

class GHX:
    def __init__(self):
        self.ID = None
        self.type = None
        self.node_in_ID = None
        self.node_out_ID = None
        self.n_rows = None
        self.n_cols = None
        self.row_spacing = None
        self.col_spacing = None
        self.beta = None
        self.height = None
        self.m_flow_ghe_design = None

        self.input = None
        self.output = None

class Zone:
    def __init__(self):
        self.name = None
        self.ID = None
        self.type = None
        self.node_in_ID = None
        self.node_out_ID = None
        self.HPmodel = None
        self.loads_file = None
        self.beta = None

        self.input = None
        self.output = None

class Building:
    def __init__(self):
        self.name = None
        self.ID = None
        self.zoneID = None

class Node:
    def __init__(self):
        self.ID = None
        self.type = None
        self.x = None
        self.y = None
        self.z = None

class Pipe:
    def __init__(self):
        self.ID = None
        self.type = None
        self.node_in_ID = None
        self.node_out_ID = None
        self.length = None

        self.input = None
        self.output = None


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

class GHEHPSystem:
    def __init__(self):
        self.title = None
        self.GHXs = []
        self.buildings = []
        self.zones = []
        self.nodes = []
        self.pipes = []
        self.HPmodels = []

    def read_GHEHPSystem_data(self, data):
        for line in data:  # loop over all the lines
            cells = [c.strip() for c in line.strip().split(',')]
            keyword = cells[0].lower()

            if keyword == 'title':
                self.title = cells[1].replace("'", "")

            if keyword == 'ghx':
                thisghx = GHX()
                thisghx.ID = str(cells[1])
                thisghx.type = str(cells[2])
                thisghx.node_in_ID = str(cells[3])
                thisghx.node_out_ID = str(cells[4])
                thisghx.n_rows = float(cells[5])
                thisghx.n_cols = float(cells[6])
                thisghx.row_spacing = float(cells[7])
                thisghx.col_spacing = float(cells[8])
                thisghx.beta = float(cells[9])
                thisghx.height = float(cells[10])
                thisghx.m_flow_ghe_design = float(cells[11])
                self.GHXs.append(thisghx)

            if keyword == 'building':
                thisbuilding = Building()
                thisbuilding.name = str(cells[1])
                thisbuilding.ID = str(cells[2])
                thisbuilding.zoneIDs = ([zones.strip() for zones in cells[3:]])
                self.buildings.append(thisbuilding)

            if keyword == 'zone':
                thiszone = Zone()
                thiszone.name = str(cells[1])
                thiszone.ID = str(cells[2])
                thiszone.type = str(cells[3])
                thiszone.node_in_ID = str(cells[4])
                thiszone.node_out_ID = str(cells[5])
                thiszone.HPmodel = str(cells[6])
                thiszone.loads_file = pd.read_csv(cells[7])
                thiszone.beta = float(cells[8])
                self.zones.append(thiszone)

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
                thispipe.node_in_ID = str(cells[3])
                thispipe.node_out_ID = str(cells[4])
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

        # end for line
        self.UpdateConnections()

    def UpdateConnections(self):

        for pipe in self.pipes:
            pipe.input = FindItemByID(pipe.node_in_ID, self.nodes)
            pipe.output = FindItemByID(pipe.node_out_ID, self.nodes)

        for GHX in self.GHXs:
            GHX.input = FindItemByID(GHX.node_in_ID, self.nodes)
            GHX.output = FindItemByID(GHX.node_out_ID, self.nodes)

        for zone in self.zones:
            zone.input = FindItemByID(zone.node_in_ID, self.nodes)
            zone.output = FindItemByID(zone.node_out_ID, self.nodes)


    def drawnetwork(self):
        pipes = self.pipes
        nodes = self.nodes
        zones = self.zones
        GHXs = self.GHXs

        # Drawing zones
        glLineWidth(5)
        glColor3f(0, 0, 0)

        for zone in zones:
            glBegin(GL_LINE_LOOP)  # begin drawing connected lines
            glVertex2f(zone.input.x, zone.input.y + 2)
            glVertex2f(zone.input.x+10, zone.input.y + 2)
            glVertex2f(zone.input.x+10, zone.input.y - 2)
            glVertex2f(zone.input.x, zone.input.y - 2)
            glEnd()

        # Drawing GHXs
        glLineWidth(5)
        glColor3f(1, 1, 1)

        for GHX in GHXs:
            glBegin(GL_LINE_LOOP)  # begin drawing connected lines
            glVertex2f(GHX.input.x, GHX.input.y - 2)
            glVertex2f(GHX.input.x - 10, GHX.input.y - 2)
            glVertex2f(GHX.input.x - 10, GHX.input.y + 2)
            glVertex2f(GHX.input.x, GHX.input.y+2)
            glEnd()

        # Drawing pipes
        glColor3f(0, 0, 1)
        glLineWidth(3)
        glBegin(GL_LINES)  # begin drawing connected lines
        for pipe in pipes:
            glVertex2f(pipe.input.x, pipe.input.y)
            glVertex2f(pipe.output.x, pipe.output.y)
        glEnd()

        # Darwing arrows
        glLineWidth(3)
        glColor3f(0, 0, 1)
        for pipe in pipes:
            xtip, ytip = pipe.output.x, pipe.output.y
            xstart, ystart = pipe.input.x, pipe.input.y
            angle = np.arctan2(ytip - ystart, xtip - xstart) * 180 / np.pi

            gl2DArrow(xtip, ytip, size=1, angleDeg=angle, widthDeg=30, toCenter=False, fill=True)

        # Drawing nodes
        #glColor3f(1, 0, 0)
        glLineWidth(3)
        radius = 0.5
        for node in nodes:
            if node.type == "splitting":
                glColor3f(1, 0, 0)
            elif node.type == "combining":
                glColor3f(1, 1, 1)
            elif node.type == "device_in":
                glColor3f(0, 0, 1)
            elif node.type == "device_out":
                glColor3f(1, 1, 0)
            else:
                glColor3f(0, 1, 0)

            gl2DCircle(node.x, node.y, radius, fill=True)

        # Writing text
        glColor3f(1, 1, 1)
        glLineWidth(3)
        hf.drawText(file_name, 35, -5, scale=2.5, slant=0, angle=0, center=True)

def FindItemByID(ID, objectlist):
    # search a list of objects to find one with a particular name
    # of course, the objects must have a "name" member
    for item in objectlist:  # all objects in the list
        if item.ID == ID:  # does it have the ID I am seeking?
            return item  # then return this one
    # next item
    return None  # couldn't find it

def main():
    f1 = open(file_path, 'r')  # open the file for reading
    data = f1.readlines()  # read the entire file as a list of strings
    f1.close()  # close the file  ... very important

    System = GHEHPSystem()
    System.read_GHEHPSystem_data(data)

    # Draw the system network, set the window width and height
    gl2d = gl2D(None, System.drawnetwork, width=2000, height=1500)
    gl2d.setViewSize(-10, 80, -10, 80, False)
    gl2d.glWait()  # wait for the user to close the window

    print("Finished drawing 1")


main()
