from ghedesigner.ghe.simulation import SimulationParameters
from ghedesigner.media import Soil, Grout, Pipe, GHEFluid
from pygfunction.boreholes import Borehole

from District_system_class import GHEHPSystem
from OpenGL_2D_class_GLFW import gl2D

import json

import time

start_time = time.time()

System = GHEHPSystem()


def read_data_from_json_file():
    with open("find_design_bi_rectangle_single_u_tube.json", 'r') as f:
        data = json.load(f)

    # Extract input values
    fluid_data = data["fluid"]
    soil_data = data["ground-heat-exchanger"]["ghe1"]["soil"]
    grout_data = data["ground-heat-exchanger"]["ghe1"]["grout"]
    pipe_data = data["ground-heat-exchanger"]["ghe1"]["pipe"]
    borehole_data = data["ground-heat-exchanger"]["ghe1"]["borehole"]
    geometric_data = data["ground-heat-exchanger"]["ghe1"]["geometric_constraints"]

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

    pipe_positions = Pipe.place_pipes(s, r_out, 1)

    pipe = Pipe(
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


def AnimationCallback(frame, nframes):
    # calculations needed to configure the picture
    # these could be done here or by calling a class method
    System.current_frame = (frame + 1) * 145


def main():
    f1 = open("1-pipe_3ghe-6hp_system_input.txt", 'r')  # open the file for reading     # "1-pipe_3ghe-6hp_system_input.txt"
    data = f1.readlines()  # read the entire file as a list of strings
    f1.close()  # close the file  ... very important

    System.read_GHEHPSystem_data(data)

    fluid, pipe, grout, soil, borehole, sim_params = read_data_from_json_file()
    System.solveSystem(fluid, pipe, grout, soil, borehole, sim_params)
    System.createOutput()
    System.current_frame = 1

    # Draw the house, set the window width and height
    gl2d = gl2D(None, System.drawnetwork, width=2000, height=1500)
    gl2d.setViewSize(-10, 50, -10, 80, False)
    gl2d.glWait()  # wait for the user to close the window

    # Draw the house, set the window width and height
    gl2d = gl2D(None, System.drawnetwork, width=2000, height=1500)
    gl2d.setViewSize(-10, 50, -10, 80, False)
    nframes = 59
    gl2d.glStartAnimation(AnimationCallback, nframes, delaytime=0.1,
                          reverse=False, repeat=False, reset=False)

    gl2d.glWait()  # wait for the user to close the window

    print("Finished drawing 1")


main()

end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")





