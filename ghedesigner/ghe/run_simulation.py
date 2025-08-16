#from District_system_class import GHEHPSystem
from checking_schematic import GHEHPSystem
from OpenGL_2D_class_GLFW import gl2D

import time

start_time = time.time()

System = GHEHPSystem()


def AnimationCallback(frame, nframes):
    # calculations needed to configure the picture
    # these could be done here or by calling a class method
    #System.current_frame = (frame + 1) * 146
    System.current_frame = (frame) * 146


def main():
    f1 = open("1-pipe_3ghe-6hp_system_input.txt", 'r')
    data = f1.readlines()  # read the entire file as a list of strings
    f1.close()  # close the file  ... very important

    System.read_GHEHPSystem_data(data)

    fluid, pipe, grout, soil, borehole, sim_params = System.read_data_from_json_file()
    System.solveSystem(fluid, pipe, grout, soil, borehole, sim_params)
    System.createOutput()
    System.current_frame = 0

    # Draw
    gl2d = gl2D(None, System.drawnetwork, width=2000, height=1500)
    gl2d.setViewSize(-10, 50, -10, 80, False)
    gl2d.glWait()  # wait for the user to close the window

    # Draw
    gl2d = gl2D(None, System.drawnetwork, width=2000, height=1500)
    gl2d.setViewSize(-10, 50, -10, 80, False)
    nframes = 60  #59
    gl2d.glStartAnimation(AnimationCallback, nframes, delaytime=0.1,
                          reverse=False, repeat=False, reset=False)

    gl2d.glWait()  # wait for the user to close the window

    print("Finished drawing 1")


if __name__ == "__main__":
    main()

end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")





