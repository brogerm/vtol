import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')  # add parent directory
import vtolParam as P
from lon_controller import LonController
from lat_controller import LatController
from signalGenerator import signalGenerator
from vtolAnimation import vtolAnimation
from plotData import plotData
from vtolDynamics import vtolDynamics
from plotObserverData import plotObserverData

# instantiate reference input classes
vtol = vtolDynamics()
lon_ctrl = LonController()
lat_ctrl = LatController()
z_reference = signalGenerator(amplitude=2, frequency=0.15)
h_reference = signalGenerator(amplitude=1, frequency=0.08)

# set disturbance input
force_disturbance = 1.0
torque_disturbance = 0.1

# instantiate the simulation plots and animation
dataPlot = plotData()
observerPlot = plotObserverData()
animation = vtolAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    z_r = z_reference.square(t)
    h_r = h_reference.square(t)

    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        u_lon = lon_ctrl.u(h_r, [vtol.states()[1], vtol.states()[4]])
        u_lat = lat_ctrl.u(z_r, [vtol.states()[0], vtol.states()[2], vtol.states()[3], vtol.states()[5]])
        sys_input = [u_lon[0] + force_disturbance, u_lat[0] + torque_disturbance]  # input to plant is control input + disturbance (formatted as a list)
        vtol.propagateDynamics(sys_input)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation
    animation.drawVtol(vtol.states())
    dataPlot.updatePlots(t, z_r, vtol.states(), [u_lon, u_lat], h_r)
    # construct combined x_hat
    x_hat = [lat_ctrl.x_hat[0], lon_ctrl.x_hat[0], lat_ctrl.x_hat[1], lat_ctrl.x_hat[2], lat_ctrl.x_hat[1], lat_ctrl.x_hat[3]]
    observerPlot.updatePlots(t, vtol.states(), x_hat)
    plt.pause(0.1)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
