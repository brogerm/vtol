import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')  # add parent directory
import vtolParam as P
from vtolController import vtolController
from signalGenerator import signalGenerator
from vtolAnimation import vtolAnimation
from plotData import plotData
from vtolDynamics import vtolDynamics

# instantiate reference input classes
vtol = vtolDynamics()
ctrl = vtolController()
z_reference = signalGenerator(amplitude=3, frequency=0.08)
h_reference = signalGenerator(amplitude=2, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = plotData()
animation = vtolAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    z_r = z_reference.square(t)
    h_r = h_reference.square(t)

    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        u = ctrl.update(z_r[0], h_r[0], vtol.states())
        vtol.propagateDynamics([u[0], u[1]])  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation
    animation.drawVtol(vtol.states())
    dataPlot.updatePlots(t, z_r, vtol.states(), u, h_r)
    plt.pause(0.1)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
