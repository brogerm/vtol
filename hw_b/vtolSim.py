import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')  # add parent directory
import vtolParam as P
from signalGenerator import signalGenerator
from vtolAnimation import vtolAnimation
from plotData import plotData
from vtolDynamics import vtolDynamics

# instantiate reference input classes
vtol = vtolDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.1)
zRef = signalGenerator(amplitude=0.5, frequency=0.1)
thetaRef = signalGenerator(amplitude=1/4*np.pi, frequency=0.1)
flRef = signalGenerator(amplitude=7.3575, frequency=.05)
frRef = signalGenerator(amplitude=7.2, frequency=.05)

# instantiate the simulation plots and animation
dataPlot = plotData()
animation = vtolAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    r = reference.square(t)

    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    fr = []
    fl = []
    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        fl = [(P.mc + 2 * P.ml) * P.g / 2]
        fr = fl
        vtol.propagateDynamics([fl[0], fr[0]])  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation
    print(fr[0])
    animation.drawVtol(vtol.states())
    dataPlot.updatePlots(t, r, vtol.states(), fr)
    plt.pause(0.1)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
