import numpy as np 
import random
import vtolParam as P


class vtolDynamics:
    '''
        Model the physical system
    '''

    def __init__(self):
        # Initial state conditions
        self.state = np.matrix([[P.z0],          # z initial horizontal position
                                [P.h0],             # h initial vertical position
                                [P.theta0],      # Theta initial orientation
                                [P.zdot0],       # zdot initial horizontal velocity
                                [P.hdot0],          # hdot initial vertical velocity
                                [P.thetadot0]])  # Thetadot initial velocity
        #################################################
        # The parameters for any physical system are never known exactly.  Feedback
        # systems need to be designed to be robust to this uncertainty.  In the simulation
        # we model uncertainty by changing the physical parameters by a uniform random variable
        # that represents alpha*100 % of the parameter, i.e., alpha = 0.2, means that the parameter
        # may change by up to 20%.  A different parameter value is chosen every time the simulation
        # is run.
        alpha = 0.2  # Uncertainty parameter
        self.ml = P.ml * (1+2*alpha*np.random.rand()-alpha)  # Mass of the left rotor, kg
        self.mc = P.mc * (1+2*alpha*np.random.rand()-alpha)  # Mass of the center pod, kg
        self.mr = P.mr * (1+2*alpha*np.random.rand()-alpha)  # Mass of the right rotor, kg
        self.d = P.d * (1+2*alpha*np.random.rand()-alpha)  # distance between center pod and rotors, m
        self.b = P.mu * (1+2*alpha*np.random.rand()-alpha)  # Damping coefficient, Ns
        self.g = P.g  # the gravity constant is well known and so we don't change it.
        self.Jc = P.Jc * (1+2*alpha*np.random.rand()-alpha)  # given to us, but assume some uncertainty

    def propagateDynamics(self, u):
        '''
            Integrate the differential equations defining dynamics
            P.Ts is the time step between function calls.
            u contains the system input(s).
        '''
        # Integrate ODE using Runge-Kutta RK4 algorithm
        k1 = self.derivatives(self.state, u)
        k2 = self.derivatives(self.state + P.Ts/2*k1, u)
        k3 = self.derivatives(self.state + P.Ts/2*k2, u)
        k4 = self.derivatives(self.state + P.Ts*k3, u)
        self.state += P.Ts/6 * (k1 + 2*k2 + 2*k3 + k4)

    def derivatives(self, state, u):
        '''
            Return xdot = f(x,u), the derivatives of the continuous states, as a matrix
        '''
        # re-label states and inputs for readability
        z = state.item(0)
        h = state.item(1)
        theta = state.item(2)
        zdot = state.item(3)
        hdot = state.item(4)
        thetadot = state.item(5)
        F = u[0]
        tau = u[1]

        fr = 1/2 * F + 1/(2 * self.d) * tau
        fl = 1/2 * F - 1/(2 * self.d) * tau
        # The equations of motion.
        zddot = (-(fl + fr)*np.sin(theta)-self.b * zdot) / (self.mc + 2 * self.ml)
        hddot = 1/(self.mc + 2 * self.ml) * ((fl + fr)*np.cos(theta) - (self.mc + 2 * self.mr) * self.g)
        thetaddot = 1/(self.Jc + 2 * self.ml * self.d**2) * ((fr - fl) * self.d)
        x1dot = zdot
        x2dot = hdot
        x3dot = thetadot
        x4dot = zddot
        x5dot = hddot
        x6dot = thetaddot
        # build xdot and return
        xdot = np.matrix([[x1dot], [x2dot], [x3dot], [x4dot], [x5dot], [x6dot]])
        return xdot

    def outputs(self):
        '''
            Returns the measured outputs as a list
            [z, theta] with added Gaussian noise
        '''
        # re-label states for readability
        z = self.state.item(0)
        h = self.state.item(1)
        theta = self.state.item(2)
        # add Gaussian noise to outputs
        z_m = z + random.gauss(0, 0.01)
        h_m = h + random.gauss(0, 0.01)
        theta_m = theta + random.gauss(0, 0.001)
        # return measured outputs
        return [z_m, h_m, theta_m]

    def states(self):
        '''
            Returns all current states as a list
        '''
        return self.state.T.tolist()[0]