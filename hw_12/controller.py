import numpy as np
import paramHw12 as P

class Controller:
    # state feedback control using dirty derivatives to estimate zdot and thetadot
    def __init__(self):
        self.z_dot = 0.0
        self.h_dot = 0.0
        self.theta_dot = 0.0
        self.z_d1 = 0.
        self.h_d1 = 0.0
        self.theta_d1 = 0.0
        self.K_lon = P.K_lon                 # state feedback gain
        self.ki_lon = P.ki_lon               # Input gain
        self.K_lat = P.K_lat               # state feedback gain
        self.ki_lat = P.ki_lat               # Input gain
        self.Flimit = P.F_max       # Maxiumum Force
        self.Taulimit = P.tau_max       # Maxiumum torque
        self.beta = P.beta           # dirty derivative gain
        self.Ts = P.Ts               # sample rate of controller
        self.h_integrator = 0.0
        self.z_integrator = 0.0  # integrator
        self.h_error_d1 = 0.0  # error signal delayed by 1 sample
        self.z_error_d1 = 0.0  # error signal delayed by 1 sample

    def u(self, y_r, y):
        # y_r is the referenced input
        # y is the current state
        z_r = y_r[0]
        h_r = y_r[1]
        z = y[0]
        h = y[1]
        theta = y[2]

        # differentiate z and theta
        self.differentiateZ(z)
        self.differentiateTheta(theta)
        self.differentiateH(h)

        # integrate longitudinal error
        h_error = h_r - h
        self.h_integrateError(h_error)

        # integrate lateral error
        z_error = z_r - z
        self.z_integrateError(z_error)

        # Construct the longitudinal and lateral states
        xlon = np.matrix([[h], [self.h_dot]])
        xlat = np.matrix([[z], [theta], [self.z_dot], [self.theta_dot]])

        # Compute the state feedback controller
        f_unsat = -self.K_lon * xlon - self.ki_lon * self.h_integrator
        f_sat = self.saturate(f_unsat, self.Flimit)

        if self.ki_lon != 0.0:
            self.h_integrator = self.h_integrator + self.Ts / self.ki_lon * (f_sat - f_unsat)

        # Compute the state feedback controller
        tau_unsat = -self.K_lat * xlat - self.ki_lat * self.z_integrator
        tau_sat = self.saturate(tau_unsat, self.Taulimit)

        if self.ki_lat != 0.0:
            self.z_integrator = self.z_integrator + self.Ts / self.ki_lat * (tau_sat - tau_unsat)

        return [f_sat.item(0), tau_sat.item(0)]

    def differentiateZ(self, z):
        '''
            differentiate z
        '''
        self.z_dot = self.beta*self.z_dot + (1-self.beta)*((z - self.z_d1) / self.Ts)
        self.z_d1 = z

    def differentiateTheta(self, theta):
        '''
            differentiate theta
        '''
        self.theta_dot = self.beta*self.theta_dot + (1-self.beta)*((theta - self.theta_d1) / self.Ts)
        self.theta_d1 = theta

    def differentiateH(self, h):
        '''
            differentiate h
        '''
        self.h_dot = self.beta*self.h_dot + (1-self.beta)*((h - self.h_d1) / self.Ts)
        self.h_d1 = h

    def h_integrateError(self, error):
        self.h_integrator = self.h_integrator + (self.Ts / 2.0) * (error + self.h_error_d1)
        self.h_error_d1 = error

    def z_integrateError(self, error):
        self.z_integrator = self.z_integrator + (self.Ts / 2.0) * (error + self.z_error_d1)
        self.z_error_d1 = error

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u

