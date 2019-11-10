import numpy as np
import lat_param_hw13 as P

class LatController:
    # state feedback control using dirty derivatives to estimate zdot and thetadot
    def __init__(self):
        self.x_hat = np.matrix([
            [0.0],  # initial estimate for z_hat
            [0.0],  # initial estimate for th_hat
            [0.0],  # initial estimate for z_hat_dot
            [0.0]])  # initial estimate for th_hat_dot
        self.tau_d1 = 0.0  # Computed torque, delayed by one sample
        self.z_dot = 0.0
        self.theta_dot = 0.0
        self.z_d1 = 0.
        self.theta_d1 = 0.0
        self.K = P.K               # state feedback gain
        self.ki = P.ki               # Input gain
        self.limit = P.u_max       # Maxiumum torque
        self.L = P.L  # observer gain
        self.A = P.A  # system model
        self.B = P.B
        self.C = P.C
        self.z_integrator = 0.0  # integrator
        self.z_error_d1 = 0.0  # error signal delayed by 1 sample
        self.beta = P.beta  # dirty derivative gain
        self.Ts = P.Ts  # sample rate of controller

    def u(self, y_r, y):
        # y_r is the referenced input
        # y is the current state
        z_r = y_r[0]

        # update the observer and extract z_hat
        self.updateObserver(y)
        z_hat = self.x_hat[0]

        # integrate error
        z_error = z_r - z_hat
        self.z_integrateError(z_error)

        # Compute the lateral state feedback controller
        tau_unsat = -self.K * self.x_hat - self.ki * self.z_integrator
        tau = self.saturate(tau_unsat)
        self.updateTorque(tau)

        return [tau.item(0)]

    def updateObserver(self, y_m):
        N = 10
        y = np.matrix([
            [y_m[0]],
            [y_m[1]]])
        for i in range(0, N):
            self.x_hat = self.x_hat + self.Ts/float(N)*(
                self.A*self.x_hat + self.B*self.tau_d1 + self.L*(y-self.C*self.x_hat)
            )

    def updateTorque(self, tau):
        self.tau_d1 = tau

    def z_integrateError(self, error):
        self.z_integrator = self.z_integrator + (self.Ts / 2.0) * (error + self.z_error_d1)
        self.z_error_d1 = error

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u

