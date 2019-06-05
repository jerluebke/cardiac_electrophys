# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from alpha import AlphaBase


def Lap1d(a, o, dx, boundary='neumann'):
    if boundary == 'neumann':
        o[1:-1] = (a[2:] - 2*a[1:-1] + a[:-2]) / dx**2
        o[0], o[-1] = o[1], o[-2]
    elif boundary == 'periodic':
        o[:] = (np.roll(a, 1) + np.roll(a, -1) - 2*a) / dx**2
    return o



class Pulse1d(AlphaBase):
    def __init__(self, xmax, dx, tmax, eta, plot_interval=1,
                 **alpha_params):
        super().__init__(**alpha_params)
        self.dx     = dx
        self.eta    = eta
        self.dt     = .1 * dx**2 / (2. * eta)
        self.t      = np.arange(0, tmax, self.dt)
        self.x      = np.arange(0, xmax, dx)
        self.steps  = self.t.size
        self.V      = np.zeros_like(self.x)
        self.W      = np.zeros_like(self.x)
        self.LV     = np.zeros_like(self.x)
        self.plot_interval = plot_interval
        self.measurements  = None


    def _log(self, arg):
        if self.plot_interval != 0:
            print(arg)


    def integrate(self):
        # pulling instance variables to local scope
        V, W, LV = self.V, self.W, self.LV
        xmax, t, dx, dt, eta = self.x[-1], self.t, self.dx, self.dt, self.eta

        # initial excitation of action potential
        V[0] = 1.
        # initial boundary condition
        boundary = 'neumann'

        # variables for measuring CV, BCL, APD
        threshold_low   = .05
        threshold_high  = .95
        previous_V      = 0.
        uptime          = -10.

        self.measurements = dict(cv=[], bcl=[], apd=[], rt=[])

        for i in range(self.steps-1):
            # compute laplacian and alpha step
            Lap1d(V, LV, dx, boundary)
            dV, dW = self.G(V, W)

            # euler step
            V += dt * dV + dt * eta * LV
            W += dt * dW

            # current measurement
            now         = 12.9 * t[i]
            current_V   = V[10]

            if current_V > threshold_low and previous_V < threshold_low:
                # detected upstroke
                CV      = xmax / (now - uptime)
                BCL     = now - uptime
                uptime  = now
                self.measurements['cv'].append(CV)
                self.measurements['bcl'].append(BCL)
                self._log('UPSTROKE:\n\tCV = %e\n\tBCL = %e\n' % (CV, BCL))

            elif current_V > threshold_high and previous_V < threshold_high:
                rise_time = now - uptime
                self.measurements['rt'].append(rise_time)
                self._log('rise time = %e\n' % rise_time)

            elif current_V < threshold_low and previous_V > threshold_low:
                # detected downstroke
                APD = now - uptime
                self.measurements['apd'].append(APD)
                self._log('DOWNSTROKE:\n\tAPD = %e\n\n' % APD)

            previous_V = current_V

            # switching boundary condition to periodic, after enough time for
            # the pulse to develop
            if i == 4_000:
                boundary = 'periodic'
                self._log('switched to peridodic boundary conditions.\n')

            if self.plot_interval != 0 and i % self.plot_interval == 0:
                yield V



def collect_measurements(name, values):
    params = dict(xmax=100, dx=.2, tmax=1000., eta=.3, plot_interval=0,
                  a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
    result_dict = dict(cv=[[], []], bcl=[[], []], apd=[[], []], rt=[[], []])
    for value in values:
        params[name] = value
        p = Pulse1d(**params)
        p.integrate()
        for key, val in p.measurements.items():
            result_dict[key][0].append(np.mean(val[:2]))
            result_dict[key][1].append(np.std(val[:2]))

    return result_dict



if __name__ == '__main__':
    params = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
    p = Pulse1d(25, .2, 1000., .3, 50, **params)
    pulse_gen = p.integrate()

    fig, ax = plt.subplots()
    ax.set(title='action potential propagation', xlabel='x', ylabel='V',
           xlim=(0, p.x[-1]), ylim=(-.1, 1.1))
    V_line, = ax.plot([], [])

    def init():
        V_line.set_data([], [])
        return V_line,

    def step(V):
        V_line.set_data(p.x, V)
        return V_line,

    anim = animation.FuncAnimation(fig, step, frames=pulse_gen, interval=20,
                                   init_func=init, blit=True, repeat=False)


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
