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
    def __init__(self, xmax, dx, tmax, eta, plot_interval=50,
                 **alpha_params):

        # set aliev-panfilov params
        super().__init__(**alpha_params)

        # set step size in space and time
        self.dx     = dx
        self.eta    = eta
        # for time step, consider cfl condition
        self.dt     = .1 * dx**2 / (2. * eta)

        # time array for measurements
        self.t      = np.arange(0, tmax, self.dt)
        # x array for plotting
        self.x      = np.arange(0, xmax, dx)
        self.xmax   = xmax
        # get number of integration steps from time array
        self.steps  = self.t.size

        # arrays for V, W and Lap(V)
        self.V      = np.zeros_like(self.x)
        self.W      = np.zeros_like(self.x)
        self.LV     = np.zeros_like(self.x)

        # yield data for plotting after every plot_interval steps
        # set plot_interval = 0 to disable plotting
        self.plot_interval = plot_interval

        # variables for measuring conduction velocity (CV), periode (base cycle
        # length BCL), action potential peak width (APD) and rise time of peak
        self.measurements   = None
        self.nom            = 0     # number of measurements
        self.threshold_low  = .05
        self.threshold_high = .95
        self.previous_V     = 0.
        self.uptime         = -10.


    def _log(self, arg):
        if self.plot_interval != 0:
            print(arg)


    def _take_measurement(self, current_V, current_t):
        now = 12.9 * current_t

        if current_V > self.threshold_low \
                and self.previous_V < self.threshold_low:
            # detected upstroke
            CV  = self.xmax / (now - self.uptime)
            BCL = now - self.uptime
            self.uptime = now
            self.measurements['cv'].append(CV)
            self.measurements['bcl'].append(BCL)
            self._log('UPSTROKE:\n\tCV = %e\n\tBCL = %e\n' % (CV, BCL))

        elif current_V > self.threshold_high \
                and self.previous_V < self.threshold_high:
            rise_time = now - self.uptime
            self.measurements['rt'].append(rise_time)
            self._log('rise time = %e\n' % rise_time)

            # measure also maximum potential
            self.measurements['vmax'].append(self.V.max())
            # and increment number of measurements
            self.nom += 1

        elif current_V < self.threshold_low \
                and self.previous_V > self.threshold_low:
            # detected downstroke
            APD = now - self.uptime
            self.measurements['apd'].append(APD)
            self._log('DOWNSTROKE:\n\tAPD = %e\n\n' % APD)

        self.previous_V = current_V


    def integrate(self, max_nom=0x7fff):
        # pulling instance variables to local scope
        V, W, LV = self.V, self.W, self.LV
        t, dx, dt, eta = self.t, self.dx, self.dt, self.eta

        # initial excitation of action potential
        V[0] = 1.
        # initial boundary condition
        boundary = 'neumann'

        # (re-)set measurement dict
        self.measurements = dict(cv=[], bcl=[], apd=[], rt=[], vmax=[])

        for i in range(self.steps-1):
            # compute laplacian and alpha step
            Lap1d(V, LV, dx, boundary)
            dV, dW = self.G(V, W)

            # euler step
            V += dt * dV + dt * eta * LV
            W += dt * dW

            self._take_measurement(V[10], t[i])

            if self.nom >= max_nom:
                print('recorded maximal number of measurements. Terminating...')
                break

            # switching to periodic boundary condition, after enough time for
            # the pulse to develop
            if i == 4_000:
                boundary = 'periodic'
                self._log('switched to peridodic boundary conditions.\n')

            # yield current data for plotting
            if self.plot_interval != 0 and i % self.plot_interval == 0:
                yield V



def collect_measurements(name, values, tmax_list, discard_list, nom):
    # default params for Pulse1d
    params = dict(xmax=100, dx=.2, tmax=tmax_list[0], eta=.3, plot_interval=0,
                  a=.15, k=8., e0=2e-3, m1=.2, m2=.3)

    # record for each value mean and std of each quantity
    result_dict = dict(cv=[[], []], bcl=[[], []], apd=[[], []],
                       rt=[[], []], vmax=[[], []])

    for value, tmax, discard in zip(values, tmax_list, discard_list):
        params[name] = value
        params['tmax'] = tmax
        p = Pulse1d(**params)
        list(p.integrate(nom+discard))
        print('recorded %d measurements after %d steps' \
              % (len(p.measurements['rt']), p.steps))
        for key, val in p.measurements.items():
            result_dict[key][0].append(np.mean(val[discard:]))
            result_dict[key][1].append(np.std(val[discard:]))

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
