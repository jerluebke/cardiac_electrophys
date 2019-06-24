# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from alpha import AlphaBase


def Lap1d(a, o, dx, boundary='neumann'):
    """Compute Laplacian (with lattice constant dx) of 1d array a with desired
    boundary condition, store output in o.

    Available boundary conditions: von Neumann, periodic
    """
    o[:] = (np.roll(a, 1) + np.roll(a, -1) - 2*a) / dx**2
    if boundary == 'neumann':
        o[0], o[-1] = o[1], o[-2]
    elif boundary == 'periodic':
        # do nothing, already imposed by `np.roll`
        pass
    return o



class Pulse1d(AlphaBase):
    """class for integrating 1d tissue with aliev-panfilov model

    In the case of continuous tissue, one adds a diffusion term to the
    action potential equation:
        dV/dt = G(V) + eta * Lap(V)
    which causes the peak to travel spatially (see `self.integrate`)
    """

    def __init__(self, xmax, dx, tmax, eta, plot_interval=50,
                 **alpha_params):
        """majority of params should self-explanatory

        notes
        =====
        plot_interval   :   yield array for plotting every `plot_interval`
                            steps.
                            setting `plot_interval = 0` disables this yielding
                            and printing of intermediate messages.
        """
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

        # variables for measuring (see `self._take_measurement` for details)
        self.measurements   = None
        self.nom            = 0     # number of recorded measurements
        self.threshold_low  = .05
        self.threshold_high = .95
        self.previous_V     = 0.
        self.uptime         = -10.


    def _log(self, arg):
        if self.plot_interval != 0:
            print(arg)


    def _take_measurement(self, current_V, current_t):
        """is called from inside of `self.integrate` each step to check if at
        the given reference point `current_V` a rising or falling pulse
        occured.
        in such a case write quantities to dict `self.measurements`.

        quantites to be measured
        ========================
        conduction velocity (cv)
        base cycle length (bcl, period)
        action potential peak width (apd)
        rise time (rt)
        maximal peak height (vmax)
        """
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
        """integrate the system using a simple euler step

        params
        ======
        max_nom :   maximal number of measurements to record.
                    the default value is 0x7fff (=32767) which should be large
                    enough to not interfere with regular plotting usage

        notes
        =====
        every `plot_interval` steps the current potential array is yielded for
        plotting, setting `plot_interval = 0` disables this behaviour.

        because of the `yield` keyword, calling `gen = self.integrate()`
        returns a so called `generator object`. one can request each next value
        by calling e.g. `next(gen)` or `list(gen)` (the last command receives
        all values).

        this is also necessary for `plot_interval = 0` (of course no value will
        be yielded in that case)
        """
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

            # check if required number of measurements was reached
            if self.nom >= max_nom:
                print('recorded maximal number of measurements after %d steps.'
                      ' Terminating...' % i)
                break

            # switching to periodic boundary condition, after enough time for
            # the pulse to develop
            if i == 4_000:
                boundary = 'periodic'
                self._log('switched to peridodic boundary conditions.\n')

            # yield current data for plotting
            if self.plot_interval != 0 and i % self.plot_interval == 0:
                yield V



def collect_measurements(name, values, tmax, discard_list, nom):
    """Given a parameter name and an array of values, integrate the system for
    each given value and record conduction velocity (cv), period (base cycle
    length, bcl), action potential peak width (adp), peak rise time (rt) and
    maximal peak height (vmax).
    For each quantity compute and mean and std deviation.

    params
    ======
    name            :   name of param to vary
    values          :   array of values for given param
    tmax            :   max integration time
    discard_list    :   especially for short domain length (small xmax) the
                        peak needs to pass through the domain several times
                        before stabilizing.
                        if only interested in equilibrium, the values in this
                        array specify how many of the measured values to
                        exclued from mean and std, i.e.:
                            measured_mean = mean(measured_array[discard:])
    nom             :   number of measurements to record per value

    returns
    =======
    dict: {quantity : [mean_list, std_list]}
    """
    # params for Pulse1d
    params = dict(xmax=100, dx=.2, tmax=tmax, eta=.3, plot_interval=0,
                  a=.15, k=8., e0=2e-3, m1=.2, m2=.3)

    # record for each value mean and std of each quantity
    result_dict = dict(cv=[[], []], bcl=[[], []], apd=[[], []],
                       rt=[[], []], vmax=[[], []])

    for value, discard in zip(values, discard_list):
        print('%s = %s' % (name, value))

        # create Pulse1d instance for current value
        params[name] = value
        p = Pulse1d(**params)

        # integrate system until desired number of measurements is reached.
        # p.integrate returns a `generator` object (because of the `yield`
        # keyword). calling `list` on it, causes it to be executed until
        # depletion
        list(p.integrate(nom+discard))

        print('recorded %d measurements\n' \
              % len(p.measurements['rt']))

        # calculate and save mean and std for the measured quantities
        for key, val in p.measurements.items():
            result_dict[key][0].append(np.mean(val[discard:]))
            result_dict[key][1].append(np.std(val[discard:]))

    return result_dict



def measurements_for_varying_xmax():
    print('\nTAKING MEASUREMENTS FOR VARYING DOMAIN LENGTH\n')

    xmax = np.array([25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200])
    iterations_until_equilibrium = \
            np.array([5, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1])
    measurement_dict = \
            collect_measurements('xmax', xmax, tmax=3000, nom=10,
                                 discard_list=iterations_until_equilibrium)

    errorbar_config = \
            dict(ecolor='black', capsize=3, ls='none',
                 marker='o', mec='red', mfc='red')

    fig1 = plt.figure()
    fig1.suptitle('ALL MEASUREMENTS')
    for i, (key, val) in enumerate(measurement_dict.items()):
        ax = fig1.add_subplot(2, 3, i+1, xlabel='length', ylabel=key)
        ax.errorbar(xmax, val[0], yerr=val[1], **errorbar_config)

    rt = measurement_dict['rt']
    cv = measurement_dict['cv']
    rt_lin = np.linspace(28, 42)

    def cv_th(rt, a=.17):
        return np.sqrt(.3 / (2.*rt) * (1. - 2.*a))

    fig2, ax = plt.subplots()
    ax.set(xlabel='rise time', ylabel='conduction velocity')
    ax.errorbar(rt[0], cv[0], xerr=rt[1], yerr=cv[1],
                label='Measurements', **errorbar_config)
    ax.plot(rt_lin, cv_th(rt_lin), label='Theoretical Estimate')
    ax.legend()

    return measurement_dict



def animation_of_propagating_action_potential():
    print('\n\n\nPLOTTING ACTION POTENTIAL PROPAGATION\n')

    params = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
    p = Pulse1d(25, .2, 1000, .3, 50, **params)
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

    return fig, anim



if __name__ == '__main__':
    md = measurements_for_varying_xmax()
    #  fig, anim = animation_of_propagating_action_potential()

    plt.show()


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
