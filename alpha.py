# -*- coding: utf-8 -*-
"""Aliev-Panfilov model: a phenomenological model of a myocytes action
potential depending on two variables (given by two coupled ODEs) and five
parameters.

Compared to Hodgkin-Huxley (or similar), this model is much less complex as all
the biophysical details of ion pump dynamics are summarized in a single
phenomenological relaxation variable W. Thus it is more suitable for large
scale simulations of tissue.
"""

import numpy as np
import matplotlib.pyplot as plt



class AlphaBase:
    def __init__(self, a, k, e0, m1, m2):
        self.__dict__.update(dict(a=a, k=k, e0=e0, m1=m1, m2=m2))


    def G(self, V, W):
        """equations for V, w:
            [dV/dt, dW/dt] = G(V, W)
        """
        return np.array([
            -self.k * V * (V - self.a) * (V - 1) - V * W,
            self.e(V, W) * (-self.k * V * (V - self.a - 1) - W)
        ])


    def e(self, V, W):
        return self.e0 + self.m1 * W / (V + self.m2)


    def integrate(self):
        pass



class Alpha(AlphaBase):
    """Aliev-Panfilov simulation class"""

    def __init__(self, V0, W0, tmax, dt,
                 _prepare_plot=True, **params):
        # set simulation parameters
        super().__init__(**params)
        self.t      = np.arange(0, tmax, dt)
        self.dt     = dt
        self.steps  = self.t.size

        # result array
        self.VW     = np.zeros((2, self.steps))
        self.VW[0,0]= V0
        self.VW[1,0]= W0

        # plotting
        if _prepare_plot:
            self.fig        = plt.figure()
            self.phaseplot  = self.fig.add_subplot(121,
                                                   title='phase plot',
                                                   xlabel='V', ylabel='W')
            self.timeplot   = self.fig.add_subplot(122,
                                                   title='development of V, W',
                                                   xlabel='t/ms', ylabel='V/mV')


    def integrate(self):
        """integrate the ODE system with a simple Euler-step
        """
        for i in range(self.steps-1):
            VW_i = self.VW[:,i]
            self.VW[:,i+1] = VW_i + self.dt * self.G(*VW_i)


    def portrait(self):
        """create a quiver plot of [dV/dt, dW/dt] to visualize the phase space
        dynamics of the ODEs
        """
        V_mesh, W_mesh = np.mgrid[0:1.1:.1,0:2.5:.2]
        V_dot, W_dot = self.G(V_mesh, W_mesh)
        self.phaseplot.quiver(V_mesh, W_mesh, V_dot, W_dot,
                              color='grey', scale=15, width=.004)


    def plot(self):
        """integrate the ODEs for the paarameters of the instance and plot the
        results, showing the development of action potential V and relaxation W
        over time.

        the obtained results need to be scaled in order to correspond to the
        'true' physical quantities:
            V_phys = (100*V-80) mV
            t_phys = 12.9*t ms
        """
        self.integrate()
        self.phaseplot.plot(self.VW[0], self.VW[1])
        self.timeplot.plot(12.9*self.t, 100*self.VW[0]-80,
                           12.9*self.t, 100*self.VW[1]-80)
        self.timeplot.legend(['action potential V', 'relaxation variable W'],
                             loc='upper right')


    def vary_param(self, key, val_list):
        """given the name of a parameter of the simulation and a list of its
        possible values, integrate the system for each value and plot the
        results, in order to get an understanding of how the given parameter
        influences the dynamics of the system
        """
        colors = plt.get_cmap('Dark2').colors
        for val, color in zip(val_list, colors):
            setattr(self, key, val)
            self.integrate()
            self.phaseplot.plot(self.VW[0], self.VW[1],
                                color=color,
                                label='%s=%s' % (key, val))
            self.timeplot.plot(12.9*self.t, 100*self.VW[0]-80,
                               color=color)
            self.timeplot.plot(12.9*self.t, 100*self.VW[1]-80,
                               color=color,
                               label='%s=%s' % (key, val))
        self.phaseplot.legend(loc='upper right')
        self.timeplot.legend(loc='upper right')



PARAMS = dict(
    V0=.2, W0=.0, tmax=60., dt=.01,
    a=.15, k=8., e0=2e-3, m1=.2, m2=.3
)


if __name__ == '__main__':
    alpha = Alpha(**PARAMS)
    alpha.portrait()
    alpha.plot()

    alpha = Alpha(**PARAMS)
    alpha.vary_param('a', [.1, .125, .15, .175, .199, .2])

    alpha = Alpha(**PARAMS)
    alpha.vary_param('k', [4., 6., 8., 10., 12., 14.])

    alpha = Alpha(**PARAMS)
    alpha.vary_param('e0', [1e-5, 1e-4, 2e-3, 1e-1, 1.])

    alpha = Alpha(**PARAMS)
    alpha.vary_param('m1', [.05, .1, .2, .5, 1., 2.])

    alpha = Alpha(**PARAMS)
    alpha.vary_param('m2', [.1, .3, 1., 2., 5.])


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
