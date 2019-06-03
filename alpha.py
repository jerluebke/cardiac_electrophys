# -*- coding: utf-8 -*-
"""Aliev-Panfilov model"""

import numpy as np
import matplotlib.pyplot as plt


class Alpha:
    def __init__(self, V0, W0, tmax, dt, a, k, e0, m1, m2):
        # set simulation parameters
        self.__dict__.update(dict(a=a, k=k, e0=e0, m1=m1, m2=m2))
        self.t      = np.arange(0, tmax, dt)
        self.dt     = dt
        self.steps  = self.t.size

        # result arrays
        self.V      = np.zeros(self.steps)
        self.W      = np.zeros(self.steps)
        self.V[0]   = V0
        self.W[0]   = W0

        # plotting
        self.fig        = plt.figure()
        self.phaseplot  = self.fig.add_subplot(121,
                                               title='phase plot',
                                               xlabel='V', ylabel='W')
        self.timeplot   = self.fig.add_subplot(122,
                                               title='development of V, W',
                                               xlabel='t/ms', ylabel='V/mV')


    def G(self, V, W):
        return (
            -self.k * V * (V - self.a) * (V - 1) - V * W,
            self.e(V, W) * (-self.k * V * (V - self.a - 1) - W)
        )


    def e(self, V, W):
        return self.e0 + self.m1 * W / (V + self.m2)


    def integrate(self):
        for i in range(self.steps-1):
            V_i, W_i    = self.V[i], self.W[i]
            V_dot, W_dot= self.G(V_i, W_i)
            self.V[i+1] = V_i + self.dt * V_dot
            self.W[i+1] = W_i + self.dt * W_dot


    def portrait(self):
        V_mesh, W_mesh = np.mgrid[0:1.1:.1,0:2.5:.2]
        V_dot, W_dot = self.G(V_mesh, W_mesh)
        self.phaseplot.quiver(V_mesh, W_mesh, V_dot, W_dot,
                              color='grey', scale=15, width=.004)


    def plot(self, legend_entry=None):
        self.integrate()
        self.phaseplot.plot(self.V, self.W)
        self.timeplot.plot(12.9*self.t, 100*self.V-80,
                           12.9*self.t, 100*self.W-80)
        self.timeplot.legend(['action potential V', 'relaxation variable W'],
                             loc='upper right')


    def vary_param(self, key, val_list):
        colors = plt.get_cmap('Dark2').colors
        for val, color in zip(val_list, colors):
            setattr(self, key, val)
            self.integrate()
            self.phaseplot.plot(self.V, self.W,
                                color=color,
                                label='%s=%s' % (key, val))
            self.timeplot.plot(12.9*self.t, 100*self.V-80,
                               color=color)
            self.timeplot.plot(12.9*self.t, 100*self.W-80,
                               color=color,
                               label='%s=%s' % (key, val))
        self.phaseplot.legend(loc='upper right')
        self.timeplot.legend(loc='upper right')



params  = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
alpha   = Alpha(.2, .0, 60, .01, **params)
alpha.portrait()
alpha.plot()

alpha   = Alpha(.2, .0, 60, .01, **params)
alpha.vary_param('a', [.1, .125, .15, .175, .199, .2])

alpha   = Alpha(.2, .0, 60, .01, **params)
alpha.vary_param('k', [4., 6., 8., 10., 12., 14.])

alpha   = Alpha(.2, .0, 60, .01, **params)
alpha.vary_param('e0', [1e-5, 1e-4, 2e-3, 1e-1, 1.])

alpha   = Alpha(.2, .0, 60, .01, **params)
alpha.vary_param('m1', [.05, .1, .2, .5, 1., 2.])

alpha   = Alpha(.2, .0, 60, .01, **params)
alpha.vary_param('m2', [.1, .3, 1., 2., 5.])


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
