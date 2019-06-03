# -*- coding: utf-8 -*-
"""Aliev-Panfilov model"""

import numpy as np
import matplotlib.pyplot as plt


class Alpha:
    def __init__(self, V0, W0, tmax, dt, a, k, m1, m2, e0):
        self.__dict__.update(dict(a=a, k=k, m1=m1, m2=m2, e0=e0))
        self.t      = np.arange(0, tmax, dt)
        self.steps  = self.t.size
        self.V      = np.zeros(self.steps)
        self.W      = np.zeros(self.steps)
        self.V[0]   = V0
        self.W[0]   = W0
        self.fig    = None

    def GV(self, V, W):
        return -self.k * V * (V - self.a) * (V - 1) - V * W

    def GW(self, V, W):
        return self.e(V, W) * (-self.k * V * (V - self.a - 1) - W)

    def e(self, V, W):
        return self.e0 + self.m1 * W / (V + self.m2)

    def integrate(self):
        for i in range(self.steps-1):
            V_i, W_i = self.V[i], self.W[i]
            self.V[i+1] = V_i + self.dt * self.GV(V_i, W_i)
            self.W[i+1] = W_i + self.dt * self.GW(V_i, W_i)

    def plot(self, use_existing=True):
        if self.fignum is None or not use_existing:
            self.fignum = plt.figure().number
        else:
            plt.figure(self.fignum)
        self.integrate()
        plt.subplot(121)
        plt.plot(self.V, self.W)
        plt.subplot(122)
        plt.plot(self.t, self.V, self.t, self.W)




#  vim: set ff=dos tw=79 sw=4 ts=8 et ic ai :
