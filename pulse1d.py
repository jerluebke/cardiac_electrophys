# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from alpha import Alpha, PARAMS


def Lap1d(a, o, dx):
    o[:,1:-1] = (a[:,1:-1] + a[:,:-2]) / dx**2
    o[:,0], o[:,-1] = o[:,1], o[:,-2]
    return o



class Pulse1d(Alpha):
    def __init__(self, xmax, dx, eta, **params):
        self.x      = np.arange(0, xmax, dx)
        self.dx     = dx
        self.xsteps = self.x.size
        params['dt']= 0.5 * dx**2 / (2 * eta)
        self.eta    = eta
        self.L_VW   = np.zeros((2, self.xsteps))
        super().__init__(_extra_dim=(self.xsteps,),
                         _prepare_plot=False,
                         **params)


    def integrate(self):
        VW, L_VW    = self.VW, self.L_VW
        dx, dt, eta = self.dx, self.dt, self.eta

        for i in range(self.steps-1):
            Lap1d(VW[:,i,:], L_VW, dx)
            VW[:,i+1,:] =                       \
                    VW[:,i,:]                   \
                    + dt * self.G(*VW[:,i,:])   \
                    + dt * eta * L_VW



#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
