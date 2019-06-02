# -*- coding: utf-8 -*-

from math import exp
import numpy as np
import matplotlib.pyplot as plt


ALPHA = {
    'n' :   lambda V: .01 * (V - 10.) / (1. - exp((10. - V) / 10.)),
    'm' :   lambda V: .1 * (V - 25.) / (1. - exp((25. - V) / 10.)),
    'h' :   lambda V: .07 * exp(-V / 20.)
}

BETA = {
    'n' :   lambda V: .125 * exp(-V / 80.),
    'm' :   lambda V: 4. * exp(-V / 18.),
    'h' :   lambda V: 1. / (1. + exp((30. - V) / 10.))
}

V_REST = {
    'nul'       :   0.,
    'muscle'    :   -80.,
    'purkinje'  :   -90.,
    'AV'        :   -65.,
    'SA'        :   -55.
}


def I_K(V, n):
    return 36. * n**4 * (V + 12),

def I_Na(V, m, h):
    return 120. * m**3 * h * (V - 115),

def I_l(V):
    return .3 * (V - 10)


def gate(f, name):
    """gating variable according to
        f_n+1 = f_n + dt * gate(f_n)
    where f in {n, m, h}

    params
    ======
    f   :   current value of gating variable
    name:   name of gating variable, used as key for ALPHA and BETA

    returns
    =======
    gate(f) = alpha * (1 - f) - beta * f
    """
    return ALPHA[name] * (1. - f) - BETA[name] * f


def VG(V, n, m, h):
    """action potential according to
        V_n+1 = V_n + dt * VG(V)

    params
    ======
    V       :   current value of action potential
    n, m, h :   current values of gating variables (to compute currents)

    returns
    =======
    VG(V, n, m, h) = -1/c * (sum I_s(V, n, m, h))
    """
    c = 1.
    return -(I_K(V, n) + I_Na(V, m, h) + I_l(V)) / c


# simulation parameters
dt      = 0.1
tmax    = 100.
t       = np.arange(0., tmax, dt)
steps   = t.size
V_0     = V_REST['nul']

# initialize result arrays
V = np.zeros(steps)
n = np.zeros(steps)
m = np.zeros(steps)
h = np.zeros(steps)

# set initial values
# assume gating variables are in initial steady state
V[0] = V_0
for g, key in zip((n, m, h), ('n', 'm', 'h')):
    g[0] = ALPHA[key](V_0) / (ALPHA[key](V_0) + BETA[key](V_0))

for i in range(steps):
    V[i+1] = V[i] + dt * VG(V[i], n[i], m[i], h[i])
    for g, key in zip((n, m, h), ('n', 'm', 'h')):
        g[i+1] = g[i] + dt * gate(g[i], key)


plt.figure()
plt.subplot(131)
plt.plot(t, V)
plt.subplot(132)
plt.plot(t, I_K(V, n),
         t, I_Na(V, m, h),
         t, I_l(V))
plt.subplot(133)
plt.plot(t, n**4, t, m**3, t, h)


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
