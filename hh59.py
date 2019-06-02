# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


ALPHA = {
    'n' :   lambda V: .01 * (V - 10.) / (1. - np.exp((10. - V) / 10.)),
    'm' :   lambda V: .1 * (V - 25.) / (1. - np.exp((25. - V) / 10.)),
    'h' :   lambda V: .07 * np.exp(-V / 20.)
}

BETA = {
    'n' :   lambda V: .125 * np.exp(-V / 80.),
    'm' :   lambda V: 4. * np.exp(-V / 18.),
    'h' :   lambda V: 1. / (1. + np.exp((30. - V) / 10.))
}


def I_K(V, n):
    return 36. * n**4 * (V + 12)

def I_Na(V, m, h):
    return 120. * m**3 * h * (V - 115)

def I_l(V):
    return .3 * (V - 10)


def gate(f, V, name):
    """gating variable according to
        f_n+1 = f_n + dt * gate(f_n)
    where f in {n, m, h}

    params
    ======
    f   :   current value of gating variable
    V   :   current value of potential
    name:   name of gating variable, used as key for ALPHA and BETA

    returns
    =======
    gate(f) = alpha * (1 - f) - beta * f
    """
    return ALPHA[name](V) * (1. - f) - BETA[name](V) * f


def VG(V, n, m, h, Ip=0):
    """action potential according to
        V_n+1 = V_n + dt * VG(V)

    params
    ======
    V       :   current value of action potential
    n, m, h :   current values of gating variables (to compute currents)
    Ip      :   source (ion pumps) current at current time

    returns
    =======
    VG(V, n, m, h) = -1/c * (sum I_s(V, n, m, h))
    """
    c = 1.
    return -(I_K(V, n) + I_Na(V, m, h) + I_l(V) - Ip) / c


# simulation parameters
#  dt      = 0.055  # max dt for reasonable results
dt      = 0.01
tmax    = 25.
t       = np.arange(0., tmax, dt)
steps   = t.size
Ip      = np.zeros_like(t)
#  Ip[int(20./dt):int(80./dt)] += 10.
V_0 = -7.

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

# integrate ODEs with Euler steps
for i in range(steps-1):
    V[i+1] = V[i] + dt * VG(V[i], n[i], m[i], h[i], Ip[i])
    for g, key in zip((n, m, h), ('n', 'm', 'h')):
        g[i+1] = g[i] + dt * gate(g[i], V[i], key)


# plotting
plt.figure()
plt.subplot(131,
            title='Potential',
            xlabel=r'$t/ms$', ylabel=r'$-V/mV$')
plt.plot(t, V)

plt.subplot(132,
            title='Currents',
            xlabel=r'$t/ms$', ylabel=r'$I/Am^{-2}$')
plt.plot(t, I_K(V, n),
         t, I_Na(V, m, h),
         t, I_l(V))
plt.legend([r'$I Ka^+$', r'$I Na^+$', r'$I$ leak'],
           loc='upper right')

plt.subplot(133,
            title='Gating Variables',
            xlabel=r'$t/ms$', ylabel='Gating Variables')
plt.plot(t, n**4,
         t, m**3,
         t, h)
plt.legend([r'Gate $Ka^+$ $n^4$', r'Gate $Na^+$ $m^3$', r'Gate $Na^+$ $h$'],
           loc='upper right')


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
