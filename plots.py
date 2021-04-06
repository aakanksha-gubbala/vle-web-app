import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
# import scipy.optimize as opt
from antoine import get_psat


def isothermalPlots(x1, y1, P, p1_s, p2_s):
    style.use('classic')

    x = np.linspace(0, 1, 50)

    P_raoult = x * p1_s + (1 - x) * p2_s
    y_raoult = x * p1_s / P_raoult

    fig1 = plt.figure(facecolor='white')
    plt.title(r"$P-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1.2 * max(P))
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$P\ (kPa)$')
    plt.scatter(x1, P)
    plt.plot(x, P_raoult, label=r"$Raoult's\ law$", color='black')
    plt.legend(loc='best', frameon=False)

    fig2 = plt.figure(facecolor='white')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(r"$y-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$y_1$')
    plt.scatter(x1, y1)
    plt.plot(x, y_raoult, label=r"$Raoult's\ law$", color='black')
    plt.plot(x, x, color='black')
    plt.legend(loc='best', frameon=False)

    return fig2, fig1


def isobaricPlots(x1, y1, T):
    style.use('classic')

    x = np.linspace(0, 1, 10)

    fig1 = plt.figure(facecolor='white')
    plt.title(r"$T-x-y$")
    plt.xlim(0, 1)
    plt.ylim(0.98 * min(T), 1.02 * max(T))
    plt.xlabel(r'$x_1, y_1$')
    plt.ylabel(r'$T\ (K)$')
    plt.scatter(x1, T, label=r'$x_1$', color='blue')
    plt.scatter(y1, T, label=r'$y_1$', color='orange')
    plt.legend(loc='best', fontsize=10, frameon=False)

    fig2 = plt.figure(facecolor='white')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(r"$y-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$y_1$')
    plt.scatter(x1, y1)
    plt.plot(x, x, color='black')

    return fig2, fig1
