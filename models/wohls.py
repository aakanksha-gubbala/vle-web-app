import scipy.constants as constants
import numpy as np
import scipy.optimize as opt
from scipy.special import xlogy
from sklearn import metrics
import matplotlib.pyplot as plt
from matplotlib import style
from volume import get_volume


def main(x1, y1, P, gamma1, G_e, x, p1_s, p2_s, T, P_raoult, s1, s2):
    class Wohls:
        def gamma1(z, A):
            return np.exp(2 * A * q1 * (1 - z) ** 2)

        def gamma2(z, A):
            return np.exp(2 * A * q2 * z ** 2)

    def get_parameter(z, gamma1):
        A, params_cov = opt.curve_fit(Wohls.gamma1, z, gamma1, p0=1000, maxfev=10000)
        return A


    def get_accuracy(gamma1, wohls_gamma1):
        return metrics.r2_score(gamma1, wohls_gamma1)


    style.use('classic')

    q1, q2 = get_volume(s1, T), get_volume(s2, T)
    z1 = x1 * q1 / (x1 * q1 + (1 - x1) * q2)
    z = x * q1 / (x * q1 + (1 - x) * q2)

    A = get_parameter(z1, gamma1)
    acc = get_accuracy(gamma1, Wohls.gamma1(z1, A))

    Ge = constants.R*T*(xlogy(x, Wohls.gamma1(z, A)) + xlogy(1-x, Wohls.gamma2(z, A)))

    fig4 = plt.figure(facecolor='white')
    plt.title(r"$G^E-x$")
    plt.xlim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$G^E\ (J/mol)$')
    plt.scatter(x1, G_e)
    plt.plot(x, Ge, label=r"$Wohls\ model$", color='red')
    plt.axhline(0, color='black')

    z = x * q1 / (x * q1 + (1 - x) * q2)
    P_Wohls = x * p1_s * Wohls.gamma1(z, A) + (1 - x) * p2_s * Wohls.gamma2(z, A)
    y_Wohls = x * p1_s * Wohls.gamma1(z, A) / P_Wohls

    fig5 = plt.figure(facecolor='white')
    plt.title(r"$P-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1.2 * max(P))
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$P\ (kPa)$')
    plt.scatter(x1, P)
    plt.plot(x, P_Wohls, label=r"$Wohls\ model$", color='red')
    plt.plot(x, P_raoult, color='black', label=r"$Raoult's\ Law$")
    plt.legend(loc='best', fontsize=10)

    fig6 = plt.figure(facecolor='white')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(r"$y-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$y_1$')
    plt.scatter(x1, y1)
    plt.plot(x, y_Wohls, label=r"$Wohls\ model$", color='red')
    plt.plot(x, x, color='black')
    plt.legend(loc='best', fontsize=10)

    return A, acc, fig4, fig5, fig6
