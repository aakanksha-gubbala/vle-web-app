import streamlit as st
import scipy.constants as constants
import numpy as np
import scipy.optimize as opt
from scipy.special import xlogy
from sklearn import metrics
import matplotlib.pyplot as plt
from matplotlib import style
from volume import get_volume
from antoine import get_psat


class Wohls:
    def __init__(self, s1, s2, T):
        self.q1 = get_volume(s1, T)
        self.q2 = get_volume(s2, T)
        self.T = T

    def Ge(self, x1, A):
        x1 = np.asfarray(x1, float)
        z1 = x1 * self.q1 / (x1 * self.q1 + (1 - x1) * self.q2)
        z2 = (1 - x1) * self.q2 / (x1 * self.q1 + (1 - x1) * self.q2)
        return constants.R * self.T * (x1 * self.q1 + (1 - x1) * self.q2) * (2 * A) * (z1 * z2)

    def gamma1(self, z, A):
        return np.exp(2 * A * self.q1 * (1 - z) ** 2)

    def gamma2(self, z, A):
        return np.exp(2 * A * self.q2 * z ** 2)

    # @st.cache(suppress_st_warning=True)
    def get_parameter(self, x, G_e):
        A, params_cov = opt.curve_fit(self.Ge, x, G_e, p0=1000, maxfev=10000)
        return A

    # @st.cache(suppress_st_warning=True)
    def get_accuracy(self, G_e, x1):
        A = self.get_parameter(x1, G_e)
        Ge = self.Ge(x1, A)
        return metrics.r2_score(G_e, Ge)


def main(x1, y1, P, G_e, T, s1, s2):
    style.use('classic')

    w = Wohls(s1, s2, T)
    A = w.get_parameter(x1, G_e)
    acc = w.get_accuracy(G_e, x1)

    x = np.linspace(0, 1, 50)

    fig4 = plt.figure(facecolor='white')
    plt.title(r"$G^E-x$")
    plt.xlim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$G^E\ (J/mol)$')
    plt.scatter(x1, G_e)
    plt.plot(x, w.Ge(x, A), label=r"$Wohls\ model$", color='red')
    plt.axhline(0, color='black')

    q1 = get_volume(s1, T)
    q2 = get_volume(s2, T)
    z = x * q1 / (x * q1 + (1 - x) * q2)

    p1_s = get_psat(s1, T)
    p2_s = get_psat(s2, T)

    P_Wohls = x * p1_s * w.gamma1(z, A) + (1 - x) * p2_s * w.gamma2(z, A)
    y_Wohls = x * p1_s * w.gamma1(z, A) / P_Wohls

    P_raoult = x * p1_s + (1 - x) * p2_s

    fig5 = plt.figure(facecolor='white')
    plt.title(r"$P-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1.2 * max(P_Wohls))
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
