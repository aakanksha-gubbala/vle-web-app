import streamlit as st
import scipy.constants as constants
import numpy as np
import scipy.optimize as opt
from sklearn import metrics
import matplotlib.pyplot as plt
from matplotlib import style


class RK2:
    def Ge(x, A, B):
        return x * (1 - x) * (A + B * (x - (1 - x)))

    def gamma1(x, A, B, T):
        return np.exp((A - B + 4 * B * x) * (1 - x) ** 2 / (constants.R * T))

    def gamma2(x, A, B, T):
        return np.exp((A + B - 4 * B * (1 - x)) * (x) ** 2 / (constants.R * T))


@st.cache(suppress_st_warning=True)
def get_parameter(x, G_e):
    [A, B], params_cov = opt.curve_fit(RK2.Ge, x, G_e, p0=[1000,1000], maxfev=10000)
    return [A, B]


@st.cache(suppress_st_warning=True)
def get_accuracy(G_e, Ge):
    return metrics.r2_score(G_e, Ge)


def main(x1, y1, P, G_e, x, p1_s, p2_s, T, P_raoult):
    style.use('classic')

    [A, B] = get_parameter(x1, G_e)
    acc = get_accuracy(G_e, RK2.Ge(x1, A, B))

    fig4 = plt.figure(facecolor='white')
    plt.title(r"$G^E-x$")
    plt.xlim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$G^E\ (J/mol)$')
    plt.scatter(x1, G_e)
    plt.plot(x, RK2.Ge(x, A, B), label=r"$RK2\ model$", color='red')
    plt.axhline(0, color='black')

    P_rk = x * p1_s * RK2.gamma1(x, A, B, T) + (1 - x) * p2_s * RK2.gamma2(x, A, B, T)
    y_rk = x * p1_s * RK2.gamma1(x, A, B, T) / P_rk

    fig5 = plt.figure(facecolor='white')
    plt.title(r"$P-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1.2 * max(P))
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$P\ (kPa)$')
    plt.scatter(x1, P)
    plt.plot(x, P_rk, label=r"$RK2\ model$", color='red')
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
    plt.plot(x, y_rk, label=r"$RK2\ model$", color='red')
    plt.plot(x, x, color='black')
    plt.legend(loc='best', fontsize=10)

    return [A, B], acc, fig4, fig5, fig6



