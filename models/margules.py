import streamlit as st
import scipy.constants as constants
import numpy as np
import scipy.optimize as opt
from sklearn import metrics
import matplotlib.pyplot as plt
from matplotlib import style


class Margules:
    def Ge(x, A):
        return A * x * (1 - x)

    def gamma1(x, A, T):
        return np.exp(A * (1 - x) ** 2 / (constants.R * T))

    def gamma2(x, A, T):
        return np.exp(A * x ** 2 / (constants.R * T))


@st.cache(suppress_st_warning=True)
def get_parameter(x, G_e):
    A, params_cov = opt.curve_fit(Margules.Ge, x, G_e, p0=1000, maxfev=10000)
    return A

@st.cache(suppress_st_warning=True)
def get_accuracy(G_e, Ge):
    return metrics.r2_score(G_e, Ge)



def main(x1, y1, P, G_e, x, p1_s, p2_s, T, P_raoult):
    style.use('classic')

    A = get_parameter(x1, G_e)
    acc = get_accuracy(G_e, Margules.Ge(x1, A))


    fig4 = plt.figure(facecolor='white')
    plt.title(r"$G^E-x$")
    plt.xlim(0, 1)
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$G^E\ (J/mol)$')
    plt.scatter(x1, G_e)
    plt.plot(x, Margules.Ge(x, A), label=r"$Margules\ model$", color='red')
    plt.axhline(0, color='black')


    P_margules = x * p1_s * Margules.gamma1(x, A, T) + (1 - x) * p2_s * Margules.gamma2(x, A, T)
    y_margules = x * p1_s * Margules.gamma1(x, A, T) / P_margules

    fig5 = plt.figure(facecolor='white')
    plt.title(r"$P-x$")
    plt.xlim(0, 1)
    plt.ylim(0, 1.2 * max(P_margules))
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$P\ (kPa)$')
    plt.scatter(x1, P)
    plt.plot(x, P_margules, label=r"$Margules\ model$", color='red')
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
    plt.plot(x, y_margules, label=r"$Margules\ model$", color='red')
    plt.plot(x, x, color='black')
    plt.legend(loc='best', fontsize=10)

    return A, acc, fig4, fig5, fig6


