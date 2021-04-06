import streamlit as st
import numpy as np
import scipy.constants as constants
from scipy.special import xlogy
from antoine import get_psat
import scipy.optimize as opt
from volume import get_volume
from matplotlib import style
import matplotlib.pyplot as plt
import scipy.interpolate as interp


class Wilson:
    def __init__(self, s1, s2, T):
        self.s1 = s1
        self.s2 = s2
        self.p1_s = get_psat(s1, T)
        self.p2_s = get_psat(s2, T)
        self.q1 = get_volume(self.s1, T)
        self.q2 = get_volume(self.s2, T)

    def Ge(self, X, A, B):
        R = constants.R
        [x1, T] = X
        x2 = 1 - x1
        lam12 = (self.q2 / self.q1) * np.exp(-A / (R * T))
        lam21 = (self.q1 / self.q2) * np.exp(-B / (R * T))
        return (-x1 * np.log(x1 + lam12 * x2) - x2 * np.log(x2 + lam21 * x1)) * R * T

    def gamma1(self, X, A, B):
        R = constants.R
        [x, T] = X
        x1 = x
        x2 = 1 - x
        lam12 = (self.q2 / self.q1) * np.exp(-A / (R * T))
        lam21 = (self.q1 / self.q2) * np.exp(-B / (R * T))
        return np.exp(-np.log(x1 + lam12 * x2) + x2 * (lam12 / (x1 + x2 * lam12) - lam21 / (x2 + lam21 * x1)))

    def gamma2(self, X, A, B):
        R = constants.R
        [x, T] = X
        x1 = x
        x2 = 1 - x
        lam12 = (self.q2 / self.q1) * np.exp(-A / (R * T))
        lam21 = (self.q1 / self.q2) * np.exp(-B / (R * T))
        return np.exp(-np.log(x2 + lam21 * x1) - x1 * (lam12 / (x1 + x2 * lam12) - lam21 / (x2 + lam21 * x1)))

    @st.cache(suppress_st_warning=True)
    def get_parameter(self, X, G_e):
        [x1, T] = X
        [A, B], params = opt.curve_fit(self.Ge, [x1, T], G_e)
        return [A, B]


def main(s1, s2, G_e, P, T, x1, y1):
    style.use('classic')
    [A, B] = Wilson(s1, s2, S).get_parameter([x1, T], G_e)

    x = np.linspace(0, 1, 5)

    '''def modRaoultsLaw(S):
        wilson = Wilson(s1, s2, S)
        x1 = x
        x2 = 1 - x1
        return x1 * wilson.gamma1([x1, S], A, B) * get_psat(s1, S) + x2 * wilson.gamma2([x1, S], A, B) * get_psat(s2, S) - P

    T_wilson = opt.fsolve(modRaoultsLaw, np.linspace(min(T), max(T), 5))
    wilson = Wilson(s1, s2, T_wilson)

    y_wilson = x * wilson.gamma1([x, T_wilson], A, B) * get_psat(s1, T_wilson) / P
    Ge_wilson = wilson.Ge([x1, T1], A, B)'''
    wilson = Wilson(s1, s2, T)
    Ge_wilson = wilson.Ge([x1, T], A, B)
    y_wilson = x1 * wilson.gamma1([x1, T], A, B) * get_psat(s1, T) / P

    fig1 = plt.figure(facecolor='white')
    plt.title(r"$G^E-x$")
    plt.plot(x1, interp.make_interp_spline(x1, Ge_wilson)(x1))
    plt.scatter(x1, G_e)
    plt.xlim(0, 1)
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$G^E (J/mol)$")
    plt.axhline(0)

    fig2 = plt.figure(facecolor='white')
    plt.title(r"$T-x-y$")
    plt.plot(x1, interp.make_interp_spline(x1, T)(x1), color='blue')
    plt.plot(y_wilson, interp.make_interp_spline(y_wilson, T)(y1), color='green')
    plt.scatter(x1, T, color='blue')
    plt.scatter(y1, T, color='green')
    plt.xlim(0, 1)
    plt.xlabel(r"$x_1, y_1$")
    plt.ylabel(r"$T (K)$")

    fig3 = plt.figure(facecolor='white')
    plt.title(r"$y-x$")
    plt.plot(x1, interp.make_interp_spline(x1, y_wilson)(x1))
    plt.plot(x, x, color='black')
    plt.scatter(x1, y1)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$y_1$")

    return fig1, fig2, fig3
