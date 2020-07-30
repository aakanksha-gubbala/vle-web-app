import streamlit as st
import numpy as np
import scipy.constants as constants
from scipy.special import xlogy
from antoine import get_psat
from volume import get_volume


class Wilson:
    def __init__(self, s1, s2, T):
        self.s1 = s1
        self.s2 = s2
        self.p1_s = get_psat(s1, T)
        self.p2_s = get_psat(s2, T)
        self.q1 = get_volume(self.s1, T)
        self.q2 = get_volume(self.s2, T)

    def exp(self, x1, y1, P):
        x2 = 1 - x1
        self.gamma1 = P * y1 / (x1 * self.p1_s)
        self.gamma2 = P * y2 / (x2 * self.p2_s)
        self.G_e = constants.R * T * (xlogy(x1, self.gamma1) + xlogy(x2, self.gamma2))

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
def get_parameter(s1, s2, X, G_e):
    [x1, T] = X
    [A, B], params = opt.curve_fit(Wilson(s1, s2, T).Ge, [x1, T], G_e)
    return [A, B]


