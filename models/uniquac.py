# import streamlit as st
import scipy.constants as constants
import numpy as np
import scipy.optimize as opt
from sklearn import metrics
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd
from scipy.constants import R

uniquac_params = pd.read_csv("uniquac_params.txt", sep='!', names=['compound', 'r', 'q'])


def get_params(s1, s2):
    r = np.zeros(2)
    q = np.zeros(2)
    for s, i in zip([s1, s2], range(2)):
        r[i] = uniquac_params.loc[uniquac_params['compound'] == s]["r"]
        q[i] = uniquac_params.loc[uniquac_params['compound'] == s]["q"]
    return r, q


class UNIQUAC:
    def __init__(self, s1, s2):
        self.s1 = s1
        self.s2 = s2
        self.r, self.q = get_params(self.s1, self.s2)


    def gamma1(self, X, A, B):
        [x, T] = X
        self.phi1 = x * self.r[0] / (x * self.r[0] + (1 - x) * self.r[1])
        self.theta1 = x * self.q[0] / (x * self.q[0] + (1 - x) * self.q[1])
        self.theta2 = 1 - self.theta1
        self.phi2 = 1 - self.phi1
        self.l1 = 5 * (self.r[0] - self.q[0]) - (self.r[0] - 1)
        self.l2 = 5 * (self.r[1] - self.q[1]) - (self.r[1] - 1)
        ln_gam_c1 = np.log(self.phi1 / x) + 5 * self.q[0] * np.log(self.theta1 / self.phi1) + self.l1 - self.phi1 / x * (x * self.l1 + (1 - x) * self.l2)
        gam_c1 = np.exp(ln_gam_c1)

        tau12 = np.exp(-A / (R * T))
        tau21 = np.exp(-B / (R * T))

        ln_gam_r1 = self.q[0] * (1 - np.log(self.theta1 + self.theta2 * tau21) -
                                 self.theta1 / (self.theta1 + self.theta2 * tau21) -
                                 self.theta2 * tau12 / (self.theta1 * tau12 + self.theta2))
        gam_r1 = np.exp(ln_gam_r1)
        return gam_r1 * gam_c1

    def gamma2(self, X, A, B):
        [x, T] = X
        ln_gam_c2 = np.log(self.phi2 / (1 - x)) + 5 * self.q[1] * np.log(self.theta2 / self.phi2) + self.l2 - \
                    self.phi2 / (1 - x) * (x * self.l1 + (1 - x) * self.l2)
        gam_c2 = np.exp(ln_gam_c2)

        tau12 = np.exp(-A / (R * T))
        tau21 = np.exp(-B / (R * T))

        ln_gam_r2 = self.q[1] * (1 - np.log(self.theta1 * tau12 + self.theta2) -
                                 self.theta1 * tau21 / (self.theta1 + self.theta2 * tau21) -
                                 self.theta2 / (self.theta1 * tau12 + self.theta2))
        gam_r2 = np.exp(ln_gam_r2)
        return gam_r2 * gam_c2

    def costfunction(self, params, X, gamma):
        [A, B] = params
        residuals = (np.concatenate((self.gamma1(X, A, B), self.gamma2(X, A, B))) - gamma) / gamma
        return residuals

    # @st.cache(suppress_st_warning=True)
    def get_parameter(self, X, gamma):
        params = opt.least_squares(self.costfunction, [1000, 1000], args=(X, gamma))
        return params

