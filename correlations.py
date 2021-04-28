import time
import streamlit as st
import numpy as np
import pandas as pd
import requests
import scipy.constants as constants
from scipy.special import xlogy
from antoine import get_psat
from volume import get_volume
import models.margules, models.redlichkister, models.vanlaar, models.alphagm, models.wohls
import lxml


def main():
    st.title("Isothermal Binary VLE Data")
    st.write(
        """ The *Margules* model, *Redlich-Kister Expansion* truncated to two terms, *van Laar* model and the *Truncated Wohls expansion* are implemented here. """
        r"In case $\alpha$ fits the data with an accuracy of 80% or above, the $\alpha_{GM}$ value is displayed.")

    compounds = ['Acetonitrile', 'Acetone', '1,2-Ethanediol', 'Ethanol',
                 'Diethyl ether', 'Ethyl acetate', 'Benzene', '1-Butanol',
                 'Chloroform', 'Cyclohexane', 'Acetic acid butyl ester', 'Acetic acid',
                 'Hexane', '2-Propanol', '1-Hexene', 'Methanol',
                 'Tetrahydrofuran', 'Water', 'm-Xylene', 'p-Xylene', '1,3-Butadiene', 'Hexadecane']

    menu_options = compounds.copy()

    for i, compound in enumerate(compounds):
        if ' ' in compound:
            compounds[i] = compound.replace(' ', '%20')

    compound1 = st.selectbox('Select compound 1', menu_options, key='compound1')
    compound2 = st.selectbox('Select compound 2', menu_options, key='compound2')

    i1 = menu_options.index(compound1)
    i2 = menu_options.index(compound2)

    st.info("You have chosen %s and %s" % (compound1, compound2))

    # @st.cache(suppress_st_warning=True)
    def link_generator(i1, i2):
        url = 'http://www.ddbst.com/en/EED/VLE/VLE%20' + compounds[i1] + '%3B' + compounds[i2] + '.php'
        if requests.get(url).status_code == 404:
            url = 'http://www.ddbst.com/en/EED/VLE/VLE%20' + compounds[i2] + '%3B' + compounds[i1] + '.php'
        return url

    try:
        if compound1 == compound2:
            st.warning('Choose different compounds')
        else:
            url = link_generator(i1, i2)

        if requests.get(url).status_code == 404:
            st.error("VLE data for this pair of compounds doesn't exist at DDBST.")

        dataframes = pd.read_html(url)

        isothermal_vledata = []
        T = []
        for i, data in enumerate(dataframes):
            col = data.columns
            if col.dtype == object:
                if len(col) == 3 and 'P' in col[0] and 'x1' in col[1] and 'y1' in col[2]:
                    T.append(float(dataframes[i - 1][1]))
                    isothermal_vledata.append(dataframes[i])

        if isothermal_vledata == []:
            st.error('There is no isothermal data available for this pair of compounds at DDBST')
        else:
            for i in range(len(T)):
                st.write('%d)' % (i + 1), 'T = ', T[i], 'K')
                st.write(isothermal_vledata[i])
            if len(T) == 1:
                choice = 1
            else:
                choice = st.number_input('Choose a dataset', value=1, min_value=1, max_value=len(T))

            st.info('Analysing dataset %d ...' % choice)
            P = isothermal_vledata[choice - 1]['P [kPa]']
            x1 = isothermal_vledata[choice - 1]['x1 [mol/mol]']
            y1 = isothermal_vledata[choice - 1]['y1 [mol/mol]']
            T = T[choice - 1]
            st.write(r'$T = %0.2f K$' % T)

            p1sat = get_psat(compounds[i1], T)
            p2sat = get_psat(compounds[i2], T)

            if p1sat > p2sat:
                st.info('The more volatile component is %s' % menu_options[i1])
                s1, s2 = compounds[i1], compounds[i2]
            else:
                st.info('The more volatile component is %s' % menu_options[i2])
                s1, s2 = compounds[i2], compounds[i1]

            p1_s = max(p1sat, p2sat)
            p2_s = min(p1sat, p2sat)

            st.write(r'$p_1^s = %0.3f kPa$' % p1_s)
            st.write(r'$p_2^s = %0.3f kPa$' % p2_s)

            x = np.linspace(0, 1, 50)
            P_raoult = x * p1_s + (1 - x) * p2_s
            y_raoult = x * p1_s / P_raoult

            n_points = len(x1) - 1

            try:
                if x1[0] == 0 and x1[n_points] != 1:
                    x1, y1, P = x1[1:], y1[1:], P[1:]
                if x1[0] != 0 and x1[n_points] == 1:
                    x1, y1, P = x1[:n_points], y1[:n_points], P[:n_points]
                if x1[0] == 0 and x1[n_points] == 1:
                    x1, y1, P = x1[1:n_points], y1[1:n_points], P[1:n_points]
            except KeyError:
                pass

            q1, q2 = get_volume(s1, T), get_volume(s2, T)
            z1 = x1 * q1 / (x1 * q1 + (1 - x1) * q2)

            gamma1 = np.divide(P * y1, x1 * p1_s)
            gamma2 = np.divide(P * (1 - y1), ((1 - x1) * p2_s))
            G_e = constants.R * T * (xlogy(x1, gamma1) + xlogy(1 - x1, gamma2))

            alpha_gm = models.alphagm.get_alpha_gm(x1, y1)
            if alpha_gm == 0:
                pass
            else:
                st.success(r"$\alpha_{GM}=%0.3f$" % alpha_gm)
                st.text("Try using this value in the McCabe-Thiele Plotter!")

            model = st.selectbox("Choose a model",
                                 ["Select", "Margules", "Redlich Kister", "van Laar", "Truncated Wohls"], key='model')

            if model == "Select":
                st.info("Select a model")
            else:
                if model != "Truncated Wohls":
                    MODELS = {"Margules": models.margules, "Redlich Kister": models.redlichkister,
                              "van Laar": models.vanlaar}
                    latest_iteration = st.empty()
                    bar = st.progress(0)

                    for i in range(100):
                        latest_iteration.text(f'{i + 1}%')
                        bar.progress(i + 1)
                        time.sleep(0.03)

                    A, acc, fig4, fig5, fig6 = MODELS[model].main(x1, y1, P, G_e, x, p1_s, p2_s, T, P_raoult)

                    if model == "Margules":
                        st.write(r"$G^E = %0.3fx_1x_2$" % A)
                    if model == "Redlich Kister":
                        st.write(r"$G^E = x_1x_2(%0.3f + (%0.3f)(x_1-x_2))$" % (A[0], A[1]))
                    if model == "van Laar":
                        st.write(r"$\frac{x_1x_2}{G^E} = \frac{x_1}{%0.3f} + \frac{x_2}{%0.3f}$" % (A[1], A[0]))
                    st.write(r"$R^2$ score = %0.3f" % acc)
                    st.write(fig4, fig5, fig6)
                else:
                    latest_iteration = st.empty()
                    bar = st.progress(0)

                    for i in range(100):
                        latest_iteration.text(f'{i + 1}%')
                        bar.progress(i + 1)
                        time.sleep(0.03)

                    A, acc, fig4, fig5, fig6 = models.wohls.main(x1, y1, P, G_e, T, s1, s2)

                    st.write(r"Molar volumes: $q_1=%0.3e$, $q_2=%0.3e$" % (q1, q2))
                    st.write(r"$\frac{G^E/RT}{x_1q_1 + x_2q_2} = 2(%0.3f)z_1z_2$" % A)
                    st.write(r"$R^2$ score = %0.3f" % acc)
                    st.write(fig4, fig5, fig6)
    except:
        pass

    st.sidebar.title("Note")
    st.sidebar.info(""" The saturation pressures are obtained from 
            [DDBST's database](http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Ethanol).""")
    st.sidebar.info(""" The densities are obtained from 
                [DDBST's database](http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=Diethyl%20ether).""")
    st.sidebar.info("These models can only be used for binary isothermal vapor-liquid equilibrium data.")
