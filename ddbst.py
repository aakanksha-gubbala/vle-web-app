import streamlit as st
import numpy as np
from matplotlib import style
import matplotlib.pyplot as plt
import pandas as pd
import requests
import scipy.constants as constants
from scipy.special import xlogy
from antoine import get_psat
import lxml


def main():
    st.title('Isothermal Binary VLE Data')
    st.markdown("Vapor-liquid equlibrium data of 30 important components from "
                "[Dortmund Data Bank](http://www.ddbst.com/en/EED/VLE/VLEindex.php) can be accessed from here. "
                "Find out which pair of components have isothermal data available and see the $y-x$, $P-x$ and $G^E-x$ graphs.")


    compounds = ['Acetonitrile', 'Acetone', '1,2-Ethanediol', 'Ethanol',
                 'Diethyl ether', 'Ethyl acetate', 'Benzene', '1-Butanol',
                 'Chloroform', 'Cyclohexane', 'Acetic acid butyl ester', 'Acetic acid',
                 'Hexane', '2-Propanol', '1-Hexene', 'Methanol',
                 'Tetrahydrofuran', 'Water', 'm-Xylene', 'p-Xylene',
                 'N-Methyl-2-pyrrolidone', '1,3-Butadiene', 'Hexadecane']

    menu_options = compounds.copy()

    for i, compound in enumerate(compounds):
        if ' ' in compound:
            compounds[i] = compound.replace(' ', '%20')

    compound1 = st.selectbox('Select compound 1', menu_options, key='compound1')
    compound2 = st.selectbox('Select compound 2', menu_options, key='compound2')

    i1 = menu_options.index(compound1)
    i2 = menu_options.index(compound2)

    st.info("You have chosen %s and %s" % (compound1, compound2))

    @st.cache(suppress_st_warning=True)
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
            st.error('There is no isothermal data available at DDBST')
        else:
            for i in range(len(T)):
                st.write('%d)' % (i + 1), 'T = ', T[i], 'K')
                st.write(isothermal_vledata[i])
            if len(T) == 1:
                choice = 1
            else:
                choice = st.number_input('Choose a dataset', value=1, min_value=1, max_value=len(T))

            if st.button('Go!'):
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
                else:
                    st.info('The more volatile component is %s' % menu_options[i2])

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

                gamma1 = np.divide(P * y1, x1 * p1_s)
                gamma2 = np.divide(P * (1 - y1), ((1 - x1) * p2_s))
                G_e = constants.R * T * (xlogy(x1, gamma1) + xlogy(1 - x1, gamma2))

                fig1 = plt.figure(facecolor='white')
                plt.gca().set_aspect('equal', adjustable='box')
                plt.title(r"$y-x$")
                plt.xlim(0, 1)
                plt.ylim(0, 1)
                plt.xlabel(r'$x_1$')
                plt.ylabel(r'$y_1$')
                plt.scatter(x1, y1)
                plt.plot(x, y_raoult, label=r"$Raoult's\ Law$", color='black')
                plt.plot(x, x, color='grey')
                plt.legend(loc='best', fontsize=10)

                st.write(fig1)

                fig2 = plt.figure(facecolor='white')
                plt.title(r"$P-x$")
                plt.xlim(0, 1)
                plt.ylim(0, 1.2*max(P))
                plt.xlabel(r'$x_1$')
                plt.ylabel(r'$P\ (kPa)$')
                plt.scatter(x1, P)
                plt.plot(x, P_raoult, color='black', label=r"$Raoult's\ Law$")
                plt.legend(loc='best', fontsize=10)

                st.write(fig2)

                fig3 = plt.figure(facecolor='white')
                plt.title(r"$G^E-x$")
                plt.xlim(0, 1)
                plt.xlabel(r'$x_1$')
                plt.ylabel(r'$G^E\ (J/mol)$')
                plt.scatter(x1, G_e)
                plt.axhline(0, color='black')

                st.write(fig3)
    except:
        ''

    st.sidebar.title("Note")
    st.sidebar.info(""" The saturation pressures are obtained from 
    [DDBST's database](http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Ethanol).""")
