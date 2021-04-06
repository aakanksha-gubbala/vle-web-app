import streamlit as st
import numpy as np
import matplotlib
from matplotlib import style
import matplotlib.pyplot as plt
import pandas as pd
import requests
import scipy.constants as constants
from scipy.special import xlogy
from antoine import get_psat
import models.wilson
from sklearn import metrics
import lxml
import html5lib


def main():
    st.title('Isobaric Binary VLE Data')
    st.markdown("Vapor-liquid equlibrium data of 30 important components from "
                "[Dortmund Data Bank](http://www.ddbst.com/en/EED/VLE/VLEindex.php) can be accessed from here. "
                "Find out which pair of components have isobaric data available and see the $y-x$, $T-x-y$ and $G^E-x$ graphs.")

    style.use("classic")

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

        isobaric_vledata = []
        P = []
        for i, data in enumerate(dataframes):
            col = data.columns
            if col.dtype == object:
                if len(col) == 3 and 'T' in col[0] and 'x1' in col[1] and 'y1' in col[2]:
                    P.append(float(dataframes[i - 1][1]))
                    isobaric_vledata.append(dataframes[i])

        if isobaric_vledata == []:
            st.error('There is no isobaric data available at DDBST')
        else:
            for i in range(len(P)):
                st.write('%d)' % (i + 1), 'P = ', P[i], 'kPa')
                st.write(isobaric_vledata[i])
            if len(P) == 1:
                choice = 1
            else:
                choice = st.number_input('Choose a dataset', value=1, min_value=1, max_value=len(P))

            st.info('Analysing dataset %d ...' % choice)
            T = isobaric_vledata[choice - 1]['T [K]']
            x1 = np.array(isobaric_vledata[choice - 1]['x1 [mol/mol]'])
            y1 = np.array(isobaric_vledata[choice - 1]['y1 [mol/mol]'])
            P = P[choice - 1]
            st.write(r'$P = %0.3f kPa$' % P)

            n_points = len(x1) - 1

            p1sat = get_psat(compounds[i1], T[0])
            p2sat = get_psat(compounds[i2], T[0])

            if p1sat > p2sat:
                st.info('The more volatile component is %s' % menu_options[i1])
                s1, s2 = compounds[i1], compounds[i2]
            else:
                st.info('The more volatile component is %s' % menu_options[i2])
                s1, s2 = compounds[i2], compounds[i1]

            x = np.linspace(0, 1, 50)

            p1_s = get_psat(s1, T)
            p2_s = get_psat(s2, T)

            try:
                if x1[0] == 0 and x1[n_points] != 1:
                    x1, y1, T = x1[1:], y1[1:], T[1:]
                if x1[0] != 0 and x1[n_points] == 1:
                    x1, y1, T = x1[:n_points], y1[:n_points], T[:n_points]
                if x1[0] == 0 and x1[n_points] == 1:
                    x1, y1, T = x1[1:n_points], y1[1:n_points], T[1:n_points]
            except KeyError or IndexError:
                pass

            gamma1 = []
            gamma2 = []
            y_roault = []
            for i in range(len(x1)):
                gamma1.append(np.divide(P * y1[i], x1[i] * p1_s[i]))
                gamma2.append(np.divide(P * (1 - y1[i]), (1 - x1[i]) * p2_s[i]))

            G_e = constants.R * T * (xlogy(x1, gamma1) + xlogy(1 - x1, gamma2))

            fig1 = plt.figure(facecolor='white')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.title(r"$y-x$")
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.xlabel(r'$x_1$')
            plt.ylabel(r'$y_1$')
            plt.scatter(x1, y1)
            plt.plot(x, x, color='grey')

            st.write(fig1)

            fig2 = plt.figure(facecolor='white')
            plt.title(r"$T-x$")
            plt.xlim(0, 1)
            plt.ylim(0.98 * min(T), 1.02 * max(T))
            plt.xlabel(r'$x_1$')
            plt.ylabel(r'$T\ (K)$')
            plt.scatter(x1, T, label=r'$x_1$', color='blue')
            plt.scatter(y1, T, label=r'$y_1$', color='green')
            plt.legend(loc='best', fontsize=10, frameon=False)

            st.write(fig2)

            model = st.selectbox('Choose a model', ['Select', 'NRTL'])
            if model == "Select":
                st.info("Select a model")
            if model == 'NRTL':
                st.info("Under development ...")
    except:
        ''

    st.sidebar.title("Note")
    st.sidebar.info(""" The saturation pressures are obtained from 
            [DDBST's database](http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Ethanol).""")



'''[theme]
primaryColor="#0013f7"
backgroundColor="#ffffff"
secondaryBackgroundColor="#12e279"
textColor="#000000"
font="monospace"   
'''