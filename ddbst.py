import streamlit as st
import numpy as np
from matplotlib import style
import pandas as pd
import requests
import scipy.constants as constants
import scipy.optimize as opt
from antoine import get_psat
import lxml
from plots import isobaricPlots, isothermalPlots


def main():
    st.title('Binary VLE Data')
    st.markdown("Vapor-liquid equlibrium data of 30 important components from "
                "[Dortmund Data Bank](http://www.ddbst.com/en/EED/VLE/VLEindex.php) can be accessed from here. ")

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

        vledata = []
        S = []
        for i, data in enumerate(dataframes):
            col = data.columns
            if col.dtype == object:
                if (len(col) == 3 and 'P' in col[0] and 'x1' in col[1] and 'y1' in col[2]) or (
                        len(col) == 3 and 'T' in col[0] and 'x1' in col[1] and 'y1' in col[2]):
                    S.append(float(dataframes[i - 1][1]))
                    vledata.append(dataframes[i])

        if vledata == []:
            st.error('Complete VLE data is not available at DDBST')
        else:
            for i in range(len(S)):
                if vledata[i].columns[0] == 'P [kPa]':
                    st.write('%d)' % (i + 1), 'T = ', S[i], 'K')
                if vledata[i].columns[0] == 'T [K]':
                    st.write('%d)' % (i + 1), 'P = ', S[i], 'kPa')
                st.write(vledata[i])
            if len(S) == 1:
                choice = 1
            else:
                choice = st.number_input('Choose a dataset', value=1, min_value=1, max_value=len(S))

            st.info('Analysing dataset %d ...' % choice)

            try:
                P = vledata[choice - 1]['P [kPa]']
                x1 = vledata[choice - 1]['x1 [mol/mol]']
                y1 = vledata[choice - 1]['y1 [mol/mol]']
                T = S[choice - 1]

                st.write(r'$T = %0.2f K$' % T)

                p1sat = get_psat(compounds[i1], T)
                p2sat = get_psat(compounds[i2], T)

                if p1sat > p2sat:
                    st.info('The more volatile component is %s' % menu_options[i1])
                else:
                    st.info('The more volatile component is %s' % menu_options[i2])

                p1_s = max(p1sat, p2sat)
                p2_s = min(p1sat, p2sat)

                fig1, fig2 = isothermalPlots(x1, y1, P, p1_s, p2_s)
                st.write(fig1, fig2)
            except:
                pass
            try:
                T = vledata[choice - 1]['T [K]']
                x1 = vledata[choice - 1]['x1 [mol/mol]']
                y1 = vledata[choice - 1]['y1 [mol/mol]']
                P = S[choice - 1]

                st.write(r'$P = %0.2f kPa$' % P)

                p1sat = get_psat(compounds[i1], T[0])
                p2sat = get_psat(compounds[i2], T[0])

                if p1sat > p2sat:
                    st.info('The more volatile component is %s' % menu_options[i1])
                    s1, s2 = compounds[i1], compounds[i2]
                else:
                    st.info('The more volatile component is %s' % menu_options[i2])
                    s1, s2 = compounds[i2], compounds[i1]

                p1_s = get_psat(s1, T)
                p2_s = get_psat(s2, T)

                fig1, fig2 = isobaricPlots(x1, y1, T)
                st.write(fig1, fig2)
            except:
                pass
    except:
        ''
    st.sidebar.title("Note")
    st.sidebar.info(""" The saturation pressures are obtained from 
            [DDBST's database](http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Ethanol).""")


