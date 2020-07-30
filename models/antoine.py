import numpy as np
import pandas as pd


def get_psat(s, T):
    antoine_url = 'http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=' + s

    antoine = pd.read_html(antoine_url)[6]
    antoine = antoine.drop(antoine.index[0:3]).drop('No.', axis=1)

    if len(antoine) == 4:
        Tmin = np.array([float(antoine['Tmin'][3]), float(antoine['Tmin'][4])])

        Tmax = np.array([float(antoine['Tmax'][3]), float(antoine['Tmax'][4])])

        A = np.array([float(antoine['A'][3]), float(antoine['A'][3])])
        B = np.array([float(antoine['B'][3]), float(antoine['B'][3])])
        C = np.array([float(antoine['C'][3]), float(antoine['C'][3])])

        for i in range(2):
            if Tmin[i] <= (T - 273.15) <= Tmax[i]:
                psat = (101.325 / 760) * np.power(10, A[i] - B[i] / (T - 273.15 + C[i]))  # in kPa
    else:
        Tmin = antoine['Tmin'][3]

        Tmax = float(antoine['Tmax'][3])

        A = float(antoine['A'][3])
        B = float(antoine['B'][3])
        C = float(antoine['C'][3])

    def psat(T):
        return (101.325 / 760) * np.power(10, A - B / (T - 273.15 + C))  # in kPa
    return psat(T)




