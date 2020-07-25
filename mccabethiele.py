import streamlit as st
import numpy as np
from matplotlib import pyplot, style
from scipy.optimize import fsolve


np.seterr(divide='ignore', invalid='ignore')


def main():
    st.title('McCabe-Thiele Plot Generator')
    st.write('The McCabe-Thiele method is used to determine the number of equilibrium stages for a distillation column.')

    style.use('classic')

    if st.checkbox('General Conditions'):
        F = st.number_input('Feed Flow Rate', value=100.000)
        zf = st.number_input('Feed concentration', value=0.500)
        xd = st.number_input('Distillate concentration', value=0.900)
        xb = st.number_input('Bottoms concentration', value=0.100)
        R = st.number_input('Reflux Ratio', value=3.000)
        q = st.number_input('Thermal Quality', value=1.000)
        a = st.number_input('Relative Volatility', value=2.500)


        def dbf(f):
            return [xd * f[0] + xb * f[1] - zf * F, f[0] + f[1] - F]


        [D, B] = fsolve(dbf, [30, 20])

        st.write('Distillate: ', round(D, 4), 'Bottoms: ', round(B, 4))

        Lr = R * D
        Ls = Lr + q * F
        Vr = Lr + D
        Vs = Vr + (q - 1) * F


        def x_eq(x):
            # x on the equlibrium curve
            return x / (a * (1 - x) + x)


        def rec_opline(x):
            # rectifying section operating line: y = (Lr/Vr)*x + (D*xd/Vr)
            return (Lr / Vr) * x + (D * xd / Vr)


        def strip_opline(x):
            # stripping section operating line: y = (Ls/Vs)*x - (B*xb/Vs)
            return (Ls / Vs) * x - (B * xb / Vs)


        # intersection point of rectifying opline and stripping opline
        def inter_pt(p):
            return [(Lr / Vr) * p[0] + (D * xd / Vr) - p[1], (Ls / Vs) * p[0] - (B * xb / Vs) - p[1]]


        [xq, yq] = fsolve(inter_pt, [0.5, 0.5])


        def min_reflux(x):
            return [q * x[0] - zf - x[1] * (q - 1), x[1] - a * x[0] / (1 + x[0] * (a - 1))]


        [xrmin, yrmin] = fsolve(min_reflux, [0.5, 0.6])

        slope = (xd - yrmin) / (xd - xrmin)
        Rmin = slope / (1 - slope)
        st.write("The minimum reflux ratio is ", round(Rmin, 4))

        # y-x equilibrium curve
        x = np.linspace(0, 1, 10000)
        y = a * x / (1 + x * (a - 1))

        gen = pyplot.figure(figsize=(7, 7), facecolor='white')
        pyplot.suptitle("McCabe-Thiele Plot")
        pyplot.xlabel('x')
        pyplot.ylabel('y')
        pyplot.plot(x, y, color='black', linewidth=1)
        pyplot.plot(x, x, color='black', linewidth=1)
        pyplot.xlim(0, 1)
        pyplot.ylim(0, 1)
        pyplot.grid(color='grey', linewidth=0.3)

        # rectifying section operating line
        x, y = [xq, xd], [yq, xd]
        pyplot.plot(x, y, label='rectifying section', color='b', linewidth=1)

        # stripping section operating line
        x, y = [xq, xb], [yq, xb]
        pyplot.plot(x, y, label='stripping section', color='g', linewidth=1)
        pyplot.legend(loc='best')

        x0, y0 = xd, xd
        for i in range(1, 100):
            x1, y1 = x_eq(y0), y0
            pyplot.plot([x0, x1], [y0, y1], color='r', linewidth=1)
            if x1 > xq:
                x2, y2 = x1, rec_opline(x1)
            if x1 < xq:
                x2, y2 = x1, strip_opline(x1)
            if (x2, y2) < (xb, xb):
                pyplot.plot([x1, x2], [y1, x2], color='r', linewidth=1)
            else:
                pyplot.plot([x1, x2], [y1, y2], color='r', linewidth=1)
            x0, y0 = x2, y2
            if y2 < x2:
                break

        pyplot.plot([xd, xd], [0, xd], linestyle='--', linewidth=1)
        pyplot.plot([xb, xb], [0, xb], linestyle='--', linewidth=1)
        pyplot.title("Number of stages = infinity", size=10) if i == 99 \
            else pyplot.title("Number of stages = %d" % i, size=10)

        st.write(gen)

    if st.checkbox('Total Reflux Conditions'):
        xd_tr = st.number_input('Top concentration', value=0.900)
        xb_tr = st.number_input('Bottom concentration', value=0.100)
        a_tr = st.number_input('Relative Volatility (average)', value=2.500)


        def x_eq(x):
            # x on the equlibrium curve
            return x / (a_tr * (1 - x) + x)


        # y-x equilibrium curve
        x = np.linspace(0, 1, 10000)
        y = a_tr * x / (1 + x * (a_tr - 1))

        tr = pyplot.figure(figsize=(7, 7), facecolor='white')
        pyplot.suptitle("McCabe-Thiele Plot - Total Reflux")
        pyplot.xlabel('x')
        pyplot.ylabel('y')
        pyplot.plot(x, y, color='black', linewidth=1)
        pyplot.plot(x, x, color='black', linewidth=1)
        pyplot.xlim(0, 1)
        pyplot.ylim(0, 1)
        pyplot.grid(color='grey', linewidth=0.3)

        x0, y0, xb_tr = xd_tr, xd_tr, xb_tr
        for i in range(1, 100):
            x1, y1 = x_eq(y0), x0
            pyplot.plot([x0, x1], [y0, y1], color='r', linewidth=1)
            x2, y2 = x1, x1
            pyplot.plot([x1, x2], [y1, x2], color='r', linewidth=1)
            x0, y0 = x2, y2
            if x2 < xb_tr:
                break
        pyplot.plot([xd_tr, xd_tr], [0, xd_tr], linestyle='--', linewidth=1)
        pyplot.plot([xb_tr, xb_tr], [0, xb_tr], linestyle='--', linewidth=1)
        pyplot.title("Number of stages = infinity", size=10) if i == 99 \
            else pyplot.title("Number of stages = %d" % i, size=10)

        st.write(tr)
