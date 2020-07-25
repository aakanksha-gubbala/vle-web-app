import streamlit as st


def main():
    st.title("Welcome.")
    st.write("This site is mainly a collection of projects related to thermodynamics.")
    st.text("Last updated: 25/7/20")

    st.subheader("Get VLE data from DDBST's online database")
    st.write("The algorithm scours through all the vapor-liquid equilibrium data at DDBST and"
             " finds the complete isothermal data of the chosen pair of compounds")
    st.image('ddbst_site.png', width=750, format="PNG")

    st.subheader("Try different correlative models ")
    st.write("Note: These models can only be used for binary isothermal vapor-liquid equilibrium data.")
    st.markdown(r'''|$Model$|$G^E$|$RTln\gamma_1$|$RTln\gamma_2$|
|:---:|:---:|:---:|:---:|
|$Margules$|$Ax_1x_2$|$Ax_2^2$|$Ax_1^2$|
|$Redlich-Kister$|$x_1x_2(A+B(x_1-x_2))$|$x_2^2(A-B+4Bx_1)$|$x_1^2(A+B-4Bx_2)$|
|$van-Laar$|$\frac{x_1x_2}{G^E}=\frac{x_1}{B} + \frac{x_2}{A}$|$\frac{A}{\big(1+\frac{Ax_1}{Bx_2}\big)^2}$|$\frac{B}{\big(1+\frac{Bx_2}{Ax_1}\big)^2}$| ''')
    st.markdown(
        ' Different models have been illustrated for the chosen dataset. To read about these models, check out these resources:')
    st.markdown(
        '* [Chemical, Biochemical, and Engineering Thermodynamics (5th Ed.) by Stanley I. Sandler](https://www.wiley.com/en-us/Chemical%2C+Biochemical%2C+and+Engineering+Thermodynamics%2C+5th+Edition-p-9781119321286)')
    st.markdown(
        '* [Introduction to Chemical Engineering Thermodynamics (7th Ed.) by Smith, Van Ness, and Abbott](http://www.learncheme.com/screencasts/thermodynamics/textbook-SVNA-7th)')

    st.subheader("Get McCabe-Thiele Plots instantly")
    st.write(
        "The McCabe-Thiele method is used to determine the number of equilibrium stages for a distillation column. "
        "The vapor-liquid equilibrium data is a crucial part of analysis and design of a distillation column.")
    st.write(
        r''' This tool gives the flow rates of distillate and bottoms streams, minimum reflux ratio and it has an option for total reflux conditions too.''')
