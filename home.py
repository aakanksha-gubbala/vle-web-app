import streamlit as st


def main():
    st.title("Welcome.")

    st.subheader("Get VLE data from [DDBST's](http://www.ddbst.com/en/EED/VLE/VLEindex.php) online database")
    st.write("The algorithm scours through all the vapor-liquid equilibrium data at DDBST and"
             " finds the complete VLE data of the chosen pair of compounds")
    # st.image('ddbst_site.png', width=750)

    st.subheader("Try different correlative models ")
    st.write("Note: These models can only be used for binary isothermal vapor-liquid equilibrium data.")
    st.markdown(r'''|Model|$G^E$|$RTln\gamma_1$|$RTln\gamma_2$|
                    |:---:|:---:|:---:|:---:|
                    |Margules|$Ax_1x_2$|$Ax_2^2$|$Ax_1^2$|
                    |Redlich-Kister|$x_1x_2(A+B(x_1-x_2))$|$x_2^2(A-B+4Bx_1)$|$x_1^2(A+B-4Bx_2)$|
                    |van-Laar|$\frac{x_1x_2}{G^E}=\frac{x_1}{B} + \frac{x_2}{A}$|$\frac{A}{\big(1+\frac{Ax_1}{Bx_2}\big)^2}$|$\frac{B}{\big(1+\frac{Bx_2}{Ax_1}\big)^2}$| 
                    |Truncated Wohls expansion|$RT(2a_{12}z_1z_2(x_1q_1+x_2q_2))$ $z_i = \frac{x_iq_i}{\sum_j x_jq_j}$|$RT(2a_{12}q_1z_2^2)$|$RT(2a_{12}q_2z_1^2)$|''')

    st.subheader("Get McCabe-Thiele Plots instantly")
    st.write(
        "The McCabe-Thiele method is used to determine the number of equilibrium stages for a distillation column. "
        "The vapor-liquid equilibrium data is a crucial part of analysis and design of a distillation column.")
    st.write(
        r''' This tool gives the flow rates of distillate and bottoms streams, minimum reflux ratio and it has an option for total reflux conditions too.''')

    # st.subheader("Updates")
    # st.text("23/3/21: Added NRTL model")
    # st.text("29/7/20: Added Truncated Wohls model, enabled viewing of isobaric datasets.")
    # st.text("25/7/20: Added correlative models.")
