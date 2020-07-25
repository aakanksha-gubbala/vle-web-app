import time
import streamlit as st
import numpy as np
import lxml
import ddbst, home, mccabethiele, correlations


np.seterr(divide='ignore', invalid='ignore')


PAGES = {"Home": home, "VLE Data": ddbst, "Correlations": correlations, "McCabe-Thiele Plots": mccabethiele}

st.sidebar.title("Navigation")
selection = st.sidebar.radio("Go to", ["Home", "VLE Data", "Correlations", "McCabe-Thiele Plots"])

with st.spinner(f'Loading {selection} ...'):
    PAGES[selection].main()


st.sidebar.title("About")
st.sidebar.info(""" This is an open source project maintained by [kraanky](https://kraanky.github.io/). 
Here's the [link](https://github.com/kraanky/binary-vle-data) to the source code.""")

