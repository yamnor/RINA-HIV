from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
from matplotlib import collections
import numpy as np
import pandas as pd

import streamlit as st
import io

def main():

  st.title("RIN Analyzer")
  
  node = {}

  uploaded_file = st.file_uploader("Choose a PDB file")

  if uploaded_file is not None:
    stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    xyz = np.array([[float(x) for x in line.split()[5:8]] for line in stringio.readlines() if ' CA ' in line])
    
    xy = TSNE(n_components = 2, random_state = 1, perplexity = 50).fit_transform(xyz)[:, ::-1]
    xy = np.array([[x[0], -x[1]] for x in xy])
    theta = 8 * 3.14 / 180
    xy = np.array([[np.cos(theta) * x[0] - np.sin(theta) * x[1], np.sin(theta) * x[0] + np.cos(theta) * x[1]] for x in xy])
    node = {'A': xy[0:99].T, 'B': xy[99:199].T}
    
    fig, ax = plt.subplots()
    for c in ['A', 'B']:
      ax.scatter(*node[c], s = 25, alpha = 0.5, label = c)
    st.pyplot(fig)

if __name__ == "__main__":
    main()
