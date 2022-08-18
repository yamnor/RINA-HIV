from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
from matplotlib import collections
import numpy as np
import pandas as pd

import streamlit as st
import io

import altair as alt

def main():

  st.title("RIN Analyzer for HIV-1 protease")
    
  rin = ['hb', 'vdw', 'ss', 'ion', 'pp', 'pc', 'iac', 'any']

  #pdb_file = st.file_uploader("Choose a PDB file")

  #if pdb_file is not None:
  #  stringio = io.StringIO(pdb_file.getvalue().decode("utf-8"))
  #  xyz = np.array([[float(x) for x in line.split()[5:8]] for line in stringio.readlines() if ' CA ' in line])
  #  
  #  xy = TSNE(n_components = 2, random_state = 1, perplexity = 50).fit_transform(xyz)[:, ::-1]
  #  xy = np.array([[x[0], -x[1]] for x in xy])
  #  theta = 8 * 3.14 / 180
  #  xy = np.array([[np.cos(theta) * x[0] - np.sin(theta) * x[1], np.sin(theta) * x[0] + np.cos(theta) * x[1]] for x in xy])
  #
  #  np.savetxt('node.txt', xy)
  
  xy = np.loadtxt('node.txt')

  node = pd.DataFrame(xy, columns=['x', 'y'])
  node['chain'] = ['A' for n in range(99)] + ['B' for n in range(99)]
  node['resid'] = [n+1 for n in range(99)] + [n+1 for n in range(99)]
  
  frac = {k: np.zeros((198, 198)) for k in rin}

  rin_file = st.file_uploader("Choose a RIN fraction file")

  if rin_file is not None:

    n = alt.Chart(node).mark_circle().encode(x = 'x', y = 'y', color = 'chain', tooltip=['resid'])

    edge = pd.read_table(rin_file)

    for x in ['x', 'y']:
      for i in ['i', 'j']:
        xi = []
        for n in range(len(edge)):
          xi.append(node.at(edge.at[n, i]-1, x))
        edge[f'{x}{i}'] = xi


    line = []
    threshold = 0.01
    e = [x for x in edge if abs(x['any']) > threshold]
    pos = [x['pos'] for x in e]
    for k in rin:
      val = np.array([abs(x[k]) for x in e])
      col = ['b' if x[k] > 0.0 else 'r' for x in e]
      xy_ = pd.DataFrame({'x' : xy[]})
      line = alt.Chartcollections.LineCollection(pos, linewidths = val * 3.0, colors = col, alpha = 0.5)

    matome = alt.layer(*[c, *l])

    st.altair_chart(matome, use_container_width = True)

    line = {}
    threshold = 0.01
    e = [x for x in edge if abs(x['any']) > threshold]
    pos = [x['pos'] for x in e]
    for k in rin:
      val = np.array([abs(x[k]) for x in e])
      col = ['b' if x[k] > 0.0 else 'r' for x in e]
      line[k] = collections.LineCollection(pos, linewidths = val * 3.0, colors = col, alpha = 0.5)

    fig, ax = plt.subplots()
    for c in ['A', 'B']:
      ax.scatter(*node[c], s = 25, alpha = 0.5, label = c)
    for n in range(len(xy)): ax.annotate(n+1, xy = xy[n], size = 6)
    ax.add_collection(line['hb'])
    st.pyplot(fig)


if __name__ == "__main__":
    main()
