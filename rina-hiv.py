import numpy as np
import pandas as pd

import streamlit as st

from matplotlib import pyplot as plt
from matplotlib import collections

def mkfig(w, h):
  figsize = [w, h]
  subplot = {
    'left':   0.10,
    'right':  0.10,
    'bottom': 0.10,
    'top':    0.10,
    'wspace': 1.50,
    'hspace': 2.00,
    'grid': True,
  }
  with plt.style.context('matplotlibrc'):
    plt.rcParams["figure.figsize"]        = figsize
    plt.rcParams["figure.subplot.left"]   = subplot['left'] / figsize[0]
    plt.rcParams["figure.subplot.right"]  = 1.00 - subplot['right'] / figsize[0]
    plt.rcParams["figure.subplot.bottom"] = subplot['bottom'] / figsize[1]
    plt.rcParams["figure.subplot.top"]    = 1.00 - subplot['top'] / figsize[1]
    plt.rcParams["figure.subplot.wspace"] = subplot['wspace'] / figsize[0]
    plt.rcParams["figure.subplot.hspace"] = subplot['hspace'] / figsize[1]
    plt.rcParams["axes.grid"]             = subplot['grid']
    fig, ax = plt.subplots()
  return fig, ax

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
  #  np.savetxt('node.txt', xy)

  xy = np.loadtxt('data/node.txt')
  
  #angle = st.slider('Rotate', -180, 180, 175)
  angle = 175
  theta = angle * 3.14 / 180
  xy = np.array([[np.cos(theta) * x[0] - np.sin(theta) * x[1], np.sin(theta) * x[0] + np.cos(theta) * x[1]] for x in xy])

  node = pd.DataFrame(xy, columns=['x', 'y'])
  node['chain'] = ['A' for n in range(99)] + ['B' for n in range(99)] + ['C']
  node['resid'] = [n+1 for n in range(99)] + [n+1 for n in range(99)] + [1]

  rin_file = st.file_uploader("Choose a RIN fraction file")

  if rin_file is not None:

    rin_selected = st.radio(
      "Which interactions would you like to examine?", ('hb', 'vdw', 'ss', 'ion', 'pp', 'pc', 'iac', 'any'),
      horizontal = True)

    rin_threshold = st.slider('Threshold', 0.1, 1.0, 0.1)

    showresid = st.checkbox('ResID')

    edge = pd.read_table(rin_file)
    for x in ['x', 'y']:
      for i in ['i', 'j']:
        xi = []
        for n in range(len(edge)):
          xi.append(node.at[edge.at[n, i]-1, x])
        edge[f'{x}{i}'] = xi
    
    lines = []
    fracs = []
    for idx, row in edge[edge[rin_selected] > rin_threshold].iterrows():
      lines.append([(row['xi'], row['yi']), (row['xj'], row['yj'])])
      fracs.append(row[rin_selected] * 2)

    lc = collections.LineCollection(lines, linewidth = fracs, alpha = 0.5)

    fig, ax = mkfig(6.0, 4.0)

    for c in ['A', 'B', 'C']:
      ax.scatter(node[node['chain'] == c].x, node[node['chain'] == c].y, s = 25, alpha = 0.5, label = c)
    ax.add_collection(lc)
    if showresid:
      for n in range(len(xy)): ax.annotate(node.at[n, 'resid'], xy = xy[n], size = 6)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for axis in ['top', 'left', 'bottom', 'right']:
      ax.spines[axis].set_linewidth(0.1)

    st.pyplot(fig)

    rin_file2 = st.file_uploader("Choose another RIN fraction file")

    if rin_file2 is not None:

      rin_selected2 = st.radio(
        "Which interactions would you like to examine?", ('hb', 'vdw', 'ss', 'ion', 'pp', 'pc', 'iac', 'any'),
        key = 'rin_selected2',
        horizontal = True)

      rin_threshold2 = st.slider('Threshold', 0.1, 1.0, 0.1, key = 'rin_threshold2')

      showresid2 = st.checkbox('ResID', key = 'showresid2')

      edge2 = pd.read_table(rin_file2)
      for x in ['x', 'y']:
        for i in ['i', 'j']:
          xi = []
          for n in range(len(edge2)):
            xi.append(node.at[edge2.at[n, i]-1, x])
          edge2[f'{x}{i}'] = xi

      for idx, row in edge.iterrows():
        row2 = edge2[(edge2['i'] == row.i) & (edge2['j'] == row.j)]
        st.write(row2)
        edge2.loc[(edge2['i'] == row.i) & (edge2['j'] == row.j), rin_selected2] = row[rin_selected2] - row2[rin_selected2]
      
      lines2 = []
      fracs2 = []
      color2 = []
      for idx, row in edge2[abs(edge2[rin_selected2]) > rin_threshold2].iterrows():
        lines2.append([(row['xi'], row['yi']), (row['xj'], row['yj'])])
        fracs2.append(abs(row[rin_selected2]) * 4)
        if row[rin_selected2] > 0.0:
          color2.append('b')
        else:
          color2.append('r')

      lc2 = collections.LineCollection(lines2, linewidth = fracs2, colors = color2, alpha = 0.5)

      fig, ax = mkfig(6.0, 4.0)

      for c in ['A', 'B', 'C']:
        ax.scatter(node[node['chain'] == c].x, node[node['chain'] == c].y, s = 25, alpha = 0.5, label = c)
      ax.add_collection(lc2)
      if showresid2:
        for n in range(len(xy)): ax.annotate(node.at[n, 'resid'], xy = xy[n], size = 6)
      ax.xaxis.set_visible(False)
      ax.yaxis.set_visible(False)
      for axis in ['top', 'left', 'bottom', 'right']:
        ax.spines[axis].set_linewidth(0.1)

      st.pyplot(fig)

if __name__ == "__main__":
    main()
