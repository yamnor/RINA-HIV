import io

from sklearn.manifold import TSNE

import numpy as np
import pandas as pd

import streamlit as st

from matplotlib import pyplot as plt
from matplotlib import collections


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
  
  xy = np.loadtxt('data/node.txt')

  node = pd.DataFrame(xy, columns=['x', 'y'])
  node['chain'] = ['A' for n in range(99)] + ['B' for n in range(99)]
  node['resid'] = [n+1 for n in range(99)] + [n+1 for n in range(99)]
  
  edge = pd.read_table('data/a_rin.fraction')
  for x in ['x', 'y']:
    for i in ['i', 'j']:
      xi = []
      for n in range(len(edge)):
        xi.append(node.at[edge.at[n, i]-1, x])
      edge[f'{x}{i}'] = xi

  lines = []
  for row in edge.itertuples():
    lines.append([(row.xi, row.yi), (row.xj, row.yj)])
  
  lc = collections.LineCollection(lines, linewidth = 2 * edge['vdw'], alpha = 0.5)

  figsize = [6.0, 4.0]
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

  ax.scatter(node[node['chain'] == 'A'].x, node[node['chain'] == 'A'].y, s = 25, alpha = 0.5)
  ax.scatter(node[node['chain'] == 'B'].x, node[node['chain'] == 'B'].y, s = 25, alpha = 0.5)
  ax.add_collection(lc)
  st.pyplot(fig)

if __name__ == "__main__":
    main()
