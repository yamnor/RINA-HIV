import streamlit as st
import streamlit.components.v1 as components

from pyvis.network import Network
import networkx as nx

st.title('Hello Pyvis')

nx_graph = nx.cycle_graph(10)
nx_graph.nodes[1]['title'] = 'Number 1'
nx_graph.nodes[1]['group'] = 1
nx_graph.nodes[3]['title'] = 'I belong to a different group!'
nx_graph.nodes[3]['group'] = 10
nx_graph.add_node(20, size=20, title='couple', group=2)
nx_graph.add_node(21, size=10, title='couple', group=2)
nx_graph.add_edge(20, 21, weight=1)
nx_graph.add_node(25, size=25, label='lonely', title='lonely node', group=3)

nt = Network("500px", "500px", notebook = True, heading='')
nt.from_nx(nx_graph)
nt.show('test.html')

HtmlFile = open("test.html", 'r', encoding='utf-8')
source_code = HtmlFile.read() 
components.html(source_code, height = 900,width=900)
