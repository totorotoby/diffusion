import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import csr_array
import sys

row_ptr = np.genfromtxt(sys.argv[1], delimiter='\n')
cols = np.genfromtxt(sys.argv[2], delimiter='\n')
perm = np.genfromtxt(sys.argv[3], delimiter='\n')
vx = np.genfromtxt(sys.argv[4], delimiter='\n')
vy = np.genfromtxt(sys.argv[5], delimiter='\n')

locs = {}
for i in range(len(row_ptr) - 1):
    locs[i] = (vx[i], vy[i])

new_label = {i:int(perm[i]) for i in range(len(row_ptr)-1)}
new_locs = {}
for i in range(len(row_ptr) - 1):
    new_locs[perm[i]] = (vx[i], vy[i])



adj_matrix = csr_array((np.ones(len(cols)), cols, row_ptr), shape=(len(row_ptr) - 1, len(row_ptr) - 1))

mesh = nx.from_scipy_sparse_array(adj_matrix)
fig, (ax1, ax2) = plt.subplots(1, 2)
nx.draw(mesh, pos=locs, with_labels=True, ax=ax1)
mesh = nx.relabel_nodes(mesh, new_label)
nx.draw(mesh, pos=new_locs, with_labels=True, ax=ax2)
plt.show()
