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
EToV = np.loadtxt(sys.argv[6], usecols = range(3), dtype=int)

locs = {}
for i in range(len(row_ptr) - 1):
    locs[i] = (vx[i], vy[i])

new_label = {i:int(perm[i]) for i in range(len(row_ptr)-1)}
new_locs = {}
for i in range(len(row_ptr) - 1):
    new_locs[perm[i]] = (vx[i], vy[i])



adj_matrix = csr_array((np.ones(len(cols)), cols, row_ptr), shape=(len(row_ptr) - 1, len(row_ptr) - 1))

mesh = nx.from_scipy_sparse_array(adj_matrix)
#fig, (ax1, ax2) = plt.subplots(1, 2)
nx.draw(mesh, pos=locs, with_labels=True)#, ax=ax1)
mesh = nx.relabel_nodes(mesh, new_label)
#nx.draw(mesh, pos=new_locs, with_labels=True, ax=ax2)

for (i, e) in enumerate(EToV):

    vx0 = vx[e[0]]
    vy0 = vy[e[0]]
    vx1 = vx[e[1]]
    vy1 = vy[e[1]]
    vx2 = vx[e[2]]
    vy2 = vy[e[2]]

    x = (vx0 + vx1 + vx2) / 3
    y = (vy0 + vy1 + vy2) / 3
    
    plt.text(x, y, str(i), fontsize=12, color='red', ha='left', va='bottom')

#ax2.set_visible(False)
plt.show()

