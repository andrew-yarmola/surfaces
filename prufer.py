"""
Prufer trees used to encode pants decompositonsb
"""

import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_pylab import draw_networkx_labels

def trivalent_prufer_codes(nodes):
  """ Generates all trivalent Prufer codes on given internal nodes. """
  if len(nodes) == 0:
    yield []
  else :
    n = nodes[0]
    for p in trivalent_prufer_codes(nodes[1:]):
      for i in range(len(p) + 1):
        yield [n] + p[:i] + [n] + p[i:]

def surface_prufer_codes(genus):
  """ We have internal nodes 2g through 4g-2 and outer nodes 0 through 2g - 1.
    They are paired as 0 to 1, 2 to 3, etc """
  nodes = range(2 * genus, 4 * genus - 2)
  for c in trivalent_prufer_codes(nodes):
    yield c

def prufer_tree(code):
    """ Returns graph associated to prufer code. """
    n = len(code) + 2
    valence = [1] * n
    for i in code:
      valence[i] +=1
  
    ptr = valence.index(1)
    leaf = ptr
    edges = []
    for i in code:
      edges.append((leaf, i))
      valence[i] -= 1
      if valence[i] == 1 and i < ptr:
        leaf = i
      else:
        ptr = valence.index(1, ptr + 1)
        leaf = ptr
    edges.append((leaf, n - 1)) 

    tree = nx.Graph()
    tree.add_edges_from(edges)
    return tree

def glued_graph(tree, genus):
  """ Takes a trivalent Prufer tree and closes it up to a trivalent graph. """
  graph = nx.MultiGraph(tree)
  for i in range(genus):
    n1 = list(graph.neighbors(2*i))[0]
    n2 = list(graph.neighbors(2*i+1))[0]
    graph.remove_nodes_from((2*i, 2*i + 1))
    graph.add_edge(n1, n2)
  return graph

def unique_surface_prufer_codes(genus):
  reps = []
  iso_reps = []
  for c in surface_prufer_codes(genus):
    g = glued_graph(prufer_tree(c), genus)
    seen = False
    for h in iso_reps:
      if nx.is_isomorphic(g, h):
        seen = True
        break
    if not seen:
      iso_reps.append(g)
      reps.append(c)
  return reps

def multigraph_plot(G):
  pos = nx.kamada_kawai_layout(G)
  nx.draw_networkx_nodes(G, pos, node_color = 'r', node_size = 100, alpha = 1)
  draw_networkx_labels(G, pos)
  ax = plt.gca()
  for e in G.edges:
    if len(e) == 2:
      e = (e[0], e[1], 0)
    if e[0] != e[1]:
      ax.annotate("",
        xy=pos[e[0]], xycoords='data',
        xytext=pos[e[1]], textcoords='data',
        arrowprops=dict(arrowstyle="-", color="0.5",
            shrinkA=5, shrinkB=5,
            patchA=None, patchB=None,
            connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
            ),
            ),
        )
    else:
      c = pos[e[0]]
      r = 0.03
      ax.add_artist(plt.Circle((c[0]-r, c[1]), r, color = "b", fill = False))
  plt.axis('off')
  plt.show()

