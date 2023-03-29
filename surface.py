""" Surfaces"""

from prufer import *
from hyperbolic import * 

from itertools import permutations
import functools
import numpy as np
import random

def sorted_tuple(x):
  return tuple(sorted(x))

def first(x):
  return next(iter(x))

class Cuff:
  """ Unoriented cuff. """
  def __init__(self, e_or_c):
    if isinstance(e_or_c, Cuff):
      self._edge = e_or_c._edge
    else:
      self._edge = sorted_tuple(e_or_c)

  def __hash__(self):
    return hash(self._edge)

  def __lt__(self, other):
    return self._edge < other._edge

  def __eq__(self, other):
    return self._edge == other._edge

  def __repr__(self):
    return f"Cuff({self._edge})"


class Ortho:
  """ Oriented orthosegment. """
  def __init__(self, pair):
    self._cuffs = tuple(map(Cuff, pair))

  def reversed(self):
    return Ortho(self._cuffs[::-1])

  @property
  def cuffs(self):
    return self._cuffs

  def __hash__(self):
    return hash(self._cuffs)

  def __lt__(self, other):
    return self.cuffs < other.cuffs

  def __eq__(self, other):
    return self._cuffs == other._cuffs

  def __repr__(self):
    return f"Ortho({self._cuffs})"


class Pant:
  def __init__(self, edges, orientation=1, lengths=None, twists=None):
    self._orientation = orientation
    self._lengths = {Cuff(e):None for e in edges}
    self._twists = {Cuff(e):None for e in edges}
    self._ortho_lengths, self._opps, self._traversal = {}, {}, {}
    for c1, c2, c3 in permutations(edges):
      c = Cuff(c3)
      o = Ortho((c1, c2))
      self._opps[c] = o
      # traversal order follows the order of canonical points around hexagon
      def traversal_order(o, c, ornt):
        cs, ce = o.cuffs
        if cs < ce < c or ce < cs < c:
          return (0, 0)
        if cs < c < ce:
          return (ornt, 0)
        if ce < c < cs:
          return (0, ornt)
        if c < cs < ce:
          return (-ornt, ornt)
        if c < ce < cs:
          return (ornt, -ornt)
      self._traversal[o] = traversal_order(o, c, self._orientation)
      self._ortho_lengths[o] = None

    if lengths is not None and twists is not None:
      self.set_coords(lengths, twists)

  def __repr__(self):
    return (f"Pant({tuple(self._lengths)},\n" +
    f"   orientation={self.orientation},\n" +
    f"   lengths={self._lengths}\n" +
    f"   twists={self._twists})\n" +
    f"   orthos={self._ortho_lengths}")

  @property
  def orientation(self):
    return self._orientation

  def traversal(self, o):
    return self._traversal[o]

  def length(self, c_or_o):
    if isinstance(c_or_o, Cuff):
      return self._lengths[c_or_o]
    else:
      return self._ortho_lengths[c_or_o]

  def twist(self, c):
    return self._twists[c]

  def set_coords(self, lengths, twists):
    for e, l in lengths.items():
      c = Cuff(e)
      self._lengths[c] = l
      self._twists[c] = twists[e]
    for c in self._lengths:
      op = self._opps[c]
      left, right = op.cuffs
      lenO, lenL, lenR = self._lengths[c], self._lengths[left], self._lengths[right]
      numerator = np.cosh(lenO / 2) + np.cosh(lenL / 2) * np.cosh(lenR / 2)
      denominator = np.sinh(lenL / 2) * np.sinh(lenR / 2)
      l = np.arccosh(numerator/denominator)
      self._ortho_lengths[op] = l
      self._ortho_lengths[op.reversed()] = l

class Surface:
  def __init__(self, prufer_code):
    self._code = prufer_code
    self._genus = 1 + len(self._code) // 4
    self._tree = prufer_tree(self._code)
    self._pants, self._gluings = {}, {}
    self._seen_nodes = set()
    self._add_recursive(0)

  def _add_recursive(self, n):
    self._seen_nodes.add(n)
    nbrs = tuple(self._tree.neighbors(n))
    # debug: print(f"Adding node {n} with neighbors {nbrs}")
    if len(nbrs) > 1:
      # setting orietnation of the new pair of pants hexagon based on neighbor
      ornt = 1
      for m in nbrs:
        if m in self._pants: # should only trigger once since we are in a tree
          m_p = self._pants[m]
          m_orient = m_p.orientation
          p_side = pow(-1, 1 + sorted(self._tree.neighbors(m)).index(n))
          ornt = p_side * m_orient
      self._pants[n] = Pant(tuple((n, m) for m in nbrs), ornt)
    else:
      n_cuff = Cuff((n, nbrs[0]))
      m_cuff = Cuff(first(self._tree.edges(n + (1 - 2 * (n % 2)))))
      self._gluings[n_cuff], self._gluings[m_cuff] = m_cuff, n_cuff
    for m in nbrs:
      if m not in self._seen_nodes:
        self._add_recursive(m)

  def __repr__(self):
    return f"Surface({self._tree})"

  @property
  def genus(self):
    return self._genus

  @property
  def code(self):
    return self._code

  def paths(self, source = 0): 
    """ Returns paths from the source node to all leaf nodes of the Prufer tree. """ 
    number_of_leaf_nodes = 2 * self.genus
    short_paths = nx.single_source_shortest_path(self._tree, source)
    paths = {}
    for i in range(number_of_leaf_nodes):
      paths[i] = short_paths[i]
    return paths

  def set_coords(self, fn_coords):
    """ Add length and twist for each cuff from fn_coords. """
    """ We provide fn_coords as a list of tuples: [(l_i, theta_i)]_{i \in I}. """
    """ The order is that of the lexographically sorted edges of the Prufer tree. """
    cuff_data = {}
    i = 0
    for e in sorted(self._tree.edges()):
      c = Cuff(e)
      if c not in cuff_data:
        cuff_data[c] = fn_coords[i]
        if c in self._gluings:
          cuff_data[self._gluings[c]] = fn_coords[i]
        i += 1
    for _, p in self._pants.items():
      lengths = {c: cuff_data[c][0] for c in p._lengths}
      twists = {c: cuff_data[c][1] for c in p._twists}
      p.set_coords(lengths, twists)
    
  def set_random_coords(self, min_len=0.0, max_len=100.00):
    lengths = [random.uniform(min_len, max_len) for _ in range(3 * self.genus - 3)]
    twists = [random.uniform(-l/2, l/2) for l in lengths]
    self.set_coords(list(zip(lengths, twists)))
    
  def _matrices_from_path(self, path):
    """ Returns sequence of translation and rotation matrices. """
    """ We assume the starting cuff lifts to (0,infty) with the canonical """
    """ basepoint at 0 + 1j. Then we procude a sequence of matrices that map """
    """ such that they follow the ortholines and twists along the path. """
    """ Taking produce ... matrices[2] @ matrices[1] @ matriices[0] maps the """
    """ canonical point of the last cuff to the staring one (i.e. 0 + 1j). """

    if len(path) <= 1:
      return [np.eye(2)], ['identity']
    pant_path = list(self._pants[i] for i in path[1:-1])
    cuff_path = list(map(Cuff, zip(path, path[1:])))
    ortho_path = list(map(Ortho, zip(cuff_path, cuff_path[1:])))
    types = []
    matrices = []
    for i, o in enumerate(ortho_path):
      p = pant_path[i]
      c_in, c_out = cuff_path[i], cuff_path[i+1]
      ts, te = p.traversal(o)
      if c_in != cuff_path[0]: # twist only at internal cuffs
        matrices.append(translation(p.twist(c_in)))
        types.append(f"twist of {p.twist(c_in)} along {c_in}")
      # if our ortho starts at a non-canonical point, we need to twist halfway
      if ts != 0:
        matrices.append(translation(ts * p.length(c_in)/2)) # check sign
        types.append(f"start fix of {ts * p.length(c_in)/2} for {o} at {c_in} with {(ts, te)}.")
      # travel the ortholine
      matrices.append(rotation(-p.length(o)))
      types.append(f"ortho traversal of {-p.length(o)} along {o}")
      # if our ortho ends at a non-canonical point, we need to twist halfway
      if te != 0:
        matrices.append(translation(te * p.length(c_out)/2)) # check sign
        types.append(f"end fix of {te * p.length(c_out)/2} for {o} at {c_out} with {(ts, te)}.")
      # TODO can optimze out canonical twisting in special cases
    return matrices, types

  def matrix_from_path(self, path):
    """ Return the matrix correspodning to mapping canonical endpoint of path """
    """ to the point 0 + 1j (i.e. the basebpoint/starting canonical point) """
    matrices, _ = self._matrices_from_path(path)
    return functools.reduce(np.matmul, matrices[::-1])

  def generators(self):
    """ Returns matrices that generate the fundamental group. """
    """ There are two types of generators: """
    """    1. going around cuff: path_to * around * path_back """
    """    2. gluing ones: path to * glue twist * path_from_glued_back """
    rot = np.array([[0,1],[-1,0]])
    paths = self.paths() # source node is 0
    # by convention, we glue leaf nodes in the pairs: 0 and 1, 2 and 3, etc.
    path_pairs = [(paths[i], paths[i+1]) for i in range(0, len(paths), 2)]
    generators = []
    for p1, p2 in path_pairs:
        #print(p1, p2)
        m1 = self.matrix_from_path(p1)
        p = self._pants[p2[-2]] # last internal node gives last pant
        c = Cuff((p2[-1], p2[-2])) # we use p2 so that p[-2] exists
        l = p.length(c)
        # Note: the twist here is relative to gluing of canonical points.
        # In particular, even with t = 0 we could be gluing at corners if
        # the canonical points + orientations match up this way.
        # This should not be an issue for us since we care about all parameters
        # but it is worthwhile to keep this in mind. When generating paths
        # we had to be consistent since we cross the same cuff for multiple
        # paths. Here, however, we only cross a cuff for one generator, so
        # the twist parameter is simply relative to 0 or +/- l/2.
        t = p.twist(c)
        generators.append(inverse(m1) @ translation(l) @ m1)
        # we need to reverse p2 for generating m2 for gluing paths
        if len(p1) > 1:
          m2 = self.matrix_from_path(p2[::-1])
          # we need a 180 rotation since we return from same side
          generators.append(rot @ m2 @ translation(t) @ m1)
        else:
          m2 = self.matrix_from_path(p2)
          # m1 should be identity
          generators.append(translation(t) @ m2)

    return generators

  def approximate_length_spectrum(self, depth=10, eps=0.0001):
    gens = self.generators()
    # note, gens_and_inverses[k] and gens_and_inverses[-k] are inverses
    gens_and_inverses = gens + [inverse(g) for g in gens[::-1]]
    n = len(gens_and_inverses)
    lengths = []
    prev_layer = [(np.array([[1,0],[0,1]]), None)]
    prev_word_layer = [()]
    for _ in range(depth):
      new_layer = []
      new_word_layer = []
      for k, (h, j) in enumerate(prev_layer):
        for i, g in enumerate(gens_and_inverses):
          if j is None or 0 != (i + j + 1) % n:
            w = (i,) + prev_word_layer[k]
            z = g @ h
            new_layer.append((z, i))
            new_word_layer.append(w)
            x = translation_length(z, eps, word=w)
            if len(lengths) == 0 or min(abs(l - x) for l in lengths) > eps:
              lengths.append(x)
            if abs(x) < eps:
              print(w)
      prev_layer = new_layer
      prev_word_layer = new_word_layer
    return lengths

  def can_exclude(self, know_sys, depth=10, eps=0.001):
    gens = self.generators()
    # note, gens_and_inverses[k] and gens_and_inverses[-k] are inverses
    gens_and_inverses = gens + [inverse(g) for g in gens[::-1]]
    n = len(gens_and_inverses)
    lengths = []
    prev_layer = [(np.array([[1,0],[0,1]]), None)]
    prev_word_layer = [()]
    for _ in range(depth):
      new_layer = []
      new_word_layer = []
      for k, (h, j) in enumerate(prev_layer):
        for i, g in enumerate(gens_and_inverses):
          if j is None or 0 != (i + j + 1) % n:
            w = (i,) + prev_word_layer[k]
            z = g @ h
            new_layer.append((z, i))
            new_word_layer.append(w)
            x = translation_length(z, eps, word=w)
            if abs(x) > eps and abs(x) < know_sys:
              return True
      prev_layer = new_layer
      prev_word_layer = new_word_layer
    return lengths

