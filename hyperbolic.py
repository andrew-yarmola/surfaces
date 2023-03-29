""" Hyperbolic equations"""

import numpy as np

def rotation(d):
  """ Returns a Mobius matrix that translates distance d along (-1, 1) [from -1 to 1]."""
  R = np.zeros((2,2))
  R[0][0] = np.cosh(d/2)
  R[0][1] = np.sinh(d/2)
  R[1][0] = np.sinh(d/2)
  R[1][1] = np.cosh(d/2)
  return R

def translation(t):
  """ Returns a Mobius matrix that translates distance d along (0, infty) [from 0 to infty]."""
  L = np.zeros((2,2))
  L[0][0] = np.exp(t/2)
  L[0][1] = 0
  L[1][0] = 0
  L[1][1] = np.exp(-t/2)
  return L

def translation_length(M, eps=0.001, word=None):
  """ Returns translation length of a hyerpbolic element M. """
  t = M.trace()
  if np.real(t) < 0:
    t = -t
  if abs(t/2) < 1:
    if abs(t/2) > 1 - eps:
      t = 2
    else:
      print(f"Bad trace: {t} for word {word} and matrix.\n{M} with det {np.linalg.det(M)}")
  return  2 * np.arccosh(t / 2)

def inverse(m):
  """ Inverse of det 1 matrix """
  return np.array([[m[1,1], -m[0,1]],[-m[1,0], m[0,0]]])

def bers_constant_bound(g):
  """ See paper 'A Short Note on Short Pants' by Hugo Parlier."""
  if g == 2:
    return 4.68
  R = np.arccosh(1/(np.sqrt(2) * np.sin(np.pi/(12 * g - 6))))
  return 4 * np.pi * (g-1) + 4 * R

