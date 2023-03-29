""" Bolza Surface Test"""

from prufer import *
from surface import *

import numpy as np

# Bolza Surface
# According to: https://our.unc.edu/wp-content/uploads/sites/1148/2022/02/Systoles-and-Parameterization-.pdf, the FN coords are given by 
#  ((l_1, theta_1), (l_2, theta_2), (l_3, theta_3)) = ((a, 0), (a, 0), (b, b/2))
# where a = 2arccosh(1+sqrt(2)) = sys(\mathbb{M}_2), b = 2arccosh(3+2sqrt2)

unique = unique_surface_prufer_codes(2)
bolza = Surface(unique[1])
multigraph_plot(bolza._tree)

### From the previous, see yellow and red highlighting, it looks like the pants-decomp they are adding these FN coords onto is unique[1]
### The leaf nodes must have the a-lengths, while the core connecting the pants must have the 

s = 2 * np.arccosh(1 + np.sqrt(2))
t = 2 * np.arccosh(3 + 2 * np.sqrt(2))
# twist gluings are the (0,1), (2,3), and (4,5) edges
#fn_coords = [(s, 0),(s, 0),(t, t/2)]
fn_coords = [(s, 0),(s, 0),(t, t/2)]

bolza.set_coords(fn_coords)

print("Pants Decompositon")
for n, p in bolza._pants.items():
  print(n, p)

print("Paths")
for i, path in bolza.paths().items():
  print(f"Path to {i}: {path}")
  for m, tp in zip(*bolza._matrices_from_path(path)):
    print(f"Type of matrix {tp}")
    print(f"{m}")


print(f"Path from 3 to 0: {[3, 5, 4, 0]}")
for m, tp in zip(*bolza._matrices_from_path([3, 5, 4, 0])):
  print(f"Type of matrix {tp}")
  print(f"{m}")


def special_m(n):
  vals = [abs((2*k+1) - n * np.sqrt(2)) for k in range(0, n)]
  k = vals.index(min(vals))
  return 2 * k + 1

spec = [2 * np.arccosh(special_m(n) + n * np.sqrt(2)) for n in range(1, 11)]

print(f"Bolza length spectrum: {spec}.")

print("Generators")
for m in bolza.generators():
  print(m, np.linalg.det(m))
  print(f"has translation length {translation_length(m)}")

lengths = bolza.approximate_length_spectrum(depth=4)
print(sorted(lengths))

""" Random test"""

g = 2
bers = bers_constant_bound(g)
unique = unique_surface_prufer_codes(g)
min_len=2.910678338574799
# for g = 3 min_len = 2.2856287974516105
for code in unique:
  s = Surface(code)
  multigraph_plot(s._tree)
  for _ in range(10):
    s.set_random_coords(min_len=min_len, max_len=bers)
    # print("Pants Decompositon")
    # for n, p in s._pants.items():
    #  print(n, p)
    # print("Paths")
    # for i, path in s.paths().items():
    #  print(f"Path to {i}: {path}")
    #  for m, tp in zip(*s._matrices_from_path(path)):
    #    print(f"Type of matrix {tp}")
    #    print(f"{m} with det {np.linalg.det(m)}")
    # print("Generators")
    # for m in s.generators():
    #  print(m, np.linalg.det(m))
    #  print(f"has translation length {translation_length(m)}")
    lengths = sorted(s.approximate_length_spectrum(depth=3))
    min_len = max(min_len, lengths[0])
    print(min_len, lengths[:5])
print(min_len)

