# MPlusA
**MPlusA** is a Python library for calculations in tropical algebra (also known as (min, +) and (max, +) algebra). For the full list of the library's capabilities refer to the [documentation](https://github.com/Hadelekw/mplusa_docs).

# Contributing
The library is open for all contributions. If you want to contribute to its development, please refer to the guidelines to verify the code standards before making a pull request. Developments from all areas of tropical mathematics are welcome.

# Installation
The easiest way to install the library is to use pip.

``` pip install mplusa ```

The only automatically installed dependency is [NumPy](https://numpy.org). Otherwise the library does make use of [Matplotlib](https://matplotlib.org/) for visualisations but it is not a required dependency and needs to be installed separately in case the user wants to use these capabilities.

# Examples
For the full list of the library's capabilities, refer to the documentation [documentation](https://github.com/Hadelekw/mplusa_docs).
This section is only meant to highlight the core features of the library.

## Basic operations
```Python
 # import maxplus instead of minplus to use the (max, +) semiring
from mplusa import minplus as mpa
import numpy as np

print(mpa.add(3, 5, 89))  # -> 3
print(mpa.mult(4, 8, 9))  # -> 21
print(mpa.power(10, 89))  # -> 890
print(mpa.modulo(10, 3))  # -> 1

A = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
])
B = np.array([
    [9, 8, 7],
    [6, 5, 4],
    [3, 2, 1],
])

print(mpa.add_matrices(A, B))  # -> [[1, 2, 3], [4, 5, 4], [3, 2, 1]]
print(mpa.mult_matrices(A, B))  # -> [[6, 5, 4], [9, 8, 7], [12, 11, 10]]
print(mpa.power_matrix(A, 3))  # -> [[3, 4, 5], [6, 7, 8], [9, 10, 11]]
print(mpa.modulo_matrices(A, np.array([[1, 2, 3]]).T))  # -> [[0, 0, 0], [0, 1, 0], [1, 2, 0]]
```
## Matrix-specific operations
```Python
 # import maxplus instead of minplus to use the (max, +) semiring
from mplusa import minplus as mpa
import numpy as np

A = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
])

print(mpa.kleene_star(A))  # -> [[0, 2, 3], [4, 0, 6], [7, 8, 0]]
print(mpa.kleene_star_series(A))  # -> all the matrices generated during kleene star
print(mpa.kleene_plus(A))  # -> [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
print(mpa.kleene_plus_series(A))  # -> all the matrices generated during kleene plus

print(mpa.tdet(A))  # -> 15
print(mpa.is_matrix_singular(A))  # -> True
print(mpa.matrix_rank_tropical(A))  # -> 1

# This function takes an argument func which accepts the following methods: 
# 1. mpa.power_algorithm (the default value)
# 2. mpa.karp_algorihtm
print(mpa.eigenvalue(A))  # -> 1
# This function takes an argument func which accepts the following methods:
# 1. mpa.power_algorithm (the default value)
# 2. mpa.kleene_star_algorithm
# 3. mpa.modified_kleene_star_algorithm
print(mpa.eigenvector(A))  # -> [[0], [3], [6]]
```
## Polynomials
```Python
 # import maxplus instead of minplus to use the (max, +) semiring
from mplusa import minplus as mpa
import numpy as np

p = mpa.MultivariatePolynomial(np.array([
    [1, np.inf, 15],
    [5, 0, np.inf],
    [5, 4, 3],
]))

print(p)  # -> (1) + (∞ * b) + (15 * b^2) + (5 * a) + (0 * a * b) + (∞ * a * b^2) + (5 * a^2) + (4 * a^2 * b) + (3 * a^2 * b^2)
print(p(8, 5))  # -> 1
print(p.get_linear_hyperplanes())  # -> [[1, 0, 0, 1], [15, 0, 2, 1], [5, 1, 0, 1], [0, 1, 1, 1], [5, 2, 0, 1], [4, 2, 1, 1], [3, 2, 2, 1]]

# A polynomial with one variable
p = mpa.Polynomial(5, 8, 9)

print(p)  # -> (5) + (8 * a) + (9 * a^2)
print(p(8))  # -> 5
print(p.get_linear_hyperplanes())  # -> [[5, 0, 1], [8, 1, 1], [9, 2, 1]]
print(p.get_line_intersections())  # -> [(-2.0, 5.0)]
print(p.get_roots())  # -> ([-2.0], [1]) where the first value is the root and the latter is its corresponding rank
```
## Geometry and visualizations
```Python
 # import maxplus instead of minplus to use the (max, +) semiring
from mplusa import minplus as mpa
import numpy as np
# For visualizations matplotlib is required and is imported inside of the submodule
# By default drawing polytopes and hyperplanes won't display the figure and that is
# why it is imported here; it's not necessary if a figure is displayed using show=True
import matplotlib.pyplot as plt

# NOTE: the library uses the first element in a point for projections; it should usually be 0
p = mpa.geometry.Polytope2D((0, 1, 2), (0, 3, 3), (0, 4, 5))

print(p.vertices)  # -> [[0, 1, 2], [0, 3, 3], [0, 4, 5]]
print(p.pseudovertices)  # -> [[0, 2, 3], [0, 4, 4]]
print(p.edges)  # -> [[[0, 1, 2], [0, 3, 3]], [[0, 3, 3], [0, 4, 5]], [[0, 4, 5], [0, 1, 2]]]

print(p.get_all_line_segments())  # -> point by point edges (including pseudovertices)

mpa.visualize.hasse_diagram(p)
mpa.visualize.draw_polytope2D(p, show=True)

for hyperplane in p.get_hyperplanes():
    mpa.visualize.draw_hyperplane2d(hyperplane)
plt.show()
```
