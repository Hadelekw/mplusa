# MPlusA
---
**MPlusA** is a small Python library for tropical algebra (also known as (min, +) and (max, +) algebra). It provides the definitions of basic operations on numbers and NumPy arrays, as well as a basic implementation of tropical polynomials.

Any improvements or fixes are always welcome.

## How to use
After having installed the library one can import one of the two modules the package consists of (`minplus` and `maxplus`) and use the full array of its capabilities. The functions are essentially the same between the modules.

**`add(*args) -> Real`**
	Tropical addition. Essentially an alias for Python's `min` function.

**`mult(*args) -> Real`**
	Tropical multiplication. Essentially an alias for Python's `sum` function.

**`add_matrices(A : np.ndarray, B : np.ndarray) -> np.ndarray`**
	Tropical addition of NumPy arrays. The summed matrices have to be of the same shape.

**`mult_matrices(A : np.ndarray, B : np.ndarray) -> np.ndarray`**
	Tropical multiplication of NumPy arrays. The multiplied matrices have to be of sizes MxN and NxP and their order matters. The result is of shape MxP.

**`modulo(a : Real, t : int) -> Real`**
	Tropical modulo operator. It can be understood as the difference between the number $a$ and $t^k$ where $k$ is the largest integer that satisfies $a \geq t^k$.

**`modulo_matrices(A : np.ndarray, b : np.ndarray) -> np.ndarray`**
	Tropical modulo operator for NumPy arrays. The input matrices should be of size MxN and Mx1. The result is an MxN matrix.

**`power(a : real, k : int) -> Real`**
	Tropical power operator.  Applies the multiplication k times.

**`power_matrix(A : np.ndarray, k : int) -> np.ndarray`**
	Tropical power operator for NumPy arrays. It multiplies the matrix k times.

**`unit_matrix(width : int, height : int) -> np.ndarray`**
	Creates a tropical unit matrix of given width and height.

**`star(A : np.ndarray) -> np.ndarray`**
	Definition of an unary operator of unique to tropical algebra, usually denoted as $\mathbf{A}^*$. It returns the value to which an infinite recursive sum of matrices converges. The input matrix has to be square and the series created in the process of calculating the value needs to be convergent.

**`Polynomial(*coefficients)`**
	This is a class that implements basic single-variable tropical polynomials. Calling an object of this class allows to take a value the polynomial takes at a given point.

For the full list of capabilities, refer to the [documentation](https://hadelekw.github.io/mplusa-docs.html)

## Example code
```
import numpy as np
from mplusa import minplus

# Basic operators
s = minplus.add(10, 6, 4, 13)  # -> 4
p = minplus.mult(10, 5, 8)  # -> 23
mod = minplus.modulo(p, s)  # -> 3

# NumPy arrays
A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
B = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])
S = minplus.add_matrices(A, B)  # -> [[1, 2, 3], [4, 5, 4], [3, 2, 1]]
P = minplus.mult_matrices(A, B)  # -> [[6, 5, 4], [9, 8, 7], [12, 11, 10]]
```

## Bibliography
1. A. Obuchowicz, K. A. D'Souza, Z. A. Banaszak. *An Algebraic Approach to Modelling and Performance Optimisation of a Traffic Route System*. International Journal of Applied Mathematics and Computer Science, vol. 8, no. 2, pp. 335-365, June 1998.
2. D. Speyer, B. Sturmfels. "*Tropical Mathematics*". Mathematics Magazine, vol. 82, no. 3, pp. 163-173, June 2009, doi: https://doi.org/10.4169%2F193009809x468760.
