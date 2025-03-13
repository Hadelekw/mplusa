import math
import numpy as np
from numbers import Real


def add(*args) -> Real:
    if -math.inf in args:
        raise ValueError(
            'Minplus.add: value out of domain.'
        )
    return min(args)


def mult(*args) -> Real:
    if -math.inf in args:
        raise ValueError(
            'Minplus.mult: value out of domain.'
        )
    return sum(args) if math.inf not in args else math.inf


def add_matrices(A : np.ndarray,
                 B : np.ndarray) -> np.ndarray:
    if A.shape != B.shape:
        raise ValueError(
            'Minplus.add_matrices: given matrices ' +\
            'are of different shape (A: {}, B: {}).'.format(A.shape, B.shape)
        )
    result = np.copy(A)
    shape = A.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            result[i, j] = add(result[i, j], B[i, j])
    return result


def mult_matrices(A : np.ndarray,
                  B : np.ndarray) -> np.ndarray:
    if A.shape[1] != B.shape[0]:
        raise ValueError(
            'Minplus.mult_matrices: given matrices ' +\
            'are of shapes not given as MxN and NxP (A: {}, B: {}).'.format(
                A.shape, B.shape
            )
        )
    result = np.zeros((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            result[i, j] = add(*[mult(A[i, k], B[k, j]) for k in range(A.shape[1])])
    return result


def modulo(a : Real,
           t : int) -> Real:
    if a < 0 or t < 0:
        raise ValueError(
            'Minplus.modulo: modulo operation is only defined for positive numbers.'
        )
    if a == math.inf:
        return math.inf
    if a == 0:
        return 0
    if t == math.inf or t == 0:
        return a
    return a - (a // t) * t


def modulo_matrices(A : np.ndarray,
                    b : np.ndarray) -> np.ndarray:
    if b.shape[1] != 1:
        raise ValueError(
            'Minplus.modulo_matrices: given matrix b ' +\
            'is not a vertical vector of shape Mx1 (has shape of {}).'.format(
                b.shape
            )
        )
    if A.shape[0] != b.shape[0]:
        raise ValueError(
            'Minplus.modulo_matrices: given matrix b ' +\
            'does not have an Mx1 shape against MxN matrix A (A: {}, b: {}).'.format(
                A.shape, b.shape
            )
        )
    if np.any(A < 0) or np.any(b < 0):
        raise ValueError(
            'Minplus.modulo_matrices: matrices contain negative values.'
        )
    result = np.zeros(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            result[i, j] = modulo(A[i, j], b[i])
    return result


def power(a : Real,
          k : int) -> Real:
    return mult(*[a for _ in range(k)])


def power_matrix(A : np.ndarray,
                 k : int) -> np.ndarray:
    if np.any(np.diagonal(A) != 0):
        raise ValueError(
            'Minplus.power_matrix: matrix contains non-zero values on the diagonal.'
        )
    if k == 0:
        result = unit_matrix(A.shape[0], A.shape[1])
    else:
        result = A.copy()
        for _ in range(k):
            result = mult_matrices(A, result)
    return result


def unit_matrix(width : int,
                height : int) -> np.ndarray:
    if width < 0 or height < 0:
        raise ValueError(
            'Minplus.unit_matrix: invalid width or height.'
        )
    result = np.eye(width, height)
    result[result == 0] = math.inf
    result[result == 1] = 0
    return result


def star(A : np.ndarray,
         iterations : int = 1000,
         eps : float = 0.001) -> np.ndarray:
    if A.shape[0] != A.shape[1]:
        raise ValueError(
            'Minplus.star: matrix is not square.'
        )
    series = [
        unit_matrix(A.shape[0], A.shape[1]),
        A.copy()
    ]
    for i in range(2, iterations):
        series.append(add_matrices(series[-1], series[-2]))
        # Very basic check if the series is convergent.
        if abs(np.max(series[-1] - series[-2])) < eps:
            break
    else:
        raise ValueError(
            'Minplus.star: the series for this matrix is not convergent ' +\
            '(within the limits of iterations and decimal places).'
        )
    return series[-1]


class Polynomial:
    """ A simple implementation of a single-variable tropical polynomial. """

    def __init__(self, *coefficients) -> None:
        for value in coefficients:
            if not isinstance(value, Real) or value == -math.inf:
                raise ValueError(
                    'Minplus.Polynomial.__init__: coefficient value out of domain.'
                )
        self.coefficients = coefficients[::-1]

    def __call__(self, x : float) -> float:
        return add(*[mult(coefficient, power(x, i)) for i, coefficient in enumerate(self.coefficients)])

    def get_hypersurface(self) -> list[float]:
        result = []
        candidates = [c2 - c1 for c1, c2 in zip(self.coefficients[1:], self.coefficients[:-1])]
        for candidate in candidates:
            if not isinstance(candidate, Real) or candidate == -math.inf:
                continue
            if abs(self(candidate) - self(candidate - 1)) != abs(self(candidate) - self(candidate + 1)):
                result.append(candidate)
        return result
