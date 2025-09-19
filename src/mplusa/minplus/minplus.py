import numpy as np

import math
import string

from .domain import validate_domain
from .. import utils



def add(*args : float) -> float:
    validate_domain(args)
    return min(args)


def mult(*args : float) -> float:
    validate_domain(args)
    return sum(args) if math.inf not in args else math.inf


def power(a : float,
          k : int) -> float:
    return mult(*[a for _ in range(k)])


def modulo(a : float,
           t : int) -> float:
    validate_domain([a, t])
    if a < 0 or t < 0:
        raise ValueError('The modulo operator is only defined for positive numbers.')
    if a == math.inf:
        return math.inf
    if a == 0:
        return 0
    if t == math.inf or t == 0:
        return a
    return a - (a // t) * t


def add_matrices(A : np.ndarray,
                 B : np.ndarray) -> np.ndarray:
    if A.shape != B.shape:
        raise ValueError('Given matrices have different shapes.')
    result = np.copy(A)
    shape = A.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            result[i, j] = add(result[i, j], B[i, j])
    return result


def mult_matrices(A : np.ndarray,
                  B : np.ndarray) -> np.ndarray:
    if A.shape[1] != B.shape[0]:
        raise ValueError('Given matrices are not of MxN and NxP shapes.')
    result = np.zeros((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            result[i, j] = add(*[mult(A[i, k], B[k, j]) for k in range(A.shape[1])])
    return result


def power_matrix(A : np.ndarray,
                 k : int) -> np.ndarray:
    if np.any(np.diagonal(A) != 0):
        raise ValueError('Matrix contains non-zero values on the diagonal.')
    if k == 0:
        result = unit_matrix(A.shape[0], A.shape[1])
    else:
        result = A.copy()
        for _ in range(k):
            result = mult_matrices(A, result)
    return result


def modulo_matrices(A : np.ndarray,
                    b : np.ndarray) -> np.ndarray:
    if b.shape[1] != 1:
        raise ValueError('Given matrix b is not a vertical vector of shape Mx1')
    if A.shape[0] != b.shape[0]:
        raise ValueError('Given matrix b does not have an Mx1 shape against the MxN matrix A.')
    if np.any(A < 0) or np.any(b < 0):
        raise ValueError('Given matrices contain negative values.')
    result = np.zeros(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            result[i, j] = modulo(A[i, j], b[i])
    return result


def unit_matrix(width : int,
                height : int) -> np.ndarray:
    result = np.eye(width, height)
    result[result == 0] = math.inf
    result[result == 1] = 0
    return result


def kleene_star(A : np.ndarray,
                iterations : int = 1000) -> np.ndarray:
    if A.shape[0] != A.shape[1]:
        raise ValueError('Matrix is not square.')
    series = [
        unit_matrix(A.shape[0], A.shape[1]),
        A.copy()
    ]
    for _ in range(iterations):
        series.append(add_matrices(series[-1], mult_matrices(series[-1], series[-2])))
    return series[-1]


class MultivariatePolynomial:
    """ An implementation of a tropical polynomial with multiple variables. """

    def __init__(self, coefficients : np.ndarray) -> None:
        validate_domain(coefficients)
        self.coefficients = coefficients
        self.dimensions = len(self.coefficients.shape) + 1
        self._symbols = string.ascii_lowercase

    def __call__(self, *variables : float) -> float:
        if len(variables) != self.dimensions - 1:
            raise ValueError('The amount of variables and coefficients differs.')
        result = [math.inf]
        for indices, coefficient in np.ndenumerate(self.coefficients):
            powers : list[float] = []
            for variable_index, i in enumerate(indices):
                powers.append(power(variables[variable_index], i))
            result.append(mult(coefficient, *powers))
        result = add(*result)
        return float(result)

    def __str__(self) -> str:
        result = ''
        for indices, coefficient in np.ndenumerate(self.coefficients):
            if coefficient.is_integer():
                result += '(' + str(int(coefficient))
            elif coefficient < math.inf:
                result += '(' + str(coefficient)
            else:
                result += '(∞'
            for variable_index, i in enumerate(indices):
                if i > 1:
                    result += ' * ' + self._symbols[variable_index] + '^' + str(i)
                elif i == 1:
                    result += ' * ' + self._symbols[variable_index]
            result += ') + '
        return result[:-3]

    def get_hyperplanes(self) -> list[list[float|int]]:
        """ Returns a list of coefficients of a linear equation for every hyperplane building the polynomial. """
        result = []
        for indices, coefficient in np.ndenumerate(self.coefficients):
            if coefficient == math.inf:
                continue
            hyperplane = [float(coefficient)]
            hyperplane.extend(indices)
            hyperplane.append(1)  # The coefficient of the last dimension (e.g. Z in 3D)
            result.append(
                list(
                    map(
                        lambda x: int(x) if x.is_integer() else float(x),
                        hyperplane
                    )
                )
            )
        return result


class Polynomial(MultivariatePolynomial):
    """ An implementation of a tropical polynomial with a single variable. """

    def __init__(self, *coefficients : float|int) -> None:
        validate_domain(coefficients)
        super().__init__(np.array(coefficients))

    def get_line_intersections(self) -> list[list[float|int]]:
        """ Returns a list of intersection points for the lines building the polynomial. """
        result = []
        lines = self.get_hyperplanes()  # Hyperplanes are lines in this case
        for line in lines:  # Change the form of the equation to a + bx from a + bx + cy
            _ = line.pop()
        lines = filter(lambda x: len(x) == 2, utils.powerset(lines))
        for line_1, line_2 in lines:
            point = [(line_2[0] - line_1[0]) / (line_1[1] - line_2[1])]
            point.append(line_1[0] + line_1[1] * point[0])
            result.append(tuple(point))
        result = list(filter(lambda point: round(point[1], 8) == round(self(point[0]), 8), result))  # Filter out the points not belonging to the polynomial
        return result

    def get_roots(self) -> tuple[list[float|int], list[int]]:
        """ Returns lists of roots of the polynomial and of their respective ranks (amount of monomials attaining the value). """
        result = {}
        points = self.get_line_intersections()
        for point in points:
            if not point[0] in result:
                result[point[0]] = 1
            else:
                result[point[0]] += 1
        return list(result.keys()), list(result.values())
