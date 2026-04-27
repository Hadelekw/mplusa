import numpy as np

from typing import Callable
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


def power(a : float, k : int) -> float:
    return mult(*[a for _ in range(k)])


def modulo(a : float, t : int) -> float:
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


def add_matrices(A : np.ndarray, B : np.ndarray) -> np.ndarray:
    validate_domain([A, B])
    return np.minimum(A, B)


def mult_matrices(A : np.ndarray, B : np.ndarray) -> np.ndarray:
    if A.shape[1] != B.shape[0]:
        raise ValueError('Given matrices are not of MxN and NxP shapes.')
    result = np.zeros((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            result[i, j] = add(*[mult(A[i, k], B[k, j]) for k in range(A.shape[1])])
    return result


def power_matrix(A : np.ndarray, k : int) -> np.ndarray:
    if k == 0:
        result = unit_matrix(A.shape[0], A.shape[1])
    else:
        result = A.copy()
        for _ in range(k):
            result = mult_matrices(A, result)
    return result


def modulo_matrices(A : np.ndarray, b : np.ndarray) -> np.ndarray:
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


def mult_arrays(A : np.ndarray, B : np.ndarray) -> np.ndarray:
    """
    Performs tropical tensor multiplication of NumPy arrays of any shape.
    The operation is defined as:

    A_{i_0, ..., i_k} * B_{j_0, ..., j_k} = C_{i_0, ..., i_k, j_0, ..., j_k}

    where each element of C is calculated by tropical multiplication of the values in the arrays.
    """
    result = np.zeros((*A.shape, *B.shape))
    for i, value in np.ndenumerate(A):
        for j, other_value in np.ndenumerate(B):
            result[*i, *j] = mult(value, other_value)
    return result


def unit_matrix(width : int, height : int) -> np.ndarray:
    result = np.eye(width, height)
    result[result == 0] = math.inf
    result[result == 1] = 0
    return result


def kleene_star(A : np.ndarray, iterations : int|None = None) -> np.ndarray:
    if A.shape[0] != A.shape[1]:
        raise ValueError('Matrix is not square.')
    series = [
        unit_matrix(A.shape[0], A.shape[1]),
        A.copy()
    ]
    result = add_matrices(series[0], series[1])
    if iterations is None:
        iterations = int(result.shape[0])
    for i in range(iterations):
        series.append(power_matrix(A, i))
        result = add_matrices(result, series[-1])
    return result


def kleene_plus(A : np.ndarray, iterations : int|None = None) -> np.ndarray:
    if A.shape[0] != A.shape[1]:
        raise ValueError('Matrix is not square.')
    series = [A.copy()]
    result = series[0]
    if iterations is None:
        iterations = int(result.shape[0])
    for i in range(1, iterations):
        series.append(power_matrix(A, i))
        result = add_matrices(result, series[-1])
    return result


def tdet(A : np.ndarray) -> float:
    used = np.zeros(A.shape[0])
    best = math.inf

    def backtrack(row : int, current_sum : float) -> None|float:
        nonlocal best
        if row == A.shape[0]:
            best = add(best, current_sum)
            return
        if current_sum >= best:
            return
        for column in range(A.shape[1]):
            if not used[column]:
                used[column] = 1
                backtrack(row + 1, mult(current_sum, A[row][column]))
                used[column] = 0

    backtrack(0, 0)
    return best


# TODO:
# def is_matrix_singular
# def matrix_rank_barvinok  <- unfortunately an NP-hard problem; approximate?
# def matrix_rank_kapranov  <- requires symbolic calculations; decide what to do
# def matrix_rank_tropical


def power_algorithm(A : np.ndarray, x_0 : np.ndarray|None = None, iterations : int = 1000) -> tuple:
    if x_0 is None:
        x_0 = np.ones((A.shape[1], 1))
    xs = [x_0]
    for i in range(iterations):
        xs.append(mult_matrices(A, xs[i]))
        for j in range(len(xs) - 1):
            if len(np.unique(xs[-1] - xs[j])) == 1:
                p = len(xs) - 1
                q = j
                c = float(np.unique(xs[-1] - xs[j])[0])
                return p, q, c, xs
    raise ValueError(f'Unable to find the values using the power algorithm within {iterations} iterations.')


def karp_algorithm(A : np.ndarray) -> float:
    """ An implementation of Karp's algorithm for finding the minimum cycle mean. """
    n = A.shape[0]
    dp = np.zeros((n + 1, n))
    for i in range(1, n + 1):
        for j in range(n):
            dp[i, j] = math.inf
            for k in range(n):
                if A[k, j] != math.inf:
                    dp[i, j] = add(dp[i, j], mult(dp[i - 1, k], A[k, j]))
    result = math.inf
    for i in range(n):
        for j in range(n):
            if dp[n, i] < math.inf and dp[j, i] < math.inf:
                if n != j:
                    avg = (dp[n, i] - dp[j, i]) / (n - j)
                    result = add(result, avg)
    return result


def kleene_star_algorithm(A : np.ndarray, column : int = 0) -> np.ndarray:
    """ Returns a non-normalized eigenvector of a tropical matrix. """
    l = karp_algorithm(A)  # the minimum tropical cycle mean -- assumed eigenvector
    B = np.add(A, -l)
    B_star = kleene_star(B)
    return B_star[:, [column]]


def modified_kleene_star_algorithm(A : np.ndarray) -> np.ndarray:
    """ A modified version of the Kleene star algorithm which uses Kleene plus to locate the critical column. """
    l = karp_algorithm(A)  # the minimum tropical cycle mean -- assumed eigenvector
    B = np.add(A, -l)
    B_star = kleene_star(B)
    B_plus_diagonal = np.diagonal(kleene_plus(B))
    return B_star[:, np.where(B_plus_diagonal == 0)[0]]


def eigenvalue(A : np.ndarray, func : Callable = power_algorithm) -> float:
    """ Calculate the eigenvalue of a matrix using one of the provided methods. """
    if func.__name__ == 'power_algorithm':
        p, q, c, _ = func(A)
        return c / (p - q)
    elif func.__name__ == 'karp_algorithm':
        return func(A)
    raise NotImplementedError('The provided function is not supported. Use one of the: power_algorithm, karp_algorihtm.')


def normalize_vector(v : np.ndarray) -> np.ndarray:
    """ Normalize a vector by substracting the minimal value. """
    validate_domain(v)
    if len(v.shape) == 1:
        return v - add(*v)
    elif len(v.shape) == 2:  # check for vertical vectors -- preferred in many calculations
        if v.shape[1] == 1:
            return v - add(*v[:, 0])
    raise ValueError('Not a valid vector.')


def eigenvector(A : np.ndarray, func : Callable = power_algorithm, func_kwargs : dict = {}, normalize : bool = True) -> np.ndarray:
    """ Calculate the eigenvector of a matrix using one of the provided methods. """
    if func.__name__ == 'power_algorithm':
        p, q, c, xs = power_algorithm(A, **func_kwargs)
        eigenvalue = c / (p - q)
        result = np.ones_like(xs[0]) * math.inf
        for i in range(1, p - q + 1):
            result = add_matrices(result, power(eigenvalue, (p - q - i)) + xs[q + i - 1])
    elif func.__name__ in ['kleene_star_algorithm', 'modified_kleene_star_algorithm']:
        result = func(A, **func_kwargs)
    else:
        raise NotImplementedError('The provided function is not supported. Use one of the: power_algorithm, kleene_star_algorithm, modified_kleene_star_algorithm.')
    if normalize:
        return normalize_vector(result)
    return result


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

    def get_linear_hyperplanes(self) -> list[list[float|int]]:
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
        lines = self.get_linear_hyperplanes()  # Hyperplanes are lines in this case
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
        """ Returns lists of roots of the polynomial and of their respective ranks (the amount of monomials attaining the value). """
        result = {}
        points = self.get_line_intersections()
        for point in points:
            if point[0] not in result:
                result[point[0]] = 1
            else:
                result[point[0]] += 1
        return list(result.keys()), list(result.values())
