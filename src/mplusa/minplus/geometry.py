import numpy as np

import math
from collections.abc import Collection

from .domain import validate_domain
from .minplus import add, mult, add_matrices


def project_point(point : tuple) -> tuple:
    validate_domain(point)
    return tuple((point[i] - point[0] for i in range(1, len(point))))


def point_type(point : Collection, vertices : Collection, indexing_start : int = 0) -> Collection:
    result = [[] for _ in range(len(vertices))]
    point = np.array(point)
    vertices = np.array(vertices)
    for i, vertex in enumerate(vertices):
        comparison = add(*(vertex - point))
        for j in range(len(vertex)):
            if vertex[j] - point[j] == comparison:
                result[j].append(indexing_start + i)
    return result


def tconv(points : Collection[Collection]) -> list:
    result = []
    for point in points:
        T = point_type(point, list(filter(lambda v: v != point, points)))
        if [] in T:
            result.append(point)
    return result


def line_segment(start_point : Collection, end_point : Collection, sort=True, unique=True) -> np.ndarray:
    x = np.array(start_point)
    y = np.array(end_point)
    result = []
    for i in range(len(x)):
        if y[i] >= x[i]:
            result.append(add_matrices((y[i] - x[i]) + x, y))
        else:
            result.append(add_matrices((x[i] - y[i]) + y, x))
    result = list(map(tuple, result))
    if unique:
        result = list(set(result))
    if sort:
        result = sort_line_segment(result, tuple(start_point))
    return np.array(result)


def sort_line_segment(points : list[tuple], start_point : tuple) -> list:
    points.remove(start_point)
    stack = [start_point]
    while points:
        # The line segments in tropical geometry have limited slopes and are always convex
        # Therefore the closest point (in Euclidean sense) is the next point in order
        distances = {
            math.sqrt(
                sum(
                    (p[i] - stack[-1][i])**2 for i in range(len(start_point))
                )
            ): p for p in points
        }
        stack.append(distances[min(distances.keys())])
        points.remove(stack[-1])
    return stack


class ConvexCone:

    def __init__(self, *vectors : Collection[float|int]) -> None:
        validate_domain(vectors)
        self.vectors = np.array(vectors)  # generators of the cone
        self.vector_count = self.vectors.shape[0]
        self.dimensions = self.vectors.shape[1]

    def get_point(self, constants : Collection[float|int]) -> np.ndarray:
        if len(constants) != self.vector_count:
            raise ValueError('The number of constants not equal the number of generators of the cone.')
        result = np.array([math.inf for _ in range(self.dimensions)])
        for constant, vector in zip(constants, self.vectors):
            result = add_matrices(result, vector + constant)  # Addition here is tropical multiplication
        return result

    def sample_points(self, constants_collection : Collection[np.ndarray]) -> np.ndarray:
        result = []
        grid = np.meshgrid(*constants_collection)
        for constants in np.nditer(grid):
            result.append(self.get_point([float(s) for s in constants]))
        result = np.array(result)
        return result


class Hyperplane:
    """ An implementation of a tropical hyperplane structure. """

    def __init__(self, *coefficients : float) -> None:
        validate_domain(coefficients)
        self.coefficients = coefficients
        self.dimension = len(coefficients)

    def get_value(self, point : Collection[float]) -> float:
        """ Calculate the value which the hyperplane achieves at a given point. """
        if len(point) != self.dimension:
            raise ValueError('The amount of the point\'s coordinates and the coefficients differs.')
        result = add(*[mult(c, p) for c, p in zip(self.coefficients, point)])
        return result

    def get_apex(self) -> tuple:
        """ Returns the coordinates of the point that serves as the apex of the hyperplane. """
        return tuple([-coefficient for coefficient in self.coefficients])


def hyperplane_from_apex(point : Collection[float|int]) -> Hyperplane:
    return Hyperplane(*[-coordinate for coordinate in point])


class AbstractPolytope:
    """ An implementation of a tropical polytope structure. """

    def __init__(self, faces_collection : list, normalize_points : bool = False) -> None:
        if normalize_points:
            faces_collection[0] = [(0, *face) for face in faces_collection[0]]
        self.dimension = len(faces_collection[0][0]) - 1  # The dimension is defined by the coordinates of the first point
        self.structure = {}
        for rank, faces in enumerate(faces_collection):
            validate_domain(faces)
            self.structure[rank] = np.array(faces)
        self._validate_vertices_dimension()

    def _validate_vertices_dimension(self) -> None:
        for vertex in self.structure[0]:
            if len(vertex) - 1 != self.dimension:
                raise ValueError('Incorrect dimension of the vertices\' coordinates.')

    @property
    def vertices(self) -> np.ndarray:
        """ Returns all points generating the tropical polytope. """
        return self.structure[0]

    @property
    def pseudovertices(self) -> np.ndarray:
        """ Returns all points that are generated by the tropical polytope. """
        result = []
        for line_segment in self.get_all_line_segments():
            for point in line_segment:
                if not np.any(np.all(self.vertices == point, axis=1)):  # if point not in vertices
                    result.append(point)
        return np.array(result)

    @property
    def edges(self) -> list:
        """ Returns all the edges (given as a pair of points) which build the tropical polytope. """
        result = []
        for vertices_identifiers in self.structure[1]:
            result.append([])
            for identifier in vertices_identifiers:
                result[-1].append(self.vertices[identifier])
        return result

    def get_line_segments(self, edge_index : int, sort : bool = True) -> np.ndarray:
        """ Returns a list of points belonging to a line segment of the polytope. """
        start_point, end_point = self.edges[edge_index]
        return np.array(line_segment(start_point, end_point, sort=sort))

    def get_all_line_segments(self) -> list:
        """ Returns all line segments building the polytope. """
        result = []
        for i in range(len(self.structure[1])):
            result.append(self.get_line_segments(i, sort=True))
        return result

    def _verify_facet_definition(self, point : Collection) -> bool:
        """ Verify if a point is facet-defining in the polynomial. """
        T = point_type(point, self.vertices)
        union = []
        for t in T:
            union.extend(t)
        union = list(set(union))
        if len(union) >= self.dimension:
            return True
        else:
            return False

    def get_apices(self) -> list:
        """ Returns an array of apices corresponding to the polynomial-defining half-spaces. """
        result = []
        for pseudovertex in self.pseudovertices:
            if self._verify_facet_definition(pseudovertex):
                result.append(pseudovertex)
        return result

    def get_hyperplanes(self) -> list:
        """ Returns a list of the facet-defining hyperplanes. """
        result = []
        for apex in self.get_apices():
            result.append(hyperplane_from_apex(apex))
        return result


class Polytope2D(AbstractPolytope):

    def __init__(self, *points : Collection) -> None:
        self.dimension = len(points[0]) - 1
        faces_collection = [list(points)]
        edges = []
        for i in range(1, len(points)):
            edges.append((i - 1, i))
        edges.append((len(points) - 1, 0))
        faces_collection.append(edges)
        faces_collection.append([tuple(range(len(edges)))])
        super().__init__(faces_collection)
