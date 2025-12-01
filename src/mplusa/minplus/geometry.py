import numpy as np

import math
from collections.abc import Collection

from .domain import validate_domain
from .minplus import add, mult, add_matrices


def point_type(point : Collection, vertices : Collection, indexing_start : int = 1) -> Collection:
    result = [[] for _ in range(len(vertices))]
    point = np.array(point); vertices = np.array(vertices)
    for i, vertex in enumerate(vertices):
        comparison = add(*(vertex - point))
        for j in range(len(vertex)):
            if vertex[j] - point[j] == comparison:
                result[j].append(indexing_start + i)
    return result


def hyperplane_from_apex(point : Collection[float|int]) -> Hyperplane:
    return Hyperplane(*[-coordinate for coordinate in point])


class Cone:

    def __init__(self, *vectors : Collection[float|int]) -> None:
        validate_domain(vectors)
        self.vectors = np.array(vectors)  # generators of the cone
        self.vector_count = self.vectors.shape[0]
        self.dimensions = self.vectors.shape[1]

    def get_point(self, constants : Collection[float|int]) -> np.ndarray:
        if len(constants) != self.vector_count:
            raise ValueError('The number of shift values not equal to the vector count.')
        result = np.array([math.inf for _ in range(self.dimensions)])
        for constant, vector in zip(constants, self.vectors):
            result = add_matrices(result, vector + constant)  # Addition here is tropical multiplication
        return result

    def sample_points(self, constants_groups : Collection[np.ndarray]) -> np.ndarray:
        result = []
        grid = np.meshgrid(*constants_groups)
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
        if len(point) != self.dimension:
            raise ValueError('The amount of the point\'s coordinates and the coefficients differs.')
        result = add(*[mult(c, p) for c, p in zip(self.coefficients, point)])
        return result

    def get_apex(self) -> tuple:
        return tuple([-coefficient for coefficient in self.coefficients])


class Polytope:
    """ An implementation of a tropical polytope structure. """

    def __init__(self, *faces_collection : Collection, adjust_coordinate_dimensions : bool = True) -> None:
        self.dimension = len(faces_collection) - 1
        self.structure = {}
        for rank, faces in enumerate(faces_collection):
            if not rank:
                if adjust_coordinate_dimensions:
                    faces = [(0, *face) for face in faces]
                self._validate_vertices_dimension(list(faces))
                validate_domain(faces)
            self.structure[rank] = np.array(faces)

    def _validate_vertices_dimension(self, vertices : list[Collection]) -> None:
        pivot = len(vertices[0])
        if pivot not in [self.dimension, self.dimension + 1]:
            raise ValueError('Incorrect dimension of the vertices\' coordinates.')
        for vertex in vertices:
            if len(vertex) != pivot:
                raise ValueError('Incorrect dimension of the vertices\' coordinates.')

    @property
    def vertices(self) -> np.ndarray:
        return self.structure[0]

    @property
    def pseudovertices(self) -> np.ndarray:
        result = []
        for line_segment in self.get_all_line_segments():
            for point in line_segment:
                if not np.any(np.all(self.vertices == point, axis=1)):  # if point not in vertices
                    result.append(point)
        return np.array(result)

    @property
    def edges(self) -> list:
        result = []
        for vertices_identifiers in self.structure[1]:
            result.append([])
            for identifier in vertices_identifiers:
                result[-1].append(self.vertices[identifier])
        return result

    def get_line_segments(self, edge_index : int) -> np.ndarray:
        """
        Returns a list of points belonging to a line segment of the polytope.
        It is not necessarily sorted.
        """
        result = []
        x, y = self.edges[edge_index]
        for i in range(len(x)):
            if y[i] >= x[i]:
                result.append(add_matrices((y[i] - x[i]) + x, y))
            else:
                result.append(add_matrices((x[i] - y[i]) + y, x))
        return np.array(sorted(result, key=lambda x: x[-1] - x[-2]))

    def get_all_line_segments(self) -> list:
        result = []
        for i in range(len(self.structure[1])):
            result.append(self.get_line_segments(i))
        return result

    def get_envelope(self) -> list:
        pass

    def get_apices(self) -> np.ndarray:
        """ Returns an array of apices corresponding to the half-spaces defined by the facet-defining hyperplanes. """
        return self.pseudovertices

    def get_hyperplanes(self) -> list:
        """ Returns a list of the facet-defining hyperplanes. """
        result = []
        for apex in self.get_apices():
            result.append(hyperplane_from_apex(apex))
        return result
