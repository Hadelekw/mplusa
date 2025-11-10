import numpy as np

import math
from collections.abc import Collection

from .domain import validate_domain
from .minplus import add_matrices


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


class Polytope:

    def __init__(self, *faces_collection : Collection, adjust_coordinate_dimensions : bool = True) -> None:
        self.structure = {}
        for rank, faces in enumerate(faces_collection):
            if not rank:
                if adjust_coordinate_dimensions:
                    faces = [(0, *face) for face in faces]
                validate_domain(faces)
            self.structure[rank] = faces
        self._lift_vertices()
        self.dimension = len(self.structure) - 1

    def _lift_vertices(self) -> None:
        pass

    @property
    def vertices(self) -> list[tuple]:
        return self.structure[0]

    @property
    def edges(self) -> Collection:
        result = []
        for vertices_identifiers in self.structure[1]:
            result.append([])
            for identifier in vertices_identifiers:
                result[-1].append(self.vertices[identifier])
        return result

    def get_hyperplanes(self) -> list:
        """ Returns a list of coefficients of the facet-defining hyperplanes. """
        pass

    def get_apices(self) -> list:
        """ Returns a list of apices corresponding to the half-spaces defined by the facet-defining hyperplanes. """
        pass
