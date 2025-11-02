import numpy as np

import math
from collections.abc import Collection

from .domain import validate_domain
from .minplus import add_matrices


class Cone:

    def __init__(self, *vectors : Collection[float|int]) -> None:
        for vector in vectors:
            validate_domain(vector)
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
