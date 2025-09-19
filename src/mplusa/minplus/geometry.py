from collections.abc import Collection

from .domain import validate_domain
from .minplus import add, mult


class Cone:

    def __init__(self, *vectors : Collection[float|int]) -> None:
        for vector in vectors:
            validate_domain(vector)
        self.vectors = vectors
