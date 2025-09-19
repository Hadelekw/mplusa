import numpy as np

import math
from typing import Any
from collections.abc import Collection


class TropicalNumberMeta(type):
    def __instancecheck__(self, instance: Any, /) -> bool:
        if isinstance(instance, float|int|np.floating|np.integer):
            if instance > -math.inf:
                return True
        return False


class TropicalNumber(metaclass=TropicalNumberMeta):
    pass


def validate_domain(values : Collection[Any]) -> None:
    if all(map(lambda x: isinstance(x, TropicalNumber), values)):
        return
    raise ValueError('Value out of domain.')
