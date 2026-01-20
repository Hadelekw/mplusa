import numpy as np

import math
from typing import Any
from collections.abc import Collection


class ArcticNumberMeta(type):
    def __instancecheck__(self, instance: Any, /) -> bool:
        if isinstance(instance, float|int|np.floating|np.integer):
            if instance < math.inf:
                return True
        return False


class ArcticNumber(metaclass=ArcticNumberMeta):
    pass


def validate_domain(value : Any) -> None:
    if isinstance(value, Collection):
        for item in value:
            validate_domain(item)
        return
    if isinstance(value, ArcticNumber):
        return
    raise ValueError('Value out of domain.')
