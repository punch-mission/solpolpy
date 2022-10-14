from typing import Dict, Any
from dataclasses import dataclass

import astropy.units as u
from astropy.units import Quantity


@dataclass(frozen=True)
class PolarizedData:
    B: Any = None
    pB: Any = None
    angles: None | Dict[Quantity[u.degree] | Quantity[u.radian], Any] = None
    M: Any = None
    Z: Any = None
    P: Any = None
    Br: Any = None
    Bt: Any = None
    alpha: Any = None

    def extend(self, *args, **kwargs):
        pass



