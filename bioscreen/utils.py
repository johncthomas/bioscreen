import pandas as pd
import numpy as np
import typing
from attrdictionary import AttrMap


class AttrMapAC(AttrMap):
    """Tab completable AttrMap"""
    def f(self):
        for thing in self:
            print(thing)
    def __dir__(self):
        super_dir = list(super().__dir__())
        string_keys = [str(key) for key in self if isinstance(key, str)]
        return super_dir + [key for key in string_keys if key not in super_dir]
AMap = AttrMapAC

is_numeric = pd.api.types.is_numeric_dtype
neglog10 = lambda p: -np.log10(p)