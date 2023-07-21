"""universal imports, or stuff I might want and dont mind having in the namespace"""


import enum
import os
import pathlib
import collections
import itertools
import typing
import logging
from copy import copy, deepcopy



import attrs
import pandas as pd

#from attrs import define, Factory, field
# has to be imported into the local module currently or the pycharm linter gets confused.

from typing import (
    Optional, Tuple, Callable, Any, Collection, Literal, Dict, List,
    Self, Hashable, Mapping
)

from jttools.data_wrangling import (
    read_csv, AttrMapAC, read_tsv, rename_columns, is_numeric,

)

from jttools.statistics import (
    neglog10, log2p1, apply_log2, normalise_abundance,
    normalise_zscore, normalise_median
)

from bioscreen.utils import(
    Pathy, AMap
)
