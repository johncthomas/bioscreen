import typing
import attrs
import pandas as pd
from attrs import define
from jttools.data_wrangling import (
    read_csv
)

from bioscreen.utils import (
    Pathy, AMap
)


@define
class Counts:
    counts:typing.Mapping[str, pd.DataFrame] = attrs.Factory(AMap)
    metadata:pd.DataFrame = None

    @classmethod
    def from_tsv(cls, filename:Pathy | dict[str, Pathy], count_type='raw',
                 sep='\t', meta_cols:list=None, ):
        metadata = None
        tables = AMap()
        if meta_cols is None:
            meta_cols = []
        if isinstance(filename, Pathy):
            filename = {count_type:filename}
        for ct, fn in filename.items():
            c = read_csv(fn, sep=sep)
            if meta_cols and (metadata is None):
                metadata = c.loc[:, meta_cols]
            c.drop(meta_cols, inplace=True, axis='columns')
            tables[ct] = c
        return cls(counts=tables, metadata=metadata)



