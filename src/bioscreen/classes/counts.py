import typing
import attrs
import pandas as pd
from attrs import define
from jttools.data_wrangling import (
    read_csv
)

from bioscreen.utils import (
    Pathy, AMap,
)

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

__all__ = ['Counts']


class Counts(AMap):
    def __init__(
            self, tables:typing.Mapping[str, pd.DataFrame],
            metadata:pd.DataFrame=None,
            as_copies=True,
    ):
        if as_copies:
            tables = AMap({k:t.copy() for k, t in tables.items()})
        super().__init__(tables)
        self.metadata = metadata
        self.validate_metadata()

    def validate_metadata(self):
        if self.metadata is not None:

            meta_idx = self.metadata.index
            for k, tab in self.items():
                idx_match = tab.index.isin(meta_idx)
                if not idx_match.all():
                    logger.warning(f"{len(~idx_match)} count index values not found in metadata, for table {k}")

    @classmethod
    def from_tsv(cls, filename:Pathy | dict[str, Pathy], count_type='raw',
                 sep='\t', meta_cols:list=None, ):
        # to load from df just use the main constructor
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
        return cls(tables=tables, metadata=metadata)

    def keys(self):
        keys = [c for c in super().keys() if c != 'metadata']
        return keys

    def items(self): # this is fine, AttrMapAC has dict2am kwarg.
        for k, v in super().items():
            if k == 'metadata':
                continue
            yield (k, v)


