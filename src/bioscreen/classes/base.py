import pandas as pd

from bioscreen._imports import *

import attrs

logger = logging.getLogger("bioscreen.classes")
logging.basicConfig()


StrOptional = Optional[str]
CollectionOptional = Optional[Collection]

CompsResultDF = pd.DataFrame
Sample = str
SampleGroup = str

__all__ = ['SigCols', 'StatColumns', 'logger', 'CompsResultDF', 'Sample', 'SampleGroup', 'StrOptional']

# # #@attrs.define(kw_only=True, )
# class UnmappedStatCol(str):
#     """Base class for StatCol, used directly for setting cross analysis
#     defaults where source column name may vary.
#
#     Args and attributes:
#         key: Key used to retrieve column from processed stat table
#         table: Col name for use in exported tables
#         label: Longer label used in figures
#         formatted: Label with MathJax style formatting if appropriate (not wrapped in $)"""
#     def __init__(
#             self,
#             *,
#             key: str,
#             table:str=None,
#             label:str=None,
#             formatted:str=None,
#     ):
#         super().__init__(key)
#         self.key = key
#         self._elsekey = elsekey = lambda s: s if s is not None else key
#         self.table = elsekey(table)
#         self.label = elsekey(label)
#         self._formatted = formatted
#
#     @property
#     def formatted(self):
#         return self._elsekey(self.formatted)
def elsekey(s, key):
    return s if s is not None else key

class StatCol(str):
    """Statistics column, providing methods for different string reps of
    the statistic. Used for mapping columns in tables written by analysis
    programs.

    Use PascalCase. Keyword args only except for key, which is passed to
    str.__new__.

    key: Key used to retrieve column from processed stat table. Should be
        useable as an attribute.
    original: Column string used in analysis output
    table: Col string for use in exported tables
    label: Longer label used in figures
    formatted: MathJax style formatting if appropriate
        (not wrapped in $)
    """

    _dictkeys = ('key', 'original', 'table', 'label', 'formatted')

    # we use new because subclassing an imutable obj.
    def __new__(
            cls,
            key: str,
            *,
            original: str = None,
            table: str = None,
            label: str = None,
            formatted: str = None,
    ):

        self = super().__new__(cls, key)
        self.key = key
        self.table = elsekey(table, key)
        self.label = elsekey(label, key)
        self._formatted = formatted
        self.original = original

        return self

    @property
    def formatted(self):
        return elsekey(self._formatted, self.key)
    @formatted.setter
    def formatted(self, v):
        self._formatted = v # linter confused because there's no init

    def to_dict(self):
        return {k:getattr(self, k) for k in StatCol._dictkeys}



class SigCols:
    """Default values for significance stat columns, and LFC.

    Use CLASSMETHOD `map` or access unmapped statcols as attributes.
    FDR, FDR10, p, p10 and LFC.

    Methods:
        map: StatCols for all defined for all defined StatCols by suppling
        key->original as kwargs.

    """
    # if you change the key, change the attr name
    FDR = StatCol(original='UNMAPPED', key='FDR', table='FDR', label='FDR', )
    FDR10 = StatCol(original='UNMAPPED', key='FDR10', table='FDR (-log10)',
                    label='-log10(FDR)',
                        formatted='-log_{10}(FDR)')
    p = StatCol(original='UNMAPPED', key='p', table='p-value', label='p', )
    p10 = StatCol(original='UNMAPPED', key='p10', table='p-value (-log10)',
                  label='-log10(p)',
                          formatted='-log_{10}(p)', )
    LFC = StatCol(original='UNMAPPED', key='LFC', table='LFC',
                  label='log2(FC)', formatted='log_2(FC)')

    @classmethod
    def map(cls, **key_to_original) -> Mapping[str, StatCol]:
        l = {}
        for k, o in key_to_original.items():
            kwargs = getattr(cls, k).to_dict()
            kwargs['key'] = k
            kwargs['original'] = o
            col = StatCol(**kwargs)
            l[k] = col
        return StatColumns(l)


class StatColumns(AMap):
    """A mapping of the column name used in my results tables to a StatCol."""
    def __init__(self, columns=Mapping[str, StatCol], order=None):
        if order is None:
            order = columns.keys()
        cols = {col:columns[col] for col in order}
        super().__init__(cols)

    @property
    def short(self) -> AMap:
        return AMap({k:col.table for k, col in self.items()})

    @property
    def long(self) -> AMap:
        return AMap({k: col.label for k, col in self.items()})

    @property
    def original(self) -> AMap:
        return AMap({k: col.original for k, col in self.items()})

    @property
    def formatted(self) -> AMap:
        return AMap({k: col.formatted for k, col in self.items()})

    def original_to_key(self) -> dict:
        return {col.original:col.key for col in self.values()}

    def get_mapping(self, askey='key', asvalue:str=StatCol) -> dict[str, StatCol|str]:
        """Return a dict mapping desired attribute to its StatCol, or
        to another attribute if asvalue is a string that names the attribute."""
        if asvalue is StatCol:
            return {getattr(c, askey):c for c in self.values()}
        return {getattr(c, askey):getattr(c, asvalue) for c in self.values()}


    @classmethod
    def from_records(cls, records:list[dict[str,str]]):
        return cls({kw['key']:StatCol(**kw) for kw in records})

    def to_records(self):
        return [c.to_dict() for c in self.values()]

    @classmethod
    def from_df(cls, df:pd.DataFrame):
        df = df.copy()
        df.columns = [c.lower() for c in df.columns]
        return StatColumns.from_records(df.to_dict('records'))

    def to_list(self):
        return list(self.values())


    def __iter__(self) -> typing.Iterator[str]:
        return super().__iter__()

    def __getitem__(self, item:str) -> StatCol:
        return super().__getitem__(key=item)

    def keys(self) -> list[str]:
        # ._mapping is where AttrMap keeps it's info

        return list(self._mapping.keys())

    def values(self) -> typing.ValuesView[StatCol]:
        return super().values()

    def items(self, dict2attrmap=False) -> typing.ItemsView[str, StatCol]:
        return super().items(dict2attrmap=dict2attrmap)


# class DFHolder:
#     def __init__(self, df:pd.DataFrame):
#         self.df = df
#
#         if type(df.columns) == pd.MultiIndex:
#             attr = df.columns.levels[0]

    #