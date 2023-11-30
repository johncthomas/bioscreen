import typing

from bioscreen._imports import *
from bioscreen.classes.base import SampleGroup


from attrs import define, Factory, field, validators

__all__ = ['Comparison', 'CompDict']

DEFAULT_COMP_JOINER = '-'

@define(kw_only=True)
class Comparison:
    """Hold test and control samples, potentially with metadata about how
    a comparison should be made.

    Methods return string representations of the comparison.

    Attributes:
        control: SampleGroup to be used as control in a comparison.
        test: SampleGroup to be compared to control.
        paired: Indicate that replicates in sample group are paired by order.
        name: A name to refer to the comparison. Defaults to "test-control".
            Is key in CompDict
        label: A longer name for figures perhaps. Defaults to name.
        groups: Membership of groups, {"groupname":bool}.

    Properties:
        differential: Is this a differential comparison?
            True if control is a Comparison.

    Methods: (note: differential comps get printed differently)
        str: "test-control". kwarg test_first controls test|control order.
            Uses global DEFAULT_COMP_JOINER
        arrow_Str: "control ➤ test"
        formula_str: "test - control"

    """
    control: SampleGroup|Self
    test: SampleGroup|Self
    paired:bool = False
    name:str = Factory(lambda self: self.joined(), takes_self=True)
    label:str = Factory(lambda self: self.name, takes_self=True)
    groups: dict = field(
        default=Factory(dict, ), validator=validators.deep_mapping(
            validators.instance_of(str), validators.instance_of(bool)
        )
    )

    other_cols:Mapping = None

    # Columns that should be used as kwargs when converting a DF
    #  note, index is lower cased before being used.
    # Groups gets handled differently, cast to separate columns
    _kwargs = pd.Index(['control', 'test', 'paired', 'name', 'label'])

    def __hash__(self):
        return hash(self.name)


    def joined(self,  joiner=DEFAULT_COMP_JOINER, test_first=True,):
        """f"{self.test}-{self.control}" by default."""

        if not self.differential:
            test = str(self.test)
            ctrl = str(self.control)
            if test_first:
                first = test
                last = ctrl
            else: # control first
                first = ctrl
                last = test
            return f"{first}{joiner}{last}"
        else:
            return f"{self.control.name} ➤ {self.test.name}"

    def __str__(self):
        return self.name

    def arrow_str(self, ):
        return self.joined(joiner=' ➤ ', test_first=False)

    def str(self,):
        return self.__str__()

    def formula_str(self):
        """For formulas used in Limma etc, just
        self.test - self.control"""
        if self.differential:
            return f"({self.test.formula_str()}) - ({self.control.formula_str()})"
        return f"{self.test} - {self.control}"

    @property
    def differential(self):
        return isinstance(self.control, Comparison)

    @classmethod
    def from_series(cls, compseries:pd.Series):
        """Convert a row from a DF into a Comparison.

        - Primary Comparison(**kwargs) are pulled from lower-case casted column
        names.
        - Boolean columns are converted to a dict stored in self.groups.
        - Other columns are placed in self.other_cols.
        """
        # lower case index for getting kwargs
        lc_series = compseries.copy()
        lc_series.index = lc_series.index.str.lower()
        kw_mask = lc_series.index.isin(cls._kwargs)

        # convert the series to kwargs, need to filter self._attr to things that
        #   exist in the DF
        kwargs = lc_series[kw_mask].to_dict()

        # # In order to keep the to/from_df mirrored, if there's a `Groups`
        # #   column it's just going into .other_cols.

        # find cols to be used as groups
        bools = [(type(v) == bool) for v in compseries.values]

        # add the columns where bool is True
        if any(bools):
            kwargs['groups'] = compseries[bools].to_dict()
        # else the default value for groups is used.

        # other columns get dumped in self.other_cols
        other_mask = ~(kw_mask|bools)
        if other_mask.any(): # if not all columns were kwargs
            kwargs['other_cols'] = compseries[other_mask].to_dict()

        return cls(**kwargs)

    def to_series(self) -> pd.Series:
        """Return series with index from attributes all with initial caps,
        self.other_cols converted to separate items and groups to boolean items

        Examples:
            cmp = Comparison(
                test='T', control='C', groups={'Good':True, 'Bad':False},
                other_cols={'Etc':'thing'}
            )
            cmp.to_series()

            [1] Control        C
                Test           T
                Paired     False
                Name         T-C
                Label        T-C
                Etc        thing
                Good        True
                Bad        False
                dtype: object
        """
        initcap = lambda s: s[0].upper() + s[1:]
        seriesdict = dict()

        for attr in self._kwargs:
            k = initcap(attr)
            seriesdict[k] = getattr(self, attr)

        if self.other_cols is not None:
            seriesdict.update(self.other_cols)


        if self.groups:
            seriesdict.update(self.groups)


        seriesdict['Differential'] = self.differential

        return pd.Series(seriesdict, index=list(seriesdict.keys()))


# @define
# class DifferentialComparison:
#     """Hold sample relationships for a differential comparison, where
#     one pair of samples (the baseline comparison) is compared to another
#     (the _other_ comparison).
#
#     Properties ending in _str (and `str`) are string representations.
#     In all except formula_str, if both internal comparisons have
#     `comp.name != None`, the names will be used.
#
#     Attributes:
#         baseline: Comparison used as baseline
#         other: Comparison to be compared to baseline.
#         paired: Indicate that replicates in sample group are paired by order.
#         name: A name to refer to the comparison. Defaults to test-control.
#             Is key in CompDict
#
#
#     Properties:
#         str: If both `comparison`s have a `name`, name will be used,
#             otherwise default is "otherTest_otherCtrl-baseTest_baseCtrl".
#         arrow_Str: "baseline ➤ other"
#         formula_str: "(otherTest - otherCtrl) - (baseTest - baseCtrl)".
#             For R formulas.
#     """
#     baseline: Comparison
#     other: Comparison
#     name: str
#     paired:bool = False
#     groups: set[str] = attrs.Factory(set)
#
#     def __attrs_post_init__(self):
#         if type(self.groups) is not set:
#             self.groups = set(self.groups)
#
#     def _compnames_bo(self) -> tuple[str, str]:
#         b = self.baseline
#         o = self.other
#         if (b.name is not None) and (o.name is not None):
#             return b.name, o.name
#         return b.str(joiner='_'), b.str(joiner='_')
#
#     def __str__(self):
#         bn, on = self._compnames_bo()
#         return f"{on}{DEFAULT_COMP_JOINER}{bn}"
#
#     def arrow_str(self) -> str:
#         bn, on = self._compnames_bo()
#         return f"{bn} ➤ {on}"
#
#     def str(self, joiner:str=None, test_first='ignored') -> str:
#         if joiner is None:
#             return self.__str__()
#         return self.__str__().replace(DEFAULT_COMP_JOINER, joiner)
#
#     def formula_str(self) -> str:
#         """For formulas used in Limma etc, just
#         self.test - self.control"""
#         return f"({self.other.test} - {self.other.control}) - ({self.baseline.test} - {self.baseline.control})"


class CompDict(AttrMapAC[str, Comparison]):
    """Mapping of Comparison with samples, paired, unpaired and filter functions.

    If invoked with a Collection[Comparison], it uses .name attr for keys."""

    # def __new__(cls, , value):
    #     if isinstance(obj, Mapping):
    #         return super().__new__(cls)
    #     return super().__new__(cls)

    def __init__(self, comps: Mapping[str, Comparison] | Collection[Comparison]):
        if not isinstance(comps, Mapping):
            comps = {c.name: c for c in comps}
        super().__init__(comps)

    def samples(self) -> list[SampleGroup]:
        """List samples used as either test or control"""
        c, s = [set([getattr(cmp, k) for cmp in self.values()]) for k in ('control', 'test')]
        return list(c | s)

    def filter_by_group(self, groups: str | Collection[str]):
        if type(groups) == str:
            groups = [groups]

        l = {c.name: c for c in self.values()
             if all([c.groups[grp] for grp in groups])}
        return CompDict(l)

    def filter_by(self, f: Callable[[Comparison], bool]) \
            -> Self:
        """Filter comparisons from this CompList using a function that
        takes a Comparison, and returns True or False.

        eg:
            comparisons.filter_by(
                lambda cmp: sample_details.SampleGroup.isin(
                    [cmp.test, cmp.control]
                ).any()
            )
        """
        return CompDict({
            c.name: c for c in self.values() if f(c)
        })

    def to_df(self) -> pd.DataFrame:
        """Return DF of Comparisons values. DifferentialComparison will have
        a Comparison as the value in Test Control columns."""
        df = pd.DataFrame({k: cmp.to_series() for k, cmp in self.items()}).T
        return df

    def to_joined(self, joiner=DEFAULT_COMP_JOINER, test_first=True) -> List[str]:
        return [c.joined(test_first=test_first, joiner=joiner) for c in self.values()]

    def to_str(self) -> List[str]:
        return [c.str() for c in self.values()]

    def to_formulas(self) -> List[str]:
        return [c.formula_str() for c in self.values()]

    def to_arrowed(self) -> List[str]:
        return [c.arrow_str() for c in self.values()]

    def names(self) -> list[str]:
        return [c.name for c in self.values()]

    @classmethod
    def from_df(cls, df: pd.DataFrame, ):
        d = {k: Comparison.from_series(row) for k, row in df.iterrows()}
        return cls(d)

    def __iter__(self) -> typing.Iterator[str]:
        return super().__iter__()

    def __getitem__(self, item:str) -> Comparison:
        return super().__getitem__(key=item)

    def keys(self) -> list[str]:
        # ._mapping is where AttrMap keeps it's info

        return list(self._mapping.keys())

    def values(self) -> typing.ValuesView[Comparison]:
        return super().values()

    def items(self, dict2attrmap=False) -> typing.ItemsView[str, Comparison]:
        return super().items(dict2attrmap=dict2attrmap)


def samples_of_comp(comparison:Comparison, sample_details):
    return sample_details.loc[sample_details.SampleGroup.isin([comparison.control, comparison.test])].index