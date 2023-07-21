from bioscreen._imports import *
from bioscreen.classes.base import SampleGroup


from attrs import define


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
        name: A name to refer to the comparison if we don't want to use
            test-control or whatever.

    Properties:
        str: "test-control" by default, can change joiner and test|control order.
        arrow_Str: "control ➤ test"
        formula_str: "test - control"

    """
    control: SampleGroup
    test: SampleGroup
    paired:bool = False
    name:str = None #Factory(lambda self: str(self), takes_self=True)

    def __str__(self, test_first=True):
        test = str(self.test)
        ctrl = str(self.control)
        if test_first:
            first = test
            last = ctrl
        else: # control first
            first = ctrl
            last = test
        return f"{first}{DEFAULT_COMP_JOINER}{last}"

    def arrow_str(self, ):
        return self.str(' ➤ ', test_first=False)

    def str(self, joiner:str=None, test_first=True):
        if joiner is None:
            return self.__str__()
        return self.__str__(test_first=test_first).replace(DEFAULT_COMP_JOINER, joiner)

    def formula_str(self):
        """For formulas used in Limma etc, just
        self.test - self.control"""
        return f"{self.test} - {self.control}"


@define
class DifferentialComparison:
    """Hold sample relationships for a differential comparison, where
    one pair of samples (the baseline comparison) is compared to another
    (the _other_ comparison).

    Properties ending in _str (and `str`) are string representations.
    In all except formula_str, if both internal comparisons have
    `comp.name != None`, the names will be used.

    Attributes:
        baseline: Comparison used as baseline
        other: Comparison to be compared to baseline.
        paired: Indicate that replicates in sample group are paired by order.
        name: A name to refer to the comparison if we don't want to use
            test-control or whatever.


    Properties:
        str: If both `comparison`s have a `name`, name will be used,
            otherwise default is "otherTest_otherCtrl-baseTest_baseCtrl".
        arrow_Str: "baseline ➤ other"
        formula_str: "(otherTest - otherCtrl) - (baseTest - baseCtrl)".
            For R formulas.


    """
    baseline: Comparison
    other: Comparison
    paired:bool = False
    name:str = None #Factory(lambda self: str(self), takes_self=True)

    def _compnames_bo(self) -> tuple[str, str]:
        b = self.baseline
        o = self.other
        if (b.name is not None) and (o.name is not None):
            return b.name, o.name
        return b.str(joiner='_'), b.str(joiner='_')

    def __str__(self):
        bn, on = self._compnames_bo()
        return f"{on}{DEFAULT_COMP_JOINER}{bn}"

    def arrow_str(self) -> str:
        bn, on = self._compnames_bo()
        return f"{bn} ➤ {on}"

    def str(self, joiner:str=None, test_first='ignored') -> str:
        if joiner is None:
            return self.__str__()
        return self.__str__().replace(DEFAULT_COMP_JOINER, joiner)

    def formula_str(self) -> str:
        """For formulas used in Limma etc, just
        self.test - self.control"""
        return f"({self.other.test} - {self.other.control}) - ({self.baseline.test} - {self.baseline.control})"


class CompList(List[Comparison|DifferentialComparison]):
    """List of Comparison with samples, paired, unpaired and filter functions.

    Can take DifferentialComparisons, but can't return a DF with them and
    things might get tricky. """
    def samples(self) -> list[SampleGroup]:
        c, s = [set([getattr(cmp, k) for cmp in self]) for k in ('control', 'test')]
        return list(c|s)

    def paired_comps(self) -> Self:
        return CompList([c for c in self if c.paired])

    def unpaired_comps(self) -> Self:
        return CompList([c for c in self if not c.paired])

    def filter(self, attribute:str, function:Callable) -> Self:
        """Filter complist by attribute value using boolean function.

        Args:
            attribute: Name of the attr to test.
            function: Function to test the attribute values.
        """
        l = []
        for comp in self:
            if function(getattr(self, attribute)):
                l.append(comp)
        return CompList(l)

    def to_df(self) -> pd.DataFrame:
        """Return DF of Comparisons values. DifferentialComparisons are ignored."""
        return pd.DataFrame([attrs.asdict(cmp) for cmp in self
                             if type(cmp) is not DifferentialComparison])

    def to_str(self, joiner=DEFAULT_COMP_JOINER, test_first=True) -> List[str]:
        return [c.str(test_first=test_first, joiner=joiner) for c in self]

    def to_formulas(self) -> List[str]:
        return [c.formula_str() for c in self]

    def to_arrowed(self) -> List[str]:
        return [c.arrow_str() for c in self]

    @classmethod
    def from_df(cls, df:pd.DataFrame,
                test_col:str='Test', ctrl_col:str='Control',
                other_cols=('Paired', 'Name')):

        comps = []

        for _, row in df.iterrows():
            t = row[test_col]
            c = row[ctrl_col]
            others = {k.lower():row[k] for k in other_cols}
            comps.append(
                Comparison(control=c, test=t, **others)
            )

        return cls(comps)
