import numpy as np

from bioscreen._imports import *
import pandas as pd

from bioscreen.classes.base import *
from bioscreen.classes.comparison import CompDict, Comparison
from bioscreen.classes.results import AnalysisResults
from attrs import define, Factory
from bioscreen.utils import ValidationError, validate_sample_details, validate_count_df



__all__ = [
    'ScreenExperiment',
    'validate_screen_input',
    'get_replicates_of_comparison'
]

SDCOLFIXES = {'treat': 'Treatment',
                 'Treat': 'Treatment',
                 'conc': 'Concentration',
                 'Conc': 'Conentration',
                 'Conc.': 'Concentration',
                 'SampleGroups': 'SampleGroup',
                 'Rep': 'Replicate',
                 'rep': 'Replicate'}



def validate_screen_input(count:pd.DataFrame, details:pd.DataFrame, comparisons:CompDict=None):
    """Check that the samples/sampleGroups match across count details and comparisons dataframes"""

    def notin(a, b):
        missing = ~a.isin(b)
        return a[missing]

    # Warn about samples omitted from either details or count
    cnotind = notin(count.columns, details.index)
    dnotinc = notin(details.index, count.columns)

    if len(cnotind) > 0:
        logger.warning(f"Count samples not in details:\n\t{','.join(cnotind)}")

    if len(dnotinc) > 0:
        logger.warning(f"Details samples not in counts:\n\t{','.join(dnotinc)}")


    #
    if comparisons is not None:
        #validate_cols(comparisons.columns, ['Test', 'Control'], 'Comparisons')
        compsampgroups = pd.Series(comparisons.samples())

        # All groups mentioned here need to be in details and the samples need
        #   to be in counts
        csg_found = compsampgroups.isin(details.SampleGroup)
        if not csg_found.all():
            raise ValidationError(f"Not all comparison sampleGroups in details: \n\t"
                               f"{', '.join(compsampgroups[~csg_found])}")

        used_samples = details.loc[details.SampleGroup.isin(compsampgroups)].index
        s_in_c = used_samples.isin(count.columns)
        if not s_in_c.all():
            raise ValidationError("Some samples required for comparisons not found in counts:\n\t"
                                  f"{list(used_samples[~s_in_c])}")


def get_replicates_of_comparison(sample_details, comparison:Comparison) -> np.ndarray[Sample]:
    """Control and test samples used in a comparison. Controls first."""
    m = sample_details.SampleGroup.isin([comparison.control, comparison.test])
    reps = sample_details.loc[m, 'Sample'].values
    return reps


@define
class ScreenExperiment:
    name: str
    version: str
    counts: typing.Mapping[str, pd.DataFrame]
    sample_details: pd.DataFrame = attrs.field(converter=lambda s: s.copy())

    #optional
    comparisons: CompDict = None
    group_details: pd.DataFrame = Factory(
        lambda self: self.group_details_from_sample(self.sample_details),
        takes_self=True,
    )
    results: typing.Mapping[str, AnalysisResults] = Factory(AMap)
    primary_counts:str='raw'
    fix_columns:bool = True

    def __attrs_post_init__(self):
        #validate_comparisons_table(self.comparisons)
        if self.fix_columns:
            df_rename_columns(self.sample_details, SDCOLFIXES,
                           inplace=True, verbose=True)
        validate_sample_details(self.sample_details)

        # Just test the first count file, assuming others are derived from it
        cnt = self.counts[self.primary_counts]
        validate_count_df(cnt)
        validate_screen_input(
            count=cnt,
            details=self.sample_details,
            comparisons=self.comparisons
        )

    def replicates_of_comparison(self, comparison:Comparison|str) -> np.ndarray[Sample]:
        """Control and test samples used in a comparison. Controls first."""
        if type(comparison) is not Comparison:
            comparison = self.comparisons[comparison]
        return get_replicates_of_comparison(self.sample_details, comparison)

    @staticmethod
    def group_details_from_sample(sample_details:pd.DataFrame):
        groupdeets = sample_details.drop_duplicates(
            'SampleGroup'
        ).set_index('SampleGroup', drop=False).drop('Sample', errors='ignore', axis=1)
        return groupdeets

    @classmethod
    def from_text_files(cls, name, version,
                        counts, sample_details, comparisons, cnt_type='raw'):
        """Generate Screen experiment from file paths pointing to CSV for
        sample_details or comparisons, and TSV for counts."""
        cnt = pd.read_csv(counts, index_col=0, sep='\t')
        deets = pd.read_csv(sample_details)
        grps = cls.group_details_from_sample(deets)

        deets.set_index(deets.columns[0], inplace=True, drop=False)

        comps =  pd.read_csv(comparisons)

        return cls(
            name=name,
            version=version,
            counts=AMap({cnt_type:cnt}),
            sample_details=deets,
            group_details=grps,
            comparisons=comps
        )




