from bioscreen.classes.base import *
from bioscreen.classes.comparison import CompList
from bioscreen.classes.results import AnalysisResults
from attrs import define, Factory
from bioscreen.utils import ValidationError, validate_sample_details, validate_count_df


def validate_screen_input(count:pd.DataFrame, details:pd.DataFrame, comparisons:CompList=None):
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


@define
class ScreenExperiment:
    name: str
    version: str
    counts: typing.Mapping[str, pd.DataFrame]
    sample_details: pd.DataFrame
    group_details: pd.DataFrame
    comparisons: CompList
    results: typing.Mapping[str, AnalysisResults] = Factory(AMap)

    def __attrs_post_init__(self):
        #validate_comparisons_table(self.comparisons)
        validate_sample_details(self.sample_details)

        # Just test the first count file, assuming others are derived from it
        cnt = list(self.counts.values())[0]
        validate_count_df(cnt)
        validate_screen_input(cnt, self.sample_details, self.comparisons)

    @staticmethod
    def grp_details_from_sample(sample_details):
        groupdeets = sample_details.dropduplicates('SampleGroup').set_index('SampleGroup')
        return groupdeets

    @classmethod
    def metadata_from_text_files(cls, name, version,
            counts, sample_details, comparisons, cnt_type='raw'):
        """Generate Screen experiment from file paths pointing to CSV for
        sample_details or comparisons, and TSV for counts."""
        cnt = pd.read_csv(counts, index_col=0, sep='\t')
        deets = pd.read_csv(sample_details)
        grps = cls.grp_details_from_sample(deets)

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


