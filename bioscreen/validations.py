import pandas as pd
import logging
from utils import (
    is_numeric
)

logging.basicConfig()
logger = logging.getLogger(__name__)
#
# logger.setLevel('WARNING')

class ValidationError(Exception):
    pass


def validate_cols(columns:pd.Index, required_cols, tablename='Some', ):

    if not all([(k in columns) for k in required_cols]):
        raise ValidationError(f"{tablename} table requires columns {required_cols}")

def validate_screen_input(count:pd.DataFrame, details:pd.DataFrame, comparisons:pd.DataFrame=None):
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
        validate_cols(comparisons.columns, ['Test', 'Control'], 'Comparisons')
        compsampgroups = pd.Series(list({*comparisons.Test.values, *comparisons.Control.values}))

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

def _test_validate_screen_input():
    logger.setLevel(logging.INFO)
    s1 = list('ABCDEF')
    s2 = list('ABCDEX')
    g1 = list('xxxyyy')
    g2 = list('zzzyyy')

    def get_tables(cnt_sample, deet_samples, deet_sampgrp, comps_sampgrp=None):
        """supply sample values for count, details, and comparison tables"""
        cnt = pd.DataFrame(1, columns=cnt_sample, index=list(range(5)), )
        deets_vals = {
            'Sample':deet_samples,
            'SampleGroup':deet_sampgrp,
            'Treatment':list(range(len(deet_samples)))}
        deets = pd.DataFrame(deets_vals).set_index('Sample', drop=False)
        if comps_sampgrp is None:
            comps = None
        else:
            comps = pd.DataFrame({'Test':comps_sampgrp, 'Control':comps_sampgrp})
        return dict(count=cnt, details=deets, comparisons=comps)

    tests = [('all same, with comps',   (s1, s1, g1, g1)),
             ('all same, no comps', (s1, s1, g1, None)),
             ('deets/count samples differ, no comps', (s2, s1, g1, None)),
             ('deets/count samples differ, with comps', (s2, s1, g1, g1)),
             ('deets/comps sampgroups differ', (s1, s1, g1, g2)),
             ('all differ', (s2, s1, g1, g2)),
             ]
    for testname, get_table_args in tests:
        logger.info('****************')

        logger.info(testname)
        tables = get_tables(*get_table_args)
        # for k, tab in tables.items():
        #     print(k, tab, sep='\n')
        try:
            validate_screen_input(**tables)
        except Exception as e:
            if type(e) is ValidationError:
                logger.info(f'Error caught: {e}',)
            else:
                raise e

    logger.setLevel(logging.WARNING)

def validate_count_df(cnt:pd.DataFrame):
    numcols = cnt.apply(is_numeric)
    if not numcols.all():
        raise ValidationError(f"Count DF contains non-numeric values in column(s): {cnt.columns[~numcols]}")

def validate_sample_details(deets:pd.DataFrame):
    validate_cols(deets.columns, ['Sample', 'SampleGroup', 'Treatment', 'Time', ])

def validate_comparisons_table(comps:pd.DataFrame):
    validate_cols(comps.columns, ['Test', 'Control', ])

if __name__ == '__main__':
    _test_validate_screen_input()