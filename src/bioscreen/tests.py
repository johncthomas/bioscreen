import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import pandas as pd

from bioscreen.utils import ValidationError
from bioscreen.classes.experiment import validate_screen_input


def test_validate_screen_input():

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