import logging
import os

import pandas as pd
from jttools.data_wrangling import AttrMapAC, is_numeric
import typing

Pathy = os.PathLike | str

AMap = AttrMapAC

DfOrPath = typing.Union[pd.DataFrame, os.PathLike, ]

validator_passthrough = lambda self, attribute, value: value



class ValidationError(Exception):
    pass


def validate_cols(columns:pd.Index, required_cols, tablename='Some', ):

    if not all([(k in columns) for k in required_cols]):
        raise ValidationError(f"{tablename} table requires columns {required_cols}")


def validate_count_df(cnt:pd.DataFrame):
    numcols = cnt.apply(is_numeric)
    if not numcols.all():
        raise ValidationError(f"Count DF contains non-numeric values in column(s): {cnt.columns[~numcols]}")


def validate_sample_details(deets:pd.DataFrame):
    validate_cols(deets.columns, ['Sample', 'SampleGroup', 'Treatment', 'Time', ])


def validate_comparisons_table(comps:pd.DataFrame):
    validate_cols(comps.columns, ['Test', 'Control', ])

