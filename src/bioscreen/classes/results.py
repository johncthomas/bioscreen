import pandas as pd

from bioscreen._imports import *
from bioscreen.classes.base import *
from bioscreen.classes.comparison import CompDict, Comparison
from attrs import define
import xlsxwriter

from jttools.excel import add_stats_worksheet

from bioscreen.utils import ValidationError

__all__ = ['AnalysisResults', 'comp_results_from_dir']

def comp_results_from_dir(
        results_dir, fn_to_comp: Callable,
        columns:StatColumns=None, sep=',',
    ) -> tuple[CompsResultDF, CompDict]:
    """Load csv in dir, return a multiindexed DF with columns parsed
    from filenames (no directory) using fn_to_col.

    args:
        results_dir: path with CSV files
        fn_to_comp: function that probably should return Comparison
        sep: passed to pd.read_csv"""
    results_dir = pathlib.Path(results_dir)
    results = {}
    comparisons = []
    for fn in os.listdir(results_dir):
        comp = fn_to_comp(fn)
        tbl = pd.read_csv(results_dir / fn,  sep=sep, index_col=0)
        if columns:
            tbl.columns = df_rename_columns(tbl, columns.original_to_key())
        results[str(comp)] = tbl
        comparisons.append( comp)
    comparisons = CompDict({c.name:c for c in comparisons})
    return (pd.concat(results, axis='columns'), comparisons)

def rename_filter_stat_cols(df, cols:StatColumns) -> pd.DataFrame:
    """Rename col.original to col.key. Drop any not included in cols."""
    df_rename_columns(df, cols.original_to_key(), inplace=True, )

    # drop unused columns
    okay = df.columns.isin(cols.keys())
    return df.loc[:, okay]

def add_log10_sig_cols(table:pd.DataFrame, sig_keys=('p', 'FDR'))\
        -> pd.DataFrame:
    """Add –log10(significance) columns to DF. New column names
    will append "10" e.g. p10, FDR10.

    Mutates the table in place, and returns it."""
    for k in sig_keys:
        table.loc[:, k+'10'] = table[k].apply(neglog10)
    return table

def convert_stats_tables(
        dfs:Mapping[str, pd.DataFrame],
        cols:StatColumns,
        log10_sig=('p', 'FDR')
) -> CompsResultDF:
    """Rename stat cols and concat dataframes into single multindex DF."""
    new_tables = {}

    for k, tab in dfs.items():
        tab = rename_filter_stat_cols(tab, cols)
        if log10_sig:
            add_log10_sig_cols(tab, log10_sig)
        new_tables[k] = tab

    table = pd.concat(new_tables, axis='columns')

    if not table.any().any():
        logger.warning("Returned table does not contain any non-null values.")
    return table

def comp_idx(resultsDF) -> pd.Index:
    return resultsDF.columns.levels[0]
def stat_idx(resultsDF) -> pd.Index:
    """resultsDF.columns.levels[1]"""
    return resultsDF.columns.levels[1]

@define(kw_only=True)
class AnalysisResults:
    """A set of comparison results."""
    table:CompsResultDF
    comparisons:CompDict
    columns:StatColumns
    scorekey:str

    def __attrs_post_init__(self):
        validate_comps_df(self.table, self.columns, self.comparisons)

    @staticmethod
    def _table_builder(tables:Mapping[str, pd.DataFrame], comparisons, columns, log10_sig=('p', 'FDR')):
        """Take dict of DF, convert stat column names and return single
        multi-indexed table.

        If it's already a properly formatted table, it's returned as is.
        (so that AnalysisResults can be built with results that have already
        been parsed or not)"""

        if isinstance(tables, pd.DataFrame):
            if validate_comps_df(tables, columns, comparisons):
                table = tables
            else:
                raise ValidationError('should have already happend.')
        elif isinstance(tables, Mapping):
            table = convert_stats_tables(tables, columns, log10_sig=log10_sig)
        else:
            raise ValueError("Pass a mapping of comparisonkey->DF")
        return table

    @classmethod
    def build(cls, tables: CompsResultDF | Mapping[str, pd.DataFrame],
              comparisons:CompDict, columns:StatColumns, scorekey=''):
        # probably this will be reimplimented by each child class
        #  with at least the scorekey value set

        table = cls._table_builder(tables, comparisons, columns)

        cls(table=table, comparisons=comparisons,
                            columns=columns, scorekey=scorekey)

    # make DF methods available
    @property
    def loc(self):
        return self.table.loc

    def __getitem__(self, k):
        return self.table[k]

    # todo access comps by attribute (with autocomplete)

    def get_stat_table(self, key) -> pd.DataFrame:
        return self.table.xs(key, level=1, axis=1)

    @property
    def score_table(self):
        return self.get_stat_table(self.scorekey)

    @property
    def fdr_table(self):
        return self.get_stat_table(SigCols.FDR.key)

    @property
    def p_table(self):
        return self.get_stat_table(SigCols.p.key)

    @property
    def fdr10_table(self):
        #return self.fdr_table.apply(neglog10)
        return self.get_stat_table(SigCols.FDR10.key)

    @property
    def p10_table(self):
        return self.get_stat_table(SigCols.p10.key)


    def result_table(self, ctrl_or_comp:typing.Union[str, Comparison],
                     test:str=None, ) -> pd.DataFrame:
        """pass a comp-string, Comparison or (samplestr, samplestr)
        and get a single table with results for that comp."""
        if isinstance(ctrl_or_comp, Comparison):
            comp = ctrl_or_comp
        elif test is not None:
            comp = str(Comparison(control=ctrl_or_comp, test=test))
        else:
            comp = str(ctrl_or_comp)

        return self.table[comp].copy()


    def write_comp_results_to_excel(
            self, filename:Pathy,
            included_comparisons:Collection[str]= 'all',
            drop_na_rows=True, xlsx_table_opts:dict=None):
        """Write comparisons table as excel workbook with one table
        per comp.

        NA values will be replaced with 0 except for 'p' and 'FDR'.

        Arguments:
            filename: where file should be saved.
            included_comparisons: set to only include certain stat columns.
            drop_na_rows: if True, rows that contain only NA will be dropped.
            xlsx_table_opts: kwargs passed to

        **xlx_table_opts passed to worksheet.add_table,
        see https://xlsxwriter.readthedocs.io/working_with_tables.html"""
        # todo formats: specified in ResultsColumns i guess, since we need to define xlsxwriter obj
        #   that will be shared across columns

        workbook = xlsxwriter.Workbook(filename)

        for compname, comp in self.comparisons.items():
            if (included_comparisons != 'all') and (compname not in included_comparisons):
                continue
            table = self.result_table(comp)

            table = out_table_formatter(
                table,
                included_columns=included_comparisons,
                score_col=self.columns.score,
                sort_table=True,
                drop_na_rows=drop_na_rows,
            )

            add_stats_worksheet(
                workbook=workbook,
                table=table,
                sheet_name=comp.arrow_str(),
                xlsx_table_opts=dict(
                    name=comp.joined('.')
                ) | xlsx_table_opts
            )

        workbook.close()


def out_table_formatter(
        table:pd.DataFrame,
        included_columns:Collection[str]='all',
        sort_table=True,
        sort_sig_col=SigCols.p,
        score_col:str=None,
        sig_cols=(SigCols.p, SigCols.FDR),
        drop_na_rows=True,
):
    """Replace NaN values with 0, or 1 for significance measures, and sort columns."""
    # can end up with missing rows
    # set the sig columns 1, and the rest to zero

    # remove dead rows
    if drop_na_rows:
        dead_rows = table.isna().all(1)
        if sum(dead_rows):
            logging.info(f"Dropping {sum(dead_rows)} empty rows")
            table = table.loc[~dead_rows]

    if table.isna().any().any():
        logger.info(f"NAs found, n={table['P'].isna().sum()}, replacing with 1 or zero.")

    for k in sig_cols:
        table.loc[table[k].isna(), k] = 1
    table.fillna(0, inplace=True)

    # # select and rename columns
    if included_columns != 'all':
        table = table.reindex(columns=list(included_columns))

    if sort_table:
        if score_col is not None:
            table.sort_values([sort_sig_col, score_col], ascending=[True, False], inplace=True)
        else:
            table.sort_values(sort_sig_col, ascending=True, inplace=True)
    return table


def validate_comps_df(df:CompsResultDF,
                      columns:StatColumns,
                      comparisons:CompDict) -> bool:
    # deal with the table, it could be an already formated results DF
    #   or mapping of comparison->DF

    # if a multi-indexed table is passed...
    if not hasattr(df.columns, 'levels'):
        raise ValidationError(f"Comparisons results table must be multiindexed, by comparison")

    # ..and it has the right columns, that's fine.
    # if not all(stat_idx(df).isin(columns.keys())):
    if not all([c in stat_idx(df) for c in columns.keys()]):
        logger.warning(f"Stat columns do not match."
                        f"\n df={stat_idx(df)}; expected={columns.keys()}")

    missing_comps = [c for c in comparisons.keys() if c not in comp_idx(df)]
    if missing_comps:
        logger.warning(f"Comparison in comparisons not found in results table:\n\t{missing_comps}")

    return True


if __name__ == '__main__':
    pass
    # import pickle
    # d = pathlib.Path('/mnt/m/tasks/NA327_Proteomics_UbPulldown/pickles2')
    # with open(d/'smolcomp.1.pickle', 'rb') as f:
    #     comp = pickle.load(f)
    # with open(d/'test_tables.1.pickle', 'rb') as f:
    #     tables = pickle.load(f)
    #
    # AnalysisResults(tables=tables, comparisons=comp, columns=get_limma_cols())
