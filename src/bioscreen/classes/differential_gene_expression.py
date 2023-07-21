from bioscreen.classes.base import logger, ComparisonsResults, StatColumns
from bioscreen._imports import *
from bioscreen.classes.results import AnalysisResults, comp_results_from_dir
from bioscreen.classes.comparison import Comparison, CompList
from attrs import define



@define
class DGEResults(AnalysisResults):
    @staticmethod
    def load_results(results_dir, columns, index_name='Gene') -> tuple[ComparisonsResults, CompList]:


        def fn_to_comp(fn):
            logger.info(f'DGEResults: fn_to_comp({fn})')
            test, ctrl, _ = fn.split('.')

            comp = Comparison(control=ctrl, test=test)
            return comp

        print("assuming filename is '{test}.{ctrl}.csv' ")
        compres, comparisons = comp_results_from_dir(
            results_dir,
            fn_to_comp,
            columns
        )
        compres.index.name = index_name
        return compres, comparisons


@define
class LimmaResults(AnalysisResults):

    _coltable = {
        'original': ['Gene', 'logFC', 'AveExpr', 't', 'F', 'P.Value', 'adj.P.Val', 'B', 'P10', 'FDR10'],
        'good': ['Gene', 'LFC', 'Expr', 'T', 'F', 'P', 'FDR', 'LogOdds', 'P10', 'FDR10'],
        'long': ['Gene', 'Log2(FC)', 'Ave. Expression', 't-statistic', 'F-statistic', 'p-value', 'FDR', 'Log odds',
                 '-log10(p)', '-log10(FDR)'],
        'short': ['Gene', 'LFC', 'Expr', 't', 'F', 'p', 'FDR', 'LogOdds', '-log10(p)', '-LOG10(FDR)']
    }

    _score = 'LFC'
    _results_cols = ['Gene', 'LFC', 'p', 'FDR',  'Expr', 'LogOdds']

    @classmethod
    def _get_columns(cls, index_name='Gene'):
        coltab = {}
        for k, v in  cls._coltable.items():
            coltab[k] = copy(v)
            coltab[k][0] = index_name
        results_cols = copy(cls._results_cols)
        results_cols[0] = index_name

        return StatColumns.from_dict(coltab, results_cols=results_cols, score='LFC')


    @classmethod
    def from_dir(cls, results_dir, index_name='Gene'):

        stat_cols =  cls._get_columns()
        (res, comparisons) = DGEResults.load_results(
            results_dir, stat_cols,
            index_name=index_name
        )
        res.index.name = index_name
        return cls(res, stat_cols, comparisons)


DfOrPath = typing.Union[pd.DataFrame, os.PathLike, ]

validator_passthrough = lambda self, attribute, value: value

