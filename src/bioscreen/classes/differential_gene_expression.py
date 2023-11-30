from bioscreen.classes.base import logger, CompsResultDF, StatColumns, SigCols
from bioscreen._imports import *
from bioscreen.classes.results import AnalysisResults, comp_results_from_dir
from bioscreen.classes.comparison import Comparison, CompDict
from attrs import define

__all__ = ['LimmaResults', 'load_results', 'LIMMACOLS']


def load_results(results_dir, columns, index_name='Gene') \
        -> tuple[CompsResultDF, CompDict]:
    def fn_to_comp(fn):
        logger.info(f'DGEResults: fn_to_comp({fn})')
        test, ctrl, _ = fn.split('.')

        comp = Comparison(control=ctrl, test=test)
        return comp

    compres, comparisons = comp_results_from_dir(
        results_dir,
        fn_to_comp,
        columns
    )
    compres.index.name = index_name
    return compres, comparisons


def get_limma_cols() -> StatColumns:
    # this ended up silly because I kept changing my mind about how to do it

    sig_cols = SigCols.map(**dict(zip(
        (SigCols.LFC.key, SigCols.p.key, SigCols.FDR.key, SigCols.p10.key, SigCols.FDR10.key),
        ('logFC', 'P.Value', 'adj.P.Val', 'P10', 'FDR10'),
    )))

    coltable = {
        'original': ['AveExpr', 't', 'F', 'B', ],
        'key': ['Expr', 't', 'F', 'LogOdds', ],
        'label': ['Ave. Expression', 't-statistic', 'F-statistic', 'Log odds'],
        'table': ['Expr', 't', 'F', 'LogOdds', ]
    }

    spec_cols = StatColumns.from_df(pd.DataFrame(coltable))
    intermediatecols: StatColumns = StatColumns(dict(spec_cols) | dict(sig_cols))

    order = ['LFC', 'p', 'FDR', 'Expr', 't', 'F', 'LogOdds', 'p10', 'FDR10']
    assert set(order) == set(intermediatecols.keys())
    limma_cols = StatColumns({k: intermediatecols[k] for k in order})

    return limma_cols
LIMMACOLS = get_limma_cols()

@define(kw_only=True)
class LimmaResults(AnalysisResults):

    @classmethod
    def build(cls, tables: CompsResultDF | Mapping[str, pd.DataFrame],
              comparisons:CompDict, columns=LIMMACOLS,
              scorekey='LFC'):

        table = cls._table_builder(tables, comparisons, columns)

        return cls(table=table, comparisons=comparisons,
                            columns=columns, scorekey=scorekey)

    @classmethod
    def from_dir(
            cls,
            results_dir,
            scorekey='LFC',
            index_name='Gene',
            comparisons:CompDict=None):


        (res, comps) = load_results(
            results_dir, LIMMACOLS,
            index_name=index_name
        )
        if comparisons is None:
            comparisons = comps
        res.index.name = index_name
        return cls(table=res, comparisons=comparisons, columns=LIMMACOLS, scorekey=scorekey)


if __name__ == '__main__':
    pass
    # import pickle
    # d = pathlib.Path('/mnt/m/tasks/NA327_Proteomics_UbPulldown/pickles2')
    # with open(d/'smolcomp.2.pickle', 'rb') as f:
    #     comp = pickle.load(f)
    # with open(d/'test_tables.1.pickle', 'rb') as f:
    #     tables = pickle.load(f)
    #
    # res = LimmaResults.build(tables=tables, comparisons=comp,)
    # print(res.comparisons[0].grimps)
