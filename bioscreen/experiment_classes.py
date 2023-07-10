import os
import pathlib
import collections
import itertools
import typing
import logging
from copy import copy, deepcopy

from bioscreen import validations

logger = logging.getLogger(__name__)
logger.setLevel('INFO')

from attrdictionary import AttrMap

import attrs
import pandas as pd

from attrs import define, Factory, field

from typing import Optional, List, Tuple, Dict, Callable, Any

from statsmodels.stats.multitest import fdrcorrection

import xlsxwriter

from bioscreen.utils import (
    is_numeric,
    neglog10,
    AMap
)

from bioscreen.validations import (
    validate_count_df,
    validate_screen_input,
    validate_comparisons_table,
    validate_sample_details,
    validate_cols
)

def rename_df_columns(tbl:pd.DataFrame, newcols:dict):
    # rename columns
    colmapper = {c: c for c in tbl.columns}
    colmapper.update(newcols)
    tbl.columns = tbl.columns.map(colmapper)
    return tbl


StrOptional = Optional[str]
CollectionOptional = Optional[typing.Collection]

# COLS = AttrMap(dict(
#     comparisons=dict(
#         test='Test',
#         ctrl='Control',
#         paired='Paired',
#     )
# ))


class StatColumns:
    """Column mappings, and primary score/p/fdr.
    Args and properties:
        original: original->good column names
        short: good->short label for interactive tables or to go in legends
        long: good->long label for Excel results or chart axis labels
        score, p, fdr: Good name of primary score/p/FDR column
        p10, fdr10: -log10(p|fdr10)
        results_cols: default written to written results tables (default: all)
        interactive_cols: default to put on interactive chart
            (default: results_cols)

        Set p, p10, fdr, fdr10 to None if not using.
        """
    def __init__(
            self,
            original:dict,
            short:dict,
            long:dict,
            score:str,
            p='P',
            p10='P10',
            fdr='FDR',
            fdr10='FDR10',
            results_cols:CollectionOptional=None,
            interactive_cols:CollectionOptional=None,
    ):

        self.columns = cols = list(original.values())
        assert all(c in cols for c in (fdr, fdr10, p, p10) if c is not None)
        self.original = original
        self.short = short
        self.long = long
        self.score = score
        self.p = p
        self.p10 = p10
        self.fdr = fdr
        self.fdr10 = fdr10

        if results_cols is None:
            results_cols = cols
        self.results_cols = results_cols
        if interactive_cols is None:
            interactive_cols = results_cols
        self.interactive_cols = interactive_cols

    @staticmethod
    def kwargs_from_dict(df: typing.Union[pd.DataFrame, dict[str, list]]):
        """KW args for StatColumns init from a  table/df with columns/keys:
        original, good, short, long"""
        kwargs = dict(
            original=dict(zip(df['original'], df['good'])),
        )

        for k in ('long', 'short'):
            kwargs[k] = dict(zip(df['good'], df[k]))

        return kwargs

    @classmethod
    def from_dict(cls, label_table:typing.Union[pd.DataFrame, dict[str, list]], **kwargs):
        """Init from a table/df with columns/keys:
        "original", "good", "short", "long"

        kwargs passed to constructor"""

        return cls(**cls.kwargs_from_dict(label_table), **kwargs)


Sample = str

# frozen so as can be hashed, also should never change
@attrs.frozen
class Comparison:
    control: Sample
    test: Sample

    def __str__(self):
        return f"{str(self.control)}-{str(self.test)}"
    @property
    def label(self  ):
        return f"{str(self.control)} ➤ {str(self.test)}"
    @property
    def str(self):
        return self.__str__()
    @property
    def strDot(self):
        # for Excel table name at least.
        return self.str.replace('-', '.')



# todo probably give up on this, just have Result.table = pd.DataFrame
class ComparisonsResults(pd.DataFrame):
    # def comparison(self, ctrl, test) -> pd.DataFrame:
    #     cmp = Comparison(ctrl, test)
    #     return self[cmp]

    @property
    def stat_names(self):
        return self.columns.levels[1]

    @property
    def comparisons(self):
        return self.columns.levels[0]


def comp_results_from_dir(
        results_dir, fn_to_comp: Callable,
        columns:StatColumns=None, sep=',',
    ) -> tuple[ComparisonsResults, list[Any]]:
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
            tbl = rename_df_columns(tbl, columns.original)
        results[str(comp)] = tbl
        comparisons.append( comp)
    return (ComparisonsResults(pd.concat(results, axis='columns')), comparisons)


@define
class AnalysisResults:
    table: ComparisonsResults
    columns: StatColumns
    comparisons: List[Comparison]
    # todo
    #   comp_table

    @property
    def score_table(self):
        return self.get_stat_table(self.columns.score)

    @property
    def fdr_table(self):
        return self.get_stat_table(self.columns.fdr)

    @property
    def p_table(self):
        return self.get_stat_table(self.columns.p)

    @property
    def fdr10_table(self):
        #return self.fdr_table.apply(neglog10)
        return self.get_stat_table(self.columns.fdr10)

    @property
    def p10_table(self):
        return self.get_stat_table(self.columns.p10)

    def get_stat_table(self, key) -> pd.DataFrame:
        return self.table.xs(key, level=1, axis=1)

    def result_table(self, ctrl_or_comp:typing.Union[str, Comparison],
                     test:str=None, ) -> pd.DataFrame:
        """pass a comp-string, Comparison or (samplestr, samplestr)
        and get a single table with results for that comp."""
        if test is not None:
            comp = str(Comparison(ctrl_or_comp, test))
        else:
            comp = str(ctrl_or_comp)

        return self.table[comp].copy()


    def write_comp_results_to_excel(self, filename, columns:Optional[Dict[str,str]]=None, **xlsx_table_opts):
        """Write comparisons table as excel workbook with one table
        per comp.

        xlx_table_opts added to dict and passed to worksheet.add_table,
        see https://xlsxwriter.readthedocs.io/working_with_tables.html"""
        workbook = xlsxwriter.Workbook(filename)
        # todo formats: specified in ResultsColumns i guess, since we need to define xlsxwriter obj
        #   that will be shared across columns

        for comp in self.comparisons:
            worksheet = workbook.add_worksheet(name=comp.label)
            tab = self.result_table(comp.str)

            if (tab.index.name is not None) and (tab.index.name not in tab.columns.values):
                tab = tab.reset_index()
            #xlsxwriter notation: (firstRow, firstCol, lastRow, lastCol)

            #can end up with missing rows
            # set the sig columns 1, and the rest to zero
            if tab.isna().any().any():
                logger.info(f"NAs found in {comp.str}, n={tab['P'].isna().sum()}, replacing with 1 or zero.")
            for k in ('P', 'FDR'):
                tab.loc[tab[k].isna(), k] = 1
            tab.fillna(0, inplace=True)


            # select and rename columns
            tab = tab.reindex(columns=list(columns.keys()))
            tab.sort_values([self.columns.p, self.columns.score], ascending=[True, False], inplace=True)
            tab.columns = tab.columns.map(columns)

            nrows, ncols = tab.shape
            opts = dict(
                data=tab.values,
                columns=[{'header':c} for c in tab.columns],
                name=comp.strDot, # dash not allowed

            )
            opts.update(xlsx_table_opts)
            worksheet.add_table(
                0, 0, nrows, ncols-1,
                opts
            )

        workbook.close()


GSCollectionName = str
GSName = str
GeneSet = set

class GeneSetCollections(Dict[GSCollectionName, Dict[GSName, GeneSet]]):

    @property
    def collections_map(self):
        coll_gs = {}
        for col, sets in self.items():
            for gs in sets:
                coll_gs[gs] = col
        return coll_gs

    def genes_by_setname(self, gset:GSName) -> GeneSet:
        coll = self.collections_map[gset]
        return self[coll][gset]

    def genes_in_collection(self, collection:GSCollectionName):
        genes = GeneSet()
        for genessets in self[collection].values():
            genes.update(genessets)
        return genes

    def to_tidy_df(self):
        """Each row a unique combination of gene and set name.
        Uses SetRank column names:
            geneID termID termName description dbName
        """
        cols = 'geneID termID termName description dbName'.split()
        dat = {c: [] for c in cols}
        for coll, gsets in self.items():
            for setid, genes in gsets.items():
                for gn in genes:
                    dat['geneID'].append(gn)
                    dat['termID'].append(setid)
                    dat['termName'].append(self.set_name_from_msig(setid))
                    dat['description'] = ''
                    dat['dbName'].append(coll)
        # return pd.DataFrame(dat)
        msig_tidy = pd.DataFrame(dat)
        return msig_tidy

    @staticmethod
    def set_name_from_msig(n):
        """e.g. 'GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE'
        to 'mitochondrial genome maintenance'
        """
        n = n.split('_', 1)[1]
        n = n.lower().replace('_', ' ')
        return n

    @staticmethod
    def sets_from_tsl(filename):
        gene_sets = {}

        with open(filename) as f:
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                setname = spline[0]
                genes = set(spline[1:])
                #gs = GeneSet(setname, genes, collection)
                gene_sets[setname] = genes
        return gene_sets

    @classmethod
    def from_tsl_dir(cls, directory):
        """Assumes filnames are '{collection name}.tsl', first value
        is gene set names, and all after are genes. """
        gsets_dir = pathlib.Path(directory)
        gcollections = {}
        for fn in os.listdir(gsets_dir):
            if fn.endswith('.tsl'):
                collname = fn.__str__().split('.')[-2]
                gcollections[collname] = GeneSetCollections.sets_from_tsl(
                    gsets_dir/fn
                )
        return cls(gcollections)



# @define
# class GeneSet:
#     name: str
#     genes: GeneList
#     collection: Optional[str]
#     # todo get info from the .info.csv files


@define
class GeneSetEnrichmentResults(AnalysisResults):
    collections: GeneSetCollections


    @staticmethod
    def compres_from_dir(results_dir, columns:StatColumns, doFDR=False) \
            -> Tuple[ComparisonsResults, List[Comparison]]:
        """Load results from directory, assuming files names as
        {geneSetCollection}.{treat}.{ctrl}.csv. Concats the collection
        results into a single table with a Collection column."""

        print('Assuming file format "{geneSetCollection}.{ctrl}.{test}.csv"')
        results_dir = pathlib.Path(results_dir)

        # load and organise individual results tables
        res_by_collctn = {}
        comparisons = set()
        for fn in os.listdir(results_dir):
            coll, ctrl, treat, _ = fn.split('.')
            comp = Comparison(control=ctrl, test=treat)
            comparisons.add(comp)
            tbl = pd.read_csv(results_dir / fn, index_col=0)
            if doFDR:
                tbl.loc[:, 'FDR'] = fdrcorrection(tbl.P)[1]
            rename_df_columns(tbl, columns.original)

            # put the table into the structure
            if coll not in res_by_collctn:
                res_by_collctn[coll] = {}
            res_by_collctn[coll][comp.str] = tbl
        comparisons = list(comparisons)
        # go through dict of dict, creating multiindex DF and adding it to list
        tables = []
        for coll, res in res_by_collctn.items():
            resdf = pd.concat(res, axis='columns')
            resdf.loc[:, 'Collection'] = coll
            tables.append(resdf)

        # concat list of tables.
        results = ComparisonsResults(
            pd.concat(tables, axis='index')
        )

        return (results, comparisons)

    def result_table(self, ctrl_or_comp, test=None, ) -> \
            pd.DataFrame:
        table = super().result_table(ctrl_or_comp, test)
        table.loc[:, 'Collection'] = self.table.Collection
        return table


@define
class PadogResults(GeneSetEnrichmentResults):

    _coltable = label_table={
        'original': ['Name', 'ID', 'Size', 'meanAbsT0', 'padog0', 'PmeanAbsT', 'Ppadog', 'P10', 'FDR', 'FDR10'],
         'good': ['Name', 'ID', 'Size', 'MeanAbsT0', 'Score', 'PUnweighted', 'P', 'P10', 'FDR', 'FDR10'],
         'long': ['Name', 'ID', 'Size', 'Unweighted T-score', 'Weighted T-score', 'p-value (unweighted)',
                  'p-value', '-log10(p)', 'FDR', '-log10(FDR)'],
         'short': ['Name', 'ID', 'Size', 'T', 'Weighted T', 'p (unweighted)', 'p', '-log10(p)', 'FDR',
                   '-log10(FDR)']
    }
    _results_cols = ['Name', 'Score', 'P', 'FDR', 'Size', 'ID',  ]
    _score='Score'


    @classmethod
    def from_dirs(cls, results_dir:str, gsets_dir):
        #todo gene sets shouldn't be required
        stat_cols = StatColumns.from_dict(
            cls._coltable,
            results_cols=cls._results_cols,
            score=cls._score,
        )

        results, comparisons = GeneSetEnrichmentResults.compres_from_dir(results_dir, stat_cols)

        nicenames = results.xs('Name', 'columns', level=1).iloc[:, 0].apply(
            GeneSetCollections.set_name_from_msig
        )
        for comp in comparisons:
            results.loc[:, (comp.str, 'Name')] = nicenames

        gcollections = GeneSetCollections.from_tsl_dir(gsets_dir)



        return cls(
            table=results,
            columns=stat_cols,
            collections=gcollections,
            comparisons=comparisons,
        )


@define
class DGEResults(AnalysisResults):
    @staticmethod
    def load_results(results_dir, columns, index_name='Gene') -> tuple[ComparisonsResults, list[Any]]:


        def fn_to_comp(fn):
            logger.info(f'DGEResults: fn_to_comp({fn})')
            test, ctrl, _ = fn.split('.')

            comp = Comparison(ctrl, test)
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

@define
class ScreenExperiment:
    name: str
    version: str
    counts: typing.Mapping[str, pd.DataFrame]
    sample_details: pd.DataFrame
    group_details: pd.DataFrame
    comparisons: pd.DataFrame
    results: typing.Mapping[str, AnalysisResults] = Factory(AttrMap)

    def __attrs_post_init__(self):
        validate_comparisons_table(self.comparisons)
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

def _test_padog(
        rdir='/mnt/m/tasks/MTX639_TAC_MoA_RNAseq/padog/ctrl_first/t3_htseq/',
        gsdir='/mnt/m/data_collections/MSigDB_parsed/mouse_symbols'):

    padres = PadogResults.from_dirs(
        results_dir=rdir,
        gsets_dir=gsdir
    )

    print(padres.table.head())

def _test_limma(dir='/mnt/m/tasks/MTX639_TAC_MoA_RNAseq/limma/t3_star_htseq'):
    r = LimmaResults.from_dir(dir)
    print(r.table.head())




if __name__ == '__main__':
    _test_padog()
    _test_limma()

