# import pathlib
# from jttools.data_wrangling import rename_columns
import typing

import pandas as pd


from bioscreen._imports import *
from bioscreen.classes.base import *
from bioscreen.classes.results import AnalysisResults
from bioscreen.classes.comparison import Comparison, CompDict
from statsmodels.stats.multitest import fdrcorrection
import numpy as np
from attrs import define

__all__ = ['GeneSetEnrichmentResults', 'GeneSetCollections', 'GeneSet', 'PadogResults']

GSCollectionName = str
GSName = str
GeneSet = set

class GeneSetCollections(dict[GSCollectionName, dict[GSName, GeneSet]]):

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


@define(kw_only=True)
class GeneSetEnrichmentResults(AnalysisResults):
    collections: GeneSetCollections

    @staticmethod
    def compres_from_dir(
            results_dir,
            columns:StatColumns,
            do_fdr=False,
            do_log10=True
    ) -> Tuple[CompsResultDF, CompDict]:
        """Load results from directory, assuming files names as
        {geneSetCollection}.{treat}.{ctrl}.csv. Concats the collection
        results into a single table with a Collection column."""

        results_dir = pathlib.Path(results_dir)

        # load and organise individual results tables
        res_by_collctn = {}
        comparisons = CompDict({})
        for fn in os.listdir(results_dir):
            coll, ctrl, treat, _ = fn.split('.')
            comp = Comparison(control=ctrl, test=treat)
            comparisons[comp.name] = comp
            tbl = pd.read_csv(results_dir / fn, index_col=0)

            columns.rename_df_columns(tbl, inplace=True)

            if do_fdr:
                tbl.loc[:, 'FDR'] = fdrcorrection(tbl.p)[1]
            if do_log10:
                tbl.loc[:, 'p10'] = tbl.p.apply(neglog10)
                tbl.loc[:, 'FDR10'] = tbl.FDR.apply(neglog10)

            # put the table into the structure
            if coll not in res_by_collctn:
                res_by_collctn[coll] = {}
            res_by_collctn[coll][comp.name] = tbl

        # go through dict of dict, creating multiindex DF and adding it to list
        tables = []
        for coll, res in res_by_collctn.items():
            resdf = pd.concat(res, axis='columns')
            resdf.loc[:, 'Collection'] = coll
            tables.append(resdf)

        # concat list of tables.
        results = CompsResultDF(
            pd.concat(tables, axis='index')
        )
        return (results, comparisons)

    def result_table(self, ctrl_or_comp, test=None, ) -> \
            pd.DataFrame:
        table = super().result_table(ctrl_or_comp, test)
        table.loc[:, 'Collection'] = self.table.Collection
        return table


@define(kw_only=True)
class PadogResults(GeneSetEnrichmentResults):

    _coltable = [{'original': 'Name', 'key': 'Name', 'label': 'Name', 'table': 'Name'},
                 {'original': 'ID', 'key': 'ID', 'label': 'ID', 'table': 'ID'},
                 {'original': 'Size', 'key': 'Size', 'label': 'Size', 'table': 'Size'},
                 {'original': 'meanAbsT0',
                  'key': 'MeanAbsT0',
                  'label': 'Unweighted T-score',
                  'table': 'T'},
                 {'original': 'padog0',
                  'key': 'Score',
                  'label': 'Weighted T-score',
                  'table': 'Weighted T'},
                 {'original': 'PmeanAbsT',
                  'key': 'PUnweighted',
                  'label': 'p-value (unweighted)',
                  'table': 'p (unweighted)'},
                 {'original': 'Ppadog', 'key': 'p', 'label': 'p-value', 'table': 'p'},
                 {'original': 'p10', 'key': 'p10', 'label': '-log10(p)', 'table': '-log10(p)'},
                 {'original': 'FDR', 'key': 'FDR', 'label': 'FDR', 'table': 'FDR'},
                 {'original': 'FDR10',
                  'key': 'FDR10',
                  'label': '-log10(FDR)',
                  'table': '-log10(FDR)'}]


    @classmethod
    def from_dirs(cls, results_dir:str, gsets_dir, ):
        #todo gene sets shouldn't be required

        stat_cols = StatColumns.from_records(
            cls._coltable,
        )

        results, comparisons = GeneSetEnrichmentResults.compres_from_dir(
            results_dir, stat_cols, do_fdr=True, do_log10=True
        )

        # lower case names for the gene set terms
        nicenames = results.xs('Name', 'columns', level=1).iloc[:, 0].apply(
            GeneSetCollections.set_name_from_msig
        )

        for comp in comparisons:
            results.loc[:, (comp, 'Name')] = nicenames

        gcollections = GeneSetCollections.from_tsl_dir(gsets_dir)



        return cls(
            table=results,
            columns=stat_cols,
            collections=gcollections,
            comparisons=comparisons,
            scorekey='Score',
        )



#from bioscreen.classes.experiment import get_replicates_of_comparison
from jttools.data_wrangling import index_of_true

def get_replicates_of_comparison(sample_details, comparison:Comparison) -> np.ndarray[Sample]:
    """Control and test samples used in a comparison. Controls first."""
    m = sample_details.SampleGroup.isin([comparison.control, comparison.test])
    reps = sample_details.loc[m, 'Sample'].values
    return reps


# take a table where the indexes represent duplicated gene symbols (or whatever),
# e.g. protein isoforms, return a table that has the "best" of the duplicates kept,
# best being determined by a stat table.
class CountRemapperWithDuplicates:
    """Get count tables of the best counts in a direction (up down) for
    passing to gene enrichment algo when the indexes of your counts don't
    map uniquely to gene symbols (e.g. with proteomics results)."""
    def __init__(self, counts:pd.DataFrame, sample_details:pd.DataFrame, symbolmapper:pd.Series,
                 result_table:pd.DataFrame, comparisons:CompDict, scorecol='LFC', sigcol='p10'):
        self.counts = counts
        self.sample_details = sample_details
        self.symbolmapper = symbolmapper
        self.result_table = result_table
        self.comparisons = comparisons
        self.scorecol = scorecol
        self.sigcol = sigcol

    def _get_bad_duplicates_by_sig(self, compk:str, direction: Literal['up', 'down'], ):
        """Get list of indexes to drop to keep the strongest hits in
        both directions for all in fdr and lfc tables.

        idmap should map from index values in fdr10/lfc to new values
        containing the duplicates that you want to select from."""

        # we'll get the order right by having things going
        # in the wrong direction be neg
        lfc = self.result_table[compk][self.scorecol].copy()
        sig = self.result_table[compk][self.sigcol].copy()

        if direction == 'up':
            wrongdir = lfc < 0
        else:
            wrongdir = lfc > 0

        sig[wrongdir] *= -1
        f = sig.sort_values()

        droppers = self.symbolmapper.loc[f.index].duplicated(keep='last')

        return index_of_true(droppers)

    def _get_count_better_duplicates(
            self, samples:Collection[Sample], worse_in_direction:Collection[str]
    ):
        """Worse in direction is a list of indexes that should be dropped for the direction we're
        looking at"""
        newcnt = self.counts.loc[:, samples].drop(worse_in_direction)
        newcnt.index = newcnt.index.map(self.symbolmapper.loc[newcnt.index])
        newcnt = newcnt.loc[~newcnt.isna().all(1)]
        newcnt = newcnt.loc[~newcnt.index.isna()]
        return newcnt

    def remap_count_index_keeping_better_duplicates(
            self,
            direction:Literal['up', 'down'],
            compk) -> pd.DataFrame:
        reps = get_replicates_of_comparison(
            self.sample_details,
            self.comparisons[compk]
        )
        worse_indexes = self._get_bad_duplicates_by_sig(
            compk, direction
        )
        return self._get_count_better_duplicates(reps, worse_indexes)

    def iter(self) -> Tuple[Literal['up', 'down'], Comparison, pd.Series]:
        """Iterate through counts, returning string indicating direction, comparison
        object, and new count."""
        for direction in ('up', 'down',):
            direction:Literal['up', 'down']
            for compk, comp in self.comparisons.items():
                t = self.remap_count_index_keeping_better_duplicates(direction, compk)
                yield direction, comp, t

    def __iter__(self):
        return self.iter()


    # used this for testing (manual check of first result)
    #   obviously needs appropriate input data. Only checks up direction on the first
    #   comparison, but I did switch them in iter and checked down.
    # k = 'MTX117575_48h_100nM'
    # xp = 'Whole'
    #
    # testres = {k: limma_tables[xp][k].dropna().head(20)}
    # testmapper = pd.Series(
    #     ['a'] * 5 + ['b'] * 5 + ['c'] * 5 + ['d'] * 5,
    #     index=t.index
    # )
    #
    # testcount = counts[xp]['med_nan'].reindex(t.index)
    #
    # countremapper = CountRemapperWithDuplicates(
    #     testcount,
    #     sample_details,
    #     testmapper,
    #     testres,
    #     comparisons
    # )
    # x = testres[k].loc[:, ['LFC', 'p10']]
    # i = x.index.map(testmapper)
    # x.index = i + '    ' + x.index
    # display(x)
    #
    # for d, k, c in countremapper.iter():
    #     print(d, k)
    #     display(c.head())
    #     display(testcount.loc[:, c.columns])
    #     break

if __name__ == '__main__':
    print('testing padog reading')
    fnp = '/mnt/m/tasks/MTX639_TAC_MoA_RNAseq/padog/test_first/t3_htseq/'
    fngs = '/mnt/m/data_collections/MSigDB_parsed/mouse_symbols'

    t = PadogResults.from_dirs(fnp, fngs)
    print(t.table.columns.levels[1])
