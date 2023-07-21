from bioscreen.classes.base import *
from bioscreen.classes.results import AnalysisResults
from bioscreen.classes.comparison import Comparison, CompList
from statsmodels.stats.multitest import fdrcorrection

from attrs import define


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
            -> Tuple[ComparisonsResults, CompList]:
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
            rename_columns(tbl, columns.original, inplace=True)

            # put the table into the structure
            if coll not in res_by_collctn:
                res_by_collctn[coll] = {}
            res_by_collctn[coll][comp.str] = tbl
        comparisons = CompList(list(comparisons))
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

    _coltable = {
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
