from bioscreen.classes.base import *
from bioscreen.classes.comparison import CompList, Comparison
from attrs import define
import xlsxwriter



def comp_results_from_dir(
        results_dir, fn_to_comp: Callable,
        columns:StatColumns=None, sep=',',
    ) -> tuple[ComparisonsResults, CompList]:
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
            tbl.columns = rename_columns(tbl, columns.original_to_good())
        results[str(comp)] = tbl
        comparisons.append( comp)
    comparisons = CompList(comparisons)
    return (pd.concat(results, axis='columns'), comparisons)


@define
class AnalysisResults:
    table: ComparisonsResults
    columns: StatColumns
    comparisons: CompList
    score = ''

    def get_stat_table(self, key) -> pd.DataFrame:
        return self.table.xs(key, level=1, axis=1)

    @property
    def score_table(self):
        return self.get_stat_table(self.score)

    @property
    def fdr_table(self):
        return self.get_stat_table('FDR')

    @property
    def p_table(self):
        return self.get_stat_table('p')

    @property
    def fdr10_table(self):
        #return self.fdr_table.apply(neglog10)
        return self.get_stat_table('FDR10')

    @property
    def p10_table(self):
        return self.get_stat_table('P10')



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


    def write_comp_results_to_excel(self, filename, columns:Optional[Dict[str,str]]=None, **xlsx_table_opts):
        """Write comparisons table as excel workbook with one table
        per comp.

        **xlx_table_opts passed to worksheet.add_table,
        see https://xlsxwriter.readthedocs.io/working_with_tables.html"""
        # todo formats: specified in ResultsColumns i guess, since we need to define xlsxwriter obj
        #   that will be shared across columns

        workbook = xlsxwriter.Workbook(filename)

        for comp in self.comparisons:
            worksheet = workbook.add_worksheet(name=comp.arrow_str())
            tab = self.result_table(comp.str())

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
                name=comp.str(joiner='.'), # dash not allowed

            )
            opts.update(xlsx_table_opts)
            worksheet.add_table(
                0, 0, nrows, ncols-1,
                opts
            )

        workbook.close()
