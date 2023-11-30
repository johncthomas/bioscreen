import gseapy
from bioscreen._imports import *
from jttools.data_wrangling import write_stats_workbook

def format_gseapy_res2d(table):
    """Rename columns, split the "tag" column into LeadingLen and Size."""
    t = table.copy()
    for col, i in  (('LeadingLen', 0), ('Size', 1)):
        t.loc[:, col] = t['Tag %'].apply(lambda x: x.split('/')[i])
    t.drop('Tag %', inplace=True, axis='columns')

    final_cols = 'Method Term ES NES pNOM FDR FWER LeadingLen Size PercLeadingAll LeadGenes'.split()
    ccols = ['Name', 'Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val', 'LeadingLen', 'Size',
           'Gene %', 'Lead_genes', ]
    t.columns = df_rename_columns(t, dict(zip(ccols, final_cols)))

    return t.reindex(columns=final_cols, )

def score_signed_p10(res:pd.DataFrame) -> pd.Series:
    """Get log10(p) with negative sign when LFC < 0."""
    neg = res.LFC < 0

    score = res.p10.copy()
    score[neg] *= -1

    score = score.sort_values()
    return score

def gsea_prerank_analysis(
        scores:pd.Series,
        gene_set_collections:dict[str, dict[str, list[str]]]
) -> dict[str, pd.DataFrame]:
    """Do pygsea.prerank analysis on scores, returning one table concating all results
    for all supplied GS collections."""


    gcr = AttrMapAC()

    for gscname, geneset in gene_set_collections.items():
        res = gseapy.prerank(rnk=scores, gene_sets=geneset, )
        tbl = format_gseapy_res2d(res.res2d)
        tbl = tbl.drop('Method', axis='columns')
        tbl.insert(0, 'GeneSetCollection', gscname)

        gcr[gscname] = tbl

    gsea_table = pd.concat(gcr.values(), axis='index').sort_values('FDR', ascending=True)
    return gsea_table

from jttools.excel import conditional_format_definitions
sigfmtxl = conditional_format_definitions.significance()
scrfmtxl = conditional_format_definitions.score()

def write_gsea_tables_to_xlsx(tables:dict[str, pd.DataFrame], outfn:Pathy):
    wrapped = dict(text_wrap=True, num_format='@')
    wb = write_stats_workbook(
        outfn, tables,
        conditional_formats=dict(FDR=sigfmtxl, FWER=sigfmtxl, NES=scrfmtxl,),
        close_workbook = False,
        other_formats={'Term': wrapped},
    )

    for sn, table in tables.items():
        sheet = wb.get_worksheet_by_name(sn)
        # plus 1 cus index
        if isinstance(table.index, pd.RangeIndex):
            i_off = 0
        else:
            i_off = 1
        term_i = list(table.columns).index('Term') + i_off
        sheet.set_column_pixels(term_i, term_i, width=500)
    wb.close()

def _test_run():
    import pickle
    with open('/mnt/m/tasks/NA327_Proteomics_UbPulldown/pickles_combined/limres_kggnorm.1.pickle', 'rb') as f:
        limres_kggnorm = pickle.load(f)
    with open('/mnt/m/tasks/NA327_Proteomics_UbPulldown/pickles_combined/prot_count_info.4.pickle', 'rb') as f:
        info = pickle.load(f)
    with open('/mnt/m/tasks/NA327_Proteomics_UbPulldown/pickles_combined/genesetcolls.2.pickle', 'rb') as f:
        genesetcolls = pickle.load(f)

    scores = {}
    for treat in limres_kggnorm.comparisons:
        res = limres_kggnorm.table[treat]
        score = score_signed_p10(res).dropna()

        # get the highest scoring per gene symbol
        ifo = info.KGG.reindex(score.index)

        dupes = ifo.GeneSymbol.duplicated(keep='last')
        s = score.loc[~dupes]
        s.index = ifo[~dupes].GeneSymbol

        scores[treat] = s
    gsea_res_kggnorm_p3 = {}
    for k, s in scores.items():
        gsea_res_kggnorm_p3[k] = gsea_prerank_analysis(s, genesetcolls)
    
    write_gsea_tables_to_xlsx(gsea_res_kggnorm_p3, '/mnt/m/tasks/NA327_Proteomics_UbPulldown/data/GSEA/kGG_wholenorm.fdrscore.P3.1.xlsx')

if __name__ == '__main__':
    pass
    # print('running _test_run')
    # _test_run()
