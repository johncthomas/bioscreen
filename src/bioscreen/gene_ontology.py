
"""The aim is to automate the functionality of the PANTHER gene enrichment
web interface. User provides a list of genes and gets a table giving the
enriched GO terms in a hierarchical format similar to the web results.

Currently on supports GOslim biological process because that's the most
interesting, but should be easy to expand."""
import typing

# API reference
# http://pantherdb.org/services/openAPISpec.jsp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import requests, pickle
from typing import Collection, List, Tuple, Dict, Any, Literal

from pkg_resources import resource_filename
#gs_fn = resource_filename(__name__, "data/PANTHERGOslim.hierarchy.pickle")


def parse_go_hierarchy(obo_fn, output_pickle_fn=None):
    """Using an OBO file, parse GO term hierarchy.

    Format as a mapping {parent_id1:[child_id1, child_id2, ...], ...]
    IDs are integers...

    Check the output too, this is some fragile code."""

    goslim_hierarchy = {}

    parent_ids = []
    daughter_id = 0
    parsing_active = False # True when '[Term]\n' encountered
    with open(obo_fn) as f:
        for line in f:
            if line == '[Term]\n':
                parsing_active = True
                # record the previous
                if daughter_id:
                    goslim_hierarchy[daughter_id] = parent_ids

                parent_ids = []
                line = next(f)

                # go
                if not line.startswith('id: '):
                    raise RuntimeError('OBO file not in the expected order, this data format is awful.')
                daughter_id = line.strip().replace('id: ', '')

            if line.startswith('is_a:') and parsing_active:
                pid = line.replace('is_a: ', '').split(' !')[0]
                parent_ids.append(pid)

            # These were all at the end when I checked, but I_guess there's no guarantee
            #    it breaks if you keep going here, as the 'is_a: ` statements stop being
            #    GO:numbers
            if line == '[Typedef]\n':
                parsing_active = False

    # # check there weren't extra
    # n = 0
    # for l in open(obo_fn):
    #     if l.startswith('[Term]'):
    #         n += 1
    #

    goslim_hierarchy[daughter_id] = parent_ids

    if output_pickle_fn is not None:

        with open(output_pickle_fn, 'wb') as f:
            pickle.dump(goslim_hierarchy, f)
    return goslim_hierarchy


def intgo(s:str) -> int:
    """convert GO:000nnnn to an integer"""
    # I don't think I had a good reason for doing this, so I've removed it's use.
    return int(s[3:])

def strgo(n:int) -> str:
    """inverse of intgo, supply int, get 'GO:<num>' with zeros padding"""
    return f"GO:{n:07}"


def load_hierarchy(gs_fn:typing.Union[dict, str]):
    if type(gs_fn) is dict:
        return gs_fn

    if gs_fn.endswith(".obo"):
        return parse_go_hierarchy(gs_fn)

    with open(gs_fn, 'rb') as f:
        goslim_hierarchy = pickle.load(f)
    return goslim_hierarchy


def query_panther(gene_list:Collection[str], reference_list:Collection[str]=None,
                  annot_type:Literal["bp", "mf", "pc", "cc", "cl", "pp", "rp"]='bp') -> requests.Response:
    """Get enrichment results from PANTHER service.

    Args:
        gene_list: List of genes to be tested for enrichment of annotations.
            Ideally as IDs, not names/symbols.
        reference_list: Reference genes list for enrichment.
        annot_type: Which PANTHER annotation set to test against?
            bp - GO Slim biological process (default)
            mf - GO Slim molecular function
            cc or cl - GO Slim cellular location
            pc - PANTHER protein class
            pp - PANTHER pathway
            rp - Reactome pathway

    Use `parse_results` to get tables etc. from returned `Response`
    """
    urlfmt = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList={geneList}&organism=9606&annotDataSet={annot}&enrichmentTestType=FISHER&correction=FDR"
    geneList = ','.join(gene_list)
    if annot_type == 'bp':
        annot= "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP"
    elif annot_type == 'mf':
        annot = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF"
    elif annot_type in ('cc', 'cl'):
        annot = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC"
    elif annot_type == 'pc':
        annot = "ANNOT_TYPE_ID_PANTHER_PC"
    elif annot_type == 'pp':
        annot = "ANNOT_TYPE_ID_PANTHER_PATHWAY"
    elif annot_type == 'rp':
        annot = "ANNOT_TYPE_ID_REACTOME_PATHWAY"
    else:
        raise RuntimeError(f"Annotation type of '{annot_type}' not recognised, see"
                           " documentation (i.e. this function's docstring) for details.")
    query_url = urlfmt.format(annot=annot, geneList=geneList)
    if reference_list is not None:
        query_url = query_url + '&refOrganism=9606&refInputList=' + ','.join(reference_list)
    res = requests.get(query_url)
    if res.status_code != 200:
        raise requests.ConnectionError(f"Query failed, status code: {res.status_code}. Query URL:\n{query_url}")
    return res

def parse_panther_json_to_df(pnthr_json, fdr_threshold=1):
    sig_res = []
    try:
        err = pnthr_json['search']['error']
        if err:
            print('Error?', err)
    except KeyError:
        pass
    for d in pnthr_json['overrepresentation']['group']:
        if d['fdr'] > fdr_threshold:
            continue
        if d['term']['label'] == 'UNCLASSIFIED':
            continue
        d['GO_ID'] = d['term']['id']
        d['GO_label'] = d['term']['label']
        del d['term']

        sig_res.append(d)
    sig_go_df = pd.DataFrame(sig_res)
    sig_go_df.index = sig_go_df.GO_ID#.apply(intgo)
    # do this once and never need to worry about it again?
    sig_go_df = sig_go_df.sort_values('pValue')

    return sig_go_df

def _term_hierarchy_from_df(go_results:pd.DataFrame,
                           hierarchy_fn:str,
                           goid_column='GO_ID') -> List[Dict[int, Any]]:
    """Backwards compat, silly to pass the table when I only need one column.

    Get the structure of the hierarchical table from results paresed with
    `parse_panther_json_df()`.

    Returns a recursive list of dicts containing lists of dicts; keys are more
    specific GOID paired with lists of less specific GOID

    [{int_GOID:[next_level_GOID, ...]}, {...}]

    Pass to `transcribe_structure()` to get the actual table with the results."""

    return term_hierarchy(go_results[goid_column], hierarchy_fn)



def term_hierarchy(go_series:pd.Series,
                   hierarchy_fn:str, ) -> List[Dict[int, Any]]:
    """Get the structure of the hierarchical from selected GO terms in go_series.

    Returns a recursive list of dicts containing lists of dicts; keys are more
    specific GOID paired with lists of less specific GOID

    [{int_GOID:[next_level_GOID, ...]}, {...}]

    Pass to `transcribe_structure()` to get the actual table with the results."""

    goslim_hierarchy = load_hierarchy(hierarchy_fn)

    present_ids = go_series#.apply(intgo)

    cnxns = go_series.apply(lambda x: [i for i in goslim_hierarchy[x] if i in present_ids.values])
    cnxns.index = go_series#.apply(intgo)

    # We need to reverse the direction of these associations, maybe should do that
    #   in the original goslim_hierachy, but haven't
    rcnx = pd.Series(index=present_ids.values)
    rcnx = rcnx.apply(lambda x: [])

    for gid, lst in cnxns.items():
        for g in lst:
            rcnx[g].append(gid)

    # Get the most specific term, these have no other term pointing to them
    top = rcnx.loc[rcnx.apply(len) == 0]

    def get_next_one(go):
        node = {go: []}
        for hit in cnxns[go]:
            next_hit = get_next_one(hit)
            node[go].append(next_hit)
        return node

    structure = [get_next_one(g) if get_next_one(g) is not None else {g:[]} for g in top.index]

    return structure


def transcribe_structure(
        structure:List[Dict[int, Any]],
        go_results:pd.DataFrame,
        label_col:str,
) -> Tuple[pd.DataFrame, List[str]]:

    """Take a structure generated by `term_hierarchy_from_df` and return TWO tables
    first a dataframe, and second a list of strings for printing.

    Label col should probably indicate the Term name."""

    # will be modified by iter_struct below
    out_table = [go_results.columns.tolist() + ['Level']]
    str_table = ['GO label hierarchy']

    # parse the string row for str_Table
    def str_level(goid, level=0):
        s = go_results.loc[goid, label_col]
        if level == 0:
            prefix = ''
        else:
            prefix = ' ' + '+ ' * level
        return prefix + s

    # go through the structure, generating rows and strings for the returned tables
    def iter_struct(d, level=0):
        for gid, next_levels in d.items():
            s = str_level(gid, level)
            if not s in str_table:
                str_table.append(s)
            row_vals = list(go_results.loc[gid].values)
            row_vals.append(level)

            out_table.append(row_vals)
            for next_d in next_levels:
                iter_struct(next_d, level=level + 1)

    for rows in structure:
        iter_struct(rows)
    return pd.DataFrame(out_table[1:], columns=out_table[0]), str_table


def parse_results(res_json, hierarchy_fn:str, fdr_threshold=0.1):
    """res_json, the JSON results from pantherdb.org.

    returns dict with keys:
        'json': full results JSON obtained from pantherdb.org,
        'full_table': the JSON in the form of a dataFrame,
        'table': a hierarchical table, excluding > fdr_threshold,
        'string': The table in the form of a string, for printing."""
    enrichment_results = parse_panther_json_to_df(res_json)
    sig_results = enrichment_results.loc[enrichment_results.fdr < fdr_threshold]
    hierarchy = _term_hierarchy_from_df(sig_results, hierarchy_fn)
    hierarch_table, string_table = transcribe_structure(hierarchy, sig_results)

    return {
        'json': res_json,
        'full_table': enrichment_results,
        'table': hierarch_table,
        'string': '\n'.join(string_table),
    }


def query_panther_and_parse_results(gene_list:Collection[str],
                                    reference_list:Collection[str]=None,
                                    fdr_threshold=0.1) -> Dict:
    """Query PANTHERdb using a gene list and return parsed results.

     Returns dict with keys:
        'json': full results JSON obtained from pantherdb.org,
        'full_table': the JSON in the form of a dataFrame,
        'table': a hierarchical table, excluding > fdr_threshold,
        'string': Terms as a multiline string, for printing.

    If you wish to refilter the results with a different fdr_threshold, pass
    the JSON to `parse_results` to avoid requerying pantherdb.org. parse_results
    Has the same return value."""
    response = query_panther(gene_list, reference_list=reference_list)
    return parse_results(response.json(), fdr_threshold=fdr_threshold)

def _test(genes = None, ):
    if genes is None:
        genes = 'HGNC:18258,HGNC:9108,HGNC:14052,HGNC:823,HGNC:1733,HGNC:29931,HGNC:20663,HGNC:1493,HGNC:15910,HGNC:3259,HGNC:16930,HGNC:9540,HGNC:2326,HGNC:7529,HGNC:848,HGNC:16684,HGNC:11957,HGNC:9534,HGNC:10411,HGNC:14048,HGNC:6215,HGNC:7530,HGNC:8648,HGNC:33998,HGNC:9788,HGNC:28762,HGNC:8982,HGNC:1020,HGNC:20220,HGNC:7910,HGNC:14044,HGNC:15934,HGNC:468,HGNC:838,HGNC:9535,HGNC:8014,HGNC:10378,HGNC:11617,HGNC:12582,HGNC:3433,HGNC:10241,HGNC:2263,HGNC:10389,HGNC:2231,HGNC:10364,HGNC:14161,HGNC:16235,HGNC:10530,HGNC:8975,HGNC:7643'.split(',')
    results = query_panther_and_parse_results(genes)
    for k, v in results.items():
        if k == 'json':

            continue
        print(k)
        print(v)
        print('********************\n')

def plot_enrichment_table(tab:pd.DataFrame, filename=None) -> (plt.Figure, plt.Axes):
    """Plot the terms in `tab`, in the order the appear.

    Write image if filename given."""
    import seaborn as sns
    fdr = tab.FDR.apply(lambda x: -np.log10(x))
    yrange = range(tab.shape[0])

    inmin, inmax = tab.InList.min(), tab.InList.max()

    def sz_convert(s, szmin=60, szmax=250):
        """Proportionally convert InList value to ranges between szmin, szmax"""
        distance = 1-(inmax-s)/(inmax-inmin)
        return szmin+(szmax-szmin)*distance

    sizemap = {n:sz_convert(n) for n in tab.InList.unique()}

    fig = plt.figure(figsize=(7, 0.3*tab.shape[0]))
    ax=plt.gca()

    sns.scatterplot(
        x=fdr,
        y=yrange,
        hue=tab.Enrichment,
        size=tab.InList,
        sizes=sizemap,
        palette='flare_r',
        ax=ax
    )

    # and now generate the legend with only 2 or 3 items for each
    #  instead of an entry for every unique value
    from matplotlib.lines import Line2D
    emin, emax = round(tab.Enrichment.min()), round(tab.Enrichment.max())
    markers = []
    markers.append(
        Line2D([0], [0], color='none', label='Enrichment')
    )
    cmap = plt.get_cmap('flare_r', )
    if (emax-emin)>5:
        ems = ( emax, emin+emax//2, emin)
    else:
        ems = (emax, emin)
    for enr in ems:

        d = 1-(emax-enr)/(emax-emin)
        l = Line2D(
            [0], [0],
            marker='o',
            color='none',
            label=str(round(enr)),
            markerfacecolor=cmap(d),
            markersize=np.sqrt(80)
        )
        markers.append(l)

    markers.append(
        Line2D([0], [0], color='none', label='In GO set')
    )
    if (inmax-inmin)>5:
        szs = (inmax, inmin+inmax//2, inmin)
    else:
        szs = (inmax, inmin)
    for sz in szs:
        l = Line2D(
            [0], [0],
            marker='o',
            color='none',
            label=str(sz),
            markerfacecolor='k',
            markersize=np.sqrt(sz_convert(sz))
        )
        markers.append(l)

    plt.legend(handles=markers, loc='lower right', alpha=1)

    plt.grid(axis='x')
    plt.yticks(yrange, tab.GoLabel)
    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=150)
    return fig, ax