import dataclasses

from typing import (
    Collection,
    Literal
)

import pandas as pd

pd.options.display.date_yearfirst = True

from bioscreen.rinterfaces.utils import *

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

from bioscreen.experiment_classes import CompDict
from bioscreen.classes.differential_gene_expression import LimmaResults, LIMMACOLS

from jttools.data_wrangling import AttrMapAC
from jttools.statistics import neglog10
from bioscreen.classes.base import SigCols


__all__ = ['Limma']


#from bioscreen.classes.experiment import ScreenExperiment


## # keeping this to stricly the Rinterface part.
# def load_limma_results(results_dir) -> pd.DataFrame:
#     """Load a set of limma comparisons structured as dir of CSV with
#     file names from contrasts made.
#     """
#     res = {}
#     for fn in os.listdir(results_dir):
#         filepath = os.path.join(results_dir, fn)
#         df = pd.read_csv(filepath, index_col=0)
#         # df.columns = df.columns.map(colmap)
#         df.loc[:, 'log10_FDR'] = df['adj.P.Val'].apply(neglog10)
#         res[fn.replace('.csv', '')] = df
#     return pd.concat(res, axis=1)


# def tabulate_star_logs(star_dir, samples=None, split_on_R1=True):
#     """Look for files following pattern f'{star_dir}/{samples[i]}/Log.final.out',
#     put contents in a DF, rows are samples.
#
#     When split_on_R1, row samples have the Illumina _R1_L001 suffix removed."""
#     if samples is None:
#         samples = os.listdir(star_dir)
#
#     star_outputs = []
#     for samp in samples:
#         row = {}
#         with open(os.path.join(star_dir, samp, 'Log.final.out')) as f:
#             for line in f:
#                 if '|' in line:
#                     k = line.split(' |')[0].strip()
#                     v = line.strip().rsplit()[-1]
#                     row[k] = v
#             star_outputs.append(row)
#     df = pd.DataFrame(star_outputs, index=[s.split('_R1_')[0] for s in samples])
#     df.index.name = 'Sample'
#     return df

R.source(os.path.join(pkgdir, "limmaFunctions.R"))

@dataclasses.dataclass
class LimmaRObjects:
    test_factors = NULL
    sample_factors = NULL
    counts = NULL
    design = NULL
    block = NULL
    correlation = NULL
    fit = NULL
    contrast_res = NULL
    sample_details = NULL

class Limma:

    def __init__(self,
                 #analysis_type:typing.Literal['rnaseq', 'proteomics'],
                 counts:pd.DataFrame,
                 sample_details:pd.DataFrame,
                 comparisons:CompDict,
                 test_groups:Collection[str],
                 block:Collection[str] = None,
                 #included_samples:Collection[str] = 'all',
                 voom_counts:bool = False,
                 analysis_version:str=None,
                 loglevel:Literal['INFO','DEBUG','WARNING']='INFO'):
        """Interface for running Limma in R.

        Implimentation not use formulas. test_groups should identify combinations of
        treatments (if multiple used) and differential contrasts can be used
        to pull out interactions. The model doesn't use an intercept, so contrasts
        are required.

        Arguments:
            counts: Appropriate abundance quantifications for limma.
            sample_details: Table indexed by Sample column, containing
                SampleGroup and columns for all model_factors.
            comparisons: List of Comparisons to form contrasts.
            test_groups: Test group names in order of appearance in counts.
                Will probably just be SampleGroups.
            block: Random factors. Variation specific to these detected with
                duplicateCorrelation and accounted for in the model.
            voom_counts: Apply voom to counts before running, e.g. they are RNAseq
                counts and need VST.
            analysis_version: A string identifying this specific analysis.
                Currently is just recorded in the class instance.

        Methods:
            run: call all functions required to run analysis and return dict of DF
                containing results.
        """

        logger.setLevel(loglevel)

        self.sample_details = sample_details
        self.comparisons = comparisons
        self.analysis_version = analysis_version #update doc if you ever do anything with this
        self.counts:pd.DataFrame = counts
        self.test_groups:Collection[str] = list(test_groups)
        self.voom_counts = voom_counts
        #self.filter_expr = filter_expr
        self.block:Collection[str] = block
        self.robj = LimmaRObjects()

    #todo Limma.from_experiment
    # does things like pulling test_groups from sample_details automatically



    def prep_rnaseq_data(self):
        """Filter low expression, optionally do voom,
         optionally calculatecorr.

        Adds updates robj if voom"""
        #todo test prep_rnaseq_data

        # cnt = self.counts.loc[:, self.included_samples]
        # sd = self.sample_details.loc[self.included_samples]

        robj = self.robj


        with pd_context():
            self.robj.counts = R.prep_rnaseq_counts(
                self.counts, robj.sample_factors
            )

        if self.voom_counts:
            # iteratively does variance corr if block
            vres = R.do_voom(
                robj.counts,
                robj.design,
                robj.block
            )
            robj.counts = vres['vcounts']
            robj.correlation = vres['correlation']

        # variance corr, already done in conjunction with voom if do_voom
        if (self.block is not None) and (not self.voom_counts):
            corr = ro.r.duplicateCorrelation(
                robj.counts, robj.design,
                block=robj.block
            )
            robj.correlation = corr



    def prep_data(self):
        """Create R objects for counts, design, test|sample_factors,
        & block. If block, get corrlelation."""
        robj = self.robj

        robj.test_factors = ro.r['as.factor'](ro.StrVector(self.test_groups))
        robj.sample_factors = ro.r['as.factor'](ro.StrVector(
            list(self.sample_details.SampleGroup)
        ))

        robj.counts = pd_convert(self.counts)
        robj.sample_details = pd_convert(self.sample_details)

        robj.design = R.get_design(
                robj.counts, robj.sample_details, robj.test_factors
            )

        if self.block is not None:
            robj.block = ro.StrVector([str(b) for b in self.block])
            robj.correlation = R.duplicateCorrelation(
                robj.counts, robj.design, block=robj.block
            )


    def fit(self):
        """Perform limma::lmFit. Adds result to self.robj.fit"""
        robj = self.robj
        if (robj.counts is NULL) or (robj.design is NULL):
            raise RuntimeError("Run prep_data first.")

        robj.fit = R.get_fit(
            robj.counts,
            robj.design,
            block=robj.block,
            correlation=robj.correlation
        )

    def fit_contrasts(self):

        comparisons = self.comparisons
        robj = self.robj

        # "testsamp - ctrlsamp"
        contrasts = ro.StrVector(
            comparisons.to_formulas()
        )

        names = ro.StrVector(comparisons.names())

        robj.contrast_res = R.fit_contrasts(robj.fit, robj.design, contrasts, names)


    def get_results(self) -> LimmaResults:

        contrast_names = list(R.colnames(self.robj.contrast_res[0]))
        tables = AttrMapAC()
        for cntrst in contrast_names:
            with pd_context():

                tables[cntrst] = table = R.get_toptable(self.robj.contrast_res, cntrst)

                table.loc[:, SigCols.p10] = table[LIMMACOLS.p.original].apply(neglog10)
                table.loc[:, SigCols.FDR10] = table[LIMMACOLS.FDR.original].apply(neglog10)
        return LimmaResults.build(
            tables,
            comparisons=self.comparisons,
            scorekey='LFC'
        )

    # def contrast_tables(self) -> dict[str, pd.DataFrame]:
    #     with pd_context():
    #         tables = {k:tab for k, tab in self.robj.contrast_res.items()}
    #
    #     return tables

    def run_rnaseq(self) -> LimmaResults:
        """Prep data (using prep_rnaseq_data), do fits and produce
        contrast tables."""
        self.prep_rnaseq_data()
        self.fit()
        self.fit_contrasts()
        return self.get_results()

    def run(self) -> LimmaResults:
        """Prep data do fits and produce
        contrast tables."""
        self.prep_data()
        self.fit()
        self.fit_contrasts()
        return self.get_results()


if __name__ == '__main__':
    pass


    import os
    from jttools.picklepot import PicklePot
    os.chdir('/mnt/m/tasks/NA327_Proteomics_UbPulldown/')
    pp = PicklePot('pickles2', )

    sd = pp.objects['smol_sd']
    comparisons = pp.objects['smol_cmp']
    counts = pp.objects['smol_cnt']

    smol = Limma(
        counts=counts,
        sample_details=sd,
        comparisons=comparisons,
        test_groups=sd.SampleGroup
    )

    res = smol.run()
    print(res.table.head())

    res.write_comp_results_to_excel('/mnt/m/fff.xlsx')

