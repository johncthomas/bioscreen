import typing

from sklearn.decomposition import PCA as skPCA
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, Collection


class CountPCA:
    def __init__(self, table: pd.DataFrame,
                 sample_details:Optional[pd.DataFrame]=None, **pcaKW):


        # check that the tables are compatable, and in same order
        if sample_details is not None:
            assert all(table.columns.isin(sample_details.index))
            sample_details = sample_details.reindex(index=table.columns)
        self.sample_details = sample_details

        pca = skPCA(**pcaKW)
        self.pca = pca
        scores = pca.fit_transform(table.T)

        pc_idx = [f"PC{i}" for i in range(1, scores.shape[1]+1)]

        self.scores = pd.DataFrame(
            scores,
            columns=pc_idx,
            index=table.columns
        )

        self.loadings = pd.DataFrame(
            pca.components_,
            columns=table.index,
            index=pc_idx
        )

        self.explained_variance_perc = pd.Series(pca.explained_variance_ratio_ * 100, index=pc_idx)

        self.samples = table.columns
        self.genes = table.index

    def scatter_plot(self, pc_x='PC1', pc_y='PC2', labels: Union[bool, Collection[str]] = False,
                     **scatterplot_kwargs):
        """Plot PCs.
        Args:
            pc_x: ONE indexed. principal components to plot the x values
            pc_y: as above but y

            labels: if True, use the table columns to label points, or
                pass a collection to be used as labels.

            scatterplot_kwargs are passed to sns.scatterplot. Use hue
            and style = pd.Series to show sample properties"""

        if type(pc_x) is int:
            pc_x = f"PC{pc_x}"
        if type(pc_y) is int:
            pc_y = f"PC{pc_y}"

        # default options, overwritten by any user supplied kwargs
        default_kw = dict(s=50)
        ax = sns.scatterplot(x=self.scores[pc_x], y=self.scores[pc_y],
                             **(default_kw|scatterplot_kwargs))

        try:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))
        # when there's no legend
        except ValueError:
            pass

        # axis labels give % variance explained
        for xylabel, pcxy in [(plt.xlabel, pc_x), (plt.ylabel, pc_y)]:
            xylabel(f"{pcxy} ({self.explained_variance_perc[pcxy]:.3}% variance explained)")

        # remove numbers from axes
        for ticks in (plt.xticks, plt.yticks):
            vals, labs = ticks()
            ticks(vals, ['' for _ in vals])

        # label points
        if labels is not False:
            if labels is True:
                labels = self.samples
            for i, lab in enumerate(labels):
                labxy = (self.scores.loc[:, pc_x][i],
                         self.scores.loc[:, pc_y][i])
                #todo label kwargs
                plt.annotate(
                    lab,
                    labxy,
                    (3, 3),
                    textcoords='offset points',
                    size=8,
                )

    def anova(self, max_pc=np.inf) -> pd.DataFrame:
        """Table of anova statistics for factor associations with PCs.

        return: DF with level 0 columns of F, p & FDR and level 1
        giving the PC results."""
        table = {'F': {}, 'p': {}}

        for pci in range(self.scores.shape[1]):
            pc = f"PC{pci + 1}"
            if pci > max_pc:
                break

            table['p'][pc] = p_res = {}
            table['F'][pc] = f_res = {}
            for factor in self.sample_details.columns:
                factor_groups = self.sample_details.groupby(factor).groups
                values = [self.scores.loc[f, pc] for f in factor_groups.values()]

                # if there's only one level to the factor
                if len(values) < 2:
                    continue

                res = stats.f_oneway(*values)

                p_res[factor] = res.pvalue  # not an error, PyCharm
                f_res[factor] = res.statistic

        ptable = pd.DataFrame(table['p'])
        fdr = sm.stats.multipletests(np.ravel(ptable), method='fdr_bh', )[1]
        fdr = pd.DataFrame(np.reshape(fdr, ptable.shape), index=ptable.index, columns=ptable.columns)
        return pd.concat({'F': pd.DataFrame(table['F']), 'p': ptable, 'FDR': fdr}, axis=1)


#todo add pca_rugplot (below)

# def pca_rugplot(self, pc:Union[int, str], hue=None):
#     if type(pc) is int:
#         pc = f"PC{pc}"
#     kw = dict(x=p.scores[pc], hue=hue)
#     sns.kdeplot(**kw)
#     sns.rugplot(**kw, height=0.1)
#     plt.xticks([], [])
#     plt.yticks([], [])

