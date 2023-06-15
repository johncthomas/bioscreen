from sklearn.decomposition import PCA as skPCA
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, Collection


class countPCA:
    def __init__(self, table:pd.DataFrame, **pcaKW):
        pca = skPCA(**pcaKW)
        scores = pca.fit_transform(table.T)
        pc_idx = [f"PC{i}" for i in range(scores.shape[1])]
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

        self.explained_variance_perc = pca.explained_variance_ratio_*100

        self.samples = table.columns
        self.genes = table.index

    def scatter_plot(self, pc_x=1, pc_y=2, labels:Union[bool, Collection[str]]=False,
             **scatterplot_kwargs):
        """Plot PCs.
        Args:
            pc_x: ONE indexed. principal components to plot the x values
            pc_y: as above but y

            labels: if True, use the table index to label points, or
                pass a collection to be used as labels.

            scatterplot_kwargs are passed to sns.scatterplot. Use hue
            and style = pd.Series to show sample properties"""
        ax = sns.scatterplot(x=self.scores.loc[:, pc_x], y=self.scores.loc[:, pc_y],
                             **scatterplot_kwargs)
        try:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))
        # when there's no legend
        except ValueError:
            pass

        for xylabel, pcxy in [(plt.xlabel, pc_x), (plt.ylabel, pc_y)]:
            xylabel(f"PC {pcxy + 1} ({self.explained_variance_perc[pcxy]:.3}% variance explained)")

        for ticks in (plt.xticks, plt.yticks):
            vals, labs = ticks()
            ticks(vals, ['' for _ in vals])

        if labels is not False:
            if labels is True:
                labels = self.samples
            for i, lab in enumerate(labels):
                labxy = (self.scores[:, pc_x][i],
                         self.scores[:, pc_y][i])
                plt.annotate(
                    lab,
                    labxy,
                    (3, 3),
                    textcoords='offset points',
                    size=8,
                )

