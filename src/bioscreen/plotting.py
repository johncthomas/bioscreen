import matplotlib.pyplot as plt
import seaborn as sns
from jttools.plotting import get_palette

import pandas as pd

def plot_abundance_violin(
        count:pd.DataFrame, color_factor:pd.Series=None
):
    """Violin plots, colored by unique values in color_factors."""

    if color_factor is not None:
        pal = get_palette(color_factor)
    else:
        pal = None

    # the figure
    plt.figure(figsize=(0.6*count.shape[1], 3))
    sns.violinplot(
        data=count,
        cut=0,
        palette=pal,
    )

    plt.ylabel('Abundance (log2)')

