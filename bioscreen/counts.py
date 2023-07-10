import pandas as pd
import numpy as np

def normalise_abundance(count:pd.DataFrame) -> pd.DataFrame:
    """Make column sums equal to the median of untransformed column sums"""
    return (count / count.sum()) * count.sum().median()

def normalise_zscore(count:pd.DataFrame) -> pd.DataFrame:
    """Make data have mean==0 and std==1"""
    return (count-count.mean())/count.std()

def normalise_median(count:pd.DataFrame) -> pd.DataFrame:
    return count/count.median()

def apply_log2(count:pd.DataFrame, pseudocount=1):
    return count.apply(lambda n: np.log2(n+pseudocount))