import pandas as pd
import rpy2
import rpy2.rinterface
import rpy2.robjects as ro
R = ro.r
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import os, logging

Pathy = os.PathLike | str

def pd_context(*args, **kwargs):
    """Return active context manager for calling R functions with, or
     directly converting, pandas objects."""
    return (ro.default_converter + pandas2ri.converter).context(*args, **kwargs)

def pd_convert(pd_obj):
    """Return an rpy2 conversion of a Pandas object."""
    with pd_context():
        return ro.conversion.get_conversion().py2rpy(pd_obj)

def r_to_pd(r_obj):
    with pd_context():
        res = ro.conversion.get_conversion().rpy2py(r_obj)
    return res


from rpy2.rinterface_lib.embedded import RRuntimeError
def rcatcher(func, verbosity=1):
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except RRuntimeError as err:
            logging.critical(f"{func.__name__} failed, R error: \n\t{str(err)}")
            raise err
    return wrapper


NULL = rpy2.rinterface.NULL

pkgdir = os.path.dirname(__file__)