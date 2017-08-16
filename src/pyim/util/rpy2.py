"""Utility functions for interacting with rpy2."""

from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.rinterface import RNULLType

pandas2ri.activate()


def pandas_to_dataframe(pd_frame, check_names=False):
    """Converts pandas dataframes to R dataframes."""

    r_frame = pandas2ri.py2ri_pandasdataframe(pd_frame)
    r_frame.colnames = pd_frame.columns

    if not check_names:
        r_frame.rownames = pd_frame.index

    return r_frame


def dataframe_to_pandas(r_frame):
    """Converts R dataframes to pandas dataframes."""

    pd_frame = pandas2ri.ri2py_dataframe(r_frame)

    # Extract column names if possible.
    col_names = robjects.r.colnames(r_frame)
    if not isinstance(col_names, RNULLType):
        pd_frame.columns = col_names

    # Extract row names if possible.
    index = robjects.r.rownames(r_frame)
    if not isinstance(index, RNULLType):
        pd_frame.index = index

    return pd_frame
