__author__ = 'Julian'

import pandas as pd


class PandasDfWrapper(object):

    def __init__(self, frame, copy=True):
        if copy:
            frame = frame.copy()

        self._frame = frame
        self._ix = IXWrapper(self._frame.ix, self.__class__)

    @property
    def ix(self):
        return self._ix

    def __getattr__(self, item):
        attr = getattr(self._frame, item)
        if callable(attr):
            def wrapper(*args, **kwargs):
                ret = attr(*args, **kwargs)
                if isinstance(ret, pd.DataFrame):
                    ret = self.__class__(ret)
                return ret
            return wrapper
        else:
            return attr

    def __getitem__(self, item):
        return self._frame[item]

    def __setitem__(self, key, value):
        self._frame[key] = value

    def __iter__(self):
        return iter(self._frame)

    def __repr__(self):
        return repr(self._frame)

    def __str__(self):
        return str(self._frame)


class IXWrapper(object):

    def __init__(self, ix, class_):
        self._ix = ix
        self._class = class_

    def __getitem__(self, *args, **kwargs):
        ret = self._ix.__getitem__(*args, **kwargs)
        if isinstance(ret, pd.DataFrame):
            ret = self._class(ret)
        return ret
