import os
import functools


class SympyCache(object):
    def __init__(self, path):
        """

        Parameters
        ----------
        path : str
           Name of file containing actual cache data

        """

        self.path = path
        self.calc = not os.path.isfile(self.path)
        if not self.calc:
            inp = open(self.path, 'r')
            self.cache = inp.readlines()
            inp.close()

    def cache(self, func):
        @functools.wraps(func)
        if self.calc:
            val = func(*args, **kwargs)
            self.cache.append(val)
        else:
            val = self.cache.pop(0)
        return val
