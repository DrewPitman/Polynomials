from collections import abc
from typing import TypeVar
from polynomials.polynomials import *

_Poly = TypeVar('_Poly', Polynomial, int, float, complex)

class PolynomialSystem(abc.MutableSequence):
    def __init__(self, polynomials=None, /, *args):
        super(PolynomialSystem, self).__init__()
        self._list = []
        if polynomials is not None:
            if isinstance(polynomials, list):
                self._list = polynomials
            else:
                self._list = [polynomials]
                if args:
                    self._list += list(args)

    def __add__(self, other):
        assert(len(self) == len(other))
        return PolynomialSystem([self[i] + other[i] for i in range(len(self))])

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        assert (len(self) == len(other))
        return PolynomialSystem([self[i] - other[i] for i in range(len(self))])

    def __rsub__(self, other):
        assert (len(self) == len(other))
        return PolynomialSystem([other[i] - self[i] for i in range(len(self))])

    def __mul__(self, other):
        pass

    def __delitem__(self, key):
        self._list.__delitem__(key)

    def __getitem__(self, item):
        return self._list[item]

    def __len__(self):
        return len(self._list)

    def __setitem__(self, key, value: _Poly):
        self._list[key] = value

    def __repr__(self):
        return "<{0} {1}>".format(self.__class__.__name__, self._list)

    def insert(self, index: int, value: _Poly) -> None:
        self._list.insert(self, int, value)

    def append(self, value: _Poly) -> None:
        self._list.append(value)

    def evaluate(self, substitutions):
        return PolynomialSystem([p.evaluate(substitutions) for p in self])

    def indeterminates(self):
        ind = list(set([i for p in self for i in p.get_indeterminates()]))
        ind.sort()
        return ind


def gradient(polynomial, variables=None):
    if isinstance(polynomial, Polynomial):
        if variables is None:
            variables = polynomial.get_indeterminates()
        return PolynomialSystem([polynomial.derivative(x) for x in variables])
    elif variables is None:
        return PolynomialSystem()
    else:
        return PolynomialSystem([0 for x in variables])
