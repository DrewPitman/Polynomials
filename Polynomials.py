import copy
import math


class Polynomial:
    def __init__(self, indeterminate=None):
        self.terms = []  # monomials and constants
        if indeterminate is not None:
            self.terms.append(Monomial(indeterminate))


# terms, coefficients, indeterminates, exponents

PLACEHOLDER = "PLACEHOLDER"


class Monomial:
    def __init__(self, indeterminate):
        self.degree = 1
        self.indeterminates = [indeterminate]
        self.exponents = {indeterminate: 1}

    def get_degree(self):
        return copy.copy(self.degree)  # calling copy might help if the class is extended to include noninteger degrees

    def get_exponents(self):
        return self.exponents.copy()

    def get_indeterminates(self):
        return self.indeterminates.copy()

    def __eq__(self, other):
        if isinstance(other, Monomial):
            if self.degree != other.degree:  # potentially time-saving
                return False
            else:
                return self.exponents == other.exponents  # note that if exponents match, so will indeterminates
        else:
            return False

    def __gt__(self, other):  # we use GRevLex order
        if isinstance(other,
                      (int, float, complex)):  # we assume that constants are never represented by Monomial objects
            return True
        elif isinstance(other, Monomial):
            if self == other:
                return False
            elif self.degree != other.degree:
                return self.degree > other.degree
            else:  # we compare monomials of the same degree by comparing the exponents of the last indeterminates lexicographically
                self_variables = self.get_indeterminates()
                other_variables = other.get_indeterminates()
                while self_variables:
                    x = self_variables.pop()
                    y = other_variables.pop()
                    if x != y:
                        return x < y
                    elif self.exponents[x] != other.exponents[x]:
                        return self.exponents[x] < other.exponents[x]
        else:
            raise TypeError("cannot compare Monomial with type " + str(type(other)))

    def __ge__(self, other):
        if self == other:
            return True
        else:
            return self > other

    def __lt__(self, other):
        if isinstance(other, Monomial):
            return other > self
        else:
            return False

    def __le__(self, other):
        if self == other:
            return True
        else:
            return self < other

    def __mul__(self, other):  # we do not handle coefficients, including 0
        if isinstance(other, Monomial):  # multiply monomials
            out_monomial = Monomial(PLACEHOLDER)
            out_monomial.indeterminates = sorted(list(set(self.indeterminates + other.indeterminates)))
            out_monomial.exponents = {x: self.exponents.get(x, 0) + other.exponents.get(x, 0) for x in
                                      out_monomial.indeterminates}
            out_monomial.degree = sum(out_monomial.exponents.values())
            return out_monomial
        elif other == 0:
            return 0
        elif isinstance(other,
                        (int, float, complex)):  # multiplying a monomial by a constant just gives you the monomial
            return copy.copy(self)
        else:
            raise TypeError("cannot multiply Monomial by type " + str(type(other)))

    def __rmul__(self, other):  # monomial multiplication is commutative
        return self * other

    def __pow__(self, power, modulo=None):  # add stuff for modulo. May not be useful for this class, though.
        if power == 0:
            return 1
        elif power == 1:
            return copy.copy(self)
        elif isinstance(power, int) and power > 1:
            out_monomial = copy.copy(self)
            out_monomial.exponents = {x: power * y for x, y in self.exponents.items()}
            if modulo is not None:
                out_monomial = divmod(out_monomial, modulo)[1]
            return out_monomial
        elif isinstance(power, int):
            raise ValueError("Cannot raise Monomial to a negative power")
        else:
            raise TypeError("Cannot raise Monomial to a power of type " + str(type(power)))

    def __divmod__(self,
                   other):  # remainder will be 0, 1, or a monomial. we do not work with coefficients in this class
        if isinstance(other, Monomial):  # divide monomials
            if self == other:
                return 1, 0  # using 1 for the remainder since it's the multiplicative identity
            elif self < other:
                return 0, copy.copy(self)
            else:
                for i in self.exponents:
                    print(i)
                out_exponents = {x: (self.exponents.get(x, 0) - other.exponents.get(x, 0)) for x in self.exponents}
                print(out_exponents)
                if [x for x in out_exponents.values() if x < 0] or [x for x in other.exponents if
                                                                    x not in self.exponents]:
                    return 0, copy.copy(self)
                else:
                    out_monomial = Monomial(PLACEHOLDER)
                    out_monomial.indeterminates = [x for x in out_exponents if out_exponents[x] > 0]
                    out_monomial.exponents = {x: y for x, y in out_exponents.items() if y > 0}
                    out_monomial.degree = sum(out_monomial.exponents.values())
                    return out_monomial, 0
        elif other == 0:
            raise ZeroDivisionError("Cannot divide by 0")
        elif isinstance(other, (int, float, complex)):  # divide by a constant
            return copy.copy(self), 0
        else:
            raise TypeError("monomials cannot be divided by type " + str(type(other)))

    def __rdivmod__(self, other):  # should only be called if we divide a number by a monomial
        if other == 0:
            return 0, 0
        elif isinstance(other, (int, float, complex)):
            return 0, 1
        else:
            raise TypeError("Cannot divide object of type" + str(type(other)) + "by a Monomial")

    def __repr__(self):
        repr_list = [x + "^" + str(self.exponents[x]) if self.exponents[x] != 1 else x for x in self.indeterminates]
        return " * ".join(repr_list)

    def evaluate(self, substitutions=None):  # returns an ordered pair with a monomial or 1 and a coefficient
        # substitutions should be a dictionary
        if substitutions is None:
            return copy.copy(self), 1
        elif not isinstance(substitutions, dict):
            raise TypeError("'substitions' should be of type " + str(dict) + " rather than " + str(type(substitutions)))
        else:
            out_coefficient = math.prod(
                [substitutions[x] ** self.exponents[x] for x in self.indeterminates if x in substitutions])
            if out_coefficient == 0:
                return 0, 0
            else:
                out_monomial = Monomial(PLACEHOLDER)
                out_monomial.indeterminates = [x for x in self.indeterminates if x not in substitutions]
                out_monomial.exponents = {x: self.exponents[x] for x in out_monomial.indeterminates}
                out_monomial.degree = sum(out_monomial.exponents.values())
                if not out_monomial.indeterminates:
                    return 1, out_coefficient
                else:
                    return out_monomial, out_coefficient


