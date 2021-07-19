import copy
import math

PLACEHOLDER = "PLACEHOLDER"


class Monomial:
    def __init__(self, indeterminate):
        self.indeterminates = [indeterminate]
        self.exponents = {indeterminate: 1}  # holds as values the power of each indeterminate

    def get_degree(self):
        return sum(self.get_exponents().values())

    def get_exponents(self):
        return self.exponents.copy()

    def get_indeterminates(self):
        return self.indeterminates.copy()

    def __eq__(self, other):  # two objects with the same attributes are equal
        if isinstance(other, Monomial):
            return self.exponents == other.exponents  # note that if exponents match, so will indeterminates
        else:
            return False

    def __gt__(self, other):  # monomials are in GRevLex order
        if isinstance(other,
                      (int, float, complex)):  # we assume that constants are never represented by Monomial objects
            return True
        elif isinstance(other, Monomial):
            if self == other:
                return False
            elif self.get_degree() != other.get_degree():
                return self.get_degree() > other.get_degree()
            else:  # compare same degree monomials by exponents of the latest indeterminates lexicographically
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
        return not self >= other

    def __le__(self, other):
        if self == other:
            return True
        else:
            return self < other

    def __mul__(self, other):  # we do not handle coefficients, except 0
        if isinstance(other, Monomial):  # multiply two monomials
            out_monomial = Monomial(PLACEHOLDER)
            out_monomial.indeterminates = sorted(list(set(self.indeterminates + other.indeterminates)))
            out_monomial.exponents = {x: self.exponents.get(x, 0) + other.exponents.get(x, 0)
                                      for x in out_monomial.indeterminates}
            return out_monomial
        elif other == 0:
            return 0
        elif isinstance(other, (int, float, complex)):  # multiplying a monomial by a constant
            return copy.deepcopy(self)
        else:
            raise TypeError("cannot multiply Monomial by type " + str(type(other)))

    def __rmul__(self, other):  # monomial multiplication is commutative
        return self * other

    def __pow__(self, power, modulo=None):
        if power == 0:
            return 1
        elif power == 1:
            return copy.deepcopy(self)
        elif isinstance(power, int) and power > 1:
            out_monomial = copy.deepcopy(self)
            # just multiply the indeterminates' exponents by the power
            out_monomial.exponents = {x: power * y for x, y in self.exponents.items()}
            if modulo is not None:
                out_monomial = divmod(out_monomial, modulo)[1]
            return out_monomial
        elif isinstance(power, int):
            raise ValueError("Cannot raise Monomial to a negative power")
        else:
            raise TypeError("Cannot raise Monomial to a power of type " + str(type(power)))

    def __divmod__(self, other):
        # quotient will be 0 if self is not divisible by other
        # remainder will be 0 or self in all cases
        if isinstance(other, Monomial):  # divide monomials
            if self == other:
                return 1, 0
            elif self < other:
                return 0, copy.deepcopy(self)
            else:
                for i in self.exponents:
                    print(i)
                out_exponents = {x: (self.exponents.get(x, 0) - other.exponents.get(x, 0)) for x in self.exponents}
                print(out_exponents)
                if [x for x in out_exponents.values() if x < 0] or [x for x in other.exponents if
                                                                    x not in self.exponents]:
                    return 0, copy.deepcopy(self)
                else:
                    out_monomial = Monomial(PLACEHOLDER)
                    out_monomial.indeterminates = [x for x in out_exponents if out_exponents[x] > 0]
                    out_monomial.exponents = {x: y for x, y in out_exponents.items() if y > 0}
                    return out_monomial, 0
        elif other == 0:
            raise ZeroDivisionError("Cannot divide by 0")
        elif isinstance(other, (int, float, complex)):  # divide by a constant
            return copy.deepcopy(self), 0
        else:
            raise TypeError("monomials cannot be divided by type " + str(type(other)))

    def __rdivmod__(self, other):  # should only be called if we divide a number by a monomial
        if other == 0:
            return 0, 0
        elif isinstance(other, (int, float, complex)):
            return 0, 1
        else:
            raise TypeError("Cannot divide object of type" + str(type(other)) + "by a Monomial")

    def __mod__(self, other):
        return self.__divmod__(other)[1]

    def __rmod__(self, other):
        return self.__rdivmod__(other)[1]

    def __hash__(self):
        return hash((tuple(self.indeterminates), tuple(self.exponents)))

    def __copy__(self):
        return copy.deepcopy(self)

    def __repr__(self):  # represent monomials in terms of common mathematical notation
        repr_list = [x + "^" + str(self.exponents[x]) if self.exponents[x] != 1 else x for x in self.indeterminates]
        return " * ".join(repr_list)

    def evaluate(self, substitutions=None):  # returns an ordered pair with a monomial or 1 and a coefficient
        # substitutions should be a dictionary
        if substitutions is None:
            return copy.deepcopy(self), 1
        elif not isinstance(substitutions, dict):
            raise TypeError(
                "'substitutions' should be of type " + str(dict) + " rather than " + str(type(substitutions)))
        else:
            out_coefficient = math.prod(
                [substitutions[x] ** self.exponents[x] for x in self.indeterminates if x in substitutions])
            if out_coefficient == 0:
                return 0, 0
            else:
                out_monomial = Monomial(PLACEHOLDER)
                out_monomial.indeterminates = [x for x in self.indeterminates if x not in substitutions]
                out_monomial.exponents = {x: self.exponents[x] for x in out_monomial.indeterminates}
                if not out_monomial.indeterminates:
                    return 1, out_coefficient
                else:
                    return out_monomial, out_coefficient


class Polynomial:
    def __init__(self, indeterminate):
        self.constant = 0
        self.terms = [(Monomial(indeterminate), 1)]  # list of terms of the form (monomial, coefficient)

    def get_degree(self):
        return self.terms[-1][0].get_degree()

    def get_constant(self):
        return copy.deepcopy(self.constant)

    def get_terms(self):
        return self.terms.copy()

    def __eq__(self, other):
        if isinstance(other, Polynomial):
            return self.constant == other.constant and self.terms == other.terms
        else:
            return False

    def __gt__(self, other):  # polynomials are in GRevLex order
        if isinstance(other, Polynomial):
            if self.terms == other.terms:
                return self.constant > other.constant
            else:
                s_terms = self.get_terms()
                o_terms = other.get_terms()
                while s_terms:
                    if s_terms.pop() > o_terms.pop():
                        return True
                return False

        elif isinstance(other, (int, float, complex)):
            return True
        else:
            raise TypeError("cannot compare Polynomial with type " + str(type(other)))

    def __ge__(self, other):
        if self == other:
            return True
        else:
            return self > other

    def __lt__(self, other):
        return not self >= other

    def __le__(self, other):
        return not self > other

    def __add__(self, other):
        if isinstance(other, Polynomial):
            out_polynomial = Polynomial(PLACEHOLDER)
            out_polynomial.constant = self.constant + other.constant
            monomial_list = sorted(list(set([x[0] for x in self.terms + other.terms])))
            self_dict = dict(self.terms)
            other_dict = dict(other.terms)
            out_polynomial.terms = [(x, self_dict.get(x, 0) + other_dict.get(x, 0)) for x in monomial_list
                                    if self_dict.get(x, 0) + other_dict.get(x, 0) != 0]
            if not out_polynomial.terms:
                return out_polynomial.constant
            else:
                out_polynomial.terms.sort()
                return out_polynomial
        elif isinstance(other, (int, float, complex)):
            out_polynomial = copy.deepcopy(self)
            out_polynomial.constant += other
            return out_polynomial
        else:
            raise TypeError("Cannot add polynomial to an object of type " + str(type(other)))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + -1 * other

    def __rsub__(self, other):
        return -1 * self + other

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            out_polynomial = self.constant * other
            iter_polynomial = Polynomial(PLACEHOLDER)
            for x in self.terms:
                iter_polynomial.constant = x[1] * other.constant
                iter_polynomial.terms = [(x[0] * a, x[1] * b) for a, b in other.get_terms()]
                out_polynomial += iter_polynomial
            return out_polynomial
        elif other == 0:
            return 0
        elif isinstance(other, (int, float, complex)):
            out_polynomial = copy.deepcopy(self)
            out_polynomial.constant *= other
            out_polynomial.terms = [(x, other * y) for x, y in out_polynomial.terms]
            return out_polynomial
        else:
            raise TypeError("cannot multiply polynomial by object of type " + str(type(other)))

    def __rmul__(self, other):
        return self * other

    def __pow__(self, power, modulo=None):
        if power == 0:
            return 1
        elif power == 1:
            return copy.deepcopy(self)
        elif isinstance(power, int) and power > 1:
            out_polynomial = 1
            for i in range(power):
                out_polynomial *= self
            if modulo is not None:
                out_polynomial = divmod(out_polynomial, modulo)[1]
            return out_polynomial
        elif isinstance(power, int):
            raise ValueError("Cannot raise Polynomial to a negative power")
        else:
            raise TypeError("Cannot raise Polynomial to a power of type " + str(type(power)))

    def __divmod__(self, other):
        pass

    def __rdivmod__(self, other):
        pass

    def __mod__(self, other):
        pass

    def __rmod__(self, other):
        pass

    def __repr__(self):
        term_list = self.get_terms()
        repr_list = [str(abs(x[1])) + " * " + x[0].__repr__() if x[1] not in [1, -1] else x[0].__repr__()
                     for x in term_list]  # get representations of monomials
        repr_str = ""
        if term_list.pop()[1] < 0:
            repr_str += "-"
        repr_str += repr_list.pop()
        while repr_list:
            if term_list.pop()[1] < 0:
                repr_str += " - "
            else:
                repr_str += " + "
            repr_str += str(repr_list.pop())
        if self.constant < 0:
            repr_str += " - " + str(abs(self.constant))
        elif self.constant > 0:
            repr_str += " + " + str(self.constant)
        return repr_str

    def __copy__(self):
        return copy.deepcopy(self)

    def evaluate(self, substitutions=None):
        pass


t = Polynomial('t')
s = Polynomial('s')
p = (t + (s**5*t)**2)**4
q = copy.copy(p)
print(p > t)
