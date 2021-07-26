import copy
import math

# some global variables that will really help us out
_PLACEHOLDER = "_PLACEHOLDER"
_NUMBER_TYPES = (int, float, complex)


# noinspection PyTypeChecker
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
        if isinstance(other, _NUMBER_TYPES):  # we assume that constants are never represented by Monomial objects
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
            out_monomial = Monomial(_PLACEHOLDER)
            out_monomial.indeterminates = sorted(list(set(self.indeterminates + other.indeterminates)))
            out_monomial.exponents = {x: self.exponents.get(x, 0) + other.exponents.get(x, 0)
                                      for x in out_monomial.indeterminates}
            return out_monomial
        elif other == 0:
            return 0
        elif isinstance(other, _NUMBER_TYPES):  # multiplying a Monomial by a constant
            return copy.deepcopy(self)
        else:
            raise TypeError("cannot multiply Monomial by type " + str(type(other)))

    def __rmul__(self, other):  # Monomial multiplication is commutative
        return self * other

    def __truediv__(self, other):
        quotient, remainder = divmod(self, other)
        if remainder == 0:
            return quotient
        else:
            raise ValueError(str(self) + " is not divisible by " + str(other))

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
                out_exponents = {x: (self.exponents.get(x, 0) - other.exponents.get(x, 0)) for x in self.exponents}
                if [x for x in out_exponents.values() if x < 0] or [x for x in other.exponents if
                                                                    x not in self.exponents]:
                    return 0, copy.deepcopy(self)
                else:
                    out_monomial = Monomial(_PLACEHOLDER)
                    out_monomial.indeterminates = [x for x in out_exponents if out_exponents[x] > 0]
                    out_monomial.exponents = {x: y for x, y in out_exponents.items() if y > 0}
                    return out_monomial, 0
        elif other == 0:
            raise ZeroDivisionError("Cannot divide by 0")
        elif isinstance(other, _NUMBER_TYPES):  # divide by a constant
            return copy.deepcopy(self), 0
        else:
            raise TypeError("monomials cannot be divided by type " + str(type(other)))

    def __rdivmod__(self, other):  # should only be called if we divide a number by a Monomial
        if other == 0:
            return 0, 0
        elif isinstance(other, _NUMBER_TYPES):
            return 0, 1
        else:
            raise TypeError("Cannot divide object of type" + str(type(other)) + "by a Monomial")

    def __mod__(self, other):
        return self.__divmod__(other)[1]

    def __rmod__(self, other):
        return self.__rdivmod__(other)[1]

    def __hash__(self):
        return hash((tuple(self.indeterminates), tuple(self.exponents)))

    def __copy__(self):  # Monomials are immutable, even if their attributes aren't
        return copy.deepcopy(self)

    def __repr__(self):  # represent monomials in terms of common mathematical notation
        repr_list = [x + "^" + str(self.exponents[x]) if self.exponents[x] != 1 else x for x in self.indeterminates]
        return " * ".join(repr_list)

    def divides(self, other):
        return divmod(other, self)[1] == 0

    def evaluate(self, substitutions=None):  # returns an ordered pair with a Monomial or 1 and a coefficient
        # substitutions should be a dictionary
        if substitutions is None:
            return copy.deepcopy(self), 1
        elif isinstance(substitutions, dict):
            out_coefficient = math.prod(
                [substitutions[x] ** self.exponents[x] for x in self.indeterminates
                 if x in substitutions and isinstance(substitutions[x], _NUMBER_TYPES)])
            if out_coefficient == 0:
                return 0, 0
            else:
                out_monomial = Monomial(_PLACEHOLDER)
                # unchanged indeterminates
                out_monomial.indeterminates = [x for x in self.indeterminates if x not in substitutions]
                extra_indeterminates = [substitutions[x] for x in self.indeterminates
                                        if x in substitutions and substitutions[x] not in out_monomial.indeterminates
                                        and isinstance(substitutions[x], str)]

                # exponents of unchanged variables that are not in the image of substitutions
                out_monomial.exponents = {x: self.exponents[x] for x in out_monomial.indeterminates
                                          if x not in substitutions.values()}
                # exponents of unchanged variables that are in the image of substitutions
                shared_exponents = {substitutions[x]: (self.exponents.get(x, 0) + self.exponents[substitutions[x]])
                                    for x in list(substitutions.keys())
                                    if substitutions.get(x) in out_monomial.indeterminates}
                # exponents of variables in the image of substitutions
                # that do not overlap with the updated list of indeterminates
                extra_exponents = {substitutions[x]: self.exponents[x] for x in self.indeterminates
                                   if substitutions.get(x) in extra_indeterminates}
                out_monomial.indeterminates += extra_indeterminates
                out_monomial.exponents = {**out_monomial.exponents, **shared_exponents, **extra_exponents}
                out_monomial.indeterminates.sort()

                if not out_monomial.indeterminates:
                    return 1, out_coefficient
                else:
                    return out_monomial, out_coefficient
        else:
            raise TypeError(
                "'substitutions' should be of type " + str(dict) + " rather than " + str(type(substitutions)))

    def derivative(self, indeterminate):
        ind = copy.deepcopy(indeterminate)
        # convert ind into a string if it is not already a string
        if isinstance(ind, Polynomial):
            ind = ind.terms[0][0]
        if isinstance(ind, Monomial):
            ind = ind.indeterminates[0]

        if ind in self.indeterminates:
            mon = copy.deepcopy(self)
            if mon.exponents[ind] == 1:
                mon.indeterminates = [x for x in mon.indeterminates if x is not ind]
                mon.exponents = {x: y for x, y in mon.exponents.items() if x is not ind}
                if mon.indeterminates:
                    return mon, 1
                else:
                    return 1, 1
            else:
                mon.exponents[ind] -= 1
                return mon, mon.exponents[ind] + 1
        else:
            return 0, 0


# noinspection PyTypeChecker
class Polynomial:
    def __init__(self, indeterminate):
        self.constant = 0
        self.terms = [(Monomial(indeterminate), 1)]  # list of terms of the form (Monomial, coefficient)

    def get_degree(self):
        return self.terms[-1][0].get_degree()

    def get_constant(self):
        return copy.deepcopy(self.constant)

    def get_terms(self):
        return self.terms.copy()

    def break_into_terms(self):
        out_list = [Polynomial(_PLACEHOLDER) for i in range(len(self.terms))]
        for i in range(len(out_list)):
            out_list[i].terms = [self.terms[i]]
        if self.constant:
            out_list.append(self.constant)
        out_list.sort()
        return out_list

    def get_indeterminates(self, string_output=False):
        indeterminate_list = [i for x, y in self.terms for i in x.get_indeterminates()]
        indeterminate_list = list(set(indeterminate_list))
        indeterminate_list.sort()
        if not string_output:
            indeterminate_list = [Polynomial(i) for i in indeterminate_list]
        return indeterminate_list

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

        elif isinstance(other, _NUMBER_TYPES):
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
            out_polynomial = Polynomial(_PLACEHOLDER)
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
        elif isinstance(other, _NUMBER_TYPES):
            out_polynomial = copy.deepcopy(self)
            out_polynomial.constant += other
            return out_polynomial
        else:
            raise TypeError("Cannot add Polynomial to an object of type " + str(type(other)))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + -1 * other

    def __rsub__(self, other):
        return -1 * self + other

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            out_polynomial = self.constant * other
            iter_polynomial = Polynomial(_PLACEHOLDER)
            for x in self.terms:
                iter_polynomial.terms = [(x[0] * a, x[1] * b) for a, b in other.get_terms()]
                out_polynomial += iter_polynomial
            return out_polynomial
        elif other == 0:
            return 0
        elif isinstance(other, _NUMBER_TYPES):
            out_polynomial = copy.deepcopy(self)
            out_polynomial.constant *= other
            out_polynomial.terms = [(x, other * y) for x, y in out_polynomial.terms]
            return out_polynomial
        else:
            raise TypeError("cannot multiply Polynomial by object of type " + str(type(other)))

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        quotient, remainder = divmod(self, other)
        if remainder == 0:
            return quotient
        else:
            raise ValueError(str(self) + " is not divisible by " + str(other))

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
        if isinstance(other, Polynomial):  # divide polynomials
            if divmod(self.terms[-1][0], other.terms[-1][0])[0] == 0:
                return 0, self
            else:
                # initialize variables
                o_term = other.get_terms()[-1]
                dividend = copy.deepcopy(self)
                quotient = 0

                # perform long division
                while isinstance(dividend, Polynomial) and 0 == divmod(dividend.terms[-1][0], o_term[0])[1]:
                    subquotient_tuple = divmod(dividend.terms[-1][0], o_term[0])
                    if isinstance(subquotient_tuple[0], Monomial):
                        subquotient = Polynomial(_PLACEHOLDER)
                        subquotient.terms = [(subquotient_tuple[0], dividend.terms[-1][1] / o_term[1])]
                    else:
                        subquotient = dividend.terms[-1][1] / o_term[1]

                    quotient += subquotient
                    subtrahend = other * subquotient
                    dividend -= subtrahend
                return quotient, dividend
        elif other == 0:
            raise ZeroDivisionError("division by zero")
        elif isinstance(other, _NUMBER_TYPES):  # divide Polynomial by a number
            return 1 / other * self, 0
        else:
            raise TypeError("Cannot divide Polynomial by object of type " + type(other))

    def __rdivmod__(self, other):
        if isinstance(other, _NUMBER_TYPES):  # divide a number by a Polynomial
            return 0, other
        else:
            raise TypeError("Cannot divide object of type " + type(other) + " by a Polynomial.")

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __rmod__(self, other):
        return divmod(other, self)[1]

    def __repr__(self):
        term_list = self.get_terms()
        repr_list = [str(abs(x[1])) + " * " + x[0].__repr__() if x[1] not in (1, -1) else x[0].__repr__()
                     for x in term_list]  # list of representations of monomials
        repr_str = ""
        if term_list.pop()[1] < 0:  # put the correct sign on the first term
            repr_str += "-"
        repr_str += repr_list.pop()

        while repr_list:  # can't use join because we have two joining characters
            if term_list.pop()[1] < 0:
                repr_str += " - "
            else:
                repr_str += " + "
            repr_str += str(repr_list.pop())
        # add constant if nonzero
        if self.constant < 0:
            repr_str += " - " + str(abs(self.constant))
        elif self.constant > 0:
            repr_str += " + " + str(self.constant)
        return repr_str

    def __hash__(self):
        return hash((tuple(self.terms), self.constant))

    def __copy__(self):
        return copy.deepcopy(self)  # polynomials are immutable, even though their attributes aren't

    def divides(self, other):
        return divmod(other, self)[1] == 0

    def evaluate(self, substitutions=None):  # plug in values for the indeterminates
        if substitutions is None:
            return copy.deepcopy(self)
        else:
            substitutions = {**{x: y for x, y in substitutions.items() if isinstance(x, str)},
                             **{x.terms[0][0].indeterminates[0]: y for x, y in substitutions.items()
                                if isinstance(x, Polynomial)}}

            out_polynomial = Polynomial(_PLACEHOLDER)
            out_polynomial.terms = [(x.evaluate(substitutions)[0], x.evaluate(substitutions)[1] * y)
                                    for x, y in self.terms]
            out_polynomial.constant = self.constant + sum([y for x, y in out_polynomial.terms if x == 1])
            out_polynomial.terms = [(x, y) for x, y in out_polynomial.terms if x not in (0, 1)]

            def coefficient_sum(x):
                return sum([j for i, j in out_polynomial.terms if i == x])

            monomial_list = list(set([x for x, y in out_polynomial.terms]))
            out_polynomial.terms = [(x, coefficient_sum(x)) for x in monomial_list]
            out_polynomial.terms.sort()
            return out_polynomial

    def homogenize(self, indeterminate):  # get a homogenous Polynomial using indeterminate
        # if 'indeterminate' is an indeterminate, save indeterminate to 't' as a Monomial
        if isinstance(indeterminate, Polynomial):
            if indeterminate.constant == 0 and len(indeterminate.terms) == 1 and indeterminate.terms[0][1] == 1 \
                    and indeterminate.terms[0][0].get_degree() == 1:
                t = indeterminate.terms[0][0]
            else:
                raise ValueError("input Polynomial is not an indeterminate")
        elif isinstance(indeterminate, Monomial):
            if indeterminate.get_degree() == 1:
                t = indeterminate
            else:
                raise ValueError("input Monomial is not an indeterminate")
        elif isinstance(indeterminate, str):
            t = Monomial(indeterminate)
        else:
            raise TypeError("indeterminate must be a Polynomial, Monomial, or string, not an object of type "
                            + type(indeterminate))

        # multiply every term by the proper power of 't' in order to homogenize
        out_polynomial = copy.deepcopy(self)
        degree = out_polynomial.get_degree()
        out_polynomial.terms = [(x * t ** (degree - x.get_degree()), y) for x, y in out_polynomial.terms]
        if out_polynomial.constant != 0:
            out_polynomial.terms.append((t ** degree, out_polynomial.constant))
        out_polynomial.constant = 0
        out_polynomial.terms.sort()
        return out_polynomial

    def normalize(self):  # make the largest term of the Polynomial have coefficient 1
        return self / self.terms[-1][1]

    def derivative(self, indeterminate):
        indet_str = indeterminate.terms[0][0].indeterminates[0]
        p = self - self.constant
        term_list = [(x.derivative(indeterminate)[0], x.derivative(indeterminate)[1] * y) for x, y in p.terms]

        p.constant = [y for x, y in term_list if x == 1]
        if p.constant:
            p.constant = p.constant[0]
        else:
            p.constant = 0
        p.terms = [(x, y) for x, y in term_list if isinstance(x, Monomial)]

        if p.terms:
            return p
        else:
            return p.constant
