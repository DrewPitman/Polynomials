# Polynomials.py

A package for polynomials defined in terms of indeterminates.

## Monomial class

This class handles individual monomials without coefficients, i.e., power products. We will let the Polynomial class
deal with coefficients.  
The class has 3 attributes, which are self-explanatory:

* **indeterminates**; a list of strings naming the monomial's indeterminates
* **exponents**; a dictionary whose keys are indeterminate names and whose values are the powers of each indeterminate

These can be accessed with *get_** methods.

### usage notes

Monomials must always be initialized with an indeterminate. users This class does not include addition or subtraction,
since the set of monomials is not closed under these operations.

*divmod* behaves strangely with Monomials. Since we do not include coefficients, dividing a constant other than 0 by a
Monomial will result in a remainder of 1. The quotient will be 0 if the dividend is not divisible by the divisor. The
remainder will be 0 or self in all cases. The modulus operators are also based on *divmod.*

## Polynomial class

This class handles polynomials in any number of variables.

This class has two attributes:

* **constant**; the polynomial's constant term
* **terms**; a list of monomial, coefficient pairs representing the rest of the polynomial's terms

### usage notes

Polynomials must always be initialized with an indeterminate. The purpose of this is so that users can make polynomials
by first defining indeterminates then operating on them with addition, multiplication, and exponentiation.

### evaluate()

use Monomial.evaluate() or Polynomial.evaluate() to

* evaluate objects for given values of some indeterminates
* substitute some indeterminates for others

