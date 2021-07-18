# Polynomials.py
User friendly polynomial package.
You can just define indeterminates and do math
with them.

## Monomial class
Not generally meant for user interaction. 
This class handles individual monomials.

Does not include addition or subtraction, 
since monomials are not closed under these operations
(the result is a polynomial).

objects from the Monomial class represent
power products, i.e. monomials without coefficients.
We will let the Polynomial class deal with coefficients

Monomials have 3 attributes:
    **degree, indeterminates,** and **exponents**

Monomials must always be initialized with an indeterminate.

**divmod** behaves oddly with Monomials.
Since we do not include coefficients, dividing
a constant other than 0 by a Monomial will result in
a remainder of 1.



