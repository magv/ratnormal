# RATNORMAL

*Ratnormal* computes a rational normal form of a sum of rational
expressions (products of integer powers of polynomials with integer
coefficients). It tries to do this fast.

# BUILDING

To build *Ratnormal* first install GiNaC [1] and FLINT [2] libraries,
and then run:

    make

[1] https://www.ginac.de/

[2] https://flintlib.org/

# MANUAL

## NAME

`ratnormal` -- compute rational normal form of a sum of rational
expressions.

## SYNOPSYS

`ratnormal` [-hC] [-j *threads*] *input* [*output*]

## DESCRIPTION

`ratnormal` reads an input file, which must contain a sum of
rational expressions (products of integer powers of polynomials with
integer coefficients), and brings the sum under a common denominator,
cancelling the fractions if needed, and preserving the factorization
of the denominator.

`ratnormal` uses the FLINT library.

## EXAMPLES

    $ echo '2/3/(x + 2*y^2)^2 - 1/(x + 6*x*y^2)' | ratnormal -
    (
     1/3 *
     1/(6*x*y^2+x) *
     1/(x+2*y^2)^2 *
     (-3*x^2+2*x-12*y^4)
    )

## OPTIONS

* -j *threads*

  Use this many worker threads, if possible (default: 1).

* -h

  Show this help message.

* -C

  Force colored output even if stderr is not a tty.

* -V

  Print version information.

## ARGUMENTS

* *input*

  Read the input expression from this file, with "-"
  meaning the standard input.

* *output*

  Write the result into this file, with "-" meaning
  the standard output (which is the default).

## AUTHORS

Vitaly Magerya
