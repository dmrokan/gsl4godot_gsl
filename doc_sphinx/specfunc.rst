*****************
Special Functions
*****************

.. index:: special functions

This chapter describes the GSL special function library.  The library
includes routines for calculating the values of Airy functions, Bessel
functions, Clausen functions, Coulomb wave functions, Coupling
coefficients, the Dawson function, Debye functions, Dilogarithms,
Elliptic integrals, Jacobi elliptic functions, Error functions,
Exponential integrals, Fermi-Dirac functions, Gamma functions,
Gegenbauer functions, Hypergeometric functions, Laguerre functions,
Legendre functions and Spherical Harmonics, the Psi (Digamma) Function,
Synchrotron functions, Transport functions, Trigonometric functions and
Zeta functions.  Each routine also computes an estimate of the numerical
error in the calculated value of the function.

The functions in this chapter are declared in individual header files,
such as :file:`gsl_sf_airy.h`, :file:`gsl_sf_bessel.h`, etc.  The complete
set of header files can be included using the file :file:`gsl_sf.h`.

Usage
=====

The special functions are available in two calling conventions, a
@dfn{natural form} which returns the numerical value of the function and
an @dfn{error-handling form} which returns an error code.  The two types
of function provide alternative ways of accessing the same underlying
code.

The @dfn{natural form} returns only the value of the function and can be
used directly in mathematical expressions.  For example, the following
function call will compute the value of the Bessel function
@math{J_0(x)},

@example
double y = gsl_sf_bessel_J0 (x);
@end example

There is no way to access an error code or to estimate the error using
this method.  To allow access to this information the alternative
error-handling form stores the value and error in a modifiable argument,

@example
gsl_sf_result result;
int status = gsl_sf_bessel_J0_e (x, &result);
@end example

The error-handling functions have the suffix @code{_e}. The returned
status value indicates error conditions such as overflow, underflow or
loss of precision.  If there are no errors the error-handling functions
return @code{GSL_SUCCESS}.

The gsl_sf_result struct
========================
.. index:: gsl_sf_result, gsl_sf_result_e10

@tindex gsl_sf_result @tindex gsl_sf_result_e10

The error handling form of the special functions always calculate an
error estimate along with the value of the result.  Therefore,
structures are provided for amalgamating a value and error estimate.
These structures are declared in the header file :file:`gsl_sf_result.h`.

The @code{gsl_sf_result} struct contains value and error fields.

@example
typedef struct
@{
  double val;
  double err;
@} gsl_sf_result;
@end example

The field @var{val} contains the value and the field @var{err} contains
an estimate of the absolute error in the value.

In some cases, an overflow or underflow can be detected and handled by a
function.  In this case, it may be possible to return a scaling exponent
as well as an error/value pair in order to save the result from
exceeding the dynamic range of the built-in types.  The
@code{gsl_sf_result_e10} struct contains value and error fields as well
as an exponent field such that the actual result is obtained as
@code{result * 10^(e10)}.

@example
typedef struct
@{
  double val;
  double err;
  int    e10;
@} gsl_sf_result_e10;
@end example

Modes
=====

The goal of the library is to achieve double precision accuracy wherever
possible.  However the cost of evaluating some special functions to
double precision can be significant, particularly where very high order
terms are required.  In these cases a @code{mode} argument allows the
accuracy of the function to be reduced in order to improve performance.
The following precision levels are available for the mode argument,

@table @code
@item GSL_PREC_DOUBLE
Double-precision, a relative accuracy of approximately @c{$2 \times 10^{-16}$}
@math{2 * 10^-16}.
@item GSL_PREC_SINGLE
Single-precision, a relative accuracy of approximately @c{$1 \times 10^{-7}$}
@math{10^-7}.
@item GSL_PREC_APPROX
Approximate values, a relative accuracy of approximately @c{$5 \times 10^{-4}$}
@math{5 * 10^-4}.
@end table

The approximate mode provides the fastest evaluation at the lowest
accuracy.

Airy Functions and Derivatives
==============================
@include specfunc-airy.texi

Bessel Functions
================
@include specfunc-bessel.texi

Clausen Functions
=================
@include specfunc-clausen.texi

Coulomb Functions
=================
@include specfunc-coulomb.texi

Coupling Coefficients
=====================
@include specfunc-coupling.texi

Dawson Function
===============
@include specfunc-dawson.texi

Debye Functions
===============
@include specfunc-debye.texi

Dilogarithm
===========
@include specfunc-dilog.texi

Elementary Operations
=====================
@include specfunc-elementary.texi

Elliptic Integrals
==================
@include specfunc-ellint.texi

Elliptic Functions (Jacobi)
===========================
@include specfunc-elljac.texi

Error Functions
===============
@include specfunc-erf.texi

Exponential Functions
=====================
@include specfunc-exp.texi

Exponential Integrals
=====================
@include specfunc-expint.texi

Fermi-Dirac Function
====================
@include specfunc-fermi-dirac.texi

Gamma and Beta Functions
========================
@include specfunc-gamma.texi

Gegenbauer Functions
====================
@include specfunc-gegenbauer.texi

Hypergeometric Functions
========================
@include specfunc-hyperg.texi

Laguerre Functions
==================
@include specfunc-laguerre.texi

Lambert W Functions
===================
@include specfunc-lambert.texi

Legendre Functions and Spherical Harmonics
==========================================
@include specfunc-legendre.texi

Logarithm and Related Functions
===============================
@include specfunc-log.texi

Mathieu Functions
=================
@include specfunc-mathieu.texi

Power Function
==============
@include specfunc-pow-int.texi

Psi (Digamma) Function
======================
@include specfunc-psi.texi

Synchrotron Functions
=====================
@include specfunc-synchrotron.texi

Transport Functions
===================
@include specfunc-transport.texi

Trigonometric Functions
=======================
@include specfunc-trig.texi

Zeta Functions
==============
@include specfunc-zeta.texi

Examples
========

The following example demonstrates the use of the error handling form of
the special functions, in this case to compute the Bessel function
@math{J_0(5.0)},

@example
@verbatiminclude examples/specfun_e.c
@end example

Here are the results of running the program,

@example
$ ./a.out 
@verbatiminclude examples/specfun_e.txt
@end example

The next program computes the same quantity using the natural form of
the function. In this case the error term @var{result.err} and return
status are not accessible.

@example
@verbatiminclude examples/specfun.c
@end example

The results of the function are the same,

@example
$ ./a.out 
@verbatiminclude examples/specfun.txt
@end example



References and Further Reading
==============================

The library follows the conventions of @cite{Abramowitz & Stegun} where
possible,
@itemize @w{}
@item
Abramowitz & Stegun (eds.), @cite{Handbook of Mathematical Functions}
@end itemize

The following papers contain information on the algorithms used 
to compute the special functions,
@cindex MISCFUN
@itemize @w{}
@item
Allan J. MacLeod, MISCFUN: A software package to compute uncommon
special functions.  @cite{ACM Trans.@: Math.@: Soft.}, vol.@: 22,
1996, 288--301

@item
G.N. Watson, A Treatise on the Theory of Bessel Functions,
2nd Edition (Cambridge University Press, 1944).

@item
G. Nemeth, Mathematical Approximations of Special Functions,
Nova Science Publishers, ISBN 1-56072-052-2

@item
B.C. Carlson, Special Functions of Applied Mathematics (1977)

@item
N. M. Temme, Special Functions: An Introduction to the Classical
Functions of Mathematical Physics (1996), ISBN 978-0471113133.

@item
W.J. Thompson, Atlas for Computing Mathematical Functions, John Wiley & Sons,
New York (1997).

@item
Y.Y. Luke, Algorithms for the Computation of Mathematical Functions, Academic
Press, New York (1977).

@item
S. A. Holmes and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and order
normalised associated Legendre functions, Journal of Geodesy, 76,
pg. 279-299, 2002.

@comment @item
@comment Fermi-Dirac functions of orders @math{-1/2}, @math{1/2}, @math{3/2}, and
@comment @math{5/2}.  @cite{ACM Trans. Math. Soft.}, vol. 24, 1998, 1-12.
@end itemize
