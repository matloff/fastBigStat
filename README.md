
# fastBigStat

Author: Norm Matloff, UC Davis

A fast introduction to classical large-sample methods in statistics.
See also [my fast introduction to
statistics](https://github.com/matloff/fastStat).

&mu;

## Prerequisites 

* Calculus-based probability theory -- random variables, cdfs, density
  functions, expected value etc.

* Linear algebra.

* Familiarity with limits.

## Notation

* Vectors will be in column form. 

* Prime (') will indicate matrix transpose.  

* Cov will represent the covariance matrix of a random vector.

# Some Preliminaries 

* Consider a random vector X = (X<sub>1</sub>,...,X<sub>p</sub>)', 
  and a constant m x p matrix A.

* E(AX) = A E(X).

* Cov(AX) = A Cov(X) A'.

* Let &Sigma; be a covariance matrix. As a real, symmetric matrix, it
  can be diagonalized, with the associated diagonal matrix having
  nonnegative elements. In this way, one can construct a matrix W such
  that W<sup>2</sup> = &Sigma;, and thus write W = D<sup>1/2</sup>.

# The Multivariate Normal Distribution Family

* Consider a random vector X = (X<sub>1</sub>,...,X<sub>p</sub>)'. Mean
  is now a vector &mu;, and standard deviation &sigma; now becomes the
  p x p covariance matrix &Sigma;.  Notation: N(&mu;,&Sigma;).

* Density: c exp[(t-&mu;)' &Sigma;<sup>-1</sup> (t-&mu;)], where the
  constant c makes the p-fold integral equal to 1.  This is for the case
  in which &Sigma; is invertible. 

* We will always assume &Sigma; is inverible  unless stated otherwise.
  It simply means that no X<sub>i</sub> a linear combination of the
  others.

* For p = 1 we have the familiar "bell-shaped curve", while
  for p = 2, we get a 3-dimensional bell.

![bivariate normal density](Bell.png)

* Every marginal distribution of X is also multivariate normal.

* For constant matrix A as above, AX is also multivariate normal. In
  particular, if m = 1, then AX is univariate normal; the converse is
  also true, i.e. if a'X is normal for all a, then X is multivariate
  normal.

* Conditional distributions, e.g. X<sub>2</sub> |
  (X<sub>1</sub>,X<sub>3</sub>)' = t, are linear in the 
  conditioning value t.  Also, the variance is constant in t.
  (Source of the standard linear model.)

* (X-&mu;)' &Sigma;<sup>-1</sup> (X-&mu;) has a chi-squared distribution
  with p degrees of freedom.

# Limit Theorems

* Multivariate Central Limit Theorem: Let X<sub>i</sub>, i = 1,2,... be
  iid random vectors with mean &mu; and covariance matrix &Sigma;.
  Define T<sub>n</sub> = X<sub>1</sub>+...+X<sub>n</sub>, which has mean
  and covariance matrix n &mu; and n &Sigma;. Then as n &rarr; &infin; the 
  (p-dimensional) cdf of 

    1/n<sup>1/2</sup> &Sigma;<sup>-1/2</sup> (T<sub>n</sub> - n &mu;)

  goes to the cdf of N(0,I), where I is the identity matrix of size p.

* A sequence of random variables V<sub>i</sub> is said to converge to a
  random variable V if for every &epsilon; > 0,

  lim<sub>n &rarr; &infin;</sub> P(|V<sub>n</sub> - V| > &epsilon;) = 0
 
* For random vectors V<sup>i</sup>, replace | | by, e.g. Euclidean
  distance.

