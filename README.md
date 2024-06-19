
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

* The identity matrix, of whatever given size, will be denoted by I.

# Some Preliminaries 

* Consider a random vector X = (X<sub>1</sub>,...,X<sub>p</sub>)', 
  and a constant m x p matrix A.

* E(AX) = A E(X).  Note that these are vector, i.e. componentwise
  means.

* Covariance

    - For random variables U and V, with means &alpha; and &beta;,
      define 

      Cov(U,V) = E[(U-&alpha;)(V-&beta;)]

      Intuitively, if often when U is above its mean then V is also
      above its mean, the the covariance will be positive, etc. It's
      like correlation, and in fact the latter is

      &rho;(U,V) = Cov(U,V) / [sqrt(Var(U) sqrt(Var(V))]

    - The covariance matrix Cov(X) of a random vector 
      X = (X<sub>1</sub>,...,X<sub>p</sub>)' is a p x p
      matrix with (i,j) element Cov(X<sub>i</sub>,X<sub>j</sub>),
      which reduces to Var(X<sub>i</sub>) when i = j.

    - In matrix terms, Cov(X) = [(X - &mu;) (X - &mu;)'], where &mu; =
      E(X).

    - Cov(AX) = A Cov(X) A'.

    - Let &Sigma; be a covariance matrix. As a real, symmetric matrix,
      it can be diagonalized, i.e. &Sigma; = P D P' for an orthogonal
      matrix O. Let W = PD<sup>1/2</sup>P', where the square root
      matrix simply takes the square roots of the diagonal elements of
      D. Then W<sup>2</sup> = &Sigma;, and we can thus write 
      W = &Sigma;<sup>1/2</sup>.

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

* M = (X-&mu;)' &Sigma;<sup>-1</sup> (X-&mu;) has a chi-squared distribution
  with p degrees of freedom.

* So if &eta; is the q-th quantitle of the &chi;<sup>2</sup>
  distribution with p degrees of freedom, then

  P(M < &eta;) = q

* In the above picture, horizontal slices, i.e. level sets, are
  ellipses. For p > 2, the regions are p-dimensional ellipsoids. This
  can be used to derive confidence regions in statistical applications,
  as will be seen below.

# Types of Convergence

* A sequence of random variables V<sub>n</sub> is said to converge 
  *in probability* to a random variable V if for every &epsilon; > 0,

  lim<sub>n &rarr; &infin;</sub> P(|V<sub>n</sub> - V| > &epsilon;) = 0
 
* For random vectors V<sup>i</sup>, replace | | by, e.g. Euclidean
  distance.

* Say we have random variables Q<sub>n</sub>, not necessarily iid. If
  for some constant c we have

  P(lim<sub>n &rarr; &infin;</sub> Q<sub>n</sub> = c) = 1

  then we say Q<sub>n</sub> converges *almost surely* to c.

* Example: the Strong Law of Large Numbers (SLLN). If the Q<sub>n</sub> are
  iid with mean &mu;, then the sample average converges to the
  distributional average: Set

  A<sub>n</sub> = (1/n) (V<sub>1</sub>+...+V<sub>n</sub>)

  Then A<sub>n</sub> converges almost surely to &mu;.

* If cdfs converge, i.e. for each t,

   lim<sub>n &rarr; &infin;</sub> P(Q<sub>n</sub> &le; t) = P(Q &le; t) 

  for some random variable Q, we say Q<sub>n</sub> converges *in
  distribution* to Q.

* Some types of convergence imply others:

  almost sure => in probability => in distribution 

* In p dimensions, we say that X<sub>i</sub> converges *in distribution* 
  to X if the p-dimensional cdfs converge.

* That definition is hard to deal with, but the Cramer-Wold device often
  makes things easy: Consider a sequence of random vectors X<sub>n</sub>
  and another random vector X. Then X<sub>n</sub> converges in
  distribution to X if and only if for any p-vector c, c'X<sub>n</sub>
  converges in distribution to c'X.

* Slutsky Theorem (due to Russian mathematician Evgeny Slutsky). Say
  X<sub>n</sub> converges in distribution to X, and Y<sub>n</sub>
  converges in probability to a constant c. The:

    - (i) X<sub>n</sub> + Y<sub>n</sub> converges in distribution to X+c.

    - (ii) X<sub>n</sub> Y<sub>n</sub> converges in distribution to Xc.

    - (iii) X<sub>n</sub> / Y<sub>n</sub> converges in distribution to X/c.
     
  The quantities can be matrix valued, in which case u/v means u
  v<sup>-1</sup>, i.e. matrix inverse.  Of course, the inverse must
  exist for all this to be valid.

# Central Limit Theorems (CLTs)

* Univariate Central Limit Theorem: Let X<sub>i</sub>, i = 1,2,... be
  iid random variables with mean &mu; and variance &sigma;<sup>2</sup>.
  Define T<sub>n</sub> = X<sub>1</sub>+...+X<sub>n</sub>, which has mean
  and variance n &mu; and n &sigma;<sup>2</sup>. Then 

    1/n<sup>1/2</sup> &sigma;<sup>-1</sup> (T<sub>n</sub> - n &mu;)

  converges in distribution to N(0,1).

* Note that while scaling by the factor n<sup>1/2</sup> is what works in
  this context of sums, other rates may apply for other quantities. For
  example, if we are interested in the maximum of
  X<sub>1</sub>,...,X<sub>n</sub> rather than their sum, then 
  we may divide by (2 log n)<sup>1/2</sup> (and the limiting 
  distribution will be something other than normal).

* Multivariate Central Limit Theorem: Let X<sub>i</sub>, i = 1,2,... be
  iid random vectors with mean vector &mu; and covariance matrix &Sigma;.
  Define T<sub>n</sub> = X<sub>1</sub>+...+X<sub>n</sub>, which has mean
  and covariance matrix n &mu; and n &Sigma;. Then 

    1/n<sup>1/2</sup> &Sigma;<sup>-1/2</sup> (T<sub>n</sub> - n &mu;)

  converges in distribution to  N(0,I). This can easily be proved using
  the Cramer-Wold device.

# Example: Linear Model with Fewer Assumptions

* Let Y be a random variable and X a p-dimensional random vector. The classic
  assumptions of the model are as follows. Let t = 
  (t<sub>1</sub>,...,t<sub>p</sub>)'. 

  - Linearity: E(Y | X = t) = 
  &beta;<sub>0</sub> + 
  &beta;<sub>1</sub> t<sub>1</sub> + ... +
  &beta;<sub>p</sub> t<sub>p</sub>. 
  for some unknown constant vector &beta; = 
  (&beta;<sub>0</sub>, &beta;<sub>1</sub>+...+ &beta;<sub>p</sub>)'.

  - Normality: The distribution of Y | X = t is normal for all t.

  - Heteroscedasticity: Var(Y | X = t) is the same for all t.

  - Independence: The random vectors (X<sub>1</sub>,Y<sub>1</sub>),...,
  (X<sub>n</sub>,Y<sub>n</sub>) are independent.

  (Here we have the X<sub>i</sub> random. In many books, they are
  fixed.)

* One estmates &beta; from the data, using the formula for the estimates
  b = (b<sub>0</sub>, b<sub>1</sub>,...,b<sub>p</sub>)':

  b = (A'A)<sup>-1</sup> A'W

  where A is the n x (p+1) matrix will all 1s in column 1 and in column
  j > 1, all n values of X<sub>j-1</sub>.  Here W is the vector of all n
  values of Y.

* The estimate vector b has an exact multivariate normal distribution
  with mean &beta; and covariance matrix (A'A)<sup>-1</sup>. This
  enables exact confidence intervals and tests. But we can use the
  Central Limit Theorem to get approximate intervals and tests, as
  follows. (Note: The notion of "exact" statistical inference is a myth.
  No distribution in real life is exactly normal etc.)

* To make things a little simpler to explain, we'll look at the
  non-intercept model, in which E(X) is 0 and   

  E(Y | X = t) = 
  &beta;<sub>1</sub> t<sub>1</sub> + ... +
  &beta;<sub>p</sub> t<sub>p</sub>. 

  (One can easily convert the full model to the interceptless one by
  subtracting means from X and Y.) So the A matrix is now n x p.

* To prepare for the CLT, rewrite b:

   b = [(1/n) A'A]<sup>-1</sup> [(1/n) A'W]



* Let's first look at the "denominator" [(1/n) A'A]<sup>-1</sup>. The
  number in row i, column j of A'A is the inner product of all n values
  of X<sub>i</sub> and X<sub>j</sub>. (Check this.) That number is a sum
  of n independent values, so the SLLN applies. Thus (1/n) A'A converges
  almostt surely (and thus in probability) to E(XX'). (Again, check
  this.) So, (1/n) A'A will play the role of "Y<sub>n</sub>" in part
  (iii) of the Slutsky Theorem.
