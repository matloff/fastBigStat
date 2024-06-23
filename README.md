
# fastBigStat

A fast introduction to asymptotoic normality of statistical estimators.
Includes the necessary material on multivariate distributions.  See also
[my fast introduction to
statistics](https://github.com/matloff/fastStat).

Author: Norm Matloff, UC Davis; 
[bio](http://heather.cs.ucdavis.edu/matloff.html)

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

* E(AX) = A E(X).  Note that these are vector/matrix quantities, i.e.
  componentwise means.

* Covariance

    - For random variables U and V, with means &alpha; and &beta;,
      define 

      Cov(U,V) = E[(U-&alpha;)(V-&beta;)]

      Intuitively, if often when U is above its mean then V is also
      above its mean, then the covariance will be positive, etc. It's
      like correlation, and in fact the latter is a scaled version of
      covariance,

      &rho;(U,V) = Cov(U,V) / [sqrt(Var(U) sqrt(Var(V))]

    - The covariance matrix Cov(X) of a random vector 
      X = (X<sub>1</sub>,...,X<sub>p</sub>)' is a p x p
      matrix with (i,j) element Cov(X<sub>i</sub>,X<sub>j</sub>),
      which reduces to Var(X<sub>i</sub>) when i = j.

    - In matrix terms, Cov(X) = E[(X - &mu;) (X - &mu;)'], where &mu; =
      E(X).

    - Cov(AX) = A Cov(X) A'.

    - Let &Sigma; be a covariance matrix. As a real, symmetric matrix,
      it can be diagonalized, i.e. &Sigma; = P D P' for an orthogonal
      matrix P and diagonal matrix D. Let W = PD<sup>1/2</sup>P', where
      the square root matrix simply takes the square roots of the
      diagonal elements of D. Then W<sup>2</sup> = &Sigma;, and we can
      thus write W = &Sigma;<sup>1/2</sup>.

* Law of Iterated Expectation

    - Let U and V be random quantities (scalars here, for convenience)
      for which the expectations below exist. Then

      E[ E(V|U ] = E(V)

      Note that E(V|U) is treated as a random variable, which it in fact
      is, since it is a function of U.

    - Intuition: Suppose we wish to find the mean height of all UC Davis
      students. We could ask each academic department to find the mean
      height of their students, then average those mean heights. (Since
      some departments are larger than others, that latter average is
      weighted.)

* Role of the Data

    - Say we have data W<sub>1</sub>,...,W<sub>n</sub>.

    - The statistical view is to treat the data as a random sample from a
      population.

    - In computer science, the data are often viewed as fixed, but some
      analyses describe the data as coming from a "data generating
      mechanism" and the like.

    - We take the statistical view here. Typically the data are assumed
      independent and identically distributed (iid). That first trait
      arises from sampling with replacement (for large n, not much of an
      issue), and the second means that each entity in the population
      has the same probability of being sampled.

      The sample mean

      (1/n) (W<sub>1</sub>+...+W<sub>n</sub>)

      is viewed as an estimate of the population mean E(W), where W has
      the distribution of the popuulation, and so on. 

# The Multivariate Normal Distribution Family

* Consider a random vector X = (X<sub>1</sub>,...,X<sub>p</sub>)'. The mean
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

* For a constant m x p matrix A, AX is also multivariate normal. In
  particular, if m = 1, then AX is univariate normal; the converse is
  also true, i.e. if a'X is normal for all a, then X is multivariate
  normal.

* Conditional distributions:  For simplicity, consider the conditional
  distribution of X<sub>i</sub> given (X<sub>1</sub>,...,X<sub>j-1</sub>,
  X<sub>j+1</sub>,...,X<sub>p</sub>) = w.  Then (inspiration for the
  familiar linear model):

  - That distribution is normal.

  - The mean of that distribution is linear in w, i.e. of the form
    c'w for some c.

  - The variance of that distribution is constant in w.

* The *quadratic form* M = (X-&mu;)' &Sigma;<sup>-1</sup> (X-&mu;) has a 
  chi-squared distribution with p degrees of freedom.

* So if &eta; is the q-th quantitle of the &chi;<sup>2</sup>
  distribution with p degrees of freedom, then

  P(M < &eta;) = q

* In the above picture, horizontal slices, i.e. level sets, are
  ellipses. For p > 2, the regions are p-dimensional ellipsoids. This
  can be used to derive confidence regions in statistical applications,
  as will be seen below.

# Example: Linear Model

* Let Y be a random variable and X =
  (X<sub>1</sub>,...,X<sub>p</sub>)' a p-dimensional random vector. E.g.
  we might wish to predict Y = college GPA from X<sub>1</sub> = high
  school GPA, X<sub>2</sub> = SAT score and so on.

  The classic assumptions of the model are as follows. Let t =
  (t<sub>1</sub>,...,t<sub>p</sub>)'. 

  - Linearity: E(Y | X = t) = 
  &beta;<sub>0</sub> + 
  &beta;<sub>1</sub> t<sub>1</sub> + ... +
  &beta;<sub>p</sub> t<sub>p</sub>. 
  for some unknown constant vector &beta; = 
  (&beta;<sub>0</sub>, &beta;<sub>1</sub>,..., &beta;<sub>p</sub>)'.

  - Normality: The distribution of Y | X = t is normal for all t.

  - Homoscedasticity: Var(Y | X = t) is the same for all t, which
    we will denote by &sigma;<sup>2</sup>.

  - Independence: The random vectors (X<sub>1</sub>,Y<sub>1</sub>),...,
  (X<sub>n</sub>,Y<sub>n</sub>) are independent.

  Here we have the X<sub>i</sub> random. In many books, they are
  fixed. We will switch back and forth, conditioning on the X<sub>i</sub>
  at times (thus treating them as fixed).

* One estmates &beta; from the data, using the formula for the estimates
  b = (b<sub>0</sub>, b<sub>1</sub>,...,b<sub>p</sub>)':

  b = (A'A)<sup>-1</sup> A'W

  where A is the n x (p+1) matrix will all 1s in column 1 and in column
  j > 1, all n values of X<sub>j-1</sub>.  Here W is the vector of all n
  values of Y.

* Write column j+1 of A as (X<sub>j1</sub>,...,X<sub>jn</sub>)'. If say
  X<sub>j</sub> is human age, then this vector consists of the ages of
  all people in our dataset. Similarly, write W as
  (Y<sub>1</sub>,...,Y<sub>n</sub>)', the vector of Y values for
  everyone in our dataset.

* By the way, since

  E(Y | X) = X<sub>extend</sub>' &beta;

  where X<sub>extend</sub> = (1,X)'

  we have that

  E(W | A) = A &beta;

* The estimate vector b has an exact multivariate normal distribution
  with mean &beta; and covariance matrix 

  &sigma;<sup>2</sup> (A'A)<sup>-1</sup>. 

  Note that A is a function of the X data, but here we take it to be a
  constant, by considering the conditional distribution of b given A.

  Thus confidence intervals (CIs) and regions we obtain below are also
  conditional. However, note that if, e.g. we have a conditional CI at
  the 95% level, then the unconditional coverage probability is also
  95%, by the Law of Iterated Expection above.

* The unknown parameter &sigma;<sup>2</sup> must be estimated from the
  data, as 

  s<sup>2</sup> = [1/(n-p-1)] &Sigma;<sub>i</sub> 
  (Y<sub>i</sub> - R<sub>i</sub> b)<sup>2</sup>

  where R<sub>i</sub> is row i of A.

* So the estimated conditional covariance matrix of b is

   s<sup>2</sup> (A'A)<sup>-1</sup>

* The stringent assumptions above enable exact statistical inference.
  This enables exact confidence intervals and tests. 

  - The quantity b<sub>i</sub> - &beta;<sub>i</sub> has a Student
    t-distribution with n-p-1 df, thus setting up a CI for
    &beta;<sub>i</sub>. For large n, Student t is basically the N(0,1)
    distribution.

  - In some cases, we may be interested in something like, say,

    &beta;<sub>2</sub> - &beta;<sub>1</sub> = c'&beta;

    where c = (0,1,-1,0,0,...,0)'. Since c'b will then have a univariate
    normal distribution, we again can easily obtain a CI for the desired
    quantity. Note that we need the estimated variance of c'b, which is
    
    s<sup>2</sup> c'(A'A)<sup>-1</sup>c

    So for large n, a 95% confidence interval for 
     &beta;<sub>2</sub> - &beta;<sub>1</sub> is 

    c'b &plusmn; 1.96 [s<sup>2</sup> c'(A'A)<sup>-1</sup>c]<sup>0.5</sup>

  - The quantity 

    (b-&beta;)'[s<sup>-2</sup> (A'A)] (b-&beta;)

    has an F-distribution with (p+1,n-p-1) df. This is the
    finite-sample version of the quadratic form discussed earlier.
    But again for large n, this is approximately &chi;<sup>2</sup> with 
    p+1 df. This sets up a confidence ellipsoid for &beta;:

    To form an approximate (1-&alpha;) 100% confidence ellipsoid for
    &beta;, let q be the upper-&alpha; quantile of the &chi;<sup>2</sup>
    distribution with p+1 df. Then the confidence ellipsoid is the set
    of all t such that 

    (b-t)'[s<sup>-2</sup> (A'A)] (b-t) &le; q

# Types of Convergence

* A sequence of random variables V<sub>n</sub> is said to converge 
  *in probability* to a random variable V if for every &epsilon; > 0,

  lim<sub>n &rarr; &infin;</sub> P(|V<sub>n</sub> - V| > &epsilon;) = 0
 
* For random vectors V<sub>i</sub>, replace | | by, e.g. Euclidean
  distance.

* Say we have random variables Q<sub>n</sub>, not necessarily iid. If
  for some constant c we have

  P(lim<sub>n &rarr; &infin;</sub> Q<sub>n</sub> = c) = 1

  then we say Q<sub>n</sub> converges *almost surely* to c. (The term
  "almost surely" is just a fancy way of saying, "With probability 1.")

* Example: the Strong Law of Large Numbers (SLLN). If the Q<sub>n</sub> are
  iid with mean &mu;, then the sample average converges to the
  distributional average: Set

  A<sub>n</sub> = (1/n) (Q<sub>1</sub>+...+Q<sub>n</sub>)

  Then A<sub>n</sub> converges almost surely to &mu;.

* If cdfs of random variables converge, i.e. for each t,

   lim<sub>n &rarr; &infin;</sub> P(Q<sub>n</sub> &le; t) = P(Q &le; t) 

  for some random variable Q, we say Q<sub>n</sub> converges *in
  distribution* to Q. Sometimes people omit mention of Q, simply
  referring to its distribution.

* Some types of convergence imply others:

  almost sure => in probability => in distribution 

* In p dimensions, we say that X<sub>i</sub> converges *in distribution* 
  to X if the p-dimensional cdfs converge.

* That definition is hard to deal with, but the Cramer-Wold device often
  makes things easy: Consider a sequence of random vectors X<sub>n</sub>
  and another random vector X. Then X<sub>n</sub> converges in
  distribution to X if and only if for all p-vectors c, c'X<sub>n</sub>
  converges in distribution to c'X.

* Slutsky Theorem (due to Russian mathematician Evgeny Slutsky). Say
  X<sub>n</sub> converges in distribution to X, and Y<sub>n</sub>
  converges in probability to a constant c. Then:

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

  Many generalizations of the CLT exist, in which they relax the iid
  assmuption.

* Scaling by the factor n<sup>1/2</sup> is "just right" to make things
  work here. If, say we were to divide instead by n, the quantity would
  go to 0, by the SLLN. That certainly would not be helpful in terms of
  finding probabilities, confidence intervals and so on.

  Other powers may apply for other quantities. For
  example, if we are interested in the maximum of
  X<sub>1</sub>,...,X<sub>n</sub> rather than their sum, then 
  we divide by (2 log n)<sup>1/2</sup> (and the limiting 
  distribution will be something other than normal).

* Multivariate Central Limit Theorem: Let X<sub>i</sub>, i = 1,2,... be
  iid random vectors with mean vector &mu; and covariance matrix &Sigma;.
  Define T<sub>n</sub> = X<sub>1</sub>+...+X<sub>n</sub>, which has mean
  and covariance matrix n &mu; and n &Sigma;. Then 

    1/n<sup>1/2</sup> (T<sub>n</sub> - n &mu;)

  converges in distribution to  N(0,&Sigma;). This can easily be proved using
  the Cramer-Wold device.

# Examples of Asymptotic Normality 

We will show two applications of the above material. Both are serious
applications, and thus rather involved. The reader's patience is
requested.

# Example: Method of Moments (MM) Estimators

* Background on MM estimators

  - Suppose we are modeling a random variable with some parametric
    distribution family, such as normal, exonential or gamma. We want to
    estimate the parameter(s) from our data
    X<sub>1</sub>,...,X<sub>n</sub>. How can we do this?

  - The most famous general method is *Maximum Likelihood Estimation*
    (MLE), but it's somewhat less-famous cousin, the *Method of Moments* 
    (MM), is sometimes the handier one. (It is also the easier one to
    explain.)

  - This is best introduced by example. Say we are modeling our
    X<sup>i</sup> as having an exponential distribution, i.e. they have
    density

    &tau; exp(-&tau;t).

    for some unknown parameter &tau;.  How do we construct the MM
    estimate T of &tau;?

    The mean of this distribution is 1/&tau;, which would be estimated by 
    1/T. On the other hand, the classical estimate of a population mean
    is the sample mean,

    X<sub>bar</sub> = X<sub>1</sub>,...,X<sub>n</sub>. So set our
    estimated mean under the exponential assumption, 1/T, to our general
    estimated mean, without that assmuption:

    1/T = X<sub>bar</sub>

    That gives us our MM estimate,

    T = 1/X<sub>bar</sub>

  - How does MM work in general, when we have more than one parameter?
    Let's use k to denote the number of parameters in the given parametric
    family. For instance, k = 2 for the univariate normal family,
    corresponding to the two parameters &mu; and &sigma;<sup>2</sup>. We
    then must bring in E(X<sup>2</sup>) and so on, as follows.

  - The r-th *moment* of a random variable X is defined to be
    m<sub>r</sub> = E(X<sup>r</sup>), r = 1,2,....  We will need to use
    the first k moments.

    Note that this is a population value;
    its intuitive sample estimate is 

    M<sub>r</sub> = 
    (1/n) (X<sub>1</sub><sup>r</sup>+...+X<sub>n</sub><sup>r</sup>)

  - Alternatively, for k > 1 one may use the r-th *central* moment, 

    E[(X-E(X))<sup>r</sup>]

    Since for r = 2 this is Var(X) and for r = 3 we get the skewness of the
    distribution, things may be more convenient this way.

    It will be clear from context whether we mean the central or noncentral
    moment. 

  - Use &tau; to denote the population value of the parameter, and T to 
    denote its sample estimate under MM.  Both are vectors if k > 1.

  - Before continuing, note again what happened when we set 1/T = A
    above..  The left-hand side is our estimate of m<sub>1</sub> under
    the exponential model, while the right-hand side is a general
    estimator of m<sub>1</sub>, not assuming that model This is how the
    Method of Moments works, by matching the two. We  will have k such
    equations, then solve them for T.

  - For an example with k = 2, consider the *gamma distribution family*:

    c &tau;<sub>1</sub><sup>&tau;<sub>2</sub></sup>
    t<sup>&tau;<sub>2 </sub> - 1</sup></sup>
    exp(-&tau;<sub>1</sub>t)

    where c is a constant to make the density integrate to 1.0. This is
    a common model used in network communications and medical survival
    analysis for example.  Here

    E(X) = &tau;<sub>2</sub> / &tau;<sub>1</sub>

    and

    Var(X) = &tau;<sub>2</sub> / &tau;<sub>1</sub><sup>2</sup>

  - So, we equate:

    M<sub>1</sub> = T<sub>2</sub> / T<sub>1</sub>

    and (using the central moment for M<sub>2</sub>)

    M<sub>2</sub> = T<sub>2</sub> / T<sub>1</sub><sup>2</sup>

  - These are nonlinear equations, but we are lucky here; they are
    solvable from simple algebra (otherwise we would need to use
    numerical approximation methods):

    T<sub>1</sub> = M<sub>1</sub> / M<sub>2</sub>

    and

    T<sub>2</sub> = M<sub>1</sub> T<sub>1</sub> =
    M<sub>1</sub><sup>2</sup> / M<sub>2</sub>


* Consistency

  - Say we have an estimator &theta;<sub>n</sub> for a parameter &theta;
    based on a sample of size n. As n goes to &infin;, we would like our
    estimator to have the property that &theta;<sub>n</sub> goes to
    &theta;. If it does almost surely, we say it is *strongly consistent*; if 
    the convergence is just in probability, it is known as *weak
    consistency*.

    Since the SLLN implies that M<sub>i</sub> is strongly consistent for
    m<sub>i</sub>, this implies the same for the T<sub>i</sub>, as long
    as the moments are continuous functions of the M<sub>i</sub>.  In
    the example above, for instance, M<sub>1</sub> and M<sub>2</sub>
    converge almost surely, hence the T<sub>i</sub> do too.

* Asymptotic normality

  - OK, we now have estimators for &tau;<sub>1</sub> and
    &tau;<sub>2</sub>, and they are consistent. But the latter is a very
    weak property. For good practical value, it would be desirable to
    obtain confidence intervals from them. For this we need asymptotic
    normality, as follows.

  - For simplicity, let's consider the case k = 1, so the density family
    has a single scalar parameter, &tau;, as we saw in the exponential
    example above. Then E(X) is some function g(&tau;). In the
    exponential example, g(&tau;) = 1/&tau;.

    Assume g is differentiable with a continuous derivative g<sub>1</sub>. Then
    Taylor's Theorem from calculus says

    g(T) = g(&tau;) + g<sub>1</sub>(T<sub>betw</sub>) (T - &tau;)

    for some T<sub>betw</sub> between &tau; and T. 

  - Now to set up using the CLT, recall that we will set

    g(T) = M<sub>1</sub>

    Then rewrite the earlier equation as

    n<sup>0.5</sup>(M<sub>1</sub> - m<sub>1</sub>) =
    n<sup>0.5</sup> [g(T) - g(&tau;)] =
    n<sup>0.5</sup> [g<sub>1</sub>(T<sub>betw</sub>) (T - &tau;)]

    and then

    n<sup>0.5</sup>(T - &tau;) =
    n<sup>0.5</sup> (M<sub>1</sub> - m<sub>1</sub>)/ 
    g<sub>1</sub>(T<sub>betw</sub>)

 -  We know that T is a strongly consistent estimator of &tau;, so by
    the continuity of g<sub>1</sub>, the denominator in the RHS
    converges almost surely to g<sub>1</sub>(&tau;). Applying the
    Slutsky Theorem and the CLT (M<sub>1</sub> is a sum of iid terms),
    we see that the RSH is asymptotically normal (though not with
    standard deviation 1).

    The variance in that asymptotically normal distribution will be

    Var(M<sub>1</sub>/[g<sub>1</sub>(&tau;)]<sup>2</sup>

  - To compute a CI, we replace the quantities by their sample analogs.
    E.g. our estimate of Var(M<sub>1</sub>) is

    s<sup>2</sup>= (1/n) &Sigma;<sub>i=1</sub><sup>n</sup> 
    [X<sub>i</sub> - M<sub>1</sub>]<sup>2</sup>

    (Divide by n-1 instead of n if you prefer, though there really is no
    reason to do so.)

    Our estimate of g<sub>1</sub>(&tau;) is 

    g<sub>1</sub>(T) = -1/T<sup>2</sup>

    Our CI, say for a 95% confidence level, is then 

    T &plusmn; 1.96 s/T<sup>2</sup>

  - For the case k = 2, we have two functions, g and h, corresponding to
    M<sub>1</sub> and M<sub>2</sub>, both having arguments T<sub>1</sub>
    and T<sub>2</sub>. We now have derivatives g<sub>1</sub> and
    h<sub>1</sub>, but they are now gradients.  The Multivariate CLT is
    applied and so on.



