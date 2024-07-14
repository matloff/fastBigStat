<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [fastBigStat](#fastbigstat)
   * [Prerequisites ](#prerequisites)
   * [Notation](#notation)
- [Some Preliminaries ](#some-preliminaries)
- ["Exact" vs. Approximate Statistical Inference: Historical Perspective](#exact-vs-approximate-statistical-inference-historical-perspective)
- [The Multivariate Normal Distribution Family](#the-multivariate-normal-distribution-family)
- [Example: Linear Models, Exact and Approximate Inference](#example-linear-models-exact-and-approximate-inference)
   * [Classical linear model](#classical-linear-model)
   * [*Fixed-X* vs. *random-X* models.](#fixed-x-vs-random-x-models)
   * [Estimation ](#estimation)
   * [Exact inference](#exact-inference)
   * [Approximate inference](#approximate-inference)
   * [Relaxing the assumptions](#relaxing-the-assumptions)
- [Multiple Inference Procedures](#multiple-inference-procedures)
   * [The Bonferroni Inequality](#the-bonferroni-inequality)
   * [The Scheffe' method](#the-scheffe-method)
- [Simulation of MV Normal Random Vectors](#simulation-of-mv-normal-random-vectors)
- [Types of Convergence](#types-of-convergence)
- [Central Limit Theorems (CLTs)](#central-limit-theorems-clts)
- [Example: Method of Moments (MM) Estimators](#example-method-of-moments-mm-estimators)
- [The Delta Method](#the-delta-method)

<!-- TOC end -->





<!-- TOC --><a name="fastbigstat"></a>
# fastBigStat

A *fast, practical introduction* to asymptotic normality of statistical
estimators.

Includes the necessary support material on random vectors and the
multivariate (MV) normal (Gaussian) distribution family, as well as some
insight into statistical inference based on normal distributions.

See also [my fast introduction to
statistics](https://github.com/matloff/fastStat).

Author: Norm Matloff, UC Davis; 
[bio](http://heather.cs.ucdavis.edu/matloff.html)

<!-- TOC --><a name="prerequisites"></a>
## Prerequisites 

* Calculus-based probability theory -- random variables, cdfs, density
  functions, expected value etc.

* Linear algebra.

* Familiarity with limits.

<!-- TOC --><a name="notation"></a>
## Notation

* Vectors will be in column form. 

* Prime (') will indicate matrix transpose.  

* The identity matrix, of whatever given size, will be denoted by I.

<!-- TOC --><a name="some-preliminaries"></a>
# Some Preliminaries 

* Consider a random vector X = (X<sub>1</sub>,...,X<sub>p</sub>)'. 

* E(X) means the vector of componentwise means E(X<sub>i</sub>).

* Covariance

    - For random variables U and V, with means &alpha; and &beta;,
      define 

      Cov(U,V) = E[(U-&alpha;)(V-&beta;)]

      Intuitively, if often when U is above its mean &alpha; then V is also
      above its mean &beta;, then the covariance will be positive, etc. It's
      like correlation, and in fact the latter is a scaled version of
      covariance,

      &rho;(U,V) = Cov(U,V) / [sqrt(Var(U) sqrt(Var(V))]

    - The covariance matrix Cov(X) of a random vector 
      X = (X<sub>1</sub>,...,X<sub>p</sub>)' is a p x p
      matrix with (i,j) element Cov(X<sub>i</sub>,X<sub>j</sub>),
      which reduces to Var(X<sub>i</sub>) when i = j.

    - In matrix terms, Cov(X) = E[(X - &mu;) (X - &mu;)'], where &mu; =
      E(X).

    - Let &Sigma; be a covariance matrix. As a real, symmetric matrix,
      it can be diagonalized, i.e. &Sigma; = P D P' for an orthogonal
      matrix P and diagonal matrix D. Let W = PD<sup>1/2</sup>P', where
      the square root matrix simply takes the square roots of the
      diagonal elements of D. Then W<sup>2</sup> = &Sigma;, and we can
      thus write W = &Sigma;<sup>1/2</sup>.

* Law of Iterated Expectation

    - Let U and V be random quantities (scalars here, for convenience)
      for which the expectations below exist. Then

      E[ E(V|U) ] = E(V)

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

    - Note that this means that each W<sub>i</sub> has the distribution
      of the population. E.g. if 29% of the units (e.g. people) in the
      population have value below 221.4, then P(W<sub>i</sub> < 221.4) =
      0.29.0.

      Similarly, if the population quantity is normally distributed,
      that means each W<sub>i</sub> has a normal distribution (with the
      same mean and variance as the population).

    - The sample mean

      W<sub>bar</sub> = (1/n) (W<sub>1</sub>+...+W<sub>n</sub>)

      being the sample analog of the population mean, is a natural
      estimate (i.e. one that we might think of first if we were 
      inventing statistics) of the population mean E(W), where W has
      the distribution of the popuulation

    - Similarly, since

      Var(W) = E[(W - E(W))<sup>2]</sup>,

      a natural estimator of Var(W) is

      s<sup>2</sup> = (1/n) &Sigma;<sub>i</sub>
      (W<sub>i</sub> - W<sub>bar</sub>)<sup>2</sup>

      Suppose now W is a vector. The similar estimator of Cov(W) is

      (1/n) &Sigma;<sub>i</sub> 
      (W<sub>i</sub> - W<sub>bar</sub>)
      (W<sub>i</sub> - W<sub>bar</sub>)'

      Many books will use n-1 rather than n as a divisor, to make
      s<sup>2</sup> *unbiased*, i.e.

      E(s<sup>2</sup>) = Var(W)

      This says the average of s<sup>2</sup>, over all possible samples,
      is equal to the population variance. This usually won't matter for
      us, as we are primarily interested in large-n settings.

<!-- TOC --><a name="exact-vs-approximate-statistical-inference-historical-perspective"></a>
# "Exact" vs. Approximate Statistical Inference: Historical Perspective

* Say we have scalars W<sub>i</sub> as above, with E(W<sub>i</sub>) = &mu;
  and Var(W<sub>i</sub>) = &sigma;<sup>2</sup>. The latter two are
  population values, and we wish to estimate &mu; in particular. 

* When statistical inference methodology -- confidence intervals (CIs)
  and hypothesis tests -- was being developed, they had no computers, so
  they just made the assumption that the sampled quantity, W above, has
  a normal distribution, to simplify things.

* The pioneers of statistics didn't have asymptotic theory available
  either, so the methodology they developed is called "exact"; the
  probabilities etc. are exactly correct (if the assumption of normal
  W<sub>i</sub> holds).

* The notion of the normal distribution goes back to the 18th and 19th
  centuries, especially due to the Central Limit Theorem (CLT).

* One might say that the normal family was there in Nature, waiting to
  be discovered.

* But related distribution familes, notably Student-t, chi-squared and
  the F-distribution are "man-made," in the following sense.

* With some clever math, they found that certain quantities of interest
  had the same distribution as that of the sum of k squares of N(0,1)
  variables. They thus decided to give this distribution family a name,
  *chi-squared*, thus inventing a new distribution family.  (And they
  laboriously tabulated this family's cdf by hand, just as they had for
  the normal.)

* With some more clever math, they found that a certain quantity had the
  distribution of a N(0,1) random variable divided by the square root of
  an independent chi-squared random variable. So they invented a new
  distribution family accordingly, the *Student-t distribution*. For
  similar reasons, they invented the *F-distribution* family.

* This gives us the test taught in all elementary statistics courses: To
  test the hypothesis H<sub>0</sub>: &mu; = c, form the expression
  T = (W<sub>bar</sub> - c) / (s/n<sup>0.5</sup>) (where s is as above, but
  with an n-1 divisor). This has a Student-t distribution with n-1
  degrees of freedom, whose probabilities can then be used to form
  cutoff levels for the test. 

* For n = 50, for instance, the lower and upper 2.5% cutoffs are -2.01
  and 2.01. So, we reject H<sub>0</sub> if (W<sub>bar</sub> - c) /
  (s/n<sup>0.5</sup>) is outside these bounds.

* This is called "exact" inference, since the probabilities are exact,
  if the W<sub>i</sub> are normal.

* Later other distribution families were developed, such as the
  exponential and gamma, and thus exact methods were developed for them
  as well. This led to development of a general method known as *maximum
  likelihood estimation* (MLE).

* Yet, people soon realized that in the expression T above,
  W<sub>bar</sub> is approximately normal by the Central Limit Theorem.
  Moreover, the quantity s is approximately &sigma;. So T is
  approximately N(0,1) distributed (this intuitive statement can be made
  rigorous; see below), even if W is not normal.

* In other words, one can compare T to the N(0,1) distribution,
  and thus achieve approximately correct inference, even if the sample
  distribution was not normal. Hence the practice of approximate
  inference. 

* This led to a branch of mathematical statistics known as *large sample
  theory*. Core among them is the fact that if a sequence of random
  variables Z<sub>n</sub> is asymptotically normal, then so is any
  smooth function of them g(Z<sub>n</sub>). Since the CLT says, roughly,
  that sums are approximately normal, this gives us a sequence
  X<sub>n</sub>, which when combined with a suitable function g can give
  us asymptotic normality in many useful situations.

* E.g. for MLEs. The likelihood is a product, so the log likelihood is a
  sum, just what we need for the Central Limit Theorem! It doesn't quite
  fit the CLT, is the terms are not independent, but it's enough, as
  we can use some theory to take care of the non-independence.

  We set the derivative of the log likelihood to 0, and solve for our
  parameter estimate. Here g becomes the inverse function that is needed
  to solve for the estimate, which is then approximately normal.  

  The context here is a little different -- we are still assuming an
  exact distribution, e.g.  exponential, for the sampled population, but
  the principle is the same.  Without large sample theory, in most cases
  we would not be able to derive the exact distribution of the MLE.

* The situation with linear regression models is similar, as will be
  seen below.

* The famous statistician George Box famously said, "All models are
  wrong, but some are useful." No parametric model is truly correct in
  practice. No one is 10 feet tall, and no one's weight is negative, so
  heights and weights cannot be exactly normal.  No actual random
  variable is continuous, as our measuring instruments are of only
  finite precision, so again models such as normal, exponential, gamma
  and so on are necessarily inexact.  Even the iid assumption is often
  problematic.  So the term *exact inference* is illusory.

* Approximate methods are very popular, including with me. I always use
  N(0,1) instead of the Student-t distribution, for instance.  but *many
  analysts today still favor "exact" statistical methods*. They note
  for example that there is always the question of "How large is 'large'?" 
  if one uses asymptotics.

<!-- TOC --><a name="the-multivariate-normal-distribution-family"></a>
# The Multivariate Normal Distribution Family

* Consider a random vector X = (X<sub>1</sub>,...,X<sub>p</sub>)'. The mean
  is now a vector &mu;, and standard deviation &sigma; now becomes the
  p x p covariance matrix &Sigma;.  Notation: N(&mu;,&Sigma;).

* Density: c exp[-(t-&mu;)' &Sigma;<sup>-1</sup> (t-&mu;)], where the
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
    c'w for some c. (A fairly simple formula for c is available but not
    shown here.)

  - The variance of that distribution is constant in w.

* For a symmetric matrix A and vector u, the scalar quantity u'Au is
  called a  *quadratic form*. For p-variate normal X, the
  quadratic form M = (X-&mu;)' &Sigma;<sup>-1</sup> (X-&mu;) turns
  out to have a chi-squared distribution with p degrees of freedom.

* So if &eta; is the q-th quantile of the &chi;<sup>2</sup>
  distribution with p degrees of freedom, then

  P(M < &eta;) = q

* In the above picture, horizontal slices, i.e. level sets, are
  ellipses. For p > 2, the regions are p-dimensional ellipsoids. This
  can be used to derive confidence regions in statistical applications,
  as will be seen later.

<!-- TOC --><a name="example-linear-models-exact-and-approximate-inference"></a>
# Example: Linear Models, Exact and Approximate Inference

Let's see how all this plays out with the linear model.

Let Y be a random variable and X = (X<sub>1</sub>,...,X<sub>p</sub>)' a
p-dimensional random vector. E.g.  we might wish to predict Y = college
GPA from X<sub>1</sub> = high school GPA, X<sub>2</sub> = SAT score and
so on. Let's call this the "GPA example."

<!-- TOC --><a name="classical-linear-model"></a>
## Classical linear model

The classic assumptions of the model are as follows. Let t =
  (t<sub>1</sub>,...,t<sub>p</sub>)'. 

* Linearity: E(Y | X = t) = 
  &beta;<sub>0</sub> + 
  &beta;<sub>1</sub> t<sub>1</sub> + ... +
  &beta;<sub>p</sub> t<sub>p</sub>. 
  for some unknown population vector &beta; = 
  (&beta;<sub>0</sub>, &beta;<sub>1</sub>,..., &beta;<sub>p</sub>)'.

* Normality: The distribution of Y | X = t is normal for all t.

* Homoscedasticity: Var(Y | X = t) is the same for all t, which
    we will denote by &sigma;<sup>2</sup>.

* Independence: The random vectors (X<sub>1</sub>,Y<sub>1</sub>),...,
    (X<sub>n</sub>,Y<sub>n</sub>) are independent.

The function r(t) = E(Y|X=t) is called the *regression function* of Y on
X. The term does NOT necessarily mean a linear model. So the first
bulleted item above says, "The regression function of Y on X is linear."

<!-- TOC --><a name="fixed-x-vs-random-x-models"></a>
## *Fixed-X* vs. *random-X* models.

  - In many if not most applications, the X<sub>i</sub> are random, so that
    the vectors (X<sub>1</sub>,Y<sub>1</sub>),...,
    (X<sub>h</sub>,Y<sub>h</sub>) are iid.  In our GPA example above, we
    have sampled n students at random, and their X values -- HS GPA, SAT
    -- are thus random as well.

  - But historically, the X values have usually been taken to be 
    fixed. A typical example would be, say, some industrial process in
    which a designed experiment is conducted, with pre-specified X
    values.

  - Today, the fixed-X view is still the typical one. In the GPA
    example, it would mean that we had sampled students with
    pre-specified HS GPA and SAT values, which would probably not be the
    case. However, we can still take the X values to be fixed, by
    conditioning on them.

  - In the fixed-X setting, we no longer have the second 'i' in 'iid',
    but it doesn't matter much. For instance, there are versions of the
    CLT that cover this situation.

<!-- TOC --><a name="estimation"></a>
## Estimation 

* One estmates &beta; from the data, using the formula for the estimates
  b = (b<sub>0</sub>, b<sub>1</sub>,...,b<sub>p</sub>)':

  b = (A'A)<sup>-1</sup> A'W

  where A is the n x (p+1) matrix will all 1s in column 1 and in column
  j > 1, all n values of X<sub>j-1</sub>.  Here W is the vector of all n
  values of Y.

* Write column j+1 of A as (X<sub>j1</sub>,...,X<sub>jn</sub>)'. If say
  X<sub>j</sub> is human age, then this vector consists of the ages of
  all people in our dataset. 

* The vector b, as a linear combination of W, has an exact MV normal
  distribution if Y|X is exactly normally distributed.

* By the way, since

  E(Y | X) = X<sub>extend</sub>' &beta;

  where X<sub>extend</sub> = (1,X)'

  we have that

  E(W | A) = A &beta;

  and thus 

  E(b | A) = &beta;

  (unbiasedness).

* What about Cov(b|A)? Start with the fact that, by assumption,

  Cov(Y|A) = &sigma;<sup>2</sup> I

  &sigma;<sup>2</sup> (A'A)<sup>-1</sup>. 

  Write B = (A'A)<sup>-1</sup> A, so that b = BW. Then after some
  algebra, we have

  Cov(B|A) = B &sigma;<sup>2</sup> I B' = &sigma;<sup>2</sup>
  (A'A)<sup>-1</sup>

  Note that A is a function of the X data, but here we take it to be a
  constant, by considering the conditional distribution of b given A.

  Thus confidence intervals (CIs) and regions we obtain below are also
  conditional. However, note that if, e.g. we have a conditional CI at
  the 95% level, then the unconditional coverage probability is also
  95%, by the Law of Iterated Expection above.

<!-- TOC --><a name="exact-inference"></a>
## Exact inference

* The unknown parameter &sigma;<sup>2</sup> must be estimated from the
  data, as 

  s<sup>2</sup> = [1/(n-p-1)] &Sigma;<sub>i</sub> 
  (Y<sub>i</sub> - R<sub>i</sub> b)<sup>2</sup>

  where R<sub>i</sub> is row i of A.

* So the estimated conditional covariance matrix of b is

   s<sup>2</sup> (A'A)<sup>-1</sup>

* The stringent assumptions above enable exact statistical inference.
  This enables exact confidence intervals and tests: 

  - The quantity b<sub>i</sub> - &beta;<sub>i</sub>, divided by the
    estimated standard deviation (square root of element i+1 in the
    estimated Cov(b|A)), has a Student t-distribution with n-p-1 df,
    thus setting up a CI for &beta;<sub>i</sub>. 

  - In some cases, we may be interested in something like, say,

    &beta;<sub>2</sub> - &beta;<sub>1</sub> = c'&beta;

    where c = (0,1,-1,0,0,...,0)'. Since c'b will then have a univariate
    normal distribution (recall the above material on MV normal
    distributions), we again can easily obtain a CI for the desired
    quantity. Note that we need the estimated variance of c'b, which is
    
    s<sup>2</sup> c'(A'A)<sup>-1</sup>c

<!-- TOC --><a name="approximate-inference"></a>
## Approximate inference

The quantity A'W consists of sums, so *b is asymptotically MV normally
distributed*.  Thus for large n, we can perform inference without
assuming a normal Y|X:

CIs for c'&beta;: An approximate 95% confidence interval for, e.g.
&beta;<sub>2</sub> - &beta;<sub>1</sub> is 

c'b &plusmn; 1.96 [s<sup>2</sup> c'(A'A)<sup>-1</sup>c]<sup>0.5</sup>

<!-- TOC --><a name="relaxing-the-assumptions"></a>
## Relaxing the assumptions

Just as Box pointed out that no model is perfect, no assumption (really,
part of a model) is perfect. Let's look at the normality assumption (Y|X
is normally distributed), and the *homoscedasticity* assumption, that
says Var(Y|X=t) is constant in t. 

*What can we say, assuming that linearity still holds?*

* Typically, the larger is E(Y|X=t), then the larger is Var(Y|X=t). Say
  we are predicting weight from height. Since taller people tend to have
  more (interperson) variability in weight, homoscedasticity would be
  seriously violated.

* If so, the estimator b will not be optimal, i.e. not Minimum Variance
  Unbiased.  For optimality, weighted least squares must then be used,
  with weights equal to the reciprocal of the conditional variance. That
  is, we should minimize

  &Sigma;<sub>i</sub> [Var(Y|X=X<sub>i</sub>)]<sup>-1 
  </sup>[Y<sub>i</sub> - b'X<sub>i</sub> ]<sup>2</sup>

  The word "should" applies here, as we usually don't know that
  conditional variance function.
   
* But in the sense of consistency, it doesn't matter what weights we use.
  If the linearity assumption holds, then b will still be a *consistent*
  estimator of &beta;, i.e. b &rightarrow; &beta; as n &rightarrow; &infin;.
  (See informal proof at the end of this section.)

* As noted earlier, lack of normality of Y|X alone does not invalidate
  statistical inference, due to the CLT. But assuming homoscedasticity
  when it is not approximately correct does invalidate inference, as
  the condtional covariance matrix of b is no longer &sigma;<sup>2</sup>
  (A'A)<sup>-1</sup>.

*The sandwich estimator*

Concerning that latter point, an important application of large-sample
theory is the *sandwich estimator*, which makes inference asymptotically
valid even when Var(Y|X=t) is not constant in t. It makes a simple
adjustment to the estimated Cov(b) matrix computed under the assumption
of homogeneous variance.

One of the CRAN packages that implements this adjustment is named, of
course, "sandwich." Here is an example:

``` r
z <- lm(mpg ~ .,mtcars)  # built-in R dataset
covz <- vcov(z)  # estimated covariance matrix of b
covzadj <- sandwich(z)  # adjusted version of covz
covz['carb','carb']  # prints 0.6868307
covzadj['carb','carb']  # prints 0.4209969
```

*Informal proof of consistency*

For convenience, assume Random-X regression. Let w(t) be our weight
function.

Note that for any nonnegative weight function w(t), the  quantity 

E[w(X)(Y - m(X))<sup>2</sup>] 

is minimized across all functions m by m(t) = E(Y|X=t), as can be seen
by conditioning on X and applying ordinary calculus. So if linearity
holds, 

E[w(X)((Y - v'X)<sup>2</sup>] 

achieves its minimum at v = &beta; for any weight function w, in
particular the w having constant value 1.  In other words, 

E[(Y - v'X)<sup>2</sup>] 

achieves its minimum at v = &beta;

Say we blindly do OLS (ordinary least squares, i.e.
unweighted), i.e. we choose v = b to minimize

(1/n) &Sigma;<sub>i</sub> [Y<sub>i</sub> - v'X<sub>i</sub> ]<sup>2</sup>

which is the sample analog of 

E[(Y - v'X)<sup>2</sup>] 

For each fixed v, that sample quantity will converge to the
corresponding population quantity. Then intuitively, the v that
minimizes 

(1/n) &Sigma;<sub>i</sub> [Y<sub>i</sub> - v'X<sub>i</sub> ]<sup>2</sup>

will converge to the v that minimizes

E[(Y - v'X)<sup>2</sup>] 

i.e. b will converge to &beta;. However, a formal proof of this requires
some real analysis.

<!-- TOC --><a name="multiple-inference-procedures"></a>
# Multiple Inference Procedures

Suppose we form two 95% confidence intervals from our dataset. Though
they are individually valid at the 95% level, this is not the case
jointly. In other words, in repeated sampling, the proportion of samples
in which both CIs contain the associated population parameter in the
same sample is less than 95%. What if we want to form many CIs, and
still have their joint confidence level at or above 95%?

This branch of statistics is called *multiple inference* or
*simultaneous inference*.

<!-- TOC --><a name="the-bonferroni-inequality"></a>
## The Bonferroni Inequality

One solution is to use the *Bonferroni Inequality*, which says that for
any two events A and B, 

P(A and B) &ge; 1 - P(not A) - P(not B)

This implies that if we want joint confidence level of two CIs to be at
least 95%, we can form each at the 97.5% level.

<!-- TOC --><a name="the-scheffe-method"></a>
## The Scheffe' method

But this conservative approach is not feasible for forming more than two
or three CIs.  Another approach is the *Scheffe' method*, based on
confidence ellipsoids. 

*How to form confidence ellipsoids*

  - Consider the linear model first. If all the assumptions hold,
    The quantity 

    (b-&beta;)'[s<sup>-2</sup> (A'A)] (b-&beta;)

    has an F-distribution with (p+1,n-p-1) df. This is the
    finite-sample version of the quadratic form discussed earlier.
    This sets up a confidence ellipsoid for &beta;:

    To form a (1-&alpha;) 100% confidence ellipsoid for &beta;, let q be
    the upper-&alpha; quantile of the F distribution with (p+1,n-p-1)
    df.  Then the confidence ellipsoid is the set of all t such that 

    (b-t)'[s<sup>-2</sup> (A'A)] (b-t) &le; q

  - For large n, and now not assuming sampling from a normal population,
    we form an approximate (1-&alpha;) 100% confidence ellipsoid for
    &beta; as follows: Let q be the upper-&alpha; quantile of the
    &chi;<sup>2</sup> distribution with p+1 df. Then the confidence
    ellipsoid is the set of all t such that 

    (b-t)'[s<sup>-2</sup> (A'A)] (b-t) &le; q

  - More generally, consider any asymptotically normal estimator
    &theta;<sub>est</sub> of a p-dimensional population parametric &theta;.
    Let Q denote its (asymptotic) covariance matrix. The confidence
    ellipsoid for the true &theta; is the set of all t such that

    (b-t)' Q<sup>-1</sup> (b-t) &le; q

*Forming simultaneous confidence intervals*

Now let's see how we can go from confidence *ellipsoid* to confidence
*intervals* for scalar quantities, and in such a way that the intervals
have a given simultaneous confidence level.

  - Again, consider any asymptotically normal estimator
    &theta;<sub>est</sub> of a p-dimensional population parametric &theta;.
    Let Q denote its (asymptotic) covariance matrix, and form the
    confidence ellipsoid as above. 

  - Consider a linear combination w'&theta;. Set b and c such that the
    lines (p=2), planes (p=3) or hyperplanes (p > 3)

    w'&theta; = b<sub>w</sub>

    and

    w'&theta; = c<sub>w</sub>

    are tangent to the ellipsoid from above and below. Say c<sub>w</sub> is the
    smaller of b<sub>w</sub> and c<sub>w</sub>.

  - Let's visualize this, using the CRAN package **ellipse**. Its
    function **ellipse.lm** draws a confidence ellipse for the output of
    a linear model.

    ``` r
    fit <- lm(mpg ~ qsec + wt , mtcars)
    plot(ellipse(fit, which = c('qsec', 'wt'), level = 0.90), type = 'l') 
    # draw a pair of parallel tangent lines, w = (2.5,-1)'
    abline(-9.05,2.5)
    abline(-5.70,2.5)  
    ```

    This is an example from the package, using a 90% confidence level.

  ![confidence ellipse](Ellipse.png){width=60%}

  - With 90% probability, &theta; is somewhere inside the ellipsoid.
    (Note that it is the ellipsoid that is random, not &theta;.)

  - And, if &theta; is in the ellipsoid, then w'&theta; will be in
    (c<sub>w</sub>,b<sub>w</sub>), though not *only* if; it is clear
    that there are some values of &theta; that are outside the ellipse
    but still between the two lines.  In other words,
    (c<sub>w</sub>,b<sub>w</sub>) is a confidence interval for w'&theta;
    of level at least 90%.

  - For each possible w, we get the above pair of tangents.
    Geometrically it is clear that &theta; will be inside the ellipsoid
    if and only if w'&theta; is in (c<sub>w</sub>,b<sub>w</sub>) 
    *for all w*.

  - Thus the probability that all the intervals
    (c<sub>w</sub>,b<sub>w</sub>) hold simultaneously is 90%.

<!-- TOC --><a name="simulation-of-mv-normal-random-vectors"></a>
# Simulation of MV Normal Random Vectors

Suppose we wish the generate random vectors V having a p-variate normal
distribution with mean &nu; and covariance matrix &Gamma;. How can we do
this?

We could find a square root matrix for &Gamma;, G. Then we generate
vectors X consisting of p independent N(0,1) random variables, and set

V = G X + &nu;

Since X has covariance matrix equal to the identity matrix I, V will
have the specified mean and covariance matrix, and it will be MV normal
by linearity.

R is famous for its CRAN package collection, so it is not surprising
that this operation has already been coded, first in an early package,
MASS and then in others. The function **MASS::mvrnorm** does this.

What about regression contexts, say linear models with the standard
assumptions? The model

E(Y | X = t) =
  &beta;<sub>0</sub> + 
  &beta;<sub>1</sub> t<sub>1</sub> + ... +
  &beta;<sub>p</sub> t<sub>p</sub>,

Var(Y|X=t) = &sigma;<sup>2</sup>

can be (and in most books, usually is) written as 

Y = &beta;<sub>0</sub> +
  &beta;<sub>1</sub> t<sub>1</sub> + ... +
  &beta;<sub>p</sub> t<sub>p</sub> + &epsilon;

where &epsilon; has a N(0,&sigma;<sup>2</sup>) distribution and is
independent of X. (In the fixed-X context, the latter assumption is
automatic.)

So, one can use **mvrnorm** to generate X, then **rnorm** to generate
&epsilon;, thus obtaining Y.

<!-- TOC --><a name="types-of-convergence"></a>
# Types of Convergence

* A sequence of random variables V<sub>n</sub> is said to converge 
  *in probability* to a random variable V if for every &epsilon; > 0,

  lim<sub>n &rarr; &infin;</sub> P(|V<sub>n</sub> - V| > &epsilon;) = 0
 
* For random vectors V<sub>i</sub>, replace | | by, e.g. Euclidean
  distance.

* Say we have random variables Q<sub>n</sub>. If
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

  That definition is hard to deal with, but the Cramer-Wold device often
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

<!-- TOC --><a name="central-limit-theorems-clts"></a>
# Central Limit Theorems (CLTs)

* Univariate Central Limit Theorem: Let X<sub>i</sub>, i = 1,2,... be
  iid random variables with mean &mu; and variance &sigma;<sup>2</sup>.
  Define T<sub>n</sub> = X<sub>1</sub>+...+X<sub>n</sub>, which has mean
  and variance n &mu; and n &sigma;<sup>2</sup>. Then 

    1/n<sup>1/2</sup> &sigma;<sup>-1</sup> (T<sub>n</sub> - n &mu;)

  converges in distribution to N(0,1).

  Many generalizations of the CLT exist, in which they relax the iid
  assmuption, especially the second 'i'.

* Scaling by the factor n<sup>1/2</sup> is "just right" to make things
  work here. If, say we were to divide instead by n, the quantity would
  go to 0, by the SLLN. That certainly would not be helpful in terms of
  finding probabilities, confidence intervals and so on.

  Other powers may apply for other quantities. For
  example, if we are interested in the maximum of
  X<sub>1</sub>,...,X<sub>n</sub> rather than their sum, then 
  we divide by (2 log n)<sup>1/2</sup> (and the limiting 
  distribution will be something other than normal).

* We will often encounter a sequence of random variables Q<sub>n</sub>
  that is asymptotically normal with mean &nu; and variance &gamma;.
  Note carefully that that does NOT mean that E(Q<sub>n</sub>) &approx;
  &nu; with a similar statement for variance. On the contrary,
  E(Q<sub>n</sub>) could be infinite or undefined. Instead, it merely
  means that cdf of 

  n<sup>0.5</sup>(Q<sub>n</sub> - &nu;)/&gamma;

  converges (pointwise) to the cdf of N(0,1).

* Multivariate Central Limit Theorem: Let X<sub>i</sub>, i = 1,2,... be
  iid random vectors with mean vector &mu; and covariance matrix &Sigma;.
  Define T<sub>n</sub> = X<sub>1</sub>+...+X<sub>n</sub>, which has mean
  and covariance matrix n &mu; and n &Sigma;. Then 

    1/n<sup>1/2</sup> (T<sub>n</sub> - n &mu;)

  converges in distribution to  N(0,&Sigma;). This can easily be proved using
  the Cramer-Wold device.

<!-- TOC --><a name="example-method-of-moments-mm-estimators"></a>
# Example: Method of Moments (MM) Estimators

This is a serious application of the above methodology, and will take
time to prepare and apply the concepts. The reader's patience is
requested.

*Background on MM estimators*

* Suppose we are modeling a random variable with some parametric
    distribution family, such as normal, exponential or gamma. We want to
    estimate the parameter(s) from our data
    X<sub>1</sub>,...,X<sub>n</sub>. How can we do this?

* The most famous general method is *Maximum Likelihood Estimation*
    (MLE), but it's somewhat less-famous cousin, the *Method of Moments* 
    (MM), is sometimes the handier one. (It is also the easier one to
    explain.)

* This is best introduced by example. Say we are modeling our
    X<sub>i</sub> as having an exponential distribution, i.e. they have
    density

    &tau; exp(-&tau;t).

    for some unknown parameter &tau;.  How do we construct the MM
    estimate T of &tau;?

    The mean of this distribution is 1/&tau;, which could be estimated
    by 1/T once we obtain T. On the other hand, the classical estimate
    of a population mean is the sample mean,

    X<sub>bar</sub> = 
    (1/n) ( X<sub>1</sub>+...+X<sub>n</sub>) 

    So set our estimated mean under the exponential assumption, 1/T, to
    our general estimated mean, without that assmuption:

    1/T = X<sub>bar</sub>

    That gives us our MM estimate,

    T = 1/X<sub>bar</sub>

    Done!

* How does MM work in general, when we have more than one parameter?
    Let's use k to denote the number of parameters in the given parametric
    family. For instance, k = 2 for the univariate normal family,
    corresponding to the two parameters &mu; and &sigma;<sup>2</sup>. We
    then must bring in E(X<sup>2</sup>) and so on, as follows.

* The r-th *moment* of a random variable X is defined to be
    m<sub>r</sub> = E(X<sup>r</sup>), r = 1,2,....  We will need to use
    the first k moments.

    Note that this is a population value;
    its intuitive sample estimate is 

    M<sub>r</sub> = 
    (1/n) (X<sub>1</sub><sup>r</sup>+...+X<sub>n</sub><sup>r</sup>)

* Alternatively, for k > 1 one may use the r-th *central* moment, 

    E[(X-E(X))<sup>r</sup>]

    Since for r = 2 this is Var(X) and for r = 3 we get the skewness of the
    distribution, things may be more convenient this way.

    It will be clear from context whether we mean the central or noncentral
    moment. 

* Use &tau; to denote the population value of the parameter, and T to 
    denote its sample estimate under MM.  Both are vectors if k > 1.

* Before continuing, let's review what happened when we set 1/T =
  X<sub>bar</sub> above..  The left-hand side is our estimate of
  m<sub>1</sub> under the exponential model, while the right-hand side
  is a general estimator of m<sub>1</sub>, not assuming that model This
  is how the Method of Moments works, by matching the two. We  will have
  k such equations, then solve them for T.

* For an example with k = 2, consider the *gamma distribution family*:

  c &tau;<sub>1</sub><sup>&tau;<sub>2</sub></sup>
  t<sup>&tau;<sub>2 </sub> - 1</sup></sup>
  exp(-&tau;<sub>1</sub>t)

  where c is a constant to make the density integrate to 1.0. This is
  a common model used in network communications and medical survival
  analysis for example.  

  Here is what the gamma family of densities looks like. (From
  *Probability and Statistics for Data Science: Math + R + Data*, N.
  Matloff, 2019. Here r and &lambda; refer to our &tau;<sub>2</sub>
  and &tau;<sub>1</sub>.)  

  ![Some densities from the gamma family](Gamma.png)

  The model is handy when one has a distribution on (0,&infin;),
  assumed unimodal.

* The population moments are

  E(X) = &tau;<sub>2</sub> / &tau;<sub>1</sub>

  and

  Var(X) = &tau;<sub>2</sub> / &tau;<sub>1</sub><sup>2</sup>

* So, we equate:

  M<sub>1</sub> = T<sub>2</sub> / T<sub>1</sub>

  and (using the central moment for M<sub>2</sub>)

  M<sub>2</sub> = T<sub>2</sub> / T<sub>1</sub><sup>2</sup>

* These are nonlinear equations, but we are lucky here; they are
  solvable from simple algebra (otherwise we would need to use
  numerical approximation methods):

  T<sub>1</sub> = M<sub>1</sub> / M<sub>2</sub>

  and

  T<sub>2</sub> = M<sub>1</sub> T<sub>1</sub> =
  M<sub>1</sub><sup>2</sup> / M<sub>2</sub>

*Consistency*

* Say we have an estimator &theta;<sub>n</sub> for a parameter &theta;
  based on a sample of size n. As n goes to &infin;, we would like our
  estimator to have the property that &theta;<sub>n</sub> goes to
  &theta;. If it does so almost surely, we say it is *strongly
  consistent*; if the convergence is just in probability, it is known
  as *weak consistency*.

  Since the SLLN implies that M<sub>i</sub> is strongly consistent for
  m<sub>i</sub>, this implies the same for the T<sub>i</sub>, as long
  as the moments are continuous functions of the M<sub>i</sub>.  In
  the example above, for instance, M<sub>1</sub> and M<sub>2</sub>
  converge almost surely, hence the T<sub>i</sub> do too.

*Asymptotic normality*

* OK, we now have estimators for &tau;<sub>1</sub> and
  &tau;<sub>2</sub>, and they are consistent. But the latter is a very
  weak property. For good practical value, it would be desirable to
  obtain confidence intervals for the &tau;<sub>i</sub> from the
  T<sub>k</sub>. For this we need asymptotic normality, as follows.

* For simplicity, let's consider the case k = 1, so the density family
  has a single scalar parameter, &tau;, as we saw in the exponential
  example above. Then E(X) is some function g(&tau;). In the
  exponential example, g(&tau;) = 1/&tau;.

  Assume g is differentiable with a continuous derivative g<sub>1</sub>. Then
  Taylor's Theorem from calculus says

  g(T) = g(&tau;) + g<sub>1</sub>(T<sub>betw</sub>) (T - &tau;)

  for some T<sub>betw</sub> between &tau; and T. 

* Now to set up using the CLT, recall that we will set

  g(T) = M<sub>1</sub>

  Then write 

  n<sup>0.5</sup>(M<sub>1</sub> - m<sub>1</sub>) =
  n<sup>0.5</sup> [g(T) - g(&tau;)] =
  n<sup>0.5</sup> [g<sub>1</sub>(T<sub>betw</sub>) (T - &tau;)]

  and then

  n<sup>0.5</sup>(T - &tau;) =
  n<sup>0.5</sup> (M<sub>1</sub> - m<sub>1</sub>)/ 
  g<sub>1</sub>(T<sub>betw</sub>)

* We know that T is a strongly consistent estimator of &tau;, so by
  the continuity of g<sub>1</sub>, the denominator in the RHS
  converges almost surely to g<sub>1</sub>(&tau;). Applying the
  Slutsky Theorem and the CLT (M<sub>1</sub> is a sum of iid terms),
  we see that the RHS is asymptotically normal (though not with
  standard deviation 1).

* In other words

  n<sup>0.5</sup>(T - &tau;) =<sub>asympt.</sub>
  n<sup>0.5</sup> (M<sub>1</sub> - m<sub>1</sub>)/ 
  g<sub>1</sub>(&tau;)

  Here =<sub>asymp.</sub> should not be viewed as "asymptotically
  equal to," but rather "having the same asymptotic distribution as."  

  The variance in that asymptotically normal distribution will be

  Var(M<sub>1</sub>)/[g<sub>1</sub>(&tau;)]<sup>2</sup>

* To compute a CI, we replace the quantities by their sample analogs.
  E.g. our estimate of Var(M<sub>1</sub>) is

  s<sup>2</sup>= (1/n) &Sigma;<sub>i=1</sub><sup>n</sup> 
  [X<sub>i</sub> - M<sub>1</sub>]<sup>2</sup>

  (Divide by n-1 instead of n if you prefer, though there really is no
  reason to do so.)

  Our estimate of g<sub>1</sub>(&tau;) is 

  g<sub>1</sub>(T) = -1/T<sup>2</sup>

  Our CI, say for a 95% confidence level, is then 

  T &plusmn; 1.96 s/T<sup>2</sup>

* For k = 2, we now have functions g and h, both having 
  arguments T<sub>1</sub> and T<sub>2</sub>. We now have 
  derivatives g<sub>&Del;</sub> and h<sub>&Del;</sub>, but 
  they are now gradients.  

  n<sup>0.5</sup>(M1 - m1) =<sub>asymp..</sub>
  n<sup>0.5</sup>
  g<sub>&Del;</sub>' (T - &tau;)

  n<sup>0.5</sup>(M2 - m2) =<sub>asymp.</sub>
  n<sup>0.5</sup>
  h<sub>&Del;</sub>' (T - &tau;)

  Write M = (M<sub>1</sub>,M<sub>2</sub>)' and similarly for m. Then

  n<sup>0.5</sup> (M - m) =<sub>asymp.</sub>
  n<sup>0.5</sup>
  A (T - &tau;)

  where the matrix A consists of g<sub>&Del;</sub>' in row 1 and
  h<sub>&Del;</sub>' in row 2.

* So, finally, we have that the asymptotic distribution of

  n<sup>0.5</sup> (T - &tau;) is MV normal with covariance matrix
 
  A<sup>-1</sup> Cov(M) A<sup>-1</sup>'
 
  This enables confidence intervals and ellipsoids as before.

<!-- TOC --><a name="the-delta-method"></a>
# The Delta Method

* Earlier we said, "if a sequence of random variables Z<sub>n</sub>
  including the case of random vectors), is asymptotically normal, then
  so is any smooth function of them g(Z<sub>n</sub>)." This fact is
  quite useful.

* Before we go on, note again that "asymptotically normal" means in the
  sense of the CLT. So if we say "g(Z<sub>n</sub>) is asymptotically
  normal," we mean that its cdf is approximately equal to a normal cdf,
  and we must scale by n<sup>0.5</sup>.

* The proof is similar to our derivation for MM estimators above. E.g.
  for dimension 2, set &mu;= (&mu;<sub>1</sub>,&mu;<sub>2</sub>) to
  denote the mean of the asymptotic normal distribution, with &Sigma;
  representing the corresponding covariance matrix.

  Now write

  g(Z<sub>n</sub>) - &mu; &approx; 
  g<sub>1</sub>(&mu;)(</sub>Z<sub>1n</sub> - &mu;<sub>1</sub>) +
  g<sub>2</sub>(&mu;)(</sub>Z<sub>2n</sub> - &mu;<sub>2</sub>) =
  A (Z<sub>n</sub> - &mu;)

  where the 1 x 2 matrix A = (g<sub>1</sub>(&mu;), g<sub>2</sub>(&mu;)).

  The right-hand side is now a linear form in an asymptotically MV
normally distributed vector, thus asymptotically MV (actually
univariate) normal, with mean 0 and variance A &Sigma; A'.
