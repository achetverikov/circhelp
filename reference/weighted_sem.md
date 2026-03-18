# Weighted standard error of the mean (SEM_w)

Computes the variance of a weighted mean following the definitions given
by Kirchner (2006).

## Usage

``` r
weighted_sem(x, w, na.rm = FALSE)
```

## Arguments

- x:

  variable to compute the SEM for

- w:

  weights

- na.rm:

  should NAs be removed

## Value

weighted standard error of the mean

## Details

James Kirchner describes two different cases when the weighted variance
is computed. The code here implements Case I where "one wants to give
more weight to some points than to others, because they are considered
to be more important" and "the weights differ but the uncertainties
associated with the individual xi are assumed to be the same" (Kirchner,
2006, p. 1). The formula used is: \\SEM_w = \sqrt{\left(\sum\_{i =
1}^{N} (w\_{i} x_i^2)-\bar{x}^2\right)\frac{\sum\_{i = 1}^{N}
w_i^2}{1-\sum\_{i = 1}^{N} w_i^2}} \\ The expected error is within 5% of
the bootstrapped SEM (at larger sample sizes).

## References

- Kirchner, J. 2006. Data Analysis Toolkit \#12: Weighted averages and
  their uncertainties.
  <https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf>.
  Retrieved on 04.07.2024.

- Bevington, P. R. 1969. Data Reduction and Error Analysis for the
  Physical Sciences. McGraw-Hill, 336 pp.

## Examples

``` r
set.seed(1)
n_obs <- 200
w <- runif(n_obs)
w <- w / sum(w)
x <- rnorm(n_obs, sd = 5)
weighted_sem(x, w)
#> [1] 0.3812216
```
