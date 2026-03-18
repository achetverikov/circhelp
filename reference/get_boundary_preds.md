# Get polynomial predictions for values at the boundaries

A helper function for
[`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md).

## Usage

``` r
get_boundary_preds(
  group,
  data,
  space,
  reassign_range,
  gam_ctrl,
  poly_deg,
  angle_diff_fun,
  weights = NULL
)
```

## Arguments

- group:

  group (bin) id

- data:

  dataset

- space:

  see
  [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md)

- reassign_range:

  see
  [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md)

- gam_ctrl:

  control object for gam models

- poly_deg:

  see
  [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md)

- angle_diff_fun:

  a function to compute difference between angles

## Value

a data.table with predicted values
