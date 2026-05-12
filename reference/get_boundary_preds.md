# Get polynomial predictions for values at the boundaries

A helper function for
[`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md).

## Usage

``` r
get_boundary_preds(
  group_i,
  group_label,
  data,
  dist_to_centers_mat,
  space,
  reassign_range,
  poly_deg,
  angle_diff_fun,
  weights = NULL,
  extract_only = FALSE,
  boundary_data = NULL
)
```

## Arguments

- group_i:

  integer group (bin) id

- group_label:

  group (bin) label

- data:

  dataset

- dist_to_centers_mat:

  precomputed distances to bin centers

- space:

  see
  [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md)

- reassign_range:

  see
  [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md)

- poly_deg:

  see
  [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md)

- angle_diff_fun:

  a function to compute difference between angles

## Value

a data.table with predicted values
