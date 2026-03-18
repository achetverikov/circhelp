# Remove cardinal biases for data with orientation (color, motion, ...) set in discrete steps

Remove cardinal biases for data with orientation (color, motion, ...)
set in discrete steps

## Usage

``` r
remove_cardinal_biases_discrete(err, x, space, init_outliers = NULL)
```

## Arguments

- err:

  a vector of errors, deviations of response from the true stimuli

- x:

  a vector of true stimuli in degrees (see space)

- space:

  circular space to use (a string: `180` or `360`)

- init_outliers:

  a vector determining which errors are initially assumed to be outliers
  (default: NULL)

## Value

returns a data.table with the following columns:

- is_outlier - 0 for outliers (defined as ±3\*predicted SD, where SD and
  mean are computed after excluding initially assumed outliers)

- be_c error corrected for biases (`be_c = observed error - pred`)
