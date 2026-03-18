# Conversion between the circular SD and kappa of von Mises

Conversion between the circular SD and kappa of von Mises

## Usage

``` r
vm_kappa_to_circ_sd(kappa)

vm_kappa_to_circ_sd_deg(kappa)

vm_circ_sd_to_kappa(sd)

vm_circ_sd_deg_to_kappa(sd_deg)
```

## Arguments

- kappa:

  von Mises kappa parameter

- sd:

  circular SD of von Mises (radians)

- sd_deg:

  circular SD of von Mises (degrees)

## Value

`vm_kappa_to_circ_sd` and `vm_kappa_to_circ_sd_deg` return circular SD
(in radians or degrees, respectively) corresponding to a given kappa.
`vm_circ_sd_to_kappa` and `vm_circ_sd_deg_to_kappa` return kappa
corresponding to a given circular SD (in radians or degrees,
respectively).

## Functions

- `vm_kappa_to_circ_sd_deg()`: get circular SD (in degrees) from kappa

- `vm_circ_sd_to_kappa()`: get kappa from circular SD (in radians)

- `vm_circ_sd_deg_to_kappa()`: get kappa from circular SD (in degrees)

## Examples

``` r
vm_kappa <- 5
vm_sd <- vm_kappa_to_circ_sd(vm_kappa)

vm_circ_sd_to_kappa(vm_sd)
#> [1] 5

x <- circular::rvonmises(10000, mu = circular::circular(0), kappa = vm_kappa)

sprintf("Expected SD: %.2f, actual SD: %.2f", vm_sd, circ_sd_rad(x))
#> [1] "Expected SD: 0.47, actual SD: 0.48"
```
