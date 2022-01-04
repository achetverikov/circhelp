# circhelp

A small helper package for circular data analyses in R, particularly useful for cognitive studies on orientation, motion direction, and other circular features. Contains functions for descriptive statistics for circular data (computing means, SD, and skewness), angular differences, and correlation. Also includes a function to correct for cardinal biases in the human estimates of circular features (e.g., orientation).

## Installation

    # install.packages("devtools")
    devtools::install_github('achetverikov/circhelp')

## Usage

Most of the functions are self-explanatory. The only (somewhat) complicated function is `remove_cardinal_biases`, see the vignette [cardinal biases](https://achetverikov.github.io/circhelp/articles/cardinal_biases.html) for an example use case.

<img src="https://achetverikov.github.io/circhelp/articles/cardinal_biases_files/figure-html/correct_biases_ex-1.png" title="Example of cardinal biases processing" width="100%"/>
