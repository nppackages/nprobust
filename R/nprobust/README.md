# Kernel Density and Local Polynomial Regression Methods

The package `nprobust` implements estimation, inference, bandwidth selection,
and graphical procedures for kernel density and local polynomial regression
methods, including robust bias-corrected confidence intervals.

- `lprobust()`: local polynomial point estimation and robust bias-corrected inference.
- `lpbwselect()`: data-driven bandwidth selection for local polynomial regression.
- `kdrobust()`: kernel density point estimation and robust bias-corrected inference.
- `kdbwselect()`: data-driven bandwidth selection for kernel density estimation.
- `nprobust.plot()`: graphical presentation of `lprobust()` and `kdrobust()` results.

See references for methodological and practical details.

Website: [https://nppackages.github.io/](https://nppackages.github.io/).

Source code: [https://github.com/nppackages/nprobust](https://github.com/nppackages/nprobust).

## Authors

Sebastian Calonico (<scalonico@ucdavis.edu>)

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Max H. Farrell (<mhfarrell@gmail.com>)

## Installation

To install/update use R:

```r
install.packages("nprobust")
```

## Usage

```r
library(nprobust)

# Cholesterol trial data used by the Python and Stata examples.
data <- read.csv("../nprobust_data.csv")
control <- data$t == 0

# Local polynomial regression with robust bias-corrected confidence intervals.
result <- lprobust(data$cholf[control], data$chol1[control])
summary(result)

# Data-driven bandwidth selection.
bw <- lpbwselect(data$cholf[control], data$chol1[control],
                 bwselect = "mse-dpi", neval = 7)
summary(bw)

# Kernel density estimation.
density <- kdrobust(data$chol1[control], neval = 30)
summary(density)

# Kernel density bandwidth selection.
summary(kdbwselect(data$chol1[control], bwselect = "imse-dpi"))

# Plot a local polynomial fit.
nprobust.plot(result, xlabel = "chol1", ylabel = "cholf")
```

- Replication: [nprobust illustration](https://github.com/nppackages/nprobust/blob/main/R/nprobust_illustration.R), [nprobust data](https://github.com/nppackages/nprobust/blob/main/R/nprobust_data.csv).

## Dependencies

- ggplot2

## References

For overviews and introductions, see [nppackages website](https://nppackages.github.io).

### Software and Implementation

- Calonico, Cattaneo and Farrell (2019): [nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf).<br>
_Journal of Statistical Software_ 91(8): 1-33.

### Technical and Methodological

- Calonico, Cattaneo and Farrell (2018): [On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf).<br>
_Journal of the American Statistical Association_ 113(522): 767-779.

- Calonico, Cattaneo and Farrell (2022): [Coverage Error Optimal Confidence Intervals for Local Polynomial Regression](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf).<br>
_Bernoulli_ 28(4): 2998-3022.
