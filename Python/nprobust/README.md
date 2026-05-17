# Kernel Density and Local Polynomial Regression Methods

The package `nprobust` implements estimation, inference, bandwidth selection,
and graphical procedures for kernel density and local polynomial regression
methods, including robust bias-corrected confidence intervals.

- `lprobust`: local polynomial point estimation and robust bias-corrected inference.
- `lpbwselect`: data-driven bandwidth selection for local polynomial regression.
- `kdrobust`: kernel density point estimation and robust bias-corrected inference.
- `kdbwselect`: data-driven bandwidth selection for kernel density estimation.

See references for methodological and practical details.

Website: [https://nppackages.github.io/](https://nppackages.github.io/).

Source code: [https://github.com/nppackages/nprobust](https://github.com/nppackages/nprobust).

## Authors

Sebastian Calonico (<scalonico@ucdavis.edu>)

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Max H. Farrell (<mhfarrell@gmail.com>)


## Installation

To install/update use pip:
```
pip install nprobust_pkg
```

## Usage

```python
from pathlib import Path

import pandas as pd
from nprobust import kdrobust, kdbwselect, lprobust, lpbwselect, plot_lprobust

# Cholesterol trial data used by the R and Stata examples.
data = pd.read_csv(Path("..") / "nprobust_data.csv")
control = data["t"] == 0

# Local polynomial regression with robust bias-corrected confidence intervals.
result = lprobust(data.loc[control, "cholf"], data.loc[control, "chol1"])
print(result.summary())

# Data-driven bandwidth selection.
bw = lpbwselect(data.loc[control, "cholf"], data.loc[control, "chol1"],
                bwselect="mse-dpi", neval=7)
print(bw.bws)

# Kernel density estimation.
density = kdrobust(data.loc[control, "chol1"], neval=30)
print(density.summary())

# Kernel density bandwidth selection.
print(kdbwselect(data.loc[control, "chol1"], bwselect="imse-dpi").bws)

# Plot a local polynomial fit.
fig = plot_lprobust(result, xlabel="chol1", ylabel="cholf")
fig.savefig("fit.png")
```

- Replication: [nprobust illustration](https://github.com/nppackages/nprobust/blob/main/Python/nprobust_illustration.py), [nprobust data](https://github.com/nppackages/nprobust/blob/main/Python/nprobust_data.csv).


## Dependencies

- numpy
- pandas
- scipy
- matplotlib (optional plotting extra)

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


<br><br>
