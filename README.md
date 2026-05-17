# Kernel Density and Local Polynomial Regression Methods

The package `nprobust` implements estimation, inference, bandwidth selection,
and graphical procedures for kernel density and local polynomial regression
methods, including robust bias-corrected confidence intervals.

- `lprobust`: local polynomial point estimation and robust bias-corrected inference.
- `lpbwselect`: data-driven bandwidth selection for local polynomial regression.
- `kdrobust`: kernel density point estimation and robust bias-corrected inference.
- `kdbwselect`: data-driven bandwidth selection for kernel density estimation.
- `plot_lprobust`/`plot_kdrobust` and `nprobust.plot`: graphical presentation of results.


## Python Implementation

To install/update in Python type:
```bash
pip install nprobust_pkg
```

- Help: [PyPI repository](https://pypi.org/project/nprobust-pkg/).

- Replication: [nprobust illustration](Python/nprobust_illustration.py), [nprobust data](Python/nprobust_data.csv).

## R Implementation

To install/update in R type:
```r
install.packages('nprobust')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/nprobust/nprobust.pdf), [CRAN repository](https://cran.r-project.org/package=nprobust).

- Examples/data: [nprobust illustration](R/nprobust_illustration.R), [nprobust data](R/nprobust_data.csv).

## Stata Implementation

To install/update in Stata type:
```
net install nprobust, from(https://raw.githubusercontent.com/nppackages/nprobust/main/stata) replace
```

- Help: [kdrobust](stata/kdrobust.pdf), [kdbwselect](stata/kdbwselect.pdf), [lprobust](stata/lprobust.pdf), [lpbwselect](stata/lpbwselect.pdf).

- Replication: [nprobust illustration](stata/nprobust_illustration.do), [nprobust data](stata/nprobust_data.dta).


## References

For overviews and introductions, see the [nppackages website](https://nppackages.github.io).

### Software and Implementation

- Calonico, Cattaneo and Farrell (2019): [nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf). _Journal of Statistical Software_ 91(8): 1-33.

### Technical and Methodological

- Calonico, Cattaneo and Farrell (2018): [On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf). _Journal of the American Statistical Association_ 113(522): 767-779. [Supplemental Appendix](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA--Supplement.pdf).
- Calonico, Cattaneo and Farrell (2022): [Coverage Error Optimal Confidence Intervals for Local Polynomial Regression](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf). _Bernoulli_ 28(4): 2998-3022. [Supplemental Appendix](https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli--Supplement.pdf).

## Funding

This work was supported by the National Science Foundation through grants [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931) and [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805).
