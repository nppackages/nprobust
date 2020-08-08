# NPROBUST

The nprobust package provides Stata and R implementations of bandwidth selection, point estimation and inference procedures for nonparametric kernel-based density and local polynomial methods.

This work was supported by the National Science Foundation through grant [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931) and [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805).

## Stata Implementation

To install/update in Stata type:
```
net install nprobust, from(https://sites.google.com/site/nppackages/nprobust/stata) replace
```

- Help: [kdrobust](stata/kdrobust.pdf), [kdbwselect](stata/kdbwselect.pdf), [lprobust](stata/lprobust.pdf), [lpbwselect](stata/lpbwselect.pdf).

- Replication: [do-file](references/nprobust_illustration.do), [nprobust_data](references/nprobust_data.do).

## R Implementation
To install/update in R type:
```
install.packages('nprobust')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/nprobust/nprobust.pdf), [CRAN repository](https://cran.r-project.org/package=nprobust).

- Replication: [R-script](R/nprobust_illustration.r), [nprobust_data](R/nprobust_data.csv).

## References

### Software and Implementation

- Calonico, Cattaneo and Farrell (2019): [nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference](references/Calonico-Cattaneo-Farrell_2019_JSS.pdf), Journal of Statistical Software 91(8): 1-33.

### Technical and Methodological

- Calonico, Cattaneo and Farrell (2018): [On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference](references/Calonico-Cattaneo-Farrell_2019_JASA.pdf), Journal of the American Statistical Association 113(522): 767-779. [Supplemental Appendix]

- Calonico, Cattaneo and Farrell (2020): [Coverage Error Optimal Confidence Intervals for Local Polynomial Regression](references/Calonico-Cattaneo-Farrell_2020_wp.pdf). [Supplemental Appendix]

<br><br>
