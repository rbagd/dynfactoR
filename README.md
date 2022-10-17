
This repository is deprecated in favor of [dfms](https://github.com/SebKrantz/dfms), a vastly improved version of whatever you can find here.

---

This is a public repository for `dynfactoR`, a package for `R` which facilitates estimation of dynamic factor models.

Current implementation of main `dfm` function supports vector auto-regressive type dynamics for factors, missing observations and some statistical identification restrictions. A dynamic factor model estimation will typically return 3 estimates, namely principal component estimator, a two-step estimator as well as quasi-maximum likelihood (QML) estimator. The two latter estimators are based on Kalman filtering and QML estimator is a particular case of EM-algorithm.

`dynfactoR` is easy to install with the help of `devtools`:
```
devtools::install_github("rbagd/dynfactoR")
```

`dynfactoR` also provides a dataset from National Bank of Belgium which contains monthly Belgian business and consumer survey data over a sufficiently long period. You can call it with:

```
data("NBBsurvey")
```

If you are not familiar with related literature, it could be useful to read a short introduction to dynamic factor models which is packaged as a vignette. Some academic references are also made available.

```
vignette("dynamic-factors")
```

In the future, the package will be extended to support more general dynamics as well as Markov-switching factor loading or transition matrices. Currently, these features are early experimental.
