# BTYDplus

[![Travis-CI Build Status](https://travis-ci.org/mplatzer/BTYDplus.svg?branch=master)](https://travis-ci.org/mplatzer/BTYDplus)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BTYDplus)](https://CRAN.R-project.org/package=BTYDplus)
![CRAN Downloads month](https://cranlogs.r-pkg.org/badges/BTYDplus)
![CRAN Downloads overall](https://cranlogs.r-pkg.org/badges/grand-total/BTYDplus)
<!-- [![Coverage Status](https://img.shields.io/codecov/c/github/mplatzer/BTYDplus/master.svg)](https://codecov.io/github/mplatzer/BTYDplus?branch=master) -->

The BTYDplus [R](https://www.r-project.org/) package provides advanced statistical methods to describe and predict customer's purchase behavior. It uses historic transaction records to fit a probabilistic model, which then allows to compute quantities of managerial interest on a cohort- as well as on a customer level (Customer Lifetime Value, Customer Equity, P(alive), etc.).

This package complements the [BTYD](https://cran.r-project.org/package=BTYD) package by providing several additional buy-till-you-die models, that have been published in the marketing literature, but whose implementation are complex and non-trivial. These models are: NBD, MBG/NBD, BG/CNBD-k, MBG/CNBD-k, Pareto/NBD (HB), Pareto/NBD (Abe) and Pareto/GGG.

## Installation

```
# install.packages("devtools")
devtools::install_github("mplatzer/BTYDplus", dependencies=TRUE)
library(BTYDplus)
```

## Getting Started

```
demo("cdnow")        # Demonstration of fitting various models to the CDNow dataset
demo("mbg-cnbd-k")   # Demonstration of MBG/CNBD-k model with grocery dataset
demo("pareto-abe")   # Demonstration of Abe's Pareto/NBD variant with CDNow dataset
demo("pareto-ggg")   # Demonstration of Pareto/NBD (HB) & Pareto/GGG model with grocery dataset
```

## Implemented Models

These R source files extend the functionality of the BTYD package by providing functions for parameter estimation and scoring for NBD, MBG/NBD, BG/CNBD-k, MBG/CNBD-k, Pareto/NBD (HB), Pareto/NBD (Abe) and Pareto/GGG.

* **NBD** Ehrenberg, Andrew SC. "The pattern of consumer purchases." Applied Statistics (1959): 26-41. \doi{10.2307/2985810}
* **MBG/NBD** Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical validation and comparison of models for customer base analysis. International Journal of Research in Marketing 24(3) 201–209. \doi{10.1016/j.ijresmar.2006.12.005}
* **(M)BG/CNBD-k** Reutterer, T., Platzer, M., & Schroeder, N. (2020). "Leveraging purchase regularity for predicting customer behavior the easy way." International Journal of Research in Marketing. \doi{10.1016/j.ijresmar.2020.09.002}
* **Pareto/NBD (HB)** Ma, Shao-Hui, and Jin-Lan Liu. "The MCMC approach for solving the Pareto/NBD model and possible extensions." Natural Computation, 2007. ICNC 2007. Third International Conference on. Vol. 2. IEEE, 2007. \doi{10.1109/ICNC.2007.728}
* **Pareto/NBD (Abe)** Abe, Makoto. "Counting your customers one by one: A hierarchical Bayes extension to the Pareto/NBD model." Marketing Science 28.3 (2009): 541-553. \doi{10.1287/mksc.1090.0502}
* **Pareto/GGG** Platzer, Michael, and Thomas Reutterer. "Ticking Away the Moments: Timing Regularity Helps to Better Predict Customer Activity." Marketing Science (2016). \doi{10.1287/mksc.2015.0963}

## Contributions

We certainly welcome all feedback and contributions to this package! Please use [GitHub Issues](https://github.com/mplatzer/BTYDplus/issues) for filing bug reports and feature requests, and provide your contributions in the form of [Pull Requests](https://help.github.com/articles/about-pull-requests/). See also [these general guidelines](https://guides.github.com/activities/contributing-to-open-source/#contributing) to contribute to Open Source projects on GitHub.
