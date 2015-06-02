BTYDplus
========

Extension to R package BTYD.

Installation
------------

```
devtools::install_github("mplatzer/BTYDplus", dependencies=TRUE)
library(BTYDplus)
demo(package="BTYDplus")
demo("cdnow")
```

BTYD Models
-----------

These R source files extend the functionality of the BTYD package by providing functions for parameter estimation and scoring for NBD, G/G/NBD, BG/NBD, CBG/NBD and CBG/CNBD-k models.

* NBD (MLE) - Ehrenberg, Asc. "The Pattern of Consumer Purchases." Quantitative techniques in marketing analysis: text and readings (1962): 355.

* Gamma/Gompertz/NBD (MLE) - Bemmaor, Albert C., and Nicolas Glady. "Modeling Purchasing Behavior with Sudden Death: A Flexible Customer Lifetime Model." Management Science 58.5 (2012): 1012-1021.

* MBG/NBD (MLE) - Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical validation and comparison of models for customer base analysis. International Journal of Research in Marketing 24(3) 201–209.

* MBG/CNBD-k (MLE) - Platzer, Michael, and Thomas Reutterer. forthcoming...

* Pareto/NBD (HB) - Ma, Shao-Hui, and Jin-Lan Liu. "The MCMC approach for solving the Pareto/NBD model and possible extensions." Natural Computation, 2007. ICNC 2007. Third International Conference on. Vol. 2. IEEE, 2007. http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404 - Abe, Makoto. "Counting your customers one by one: A hierarchical Bayes extension to the Pareto/NBD model." Marketing Science 28.3 (2009): 541-553.

* Pareto/NBD variant (HB) - Abe, Makoto. "Counting your customers one by one: A hierarchical Bayes extension to the Pareto/NBD model." Marketing Science 28.3 (2009): 541-553.

* Pareto/GGG (HB) - Platzer, Michael, and Thomas Reutterer. forthcoming...
