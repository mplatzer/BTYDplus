BTYDplus
========

Extension to R package BTYD.

Installation
------------

```
devtools::install_github("mplatzer/BTYDplus")
library(BTYDplus)
demo(package="BTYDplus")
demo("cdnow")
```

BTYD Models
-----------

These R source files extend the functionality of the BTYD package by providing functions for parameter estimation and scoring for NBD, G/G/NBD, BG/NBD, CBG/NBD and CBG/CNBD-k models.

* NBD
EHRENBERG, ASC. "The Pattern of Consumer Purchases." Quantitative techniques in marketing analysis: text and readings (1962): 355.

* G/G/NBD
Bemmaor, Albert C., and Nicolas Glady. "Modeling Purchasing Behavior with Sudden “Death”: A Flexible Customer Lifetime Model." Management Science 58.5 (2012): 1012-1021.

* BG/NBD
Fader, Peter S., Bruce GS Hardie, and Ka Lok Lee. "Counting Your Customers the Easy Way: An Alternative to the Pareto/NBD Model." Marketing Science 24.2 (2005): 275-284.

* CBG/NBD
Hoppe, Daniel, and Udo Wagner. "Customer base analysis: The case for a central variant of the Betageometric/NBD Model." Marketing Journal of Research and Management 3.2 (2007): 75-90.

* CBG/CNBD-k
Platzer, Michael. "Stochastic models of noncontractual consumer relationships." Master of Science in Business Administration thesis, Vienna University of Economics and Business Administration, Austria (2008).
https://sites.google.com/site/michaelplatzer/stochastic-models-of-noncontractual-consumer-relationships

* Pareto/NBD MCMC
* Pareto/CNBD MCMC

Note, that the Pareto/NBD model is already included as part of the BTYD-package, but we provide a helper method to generate artificial data following the Pareto/NBD assumptions.

A high-level introduction to CBG/CNBD-k and some reasoning for incorporating regularity in timing patterns into the modelling can also be found at http://www.slideshare.net/mplatzer/my-entry-to-the-dmef-clv-contest
