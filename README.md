# PredictingBlackSwans
An R-Package to reproduce the results of my master thesis.

To install this package from R, you first need to install the following package

```R
install.packages("devtools")
```

Having installed it, then you can install the package by just running

```R
devtools::install_github("PabloRMira/PredictingBlackSwans/PredictingBlackSwans")
```

Remember to load the package in your workspace before using:

```R
require(PredictingBlackSwans)
```

To replicate our first simulation study (Prediction Contest) just run

```R
predictionSim()
```

To replicate our second simulation study (High-Dimensional Inference Study) we recommend to run the script `inferenceSim` in script modus sequentially (in contrast to just run the entire function) to split the huge computational burden involved in it. Note: For a UNIX computer with many kernels and great computational power, the implemented parallelization may shorten the computation time by a great margin.

The results of our application concerning prediction (Predicting Black Swans) can be replicated by running

```R
predictingBS()
```

Finally, to replicate the results of the second part of our application (Analyzing the Imbalances prior to a Financial Crisis), this can be accomplished by running


```R
analyzingBS()
```
