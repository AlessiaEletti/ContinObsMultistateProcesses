# "A Spline-Based Framework for the Flexible Modelling of Continuously Observed Multistate Survival Processes"

### Alessia Eletti, Giampiero Marra and Rosalba Radice

This repository contains the code used to produce the analysis of primary breast cancer data present in the paper 

> Eletti, A., Marra, G. and Radice, R. A Spline-Based Framework for the Flexible Modelling of Continuously Observed Multistate Survival Processes *(currently under submission)*.

The data are openly accessible and originate from the *Rotterdam Breast Cancer Study*. Note that the analysis should not be used to make epidemiological claims and serves an illustrative purpose only. 

The models are fitted using the R package `GJRM`. Prediction on the transition probabilities scale are carried out using a modified version of the `mssample()` function from the R package `mstate`.
