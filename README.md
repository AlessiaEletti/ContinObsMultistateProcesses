# "A Spline-Based Framework for the Flexible Modelling of Continuously Observed Multistate Survival Processes"

### Alessia Eletti, Giampiero Marra and Rosalba Radice

This repository contains the code used to produce the case study present in the paper 

> Eletti, A., Marra, G. and Radice, R. A Spline-Based Framework for the Flexible Modelling of Continuously Observed Multistate Survival Processes *(currently under submission)*.

The dataset is openly accessible and involves primary breast cancer information originating from the *Rotterdam Breast Cancer Study*. It is modelled by means of a continuously observed multistate survival process to account for the multiple outcomes of interest. Note that the analysis should not be used to make epidemiological claims and serves an illustrative purpose only. 

The models are fitted using the R package `GJRM`. Prediction on the transition probabilities scale are carried out using a modified version of the `mssample()` function from the R package `mstate`.
