# Hierarchical Bayesian Survival model

This is an R implementation of the paper ["A Hierarchical Bayesian Model for Personalized Survival Predictions"](http://medianetlab.ee.ucla.edu/papers/FINAL_A_Hierarchical_Bayesian_Model_for_Personalized_Survival_Predictions.pdf). 

In this project we propose a probabilistic survival model which flexibly captures individual traits through a hierarchical latent variable formulation. We estimate survival trajectories by jointly sampling the location and shape of the individual survival distribution resulting in patient-specific curves with quantifiable uncertainty estimates. Our aim is to improve predictions by making them personalized to the patient-at-hand, to better understand diseases and their risk factors, and to provide interpretable model outputs to clinicians.

This repository includes the main model training function HBS.R and various utility functions to compute performance, variable importance, variable interaction importance, tests and other visualizations. 

*Please cite the above paper if this resource is used in any publication*

## Requirements

* R version 3.5
* Packages: "BayesTree","flexsurv" for the main model; packages "bartMachine" and "ICEbox" for the variable importance procedures; and package "pec" for performance computations.

## First steps
To get started, check demo.R which will guide you through an application of our algorithm with variables importance and prediction computations. 

If you have questions or comments about anything regarding this work, please do not hesitate to contact [Alexis](https://alexisbellot.github.io/Website/)
