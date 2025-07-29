# Bayesian Response-Adaptive Randomization for Cluster Randomized Controlled Trials

This repository contains the code associated with our paper *"Bayesian Response-Adaptive Randomization for Cluster Randomized Controlled Trials."*

## Contents

- `BARAcRCT_design/`: R code for trial design and simulation studies.
- `app.R`: Shiny application implementing the Bayesian response-adaptive randomization.
- `BRAR_cRCT_vary.stan`: Supporting Stan model file used by the Shiny app.

## Launch the Shiny App

You can launch the app directly from R using the following command:

```r
shiny::runGitHub("BRARcRCT", "Renieliu")
