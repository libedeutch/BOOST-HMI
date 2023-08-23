# BOOST-HMI

BOOST-HMI is a generalized energy-based Bayesian framework method that is proposed to identify spatial variable (SV) genes from the imaging-based spatially resolved transcriptomic (SRT) dataset. A hidden gene expression level indicator is introduced to dichotomize the gene expression counts into two levels: highly-expressed and low-expressed. It models the gene expression counts using a zero-inflated negative binomial mixture distribution and the interaction between highly-expressed and lowly-expressed cells is characterized by a hidden Bayesian mark interaction model. 

The following R packages are required to run the model <br/>
+ ```{r}
  Rcpp
  ```
+ ```{r}
  RcppArmadillo
  ```
  

