# BOOST-HMI

BOOST-HMI is a generalized energy-based Bayesian framework method that is proposed to identify spatial variable (SV) genes from the imaging-based spatially resolved transcriptomic (SRT) dataset. A hidden gene expression level indicator is introduced to dichotomize the gene expression counts into two levels: highly-expressed and low-expressed. It models the gene expression counts using a zero-inflated negative binomial mixture distribution and the interaction between highly-expressed and lowly-expressed cells is characterized by a hidden Bayesian mark interaction model. 

The following R packages are required to run the model <br/>
+ Rcpp
+ RcppArmadillo

# Usage

The main steps of BOOST-HMI are the following:

1. Preparation of gene expression count file <br/>
2. Preparation of spatial location file <br/>
3. BOOST-HMI analysis <br/>
4. Downstream analysis <br/>

## Example data
One replicate of the Mouse hippocampus seqFISH data cohort is  #Data/hippocampus_field_43.Rdata 




  

