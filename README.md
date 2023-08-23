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

One replicate of the Mouse hippocampus seqFISH data cohort is Data/hippocampus_field_43.Rdata <br/>
The gene expression count file has the following format

|  |Tal1|Dmbx1|Emx2|...|fbll1|
|-----|-----|-----|-----|-----|-----|
|Cell 1| 3|11|13|...|28|
|Cell 2|3|9|14|...|20|
|Cell 3|12|3|1|...|25|
|...|...|...|...|...|...|
|Cell 257|1|0|0|...|6|

The geostatistical profile has the following format

|  |x|y|
|-----|-----|-----|
|Cell 1| 686|352|
|Cell 2|488|572|
|Cell 3|178|624|
|...|...|...|
|Cell 257|461|234|

## BOOST-HMI analysis

```{R}
source("code/functions.R")
Rcpp::sourceCpp("code/zinbm.cpp")
load("data/hippocampus_field_43.Rdata")
# prior specification
Q <- 2
c <- 0.15
iter <- 50000
burn <- iter/2
M <- 1
rowsum = apply(count,1,sum)
  scaler = exp(-1/nrow(count) * sum(log(rowsum)))
  si = scaler * rowsum
  a_phi = 0.01
  b_phi= 100
  a_mu = 0.01; 
  b_mu= 100; 
  a_lambda <- 0.01;
  b_lambda<- 100;
  mu_theta <- 0;
  sigma_theta <- 3.5; # 95.44% theta in (-4,4)
  mu_omega <- 1;
  sigma_omega <- 0.5;
  # Prior in proposal distribution
  tau_mu <- 0.1 
  tau_theta<- 0.1
  tau_omega = 0.1
  tau_phi = 0.1   
  # parameter configuration
  theta_start <- rnorm(1, mu_theta, sigma_theta); 
  theta_start <- matrix(c(0,theta_start,theta_start,0), ncol =2)
  omega_start <- rep(1,Q)
  H = rep(0,nrow(count))
  phi_start <- c(1,1)
  mu_start <- c(0,0)
  # build edges and distance
  id_start <- 1;
  id_end <- nrow(count);
  build <- dist_list(x, y, c);
  edge <- build$edge;
  distance <- build$distance;
  duplicate <- build$duplicate;
  flag_start <- build$flag_start;
  flag_end <- build$flag_end;
  # Implement MCMC algorithm
   re <- model_estimator(H,z_start,count[,g], edge, distance, duplicate, id_start, id_end, flag_start, flag_end, theta_start, omega_start, 
                        lambda_start, mu_theta, sigma_theta, mu_omega, sigma_omega, a_lambda, b_lambda, iter, burn, M, 
                        tau_mu, phi_lambda ,mu_start, phi_start, si,tau_lambda,tau_theta,tau_omega,tau_phi,a_phi,b_phi,a_mu,b_mu);
```




  

