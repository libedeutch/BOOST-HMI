# ZINB seqfish 

#***********************************************************
#* section 1: load packages
#* Revision: 04/03/2023
#***********************************************************

#install.packages("optparse)
# Add control paramater
library(optparse)
library(ggplot2)
option_list<- list(
  make_option("--model_estimate", default =  FALSE, help = "Bayesian estimate of the model"),
  make_option("--simulation", default = TRUE, help = "Generate simulation datasets"),
  make_option("--trace_plot", default = FALSE, help = "Trace plot of posterior samples"),
  make_option("--ratio", default = "1_3", help = "Ratio to generate simulation"),
  make_option("--idx", default = 1, help = "control replicates")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


#source("~/fish/simulation/functions.R") 
source("～/fish/mark_interaction_model/functions.R")
#Rcpp::sourceCpp("/Users/gilly/Library/CloudStorage/OneDrive-Personal/research/fish/mark_interaction_model/functions_nb_ratioTest.cpp")
Rcpp::sourceCpp("～/fish/mark_interaction_model/zinbm.cpp")
load("~/hippocampus_field_43.Rdata")
#load("~/fish/simulation/hippocampus_field_43.Rdata")

index  = which(loc[,1]<203 | loc[,1]>822 | loc[,2]<203 | loc[,2]>822)
count = count[-index,]
loc= loc[-index,1:2]
dd = max(max(loc[,1] ) - min(loc[,1]),max(loc[,2] ) - min(loc[,2]) )
x <- (loc[,1] - min(loc[,1])) / (dd)
y <-(loc[,2] - min(loc[,2])) / (dd)

rre = list() #       H = rep(0,nrow(y_count))
for (g in 1:ncol(count)){
  print(g)
  Q <- 2;          # Number of marks
  c <- 0.15;        # Tuning parameter c that defines the neighborhood for each point
  lambda_start <-  20; 
  # c change to 0.18 in this case based on MCF plot
  iter <-50000;   # Number of iterations
  burn <- iter/2;  # Number of burn-in 
  M <- 1;
  #si <- exp(rnorm(dim(count)[1], mean = 0, sd = 0.2))
  #spark@lib_size <-apply(spark@counts, 2, sum)
  rowsum = apply(count,1,sum)
  scaler = exp(-1/nrow(count) * sum(log(rowsum)))
  si = scaler * rowsum
  a_phi = 0.01
  b_phi= 100
  a_mu = 0.01; # a_mu can't be too small as 0.7
  b_mu= 100; 
  a_lambda <- 0.01;
  b_lambda<- 100;
  mu_theta <- 0;
  sigma_theta <- 3.5; # 95.44% theta in (-4,4)
  mu_omega <- 1;
  sigma_omega <- 0.5;
  
  # Prior in proposal distribution
  tau_mu <- 0.1 #sqrt(0.01)
  tau_theta<- 0.1
  tau_lambda = 0.1 # specify a smaller tau for lambda
  tau_omega = 0.1
  tau_phi = 0.1 # if we want to limit the movement of tau_phi, it cant be larger than 1
  
  phi_lambda <- 30;  # if lambda is proposed from Gamma distribution
  
  # parameter configuration
  theta_start <- rnorm(1, mu_theta, sigma_theta); 
  theta_start <- matrix(c(0,theta_start,theta_start,0), ncol =2)
  omega_start <- rep(1,Q)
  H = rep(0,nrow(count))
  
  tmp =  kmeans(count[,g], centers=2)
  if(tmp$centers[1] > tmp$centers[2]){
    z_start <-2 - tmp$cluster
    mu_start <- c(tmp$centers[2],tmp$centers[1]);
    phi_start <- mu_start^2/
      (c(var(count[which(tmp$cluster==2),g]), var(count[which(tmp$cluster==1),g]))- mu_start);
  }
  else {
    z_start <-tmp$cluster -1
    mu_start <- c(tmp$centers[1],tmp$centers[2]);
    phi_start <- mu_start^2/
      (c(var(count[which(tmp$cluster==1),g]), var(count[which(tmp$cluster==2),g]))- mu_start);
    
  }
  print(mu_start)
  phi_start[phi_start<0] <- 1
  phi_start <- c(1,1)
  # build edges and distance
  id_start <- 1;
  id_end <- nrow(count);
  build <- dist_list(x, y, c);
  edge <- build$edge;
  distance <- build$distance;
  duplicate <- build$duplicate;
  flag_start <- build$flag_start;
  flag_end <- build$flag_end;
  # ========================================================================================
  # ========================================================================================
  # Implement MCMC algorithm (DO NOT EDIT THE CODE IN THIS BLOCK)
  # ========================================================================================
  # ========================================================================================
  start_time <- proc.time();
  # Y <- model_estimator(z, edge, distance, duplicate, id_start, id_end, flag_start, flag_end, Theta_start, omega_start, lambda_start, mu_theta, sigma_theta, mu_omega, sigma_omega, a, b, iter, burn, M, tau, phi);
  re <- model_estimator(H,z_start,count[,g], edge, distance, duplicate, id_start, id_end, flag_start, flag_end, theta_start, omega_start, 
                        lambda_start, mu_theta, sigma_theta, mu_omega, sigma_omega, a_lambda, b_lambda, iter, burn, M, 
                        tau_mu, phi_lambda ,mu_start, phi_start, si,tau_lambda,tau_theta,tau_omega,tau_phi,a_phi,b_phi,a_mu,b_mu);
  end_time <- proc.time();
  time <- end_time - start_time;
  if(sum(is.na(re$z_colmeans)) > 0) print(paste0("Alert", g ))
  rre[[g]] = re
}

save(rre, file = "~/OneDrive/research/fish/real/MIST/seqfish_lambda20.RData")


      
     
