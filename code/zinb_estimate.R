
#***********************************************************
#* section 1: load packages
#* Revision: 03/22/2023
#* R file to generate simulation data
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


source("~/fish/simulation/functions.R") 
Rcpp::sourceCpp("~/fish/simulation/zinbm.cpp")


if(opt$model_estimate){
#  PathP = "/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/count2/1_3/"
ratio = c("1_2","1_3","1_5","1_10")
for (rr in ratio){
  #PathP = paste0("fish/simulation/count2/",opt$ratio,"/")
  PathP = paste0("fish/simulation/count2/",rr,"/")
  file_tot = list.files(paste0(PathP,"count/"),".RData")
  file_tot2 = list.files(paste0(PathP,"BOOST_HM/"),".RData")
  idx = which(!file_tot %in% file_tot2)
  #first_idx = opt$idx
  #for (kk in first_idx:(first_idx+23)){ # 119 for Qiwei node, 59 for normal or smallmem
  for (kk in idx){
    #idx = 30
    load(paste0(PathP,"count/", file_tot[kk]))
    print(file_tot[kk])
    rre=list()
    idx_end = ncol(y_count)
    si = si
    #si = exp(log(si) - sum(log(si))/length(si))
    for(g in 1:idx_end){
      print(g)
      # prod(si) = 1
      x = loc[,1]
      y = loc[,2]
      Q <- 2;          # Number of marks
      c <- 0.1;        # Tuning parameter c that defines the neighborhood for each point
      # c change to 0.18 in this case based on MCF plot
      iter <-20000;   # Number of iterations
      burn <- iter/2;  # Number of burn-in 
      M <- 1;
      #tau <- sqrt(0.01); # small tau leads to steady result
      
      a_phi = 0.01
      b_phi= 100
      a_mu = 0.01; # a_mu can't be too small as 0.7
      b_mu= 100; 
      a_lambda <- 0.01;
      b_lambda<- 100;
      mu_theta <- 0;
      sigma_theta <- 3.5; # 95.44% theta in (-4,4)
      mu_omega <- 0;
      sigma_omega <- 0.5;
      
      # Prior in proposal distribution
      tau_mu <- 0.5 #sqrt(0.01)
      tau_theta<- 0.5
      tau_lambda = 0.5 # specify a smaller tau for lambda
      tau_omega = 0.5
      tau_phi = 0.5 # if we want to limit the movement of tau_phi, it cant be larger than 1
      
      phi_lambda <- 30;  # if lambda is proposed from Gamma distribution
      
      # parameter configuration
      # theta_start <- rnorm(1, mu_theta, sigma_theta); 
      # theta_start <- matrix(c(0,theta_start,theta_start,0), ncol =2)
      theta_start <- matrix(c(0,parameter_ref[g,'theta_01'],parameter_ref[g,'theta_01'],0), ncol =2)
      omega_start <- c(parameter_ref[g,'omega_0'],1)  #rep(1, Q);
      lambda_start <-  30; 
      mu_start = mu;
      phi_start <- c(phi[g],phi[g]);
      H = rep(0,nrow(y_count))
      
      # build edges and distance
      id_start <- 1;
      id_end <- nrow(y_count);
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
      re <- model_estimator(H, xi[,g],y_count[,g], edge, distance, duplicate, id_start, id_end, flag_start, flag_end, theta_start, omega_start, 
                            lambda_start, mu_theta, sigma_theta, mu_omega, sigma_omega, a_lambda, b_lambda, iter, burn, M, 
                            tau_mu, phi_lambda ,mu_start, phi_start, si,tau_lambda,tau_theta,tau_omega,tau_phi,a_phi,b_phi,a_mu,b_mu);
      end_time <- proc.time();
      time <- end_time - start_time;
      rre[[g]] = re
    }
    
    #Path_ree = paste0("~/fish/simulation/count2/", opt$ratio,"/BOOST_HM/")
    Path_ree = paste0("~/fish/simulation/count2/",rr,"/BOOST_HM/")    
save(rre,file = paste0(Path_ree,file_tot[kk]))
  }
}  
}
