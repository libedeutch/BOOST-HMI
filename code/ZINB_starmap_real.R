#***********************************************************
#* section 1: load packages
#* STARmap read data analysis
#* Apr 4 2023
#***********************************************************

#install.packages("optparse)
# Add control paramater
library(optparse)
option_list<- list(
  make_option("--model_estimate", default =  TRUE, help = "Bayesian estimate of the model"),
  make_option("--simulation", default = FALSE, help = "Generate simulation datasets"),
  make_option("--chainNum", default = 1, help = "Generate simulation datasets"),
  make_option("--idx", default = 1, help = "Gene Index to start"),
  make_option("--phi_start", default = TRUE, help = "Initialize phi from kmeans modelling")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#***********************************************************
#* section 2: load real data and cpp code
#***********************************************************

source("~/fish/real/functions.R") 
Rcpp::sourceCpp("~/fish/real/starmap/zinbm.cpp")
load("~/fish/real/starmap.RData")


#***********************************************************
#* section 3: Estimate
#***********************************************************

if(opt$model_estimate){
  # rre = list()
  # Quality control
  
  # QC step 1: spots with at least 100 counts
  index.qc1 = which(rowSums(count)<100)  
  count <- count[-index.qc1,]
  loc <- loc[-index.qc1,]
  
  # QC step 2: genes express on at least 10\% 
  cols = apply(count,2,function(x){sum(x!=0)})
  index.qc2 = which(cols/nrow(count)>=0.1)
  count = count[,index.qc2]
  
  # QC step 3: maximum count at least 10
  
  max_count = apply(count,2,max)
  index.qc3 = which(max_count>=10)
  count = count[,index.qc3]
   #  for(i in 1:nrow(index.qc3)){
   #    if (count[index.qc3[i,1],index.qc3[i,2]]<10){
   #      print(i)
   #    }
   #  }
    # 
    
  rowsum = rowSums(count)
  scaler = exp(-1/nrow(count) * sum(log(rowsum)))
  si = scaler * rowsum
  idx_end = ncol(count)
  rre = list()
  #for(gg in opt$idx:(opt$idx+421)){
  for(gg in 1:idx_end){
 # for(gg in 1:1){
    print(gg)
    # prod(si) = 1
    #si = exp(rnorm(nrow(xi),0,0.2)); 
    #si = exp(log(si) - sum(log(si))/length(si))
    x = loc[,1]
    y = loc[,2]
    dd = max(max(x) - min(x),max(y) - min(y) )
    x= (x-min(x))/dd
    y= (y-min(y))/dd
    
    Q <- 2;          # Number of marks
    c <- 0.05;        # Tuning parameter c that defines the neighborhood for each point
    # c change to 0.18 in this case based on MCF plot
    iter <-10000;   # Number of iterations
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
    sigma_theta <- 2; # 95.44% theta in (-4,4)
    mu_omega <- 1;
    sigma_omega <- 1;
    
    # Prior in proposal distribution
    tau_mu <- 0.1 #sqrt(0.01)
    tau_theta<- 0.1
    tau_lambda = 0.1 # specify a smaller tau for lambda
    tau_omega = 0.1
    tau_phi = 0.1 # if we want to limit the movement of tau_phi, it cant be larger than 1
    phi_lambda = 10;
    # parameter configuration
    theta_start <- rnorm(Q*(Q + 1)/2, mu_theta, sigma_theta); 
    theta_start[length(theta_start)] <- 0;
    theta_start[1] <-0
    Theta_start <- array2matrix_r(theta_start, Q);
    #Theta_start <- matrix(c(0,-1.2,-1.2,0), ncol =2)
    omega_start <- #c(1.8472979,1)#
      rep(1, Q);
    #lambda_start <- rgamma(1, 10, 1);
    lambda_start = 60;
    #mu_start<- c(10,50); #??? can't directly specify
    #phi_start <- c(1,1) 
    #mu_start <- rgamma(2,1,1)
    #mu_start = c(10,50)
    
    kmeans.mu = kmeans(count[,gg], centers = 2)
    xi = kmeans.mu$cluster-1
    if(kmeans.mu$centers[1] > kmeans.mu$centers[2]){
      xi = 1-xi
    }
    mu_start = sort(kmeans.mu$centers)
    mu_start[which(mu_start==0)] <- 1 # take care of
    H = rep(0,nrow(count))    
    
    #phi_start <- rgamma(2,20,1)
    #mu_start <-rgamma(2,1,1)
    # if(opt$phi_start){
    #   group0 = mu_start[1]*mu_start[1] / (var(count[which(xi==0),gg]) - mu_start[1])
    #   group1 = mu_start[2]*mu_start[2] / (var(count[which(xi==1),gg]) - mu_start[2])
    #   phi_start = c(group0, group1)
    # }
    # else{
    #   group0 = mu_start[1]*mu_start[1] / (var(count[,gg]) - mu_start[1])
    #   group1 = mu_start[2]*mu_start[2] / (var(count[,gg]) - mu_start[2])
    #   phi_start <- c(group0,group1)
    # }
    phi_start <- c(1,1)
    print(paste0("mu", mu_start))
    print(paste0("phi",phi_start))
    #yy = y_simu
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
    # model_estimator(arma::vec z, NumericVector y,IntegerMatrix E, NumericVector d, LogicalVector dd, IntegerVector id_start, IntegerVector id_end, 
    #                      IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta_s, NumericVector omega_s, double lambda_s, double mu_theta, double sigma_theta, 
    #                      double mu_omega, double sigma_omega, double a_lambda, double b_lambda, int iter, int burn, int M, double tau_mu, double phi_lambda,NumericVector mu_y_s, NumericVector phi_y_s,
    #                      NumericVector si,double tau_lambda, double tau_theta,double tau_omega, double tau_phi, double a_phi, double b_phi, double a_mu,double b_mu ) {
    #   
    re <- model_estimator(H,xi,count[,gg], edge, distance, duplicate, id_start, id_end, flag_start, flag_end, Theta_start, omega_start,
                          lambda_start, mu_theta, sigma_theta, mu_omega, sigma_omega, a_lambda, b_lambda, iter, burn, M,
                          tau_mu, phi_lambda ,mu_start, phi_start, si,tau_lambda,tau_theta,tau_omega,tau_phi,a_phi,b_phi,a_mu,b_mu);
    
    end_time <- proc.time();
    time <- end_time - start_time;
    print(time)
    rre[[gg]] = re
    }
    #Path_rre = "/Users/gilly/Library/CloudStorage/OneDrive-Personal/research/fish/real/STARmap/HM/"
    #Path_rre = "~/OneDrive/research/fish/real/STARmap/starmapHM.RData"
    #Path_rre = paste0("~/fish/real/starmap/chain",opt$chainNum,"/starmapHM_",gg, ".RData")
    Path_rre = paste0("~/fish/real/starmap/","starmapHM_chain_",opt$chainNum, ".RData")
    save(rre,phi_start,time,file = Path_rre)
  
  
}
