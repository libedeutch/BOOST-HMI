// rcpp for zero-inflated negative binomial mixture model
// All parameters:
// H for eta_ij 



#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
/* 
  In this document, gamma distribution is parameterized with shape and scale (not rate)
This version is for simulation data without excessive zeros 
*/
arma::vec rowSums(IntegerMatrix z_store);
arma::vec colMeans(IntegerMatrix z_store);
double rnorm_trunc(double mu, double sigma, double lower, double upper);
double logenergy(arma::vec z, NumericVector zi, IntegerMatrix E, NumericVector d, LogicalVector dd, NumericMatrix Theta, NumericVector omega, double lambda);
arma::vec model(arma::vec z, IntegerMatrix E, NumericVector d, int id_start, int id_end, IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta, NumericVector omega, double lambda);
double logenergy1(arma::vec z,  int qindex,int ii,IntegerMatrix E, NumericVector d, LogicalVector dd, NumericMatrix Theta, NumericVector omega, double lambda);
// [[Rcpp::export]]
Rcpp::List model_estimator(arma::vec H,arma::vec z, NumericVector y,IntegerMatrix E,  NumericVector d, LogicalVector dd, IntegerVector id_start, IntegerVector id_end, 
                           IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta_s, NumericVector omega_s, double lambda_s, double mu_theta, double sigma_theta, 
                           double mu_omega, double sigma_omega, double a_lambda, double b_lambda, int iter, int burn, int M, double tau_mu, double phi_lambda,NumericVector mu_y_s, NumericVector phi_y_s,
                           NumericVector si,double tau_lambda, double tau_theta,double tau_omega, double tau_phi, double a_phi, double b_phi, double a_mu,double b_mu ) {
  int N = id_start.size(); // ? what is N. basically it is one, but why use this 
  int n = z.n_rows;
  int Q = Theta_s.nrow();
  int i, ii, q, qq, qqq, qqqq, count, k;
  int count_2 = 10;
  double hastings = 0;
  double accept_theta = 0;
  double accept_lambda = 0;
  double accept_mu_y = 0;
  double accept_phi_y = 0;
  double accept_omega = 0;
  double lambda = lambda_s;
  double lambda_temp;
  arma::vec z_temp(n);
  arma::vec z_label_switch(n);
  arma::vec zi_label(2);
  NumericMatrix theta_store(iter, Q*(Q + 1)/2);
  NumericVector lambda_store(iter);
  NumericMatrix omega_store(iter, Q);
  NumericMatrix mu_store(iter, Q);
  NumericMatrix phi_store(iter, Q);
  NumericMatrix Theta(Q, Q);
  //NumericVector Theta(Q);
  NumericMatrix Theta_temp(Q, Q);
  NumericVector mu_y(Q);
  NumericVector mu_y_temp(Q);
  NumericVector sum_y(Q);
  NumericVector phi_y(Q);
  //NumericVector phi_y(n);
  NumericVector phi_y_temp(Q);
  //NumericVector Theta_temp(Q);
  NumericVector omega(Q);
  NumericVector omega_temp(Q);
  NumericVector rTemp(Q); // used to ? 
  NumericVector zi(Q); // count z =0 and z = 1
  NumericVector zi_temp(Q); 
  NumericMatrix z_store(iter,n);
  NumericVector rTempMu(Q); // used to? 
  IntegerVector H_sum(iter,0); // save H_sum per spot per gene across iterations 
  IntegerVector H_sum_temp(n,0); // 
  NumericVector pi(n,0.5);
  NumericVector prob_temp_H(2);
  double max_temp;
  double prob0;
  double prob1;
  double sum_temp;
  double H_temp;


  // Initialization
  for(q = 0; q < Q; q++)
  {
    zi(q) = 0;
    zi_temp(q) = 0;
    omega(q) = omega_s(q);
    mu_y(q) = mu_y_s(q);
    phi_y(q) = phi_y_s(q);
    for(qq = 0; qq < Q; qq++)
    {
      Theta(q, qq) = Theta_s(q, qq);
    }
    //Theta(q) = Theta_s(q);
  }

  for(int it = 0; it<n; it++){
    if(y(it)!=0){
    	H(it) = 0;
    }
    else{
    	H(it) = 0;
    	H_sum_temp(it) = H_sum_temp(it) + H(it);
    }
}

  // MCMC
  for(i = 0; i < iter; i++)
  {
    hastings = 0;
    for(int q = 0; q< Q; q++){
      sum_y(q) = 0;
      rTemp(q) = 0;
      rTempMu(q) =0;
      zi(q) =0;// what is zi
    }
    
    
    for(ii = 0; ii < n; ii++)
    {
      //phi_y(ii) = phi_y_s(ii);
      for(q = 0; q < Q; q++) {
        //   if(z(ii) == q + 1)
          if(z(ii) == q) // revise part
        {
          zi(q) = zi(q) + 1;
          sum_y(q) = sum_y(q) +y(ii); // \sum_i^n y_i I(z_i = 0)
          //break;
        }
      }
    }
    
    //Rcout<<sum_y(0);

    // Update H 
    for(int i4h =0; i4h<n; i4h++){
		if(y(i4h)==0){
			if(z(i4h)==0){
				prob_temp_H(0) = lgamma(y(i4h) + phi_y(0)) - lgamma(y(i4h)+1) - lgamma(phi_y(0));
				prob_temp_H(0) = prob_temp_H(0)+ phi_y(0)*(log(phi_y(0)) - log(si(i4h)*mu_y(0)+phi_y(0))) + log(1-pi(i4h)); 
			}
			else{
				prob_temp_H(0) = lgamma(y(i4h) + phi_y(1)) - lgamma(y(i4h)+1) - lgamma(phi_y(1));
				prob_temp_H(0) = prob_temp_H(1)+ phi_y(1)*(log(phi_y(1)) - log(si(i4h)*mu_y(1) +phi_y(1))) + log(1-pi(i4h)); 
			}
		prob_temp_H(1) = log(pi(i4h));
		max_temp = max(prob_temp_H);
		prob_temp_H(0) = prob_temp_H(0) - max_temp;
        prob_temp_H(1) = prob_temp_H(1) - max_temp;
        prob_temp_H(0) = exp(prob_temp_H(0));
        prob_temp_H(1) = exp(prob_temp_H(1));
        sum_temp = prob_temp_H(0) + prob_temp_H(1);
        prob_temp_H(0) = prob_temp_H(0)/sum_temp;
        prob_temp_H(1) = prob_temp_H(1)/sum_temp;
        H_temp = H(i4h);
        H(i4h) = rbinom(1, 1, prob_temp_H(1))(0);
        H_sum_temp(i4h) = H_sum_temp(i4h) - H_temp + H(i4h);

    }
}

    //   Rcout<<"H"; 
    // update mu_y: mu_y are means in NB model 
    
    for (q = 0; q < Q ; q++)
    {	
      rTempMu(q) =0;
      mu_y_temp(q) = mu_y(q);	
      // Restrict mu(1) > mu(0) 
      //mu_y_temp(q) = rgamma(1,mu_y(q)*0.5, 0.5)(0);
      do{
        mu_y_temp(q) = exp(rnorm(1,log(mu_y(q)),tau_mu)(0)); 
        // mu can not be 0 
      } while (mu_y_temp(q) < 0.01);
      

      hastings = (a_mu-1)*(log(mu_y_temp(q)) - log(mu_y(q))) - 1/b_mu*(mu_y_temp(q) - mu_y(q));
      //mu_y_temp(q) = rlnorm(1,log(mu_y(q)*mu_y(q)/sqrt(mu_y(q)*mu_y(q)+tau*tau)),sqrt(log(1+tau*tau/mu_y(q)/mu_y(q))))(0);

      for(int ii =0;ii<n ;ii++){
      	if(H(ii)==0){
      		  if(z(ii)==q){
          rTempMu(q) = rTempMu(q) - (phi_y(q)+y(ii))*(log(si(ii)*mu_y_temp(q)+phi_y(q)) -log(si(ii)*mu_y(q)+phi_y(q)))  + y(ii)*(log(mu_y_temp(q)) - log(mu_y(q)));
        }
      }
 
      }
      hastings = hastings + rTempMu(q);
      //hastings = -log(si*mu_y_temp(q) + phi_y(q)) * (zi(q)*phi_y(q)+sum_y(q)) + log(mu_y_temp(q)) * (sum_y(q) - zi(q)) - log(mu_y_temp(q))* log(mu_y_temp(q))/2/h/phi_y(q);
      //hastings = hastings + log(si*mu_y(q) + phi_y(q)) * (zi(q)*phi_y(q)+sum_y(q)) - log(mu_y(q)) * (sum_y(q) - zi(q)) + log(mu_y(q))* log(mu_y(q))/2/h/phi_y(q);
      //Rcout<<hastings<<"hist"<<mu_y(0)<<"ist"<<mu_y_temp(0)<<"end";
      if (hastings >= log(double(rand() % 10001) / 10000))
      {
        mu_y(q) = mu_y_temp(q); 
        if (i > burn) {
          accept_mu_y++;
        }
      }
    }
    //Rcout<<mu_y<<" mu "; 



    // Update phi_y
    // gamma distribution scale form is used
    for (q = 0; q < Q; q++)
    {
      rTemp(q)=0;
      phi_y_temp(q) = phi_y(q);
      // constrain phi to reasonable level
      do {
        phi_y_temp(q) = exp(rnorm(1,log(phi_y(q)),tau_phi)(0));
        //phi_y_temp(q) = rgamma(1,log(phi_y(q))/tau_phi,tau_phi)(0);
      } while (phi_y_temp(q) < 0.01);
      
      // z is needed
      //hastings = -lgamma(phi_y_temp(q))*n + lgamma(phi_y(q))*n  + n*phi_y_temp(q)* log(phi_y_temp(q))- n*phi_y(q)* log(phi_y(q)) ;
      
      //Rcout<< hastings<< " + "<<zi(q) <<"&" <<zi(q)*phi_y_temp(q)* log(phi_y_temp(q))- zi(q)*phi_y(q)* log(phi_y(q))  << "* "<< phi_y_temp(q) <<"#"<<phi_y(q)<< " "; 
      hastings = (a_phi-1)*log(phi_y_temp(q))  - (a_phi-1)*log(phi_y(q))- phi_y_temp(q)/b_phi + phi_y(q)/b_phi ;
      
      // hastings = hastings - count_z(q) / 2 * log(sigma_y_temp(q)) + count_z(q) / 2 * log(sigma_y(q)) - z.t() * (y - mu_y(q)) * (y - mu_y(q))/2/ sigma_y_temp(q) + z.t() * (y - mu_y(q)) * (y - mu_y(q)) / 2 / sigma_y(q);
      for( ii =0 ; ii<n ;ii++){
      	if(H(ii)==0){
      	  if(z(ii) == q){
          rTemp(q)= rTemp(q) + lgamma(y(ii)+phi_y_temp(q)) - lgamma(y(ii)+phi_y(q)) - (log(si(ii)*mu_y(q) + phi_y_temp(q)))*(phi_y_temp(q) + y(ii)) + (log(si(ii)*mu_y(q) + phi_y(q)))*(phi_y(q) + y(ii));
          //Rcout<<-(log(si(ii)*mu_y(q) + phi_y_temp(q)))*(phi_y_temp(q) + y(ii)) + (log(si(ii)*mu_y(q) + phi_y(q)))*(phi_y(q) + y(ii))<<" ";
          hastings = hastings  -lgamma(phi_y_temp(q)) + lgamma(phi_y(q))  + phi_y_temp(q)* log(phi_y_temp(q))- phi_y(q)* log(phi_y(q)) ;

        }
      		
      	}
        
      }
      hastings = hastings + rTemp(q);
      
      //Rcout<< "{" << hastings<< "}";
      if (hastings >= log(double(rand() % 10001) / 10000))
      {
        //  if (q==1) {Rcout<< "{" << hastings<<  "}" <<i<<" & " << phi_y_temp(q)<<"+"<<phi_y(q); }
        phi_y(q) = phi_y_temp(q);
        if (i > burn) {
          accept_phi_y++;
        }
      }
    }
    
    //Rcout<<phi_y<<" phi ";



    // Update omega
    // only update omega_0
    for(q = 0; q < Q - 1; q++) 
    {
      for(qq = 0; qq < Q; qq++)
      {
        omega_temp(qq) = omega(qq);
      }
      omega_temp(q) = rnorm(1, omega(q), tau_omega)(0);
      for(ii = 0; ii < n; ii++)
      {
        z_temp(ii) = z(ii);
      }
      
      for(ii = 0; ii < M; ii++) // M = 1, one round of gibbs sampler of auxiliary z 
      {
        for(k = 0; k < N; k++)
        {
          z_temp = model(z_temp, E, d, id_start(k), id_end(k), flag_start, flag_end, Theta, omega_temp, lambda);
        }
      }
      
      for(qq = 0; qq < Q; qq++)
      {
        zi_temp(qq) = 0;
      }
      
      for(ii = 0; ii < n; ii++)
      {
        for(qq = 0; qq < Q; qq++) 
        {
          if(z_temp(ii) == qq )
          {
            zi_temp(qq)++;
            //break;
          }
        }
      }
      //hastings = -logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) + logenergy(z, zi, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta, omega_temp, lambda) + logenergy(z_temp, zi_temp, E, d, dd, Theta, omega_temp, lambda);
      hastings = logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta, omega, lambda) + logenergy(z, zi, E, d, dd, Theta, omega_temp, lambda) - logenergy(z_temp, zi_temp, E, d, dd, Theta, omega_temp, lambda);
      //Rcout<< hastings<< "hist1" << logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) <<"ist" << - logenergy(z, zi, E, d, dd, Theta, omega, lambda)<<"end" ;
      hastings = hastings - (omega_temp(q) - mu_omega)*(omega_temp(q) - mu_omega)/2/sigma_omega/sigma_omega+ (omega(q) - mu_omega)*(omega(q) - mu_omega)/2/sigma_omega/sigma_omega;
      //Rcout<< hastings << "hist" << omega_temp(0) << "omega";
      if (hastings >= log(double(rand()%10001)/10000))
      {
        omega(q) = omega_temp(q); // omega becomes new omega
        if (i > burn) {
          accept_omega++;
        }
      }
    }
    
    
    // Update Theta
    // sum_temp = 0;
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = (q+1); qq < Q; qq++)
      {
        //if(q != Q - 1 || qq != Q - 1)
        {
          for(qqq = 0; qqq < Q; qqq++)
          {
            for (qqqq = 0; qqqq < Q; qqqq++)
            {
              Theta_temp(qqq, qqqq) = Theta(qqq, qqqq);
            }
          }
          // do {
            Theta_temp(q, qq) = rnorm(1, Theta(q, qq), tau_theta)(0);
            //lambda_temp = exp(rnorm(1,log(lambda), tau)(0));
            //} while (std::abs(Theta_temp(q,qq)) < 1 );
          
          Theta_temp(qq, q) = Theta_temp(q, qq);
          // Theta_temp(Q - 1, Q - 1) = Theta(Q - 1, Q - 1) + (Theta(q, qq) - Theta_temp(q, qq));
          for(ii = 0; ii < n; ii++)
          {
            z_temp(ii) = z(ii);
          }
          for(ii = 0; ii < M; ii++)
          {
            for(k = 0; k < N; k++)
            {
              z_temp = model(z_temp, E, d, id_start(k), id_end(k), flag_start, flag_end, Theta_temp, omega, lambda);
            }
          }
          for(qqq = 0; qqq < Q; qqq++)
          {
            zi_temp(qqq) = 0;
          }
          for(ii = 0; ii < n; ii++)
          {
            for(qqq = 0; qqq < Q; qqq++) 
            {
              if(z_temp(ii) == qqq )
              {
                zi_temp(qqq)++;
                break;
              }
            }
          }
          //hastings = -logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) +logenergy(z, zi, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta_temp, omega, lambda) + logenergy(z_temp, zi_temp, E, d, dd, Theta_temp, omega, lambda);
          hastings = logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) -logenergy(z, zi, E, d, dd, Theta, omega, lambda) + logenergy(z, zi, E, d, dd, Theta_temp, omega, lambda) - logenergy(z_temp, zi_temp, E, d, dd, Theta_temp, omega, lambda);
          
          hastings = hastings - (Theta_temp(q, qq) - mu_theta)*(Theta_temp(q, qq) - mu_theta)/2/sigma_theta/sigma_theta + (Theta(q, qq) - mu_theta)*(Theta(q, qq) - mu_theta)/2/sigma_theta/sigma_theta;
          // hastings = hastings - (Theta_temp(Q - 1, Q - 1) - mu)*(Theta_temp(Q - 1,  Q - 1) - mu)/2/sigma/sigma + (Theta(Q - 1, Q - 1) - mu)*(Theta(Q - 1, Q - 1) - mu)/2/sigma/sigma;
          if (hastings >= log(double(rand()%10001)/10000))
          {
            Theta(q, qq) = Theta_temp(q, qq);
            Theta(qq, q) = Theta(q, qq);
            // Theta(Q - 1, Q - 1) = Theta_temp(Q - 1, Q - 1);
            if (i > burn) {
              accept_theta++;
            }
          }
        }
        //sum_temp = sum_temp + Theta(q, qq);
      }
    }
    
    
    //Theta(Q - 1, Q - 1) = -sum_temp;
    
    // Updata lambda
    // use rate form of gamma distribution exp(-beta x)
    //mean = alpha/beta, variance = alpha / beta^2 
    // propose to mean = lambda, variance = phi
    /*
      int index_lambda =0;
      do {
        // lambda_temp = rgamma(1, lambda*lambda/phi_lambda, phi_lambda/lambda)(0);
        lambda_temp = exp(rnorm(1,log(lambda), tau_lambda)(0));
        if(index_lambda>500) {Rcout << "loop lambda"; break;}
        else{index_lambda++;}
        // truncated by 20 
      } while (lambda_temp < 15 | lambda_temp > 40);
      hastings = (a_lambda - 1)*(log(lambda_temp) - log(lambda)) - 1/b_lambda*(lambda_temp - lambda);
      for(ii = 0; ii < n; ii++)
      {
        z_temp(ii) = z(ii);
      }
      for (ii = 0; ii < M; ii++)
      {
        for(k = 0; k < N; k++)
        {
          z_temp = model(z_temp, E, d, id_start(k), id_end(k), flag_start, flag_end, Theta, omega, lambda_temp);
        }
      }
      for(q = 0; q < Q; q++)
      {
        zi_temp(q) = 0;
      }
      for(ii = 0; ii < n; ii++)
      {
        for(q = 0; q < Q; q++) 
        {
          if(z_temp(ii) == q )
          {
            zi_temp(q)++;
            break;
          }
        }
      }
      hastings = hastings + logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta, omega, lambda) +logenergy(z, zi, E, d, dd, Theta, omega, lambda_temp) - logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda_temp);
      if (hastings >= log(double(rand()%10001)/10000))
      {
        //Rcout<<lambda_temp<<"\n";
        lambda = lambda_temp;
        if (i > burn) {
          accept_lambda++;
        }
      }
      */



      // Update z  100% accepted
      // Rcout<< lambda<<"lambda\n";
      int index =0;
      //do{
        for(ii =0;ii<n ;ii++){ // ii iterate over edges pair or spots
          prob0 = 0;
          prob1 = 0;
          prob0 = logenergy1(z,0,ii,E,d,dd,Theta,omega,lambda);
          prob1 = logenergy1(z,1,ii,E,d,dd,Theta,omega,lambda);
        if(H(ii)==0){
          prob0 = prob0+ lgamma(y(ii)+phi_y(0)) - lgamma(y(ii)+1) - lgamma(phi_y(0)) +  phi_y(0)*log(phi_y(0)) +  y(ii)*log(si(ii)*mu_y(0))  -(phi_y(0)+y(ii))*log(si(ii)*mu_y(0) + phi_y(0));
          prob1 = prob1+ lgamma(y(ii)+phi_y(1)) - lgamma(y(ii)+1) - lgamma(phi_y(1)) + phi_y(1)*log(phi_y(1)) + y(ii)*log(si(ii)*mu_y(1)) -(phi_y(1)+y(ii))*log(si(ii)*mu_y(1)+ phi_y(1));
         } 
         // label switching 
          for(int q= 0; q<Q;q++){
            zi_label(q) = 0;
            for(int jj=0;jj<n;jj++){
              if(z(jj)==q) {
                zi_label(q) = zi_label(q)+1;
              }
            }
          }
          
          
          if(zi_label(1) >= n-10) {
            z(ii) = 0;
          }
          if(zi_label(1) <= 10  ){
            z(ii) = 1;
          }
          if(zi_label(1) > 10 & zi_label(1) < (n-10) ){
            z(ii) = rbinom(1,1,1.0/(exp(prob0 - prob1)+1))(0);
          }

        }

          // z_label_switch(ii) = z(ii);
          // Rcout<< "{" << 1.0/(exp(prob0 - prob1)+1) << "\n";
          // Rcout<< prob0<<"~"<< prob1<<"/n";
          // Rcout<<1/(exp(log(prob0) - log(prob1))+1)<<"/n";
          // z(ii) = rbinom(1,1,1.0/(exp(prob0 - prob1)+1))(0);
        
        
        
        if(index >500) {Rcout<< " loop in z "<< prob0 << " & " << prob1 <<  " PHI " <<phi_y(0) << "&"<< phi_y(1)<< " MU "<< mu_y(0) << "&"<<mu_y(1)<<"\n" ; break;}
        else{
        	index++;}
        //} while (zi_label(1)<10  );
      
      
	      // label switching 
	      arma::vec zzz(n,1);

	      double tmp; 
	      // if generated z interchange the abundance of  high and low expressed states 
	      if(mu_y(0) > mu_y(1)) {
	        // z = zzz - z; // Do we need to also switch z? 
	        tmp = mu_y(0);
	        mu_y(0) = mu_y(1);
	        mu_y(1) = tmp;
	      }
      
      // Monitor the process
      //
      if (i*100/iter == 100)
      {
        Rcout<<100<< "% has been done\n";
        count_2 = count_2 + 10;
      }
      count = 0;
      for(ii =0; ii<n; ii++){
        z_store(i,ii) = z(ii);
        H_sum(i) = H_sum(i) + H_sum_temp(ii); 
      }
      for(q = 0; q < Q; q++)
      {
        mu_store(i, q) = mu_y(q);
        phi_store(i, q) = phi_y(q);
        omega_store(i, q) = omega(q);
        for(qq = q; qq < Q; qq++)
        {
          theta_store(i, count) = Theta(q, qq);
          count++;
        }
      }
      lambda_store(i) = lambda;
  }
  accept_omega = accept_omega/(iter - burn)/(Q - 1);
  accept_theta = accept_theta/(iter - burn)/(Q*(Q + 1)/2 - 1);
  accept_lambda = accept_lambda/(iter - burn);
  accept_mu_y = accept_mu_y / (iter - burn) / (Q);
  accept_phi_y = accept_phi_y / (iter - burn) / (Q);
  NumericVector z_summed_by_row = rowSums(z_store(Range(iter/2,iter),_));
  NumericVector z_colmeans = colMeans(z_store(Range(iter/2,iter),_));
  
  return Rcpp::List::create(Rcpp::Named("H_sum") = H_sum,Rcpp::Named("mu") = mu_store, Rcpp::Named("phi") = phi_store, Rcpp::Named("z_summed_by_row") = z_summed_by_row, Rcpp::Named("z_colmeans") = z_colmeans,Rcpp::Named("omega") = omega_store, Rcpp::Named("Theta") = Theta,Rcpp::Named("theta") = theta_store, Rcpp::Named("lambda") = lambda_store, Rcpp::Named("accept_theta") = accept_theta,Rcpp::Named("accept_omega") = accept_omega, Rcpp::Named("accept_lambda") = accept_lambda, Rcpp::Named("accept_mu") = accept_mu_y, Rcpp::Named("accept_phi") = accept_phi_y );
//return Rcpp::List::create(Rcpp::Named("mu") = mu_store, Rcpp::Named("phi") = phi_store, Rcpp::Named("z") = z_store, Rcpp::Named("omega") = omega_store, Rcpp::Named("Theta") = Theta,Rcpp::Named("theta") = theta_store, Rcpp::Named("lambda") = lambda_store, Rcpp::Named("accept_theta") = accept_theta,Rcpp::Named("accept_omega") = accept_omega, Rcpp::Named("accept_lambda") = accept_lambda, Rcpp::Named("accept_mu") = accept_mu_y, Rcpp::Named("accept_phi") = accept_phi_y );

}

// [[Rcpp::export]]
Rcpp::List dist_list(NumericVector x, NumericVector y, double c) {
  int n = x.size();
  int i,j;
  int count = 0;
  double temp;
  IntegerMatrix net_temp(n*n, 2);
  NumericVector dist_temp(n*n);
  LogicalVector dup_temp(n*n);
  IntegerVector flag_start(n);
  IntegerVector flag_end(n);
  for(i = 0; i < n; i++)
  {
    flag_start(i) = count + 1;
    for(j = 0; j < n; j ++)
    {
      if(i != j && std::abs(x(i) - x(j)) <= c && std::abs(y(i) - y(j)) <= c)
      {
        temp = (x(i) - x(j))*(x(i) - x(j)) + (y(i) - y(j))*(y(i) - y(j));
        if(temp <= c*c)
        {
          net_temp(count, 0) = i + 1;
          net_temp(count, 1) = j + 1;
          dist_temp(count) = sqrt(temp);
          if(i > j)
          {
            dup_temp(count) = TRUE;
          }
          else
          {
            dup_temp(count) = FALSE;
          }
          count++;
        }
      }
    }
    flag_end(i) = count;
  }
  IntegerMatrix net(count, 2);
  NumericVector dist(count);
  LogicalVector dup(count);
  for (i = 0; i < count; i++)
  {
    for (j = 0; j < 2; j++)
    {
      net(i, j) = net_temp(i, j);
    }
    dist(i) = dist_temp(i);
    dup(i) = dup_temp(i);
  }
  return Rcpp::List::create(Rcpp::Named("edge") = net, Rcpp::Named("distance") = dist, Rcpp::Named("duplicate") = dup, Rcpp::Named("flag_start") = flag_start, Rcpp::Named("flag_end") = flag_end);
}

// [[Rcpp::export]] 
// This function is to generate z 
arma::vec model(arma::vec z, IntegerMatrix E, NumericVector d, int id_start, int id_end, IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta, NumericVector omega, double lambda) {
  int Q = Theta.nrow();
  NumericVector prob_temp(Q);
  IntegerVector state(Q);
  int i, j, q,qq;
  double temp = 0;
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for(i = id_start - 1; i < id_end; i++)
  {
    if(flag_end(i) >= flag_start(i))
    {
      //Rcout<< i<< " "<<flag_start(i)<< "\n";
      for(q = 0; q < Q; q++)
      {
        prob_temp(q) = 0;
      }
      for(j = flag_start(i) - 1; j < flag_end(i); j++)
      {
        for(q = 0; q < Q; q++)
        {
          //if(z(E(j,1)-1) != q){
            //prob_temp(q) = prob_temp(q) - Theta(q, z(E(j, 1) - 1))*exp(-lambda*d(j));
            // to align with boost-mi, update the way calculate energy
            prob_temp(q) = prob_temp(q) + Theta(q, z(E(j, 1) - 1))*exp(-lambda*d(j));
            //}
          //prob_temp(q) = prob_temp(q) - Theta(q, z(E(j, 1) - 1)-1 )*exp(-lambda*d(j));
          //if(E(j,0)-1== i) {Rcout<< "oksgdu";}
          //prob_temp(q) = prob_temp(q) - Theta(q, z(E(j, 1) - 1) )*exp(-lambda*d(j));
        }
      }
      for(q = 0; q < Q; q++)
      {
        prob_temp(q) = exp(prob_temp(q) + omega(q));
        //prob_temp(q) = exp(prob_temp(q) - omega(q));
      }
      temp = 0;
      for (q = 0; q < Q; q++)
      {
        temp = temp + prob_temp(q);
      }
      for (q = 0; q < Q; q++)
      {
        prob_temp(q) = prob_temp(q)/temp;
      }
      // z(i) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      z(i) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) ;
    }
    else
    {
      for(q = 0; q < Q; q++)
      {
        prob_temp(q) = exp(omega(q));
        //prob_temp(q) = exp(-omega(q));
      }
      temp = 0;
      for (q = 0; q < Q; q++)
      {
        temp = temp + prob_temp(q);
      }
      for (q = 0; q < Q; q++)
      {
        prob_temp(q) = prob_temp(q)/temp;
      }
      //z(i) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      z(i) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) ;
    }
  }
  return z;
}

// [[Rcpp::export]]
double logenergy(arma::vec z, NumericVector zi, IntegerMatrix E, NumericVector d, LogicalVector dd, NumericMatrix Theta, NumericVector omega, double lambda) {
  double energy = 0;
  int m = E.nrow();
  int Q = zi.size();
  int j, q;
  for(j = 0; j < m; j++)
  {
    if(!dd(j))
    {
      // Add one condition
      // Condition is unnecessary
      // if(z(E(j, 0) - 1) != z(E(j, 1) - 1)){
        energy = energy + Theta(z(E(j, 0) - 1) , z(E(j, 1) - 1) )*exp(-lambda*d(j));
        //}
      // energy = energy + Theta(z(E(j, 0) - 1) - 1, z(E(j, 1) - 1) - 1)*exp(-lambda*d(j));
      
    }
  }
  for(q = 0; q < Q; q++)
  {
    //energy = energy + zi(q)*omega(q);
    energy = energy + zi(q)*omega(q);
  }
  return energy;
}

double logenergy1(arma::vec z,int qindex,int ii,IntegerMatrix E, NumericVector d, LogicalVector dd, NumericMatrix Theta, NumericVector omega, double lambda) {
  double energy = 0;
  int m = E.nrow();
  // int Q = zi.size();
  int j; // iterate over edges pair 
  for( j =0; j<m;j++)
  {
    if((E(j, 0) - 1) == ii){
      // if(z(E(j, 0) - 1)== qindex && z(E(j, 1) - 1)!=qindex ){
        //if(z(E(j, 0) - 1)==qindex ){
          energy = energy + Theta(qindex,z(E(j, 1)-1))*exp(-lambda*d(j));
          //energy = energy - Theta(0,1)*exp(-lambda*d(j));
          //  }
      }
    }
    energy = energy + omega(qindex);// without zi
    
    return energy;
  }
  
  arma::vec rowSums(IntegerMatrix z_store){
    int nrow = z_store.nrow();
    int ncol = z_store.ncol();
    NumericVector rowsums(nrow);
    for (int hh = 0; hh<nrow;hh++){
      for (int kk = 0; kk<ncol;kk++){
        rowsums(hh) = rowsums(hh)+ z_store(hh,kk);
      }
    }
    Rcout<<"calling";
    return rowsums;
  }
  
  arma::vec colMeans(IntegerMatrix z_store){
    int nrow = z_store.nrow();
    int ncol = z_store.ncol();
    NumericVector colsums(ncol);
    for (int hh = 0; hh < ncol;hh++){
      for (int kk = 0; kk < nrow;kk++){
        colsums(hh) = colsums(hh)+ z_store(kk,hh);
      }
      colsums(hh) = colsums(hh)/nrow;
    }
    
    return colsums;
  }
  
  
  
  
