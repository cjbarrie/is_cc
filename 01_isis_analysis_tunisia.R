rm(list=ls())
gc()
options(scipen=999)
# library(rstan) may encounter install errors. Follow procedure below to install rstan with
# dependencies, as well as INLA (and graph/Rgraphviz if necessary).
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)

# options(timeout=1000)
# install.packages("INLA", repos=c(getOption("repos"), 
#                                  INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# #
library(R2jags)
library(rstan)
library(INLA)

#LOAD ANALYSIS DATA

load(file = "data/analysis/full_d_tunisia.RData")

model_code_tunisia= "

data {
                              int<lower = 1> n; // total number of observations
                              int<lower = 1> p; // number of covariates in design matrix

                               real log_offset; // offset
                                    real theta; // Pr(Y = 1 | r = 1, s = 1)

                                matrix[n, p] X; // design matrix
                           int<lower = 0> Y[n]; // vector of labels

                         int<lower = 1> N_dist; // number of small-areas
                          int<lower = 1> N_gov; // number of large-areas
  
                     int<lower = 1> dist_id[n]; // small-area id
                      int<lower = 1> gov_id[n]; // large-area id

                     int<lower=0> N_dist_edges; // Data for improper, efficient spatial prior:
int<lower=1, upper=N_dist> node1[N_dist_edges]; // node1[i] adjacent to node2[i]
int<lower=1, upper=N_dist> node2[N_dist_edges]; // and node1[i] < node2[i]
                           real scaling_factor; // scaling factor derived from the adjacency matrix
}

parameters {
                 vector[N_dist] phi; // small-area random effect - main
                 vector[N_dist] psi; // small-area random effect - spatial
                  vector[N_gov] eta; // gov random effect
   real<lower = 0,upper = 1> lambda; // mixing weight
                real<lower = 0> sigma; // sd of small-area random effect
            real<lower = 0> tau_eta; // precision of large-area random effect
                     vector[p] beta; // fixed effects coefficients
}

transformed parameters{
            real<lower = 0> sigma_eta; // sd of large-area random effect
            real<lower = 0> tau; // precision of small-area random effect

             vector[N_dist] gamma; // convoluted small-area effect
                     vector[n] mu; // logit-scale propensity to be a recruit - according to the non-probability sampling framework
                     vector[n] mu_star; // logit-scale propensity to be a recruit - in absence of the non-probability sampling framework

real<lower = 0, upper = 1> rho[n]; // propensity to be a recruit - according to the non-probability sampling framework
real<lower = 0, upper = 1> rho_star[n]; // propensity to be a recruit - in absence of the non-probability sampling framework

 real<lower = 0, upper = 1> xi[2]; // probability of being labeled a recruit


      sigma_eta = sqrt(1/tau_eta); // large-area scale
      tau = pow(sigma,-2); // small-area scale


// variance of each component should be approximately equal to 1
     gamma =  sqrt(1-lambda) * phi + sqrt(lambda / scaling_factor) * psi ; 
 // linear function of the logit-scale propensity to be a recruit
mu = log_offset + X * beta + gamma[dist_id]*sigma + eta[gov_id]*sigma_eta;
mu_star = X * beta + gamma[dist_id]*sigma + eta[gov_id]*sigma_eta;
   
for(i in 1:n){
       rho[i] = inv_logit(mu[i]); // propensity to be a recruit
       rho_star[i] = inv_logit(mu_star[i]); // propensity to be a recruit
}
                   xi[1] = theta; // defining the propensities of being a label for each of the candidate models in the mixture
                       xi[2] = 0;
}

model {

// // // Coefficients
        beta[1] ~ cauchy(0, 10);
for(i in 2:p){
       beta[i] ~ cauchy(0, 2.5); // prior on fixed-effects, excluding the intercept
}

    
// // // Random Effect on Small-Area
              phi ~ normal(0,1);
        
// // // ICAR priors
  target += -0.5 * dot_self(psi[node1] - psi[node2]);
// soft sum-to-zero constraint on psi,
// equivalent to mean(psi) ~ normal(0,0.01)
  sum(psi) ~ normal(0, 0.01 * N_dist);
  
  // // // Random Effect on Large-Area
  eta ~ normal(0,1);
  tau_eta ~ gamma(0.1,0.1);
// // //  Convolution Parameters
    lambda ~ beta(0.5,0.5); //dirichlet(alpha);
    sigma ~ normal(0,1);

// // // Contamination 
for (i in 1:n){
 target += log_mix(rho[i],
                   bernoulli_lpmf(Y[i] | xi[1] ),
                   bernoulli_lpmf(Y[i] | xi[2] ) ); // labels distributed as mixture of bernoulli distributions
}

}"

##################### ESTIMATE MODEL

full_fit_tunisia <- stan(model_code = model_code_tunisia, 
                         data = full_d_tunisia, 
                         iter = 10000,
                         warmup = 9500,
                         thin = 4,
                         pars = c('beta',
                                  'phi',
                                  'psi',
                                  'eta','tau_eta',
                                  'lambda',
                                  'tau',
                                  'gamma',
                                  'rho','rho_star'),
                         cores =4,
                         chains = 4,
                         control = list(max_treedepth = 10),
                         verbose = TRUE)

save(full_fit_tunisia,file = 'data/analysis/full_fit_tunisia.RData',compress = TRUE)