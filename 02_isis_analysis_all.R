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

load(file = "data/analysis/full_d_all.RData")

model_code_all= "
data {
        int<lower = 1> n;           // total number of observations
        int<lower = 1> p;           // number of covariates in design matrix
        matrix[n, p] X;           // design matrix
        int<lower = 0> Y[n];            // vector of labels

        int<lower = 1> N_dist;          // number of small-areas
        int<lower = 1> N_gov;           // number of large-areas
  
        int<lower = 1> dist_id[n];          // small-area id
        int<lower = 1> gov_id[n];           // large-area id

        int<lower=0> N_dist_edges;          // Data for improper, efficient spatial prior:
        int<lower=1, upper=N_dist> node1[N_dist_edges];         // node1[i] adjacent to node2[i]
        int<lower=1, upper=N_dist> node2[N_dist_edges];         // and node1[i] < node2[i]

        real scaling_factor;         // scaling factor derived from the adjacency matrix
                           
        vector[N_gov] log_offset;         // offset
        vector[N_gov] theta;         // Pr(Y = 1 | r = 1, s = 1)
}

parameters {
            vector[p] beta;         // fixed effects coefficients
            vector[N_dist] phi;         // small-area random effect - main
            vector[N_dist] psi;         // small-area random effect - spatial
            vector[N_gov] eta;          // large-area random effect
        
            real<lower = 0,upper = 1> lambda;           // mixing weight
        
            real<lower = 0> sigma;            // sd of small-area random effect
            real<lower = 0> tau_eta;            // precision of large-area random effect
}

transformed parameters{
                        vector[N_dist] gamma;           // convoluted small-area effect
                        vector[n] mu;           // logit-scale propensity to be a recruit
                        vector[n] mu_star;  
                        real<lower = 0, upper = 1> rho[n];          // propensity to be a recruit
                        real<lower = 0, upper = 1> rho_star[n]; 
        
                        real<lower = 0> tau;          // precision of small-area random effect
                        real<lower = 0> sigma_eta;          // sd of large-area random effect
        
                        matrix<lower = 0, upper = 1>[2,N_gov] xi;           // probability of being labeled a recruit

                        sigma_eta = sqrt(1/tau_eta);            // large-area scale
                        tau = pow(sigma,-2);             // small-area precision


                        gamma =  sqrt(1-lambda) * phi + sqrt(lambda / scaling_factor) * psi ;
                        // variance of each component should be approximately equal to 1
              

                        mu = log_offset[gov_id] + X * beta + gamma[dist_id]*sigma + eta[gov_id]*sigma_eta ;
                        mu_star = X * beta + gamma[dist_id]*sigma + eta[gov_id]*sigma_eta ;
                        // linear function of the logit-scale propensity to be a recruit
   
                        for(i in 1:n){
                        rho[i] = inv_logit(mu[i]);          // propensity to be a recruit
                        rho_star[i] =inv_logit(mu_star[i]); 
                        }


                        for(j in 1:N_gov){
                                    xi[1,j] = theta[j]; 
                                    xi[2,j] = 0;
                        }
                        // defining the propensities of being labeled a recruit for each of the candidate models in the mixture
}

model {

        beta[1] ~ cauchy(0, 10); // prior on intercept
        for(i in 2:p){
                beta[i] ~ cauchy(0, 2.5);           // prior on fixed-effects, excluding the intercept
        }

        phi ~ normal(0,1);      // unstructured random effect on small-area
        target += -0.5 * dot_self(psi[node1] - psi[node2]);         // ICAR prior
        sum(psi) ~ normal(0, 0.01 * N_dist);            // soft sum-to-zero, equivalent to mean(psi) ~ normal(0,0.01)
        lambda ~ beta(0.5,0.5);         //  mixing weight prior
        sigma ~ normal(0,1);          // half-normal prior for small-area sd
  
        eta ~ normal(0,1);            // random effect on large area
        tau_eta ~ gamma(1,1);           // conjugate gamma prior on large-area precision

        for (i in 1:n){
                target += log_mix(rho[i],
                bernoulli_lpmf(Y[i] | xi[1,gov_id[i]] ),
                bernoulli_lpmf(Y[i] | xi[2,gov_id[i]] ) );          // labels distributed as mixture of bernoulli distributions 
        } 
}"

##################### ESTIMATE MODEL

begin.time = Sys.time()
full_fit_all <- stan(model_code = model_code_all, 
                     data = full_d_all, 
                     iter = 10000,
                     warmup = 9500,
                     thin =4,
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
                     control = list(max_treedepth = 10, adapt_delta = 0.9),
                     verbose = TRUE)
time.taken = Sys.time() - begin.time 

save(full_fit_all,file = 'data/analysis/full_fit_all.RData',compress = TRUE)