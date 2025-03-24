

library(cubature)
library(MASS)
library(pracma)
library(mvtnorm)

library(LogConcDEAD)
library(logcondens) 
library(mclust)

source('lcic.r')


# Define "sweep" function

run_splitting_comparisons_over_n <- function(all_n, d, num_repeats_full, Sigma_max, eigensep, split_r, savefilename) {

    num_n <- length(all_n)
    num_exp <- num_n * num_repeats_full

    mu = rep(0,d)

    iexp <- 0
    all_results = list()
    
    # Print
    print("Num n values and experiment repeats:")
    print(num_n)
    print(num_repeats_full)
    flush.console()
    Sys.sleep(0.2)

    for (nind in 1:num_n) {
        for (irep in 1:num_repeats_full) {

            print('n index and repeat index:')
            print(nind)
            print(irep)
            flush.console()
            Sys.sleep(0.2)
            iexp <- iexp + 1
            
            n <- all_n[nind]

            # Simulate data
            SimData <- get_heteroskedastic_gaussian_data(d, n, true_mean_vec=mu, Sigma_max=Sigma_max, eigensep=eigensep)

            # Proposed estimator via logcondens: WITH sample splitting
            t_start <- Sys.time()
            my_estimator_logcondens <- generate_estimator_with_logcondens(SimData, r=split_r, plotting=FALSE)
            t_end <- Sys.time()
            time_taken_logcondens <- difftime(t_end, t_start, units="secs")

            # Proposed estimator via logcondens: WITHOUT sample splitting
            t_start <- Sys.time()
            my_estimator_no_split <- generate_estimator_with_logcondens(SimData, r=-1, plotting=FALSE)
            t_end <- Sys.time()
            time_taken_no_split <- difftime(t_end, t_start, units="secs")

            # Helper functions for computing sq hellinger errors
        
            hfun_logcondens <- function(X_samps){
                density_ratio <- evaluate_logcondens_estimator_vectorized(X_samps, my_estimator_logcondens)/heteroskedastic_gaussian_pdf_vectorized(X_samps, SimData)
                return(0.5*(sqrt(density_ratio)-1)^2)
            }
            
            
            hfun_no_split <- function(X_samps){
                density_ratio <- evaluate_logcondens_estimator_vectorized(X_samps, my_estimator_no_split)/heteroskedastic_gaussian_pdf_vectorized(X_samps, SimData)
                return(0.5*(sqrt(density_ratio)-1)^2)
            }
            
            generate_heteroskedastic_gaussian_samples_for_monte_carlo <- function(K_samps){
                return(mvrnorm(K_samps, mu=mu, Sigma=SimData$covariance_X))
            }

            # Error of proposed
            K_samps <- 10000
            # K_samps <- 1000
            num_repeats <- 50
            
            hellinger_error_estimate_statistics_logcondens <- naive_monte_carlo_integrate_repeated(hfun_logcondens, 
                                                                generate_heteroskedastic_gaussian_samples_for_monte_carlo, K_samps, num_repeats)
    
            
            hellinger_error_estimate_statistics_no_split <- naive_monte_carlo_integrate_repeated(hfun_no_split, 
                                                            generate_heteroskedastic_gaussian_samples_for_monte_carlo, K_samps, num_repeats)



            # Collect results
        
            exp_results_collect <- list("n"=n, "d"=d, "irep"=irep,
                                   "time_taken_logcondens"=time_taken_logcondens, "time_taken_no_split"=time_taken_no_split,
                                   "hell_err_mean_logcondens"=hellinger_error_estimate_statistics_logcondens$mean_val,
                                   "hell_err_mean_no_split"=hellinger_error_estimate_statistics_no_split$mean_val)

            all_results[[iexp]] <- exp_results_collect
        
            print("Done experiment num: ")
            print(iexp)
            print("################################################")
            flush.console()
            Sys.sleep(0.2)
        }
    }
    
    # Save
    saveRDS(all_results, file=savefilename)

    return(all_results)
}


all_n <- c(100, 500, 1000, 2000, 3000)
Sigma_max <- 15
eigensep <- 1
num_repeats_full <- 5
split_r = 0.4


for (d in c(5,10,15)) {
    savefilename <- sprintf('Results/no_split_comparisons_over_n_d=%s.RData', d)
    my_results <- run_splitting_comparisons_over_n(all_n, d, num_repeats_full, Sigma_max, eigensep, split_r, savefilename)
}
