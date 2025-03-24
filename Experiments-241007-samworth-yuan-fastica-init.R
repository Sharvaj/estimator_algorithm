
library(cubature)
library(MASS)
library(pracma)
library(mvtnorm)

library(LogConcDEAD)
library(logcondens) 
library(mclust)
library(ica)

source('lcic.r')


# Define "sweep" function

run_comparisons_over_n <- function(all_n, d, num_repeats_full, savefilename) {

    num_n <- length(all_n)
    num_exp <- num_n * num_repeats_full

    # mu = rep(0,d)
    true_mean_vec <- rep(0,d) 
    scalings_vec = rep(1, d)


    iexp <- 0
    all_results = list()
    
    # Print
    cat("Num n values and experiment repeats: ", num_n, num_repeats_full)

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
            W <- fix_signs_fun(randortho(d, type="orthonormal"))
            SimData <- get_uniform_data(d, n, scalings_vec=scalings_vec, true_mean_vec=true_mean_vec, W=W)


            # Yuan and Samworth estimator with fastICA init
            t_start <- Sys.time()
            # ICA init
            my_output <- icafast(SimData$pre_data, nc=d, center=FALSE, fun='kur')
            qrmat <- qr(t(my_output$W))
            W_init <- t(qr.Q(qrmat))
            # Alternating algorithm
            algo_output <- yuan_samworth_alternating(SimData, W_init, max_iter=20)
            t_end <- Sys.time()
            
            time_taken <- difftime(t_end, t_start, units="secs")

            # Helper functions for computing sq hellinger errors

            hfun_yuan_samworth <- function(X_samps){
                density_ratio <- evaluate_yuan_samworth_estimator_vectorized(X_samps, algo_output$estimator)/uniform_pdf_vectorized(X_samps, scalings_vec=scalings_vec, true_mean_vec=true_mean_vec, W=W)
                return(0.5*(sqrt(density_ratio)-1)^2)
            }

            generate_samples_for_monte_carlo <- function(K_samps) get_uniform_data(d, K_samps, scalings_vec=scalings_vec, true_mean_vec=true_mean_vec, W=W)$pre_data

            K_samps <- 10000
            num_repeats <- 50        
            
            
            hellinger_error_estimate_statistics <- naive_monte_carlo_integrate_repeated(hfun_yuan_samworth, generate_samples_for_monte_carlo, K_samps, num_repeats)
    
            
            # Collect results
        
            exp_results_collect <- list("n"=n, "d"=d, "irep"=irep,
                                   "time_taken"=time_taken,
                                   "hell_err_mean"=hellinger_error_estimate_statistics$mean_val)

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

num_repeats_full <- 5


for (d in c(2,3,4)) {
    savefilename <- sprintf('Results/yuan_samworth_fastica_init_over_n_d=%s.RData', d)
    my_results <- run_comparisons_over_n(all_n, d, num_repeats_full, savefilename)
}
