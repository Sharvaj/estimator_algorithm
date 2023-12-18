### Import libraries

library(cubature)
library(MASS)
library(pracma)
library(mvtnorm)
library(LogConcDEAD)
library(logcondens) 
library(mclust)

###

fix_signs_fun <- function(my_mat) {
    # Make the signs of the fist column positive
    signs_first_column = ifelse(my_mat[,1] < 0, -1, 1)
    return(my_mat*signs_first_column)
}

center_data_fun <- function(data_mat) {
    # The rows should be the samples, and the columns should be dimensions/features
    mean_vec <- colMeans(data_mat)
    centered_data = t(t(data_mat) - mean_vec)
    return(list("centered_data"=centered_data, "mean_vec"=mean_vec))
}

generate_weighted_samples <- function(data_mat, theta_weights, num_resample) {
    n <- NROW(data_mat)
    sampled_row_indices <- sample(x=c(1:n), size=num_resample, replace=TRUE, prob=theta_weights)
    return(data_mat[sampled_row_indices,])
}

## Generates synthetic data to be worked with. This should be modified
## according to the desired application and data. Creates
## a matrix whose rows represent each data point.


get_heteroskedastic_gaussian_data <- function(d, n, true_mean_vec=rep(0,d), Sigma_max=d, eigensep=1) {
  # print(true_mean_vec)
  covariance_Z <- matrix(0, ncol = d, nrow = d)
  Sigma_min <- Sigma_max - (d-1)*eigensep
  diag(covariance_Z) <- seq(from=Sigma_max, to=Sigma_min, by=(-eigensep))
  print("Diagonal covariance: ")
  print(covariance_Z)
  W <- fix_signs_fun(randortho(d, type="orthonormal"))
  covariance_X <- t(W) %*% covariance_Z %*% W
  pre_data <- mvrnorm(n, mu=true_mean_vec, Sigma=covariance_X)
  my_output <- list("n"=n, "d"=d, "covariance_Z"=covariance_Z, "true_mean_vec"=true_mean_vec, "covariance_X"=covariance_X, "W"=W, "pre_data"=pre_data)
  return(my_output)
}



get_gaussian_mixture_data <- function(d, n, true_mean_vec_list, true_proportions) {
    num_clusters <- length(true_proportions)
    cluster_inds <- sample(x=c(1:num_clusters), size=n, replace=TRUE, prob=true_proportions)
    all_cluster_SimData <- list()
    pre_data_full <- matrix(0, nrow=n, ncol=d)

    for (kind in 1:num_clusters){
        where_this_cluster <- (cluster_inds == kind)
        num_samples_this_cluster <- sum(where_this_cluster)
        all_cluster_SimData[[kind]] <- get_heteroskedastic_gaussian_data(d, num_samples_this_cluster, true_mean_vec_list[[kind]])
        pre_data_full[where_this_cluster,] <- all_cluster_SimData[[kind]]$pre_data
    }

    return(list("n"=n, "d"=d, "pre_data"=pre_data_full, "all_cluster_SimData"=all_cluster_SimData, "cluster_inds"=cluster_inds))
}

get_axis_aligned_heteroskedastic_gamma_data<- function(d, n, Sigma_max=d+1, eigensep=1) {
    # W is the identity here! Z = X

    init_data <- matrix(0, nrow=n, ncol=d)
    true_mean_vec <- rep(0,d)
    covariance_Z <- matrix(0, ncol = d, nrow = d)
    Sigma_min <- Sigma_max - (d-1)*eigensep
    Sigma_all <- seq(from=Sigma_max, to=Sigma_min, by=(-eigensep))

    for (dind in 1:d) {
        init_data[,dind] <- rgamma(n, shape=Sigma_all[dind]) 
        true_mean_vec[dind] <- Sigma_all[dind]
        covariance_Z[dind, dind] <- Sigma_all[dind]
    }

    covariance_X <- covariance_Z
    W <- diag(d)
    my_output <- list("n"=n, "d"=d, "covariance_Z"=covariance_Z, "true_mean_vec"=true_mean_vec, "covariance_X"=covariance_X, "W"=W, "pre_data"=init_data)
    return(my_output)
}

get_heteroskedastic_gamma_data <- function(d,n) {
    
#   diag(covariance_Z) <- c(d:1)
#   print("Diagonal covariance: ")
#   print(covariance_Z)
    W <- fix_signs_fun(randortho(d, type="orthonormal"))
    
    init_data <- matrix(0, nrow=n, ncol=d)
    true_mean_vec <- rep(0,d)
    covariance_Z <- matrix(0, ncol = d, nrow = d)
    for (dind in 1:d){
        init_data[,dind] <- rgamma(n, shape=d+2-dind) # CHANGE to reverse
        true_mean_vec[dind] <- d+2-dind
        covariance_Z[dind, dind] <- d+2-dind
    }
    covariance_X <- t(W) %*% covariance_Z %*% W
    centered_init_data <- t(t(init_data) - true_mean_vec)
    
    pre_data <- centered_init_data %*% W
    pre_data <- t(t(pre_data) + true_mean_vec)

    my_output <- list("n"=n, "d"=d, "covariance_Z"=covariance_Z, "true_mean_vec"=true_mean_vec, "covariance_X"=covariance_X, "W"=W, "pre_data"=pre_data)
    return(my_output)
}

package_real_data <- function(real_data_mat) {
    n = NROW(real_data_mat)
    d = NCOL(real_data_mat)

    return(list("n"=n, "d"=d, "pre_data"=real_data_mat))
}

## Randomizes and splits the observations into two sets, with fraction r 
## (0 <= r <= 1) of them going into the first set.
## Takes as input the output of the data-generating functions

randomize_center_and_split <- function(SimData, r) {
  num_samples_1 <- floor(r * SimData$n)
  num_samples_2 <- SimData$n - num_samples_1
  randomized_centered_data <- center_data_fun(SimData$pre_data[sample(SimData$n),])
  data_1 <- randomized_centered_data$centered_data[1:num_samples_1,]
  data_2 <- randomized_centered_data$centered_data[(num_samples_1+1):SimData$n,]
  return(list("data_cov_estimation"=data_1, "data_marginal_estimation"=data_2, "mean_vec"=randomized_centered_data$mean_vec,
              "n_cov_estimation"=num_samples_1, "n_marginal_estimation"=num_samples_2))
}

## Generates the density estimate using the algorithm above, 
## based on the two sets of data from randomize_center_and_split

generate_estimator_with_logcondens <- function(SimData, r=0.9, plotting=FALSE) {
    
    SplitData <- randomize_center_and_split(SimData, r)
    
    cov_compute <- cov.wt(x=SplitData$data_cov_estimation, center=FALSE, method="ML")
    # weights_vec_renormalized <- weights_vec/(sum(weights_vec)*SplitData$n_cov_estimation)
    # weighted_cov_compute <- cov.wt(x=SplitData$data_cov_estimation, wt=weights_vec, center=center, method="ML")
    # weighted_emp_cov_mat <- cov_compute$cov
    # weighted_emp_mean <- cov_compute$center

    W_hat <- t(princomp(covmat=cov_compute$cov, fix_sign=TRUE)$loadings)
    unmixed_obs <- SplitData$data_marginal_estimation %*% t(W_hat)
    print("PCA done!")
    
    marginals <- list()
    
    for (i in 1:SimData$d) {
        print("Marginal: ")
        print(i)
        out1 <- logConDens(unmixed_obs[,i], smoothed = FALSE)
        if (plotting){
            plot(out1)
        }
        marginals[[i]] <- out1
    }
    return(list("W_hat"=W_hat, "marginals"=marginals, "mean_vec"=SplitData$mean_vec))
}

generate_estimator_with_logconcdead <- function(SimData, r=0.9, plotting=FALSE) {
    
    SplitData <- randomize_center_and_split(SimData, r)
    
    cov_compute <- cov.wt(x=SplitData$data_cov_estimation, center=FALSE, method="ML")
    # weights_vec_renormalized <- weights_vec/(sum(weights_vec)*SplitData$n_cov_estimation)
    # weighted_cov_compute <- cov.wt(x=SplitData$data_cov_estimation, wt=weights_vec, center=center, method="ML")
    # weighted_emp_cov_mat <- cov_compute$cov
    # weighted_emp_mean <- cov_compute$center

    W_hat <- t(princomp(covmat=cov_compute$cov, fix_sign=TRUE)$loadings)
    unmixed_obs <- SplitData$data_marginal_estimation %*% t(W_hat)
    print("PCA done!")
    
    marginals <- list()
    
    for (i in 1:SimData$d) {
        print("Marginal: ")
        print(i)
        out1 <- mlelcd(unmixed_obs[,i], verbose=0)
        if (plotting){
            plot(out1)
        }
        marginals[[i]] <- out1
    }
    return(list("W_hat"=W_hat, "marginals"=marginals, "mean_vec"=SplitData$mean_vec))
}

generate_weighted_estimator_with_logcondens <- function(SimData, theta_weights, resample_factors=c(1,1), plotting=FALSE) {
    
    M_resample <- resample_factors[1]*SimData$n
    N_resample <- resample_factors[2]*SimData$n
    
    resampled_data <- generate_weighted_samples(SimData$pre_data, theta_weights, M_resample+N_resample)
    
    centered_resampled_data <- center_data_fun(resampled_data)
    
    centered_data_for_cov <- centered_resampled_data$centered_data[1:M_resample,]
    centered_data_for_dens <- centered_resampled_data$centered_data[(M_resample+1):(M_resample+N_resample),]
    

    cov_compute <- cov.wt(x=centered_data_for_cov, center=FALSE, method="ML")
    # weights_vec_renormalized <- weights_vec/(sum(weights_vec)*SplitData$n_cov_estimation)
    # weighted_cov_compute <- cov.wt(x=SplitData$data_cov_estimation, wt=weights_vec, center=center, method="ML")
    # weighted_emp_cov_mat <- cov_compute$cov
    # weighted_emp_mean <- cov_compute$center

    W_hat <- t(princomp(covmat=cov_compute$cov, fix_sign=TRUE)$loadings)
    unmixed_obs <- centered_data_for_dens %*% t(W_hat)
    print("PCA done!")
    
    marginals <- list()
    
    for (i in 1:SimData$d) {
        print("Marginal: ")
        print(i)
        out1 <- logConDens(unmixed_obs[,i], smoothed = FALSE)
        if (plotting){
            plot(out1)
        }
        marginals[[i]] <- out1
    }
    return(list("W_hat"=W_hat, "marginals"=marginals, "mean_vec"=centered_resampled_data$mean_vec))
}



generate_fulldim_log_concave_estimator <- function(SimData){

    my_fulldim_estimator <- mlelcd(SimData$pre_data, verbose=0)
    optflag <- 0
    if (my_fulldim_estimator$NumberOfEvaluations[1]<0){
        optflag <- 1
    }
    return(list("my_fulldim_estimator"=my_fulldim_estimator, "optflag"=optflag))
}


evaluate_logcondens_estimator_vectorized <- function(X_mat, my_estimator){
    # X_mat must be organized such that each row is a point in R^d
    d <- NCOL(X_mat)
    n_eval <- NROW(X_mat)
    
    X_mat_centered = t(t(X_mat) - my_estimator$mean_vec) # subtract the mean computed using the training data!

    Z_mat_hat = X_mat_centered %*% t(my_estimator$W_hat)
    
    result = rep(1.0, n_eval)
    for (dind in 1:d){
        result = result * evaluateLogConDens(Z_mat_hat[,dind], my_estimator$marginals[[dind]], which=2)[,3]
    }
    return(result)
}

evaluate_logconcdead_estimator_vectorized <- function(X_mat, my_estimator){
    # X_mat must be organized such that each row is a point in R^d
    d <- NCOL(X_mat)
    n_eval <- NROW(X_mat)
    
    X_mat_centered = t(t(X_mat) - my_estimator$mean_vec) # subtract the mean computed using the training data!

    Z_mat_hat = X_mat_centered %*% t(my_estimator$W_hat)
    
    result = rep(1.0, n_eval)
    for (dind in 1:d){
        result = result * dlcd(matrix(Z_mat_hat[,dind], ncol = 1), my_estimator$marginals[[dind]])
    }
    return(result)
}

evaluate_mixture_density_vectorized <- function(cluster_densities, pi_vec, X_mat) {
    n_eval <- NROW(X_mat)
    num_clusters <- length(pi_vec)
    evaluated_clusters <- matrix(0, nrow=n_eval, ncol=num_clusters)

    for (kind in 1:num_clusters) {
        evaluated_clusters[,kind] <- evaluate_logcondens_estimator_vectorized(X_mat, cluster_densities[[kind]])
    }
    return(evaluated_clusters %*% pi_vec)
}

evaluate_fulldim_estimator_vectorized <- function(X_mat, my_fulldim_estimator){
    return(dlcd(X_mat, my_fulldim_estimator))
}

## Evaluates the pdf of the Gaussian from get_gaussian_test_data at the point x.
gaussian_test_data_pdf <- function(x, d, covariance_Z, W) {
  Kmat_Z <- matrix(0, ncol = d, nrow = d)
  diag(Kmat_Z) <- 1/diag(covariance_Z)
  Kmat_X = t(W) %*% Kmat_Z %*% W
  return(prod(diag((2*pi)*covariance_Z))^(-1/2) * exp(-(1/2)*(t(x)%*% Kmat_X %*% x)))
}


visualize_theta_weights <- function(theta_weights, data_mat) {
    return(0)
}



initialize_EM <- function(SimData, num_clusters) {

    highclust <- hc( modelName="VVV", data=SimData$pre_data)
    class <- c( hclass( highclust, num_clusters ) )

    pi_vec = rep(0,num_clusters)
    y <- matrix( 0, nrow=SimData$n, ncol=num_clusters )
    theta_mat <- matrix(0, nrow=SimData$n, ncol=num_clusters)

    for( i in 1:num_clusters ) {
      pi_vec[ i ] <- sum( class==i ) / SimData$n
      ss <- SimData$pre_data[ class==i, ]
      y[ , i ] <- dmvnorm( SimData$pre_data, mean=colMeans(ss), sigma=var( ss ))
      theta_mat[,i] <- pi_vec[i] * y[ , i ] # normalized below!
    }
    theta_mat <- theta_mat / rowSums(theta_mat) # are row sums zero?
    theta_mat[is.na(theta_mat)] <- 0
    
    return(theta_mat)
}


EM_with_lcic <- function(SimData, theta_mat_init, num_clusters, resample_factors=c(1,1), num_iter=10) {
    # theta_mat_init has n rows and K columns theta_ik = f(k|Xi)

    cluster_densities = list()
    theta_mat <- theta_mat_init

    all_likelihoods <- rep(NA, num_iter)

    
    for (itind in 1:num_iter){
        print("Iter: ")
        print(itind)
        # M-step (update densities and proportions) (typically this follows the E-step?)
        for (kind in 1:num_clusters){
            theta_weights <- theta_mat[,kind]
            theta_weights <- theta_weights/sum(theta_weights)
            cluster_densities[[kind]] <- generate_weighted_estimator_with_logcondens(SimData, theta_weights, resample_factors=resample_factors)
        }
        pi_vec <- colSums(theta_mat)
        # print('Pi vec')
        # print(pi_vec)
        pi_vec <- pi_vec/sum(pi_vec)
        # print(pi_vec)

        # E-step (update theta_mat)
        for (kind in 1:num_clusters) {
            theta_mat[,kind] <- pi_vec[kind] * evaluate_logcondens_estimator_vectorized(SimData$pre_data, cluster_densities[[kind]])
        }
        # print("Theta mat")
        # print(theta_mat)
        theta_mat <- theta_mat / rowSums(theta_mat) # are row sums zero?
        theta_mat[is.na(theta_mat)] <- 0
        # print("Row sums:")
        # print(rowSums(theta_mat))
        # print("Theta mat")
        # print(theta_mat)

        # Collect densities and evaluate likelihood
        my_evaluations <- evaluate_mixture_density_vectorized(cluster_densities, pi_vec, SimData$pre_data)
        all_likelihoods[itind] <- (1/SimData$n) * sum(log(my_evaluations))
        print("Log-likelihood: ")
        print(all_likelihoods[itind])
    }
    return(list("cluster_densities"=cluster_densities, "pi_vec"=pi_vec, "theta_mat"=theta_mat, "all_likelihoods"=all_likelihoods))
}



heteroskedastic_gaussian_pdf_vectorized <- function(X_mat, SimData){
    # X_mat must be organized such that each row is a point in R^d
    my_mean <- SimData$true_mean_vec
    return(mvtnorm::dmvnorm(X_mat, my_mean, SimData$covariance_X))
}

axis_aligned_heteroskedastic_gamma_pdf_vectorized <- function(X_mat, SimData) {
    n_eval <- NROW(X_mat)
    density_data <- rep(1,n_eval)
    for (dind in 1:d){
        density_data <- density_data * dgamma(X_mat[,dind], shape=SimData$covariance_Z[dind,dind]) 
    }
    return(density_data)
}

naive_monte_carlo_integrate <- function(F_fun, random_samples_generate, K_samps){
    # F_fun shoud take rows of a matrix as input points
    # Estimates int Fdp
    
    X_samps <- random_samples_generate(K_samps)
    F_samps <- F_fun(X_samps)
    return(mean(F_samps))
}

naive_monte_carlo_integrate_with_convergence <- function(F_fun, random_samples_generate, K_samps, step=100, ylim=c(0,0.3)){
    K_samp_range <- seq(step, K_samps, by=step)
    num_K_samps <- length(K_samp_range)
    all_F_means <- rep(NULL, num_K_samps)
    
    for (kind in 1:num_K_samps) {
        all_F_means[kind] <- naive_monte_carlo_integrate(F_fun, random_samples_generate, K_samp_range[kind])
    }
    plot(x=K_samp_range, y=all_F_means, type="b", ylim=ylim)
}

naive_monte_carlo_integrate_repeated <- function(F_fun, random_samples_generate, K_samps, num_repeats=100){
    all_F_estimates <- rep(NULL, num_repeats)
    for (rind in 1:num_repeats){
        all_F_estimates[rind] <- naive_monte_carlo_integrate(F_fun, random_samples_generate, K_samps)
    }
    return(list("mean_val"=mean(all_F_estimates), "median_val"=median(all_F_estimates), 
                "sd_val"=sd(all_F_estimates), "all_vals"=all_F_estimates))
}

# hellinger_fun_for_integration <- function(X_samps, estimate_vectorized_eval_fun, my_estimator){

#     density_ratio <- estimate_vectorized_eval_fun(X_samps, my_estimator)/test_pdf_vectorized(X_samps, SimData)
#     return(0.5*(sqrt(density_ratio)-1)^2)

# }

# generate_heteroskedastic_gaussian_samples_for_monte_carlo <- function(K_samps, SimData){
#     return(mvrnorm(K_samps, mu=rep_len(0,d), Sigma=SimData$covariance_X))
# }

visualize_independent_directions <- function(SimData, my_estimator, xlim=c(-1,1), ylim=c(-1,1)){
    d = SimData$d
    W = SimData$W
    W_hat = my_estimator$W_hat
    plot(NULL, xlim=xlim, ylim=ylim, xlab=expression(x[1]), ylab=expression(x[2]), cex.lab=1.5, cex.axis=1.5)
    arrows(rep(0,d), rep(0,d), W[,1], W[,2], col="blue")
    arrows(rep(0,d), rep(0,d), W_hat[,1], W_hat[,2], col="red")
    return(W_hat %*% t(W))
}

direction_inner_products <- function(SimData, my_estimator){
    W = SimData$W
    W_hat = my_estimator$W_hat
    return(W_hat %*% t(W))
}



################
# TESTS AND PLOTS


tests_for_split_ratio_r_gaussian <- function(d, n, all_split_r_vals, Sigma_max, eigensep, num_repeats_full, savefilename) {

    # print('Here')
    mu = rep(0,d)

    num_r_vals <- length(all_split_r_vals)
    num_exp <- num_r_vals*num_repeats_full

    iexp <- 0
    all_results = list()

    # Print
    print("Num r values and experiment repeats:")
    print(num_r_vals)
    print(num_repeats_full)
    flush.console()
    Sys.sleep(0.2)

    for (rind in 1:num_r_vals) {
        for (irep in 1:num_repeats_full) {
            print('r index and repeat index:')
            print(rind)
            print(irep)
            flush.console()
            Sys.sleep(0.2)
            iexp <- iexp + 1

            # Simulate data
            SimData <- get_heteroskedastic_gaussian_data(d, n, true_mean_vec=mu, Sigma_max=Sigma_max, eigensep=eigensep)
            
            # Compute estimate
            my_estimator_logcondens <- generate_estimator_with_logcondens(SimData, r=all_split_r_vals[rind], plotting=FALSE)

            ip_mat <- direction_inner_products(SimData, my_estimator_logcondens)

            # Assess error
            hfun_logcondens <- function(X_samps){
                density_ratio <- evaluate_logcondens_estimator_vectorized(X_samps, my_estimator_logcondens)/heteroskedastic_gaussian_pdf_vectorized(X_samps, SimData)
                return(0.5*(sqrt(density_ratio)-1)^2)
            }
            
            generate_heteroskedastic_gaussian_samples_for_monte_carlo <- function(K_samps){
                return(mvrnorm(K_samps, mu=mu, Sigma=SimData$covariance_X))
            }

            K_samps <- 10000
            num_repeats_mc <- 50
            
            hellinger_error_estimate_statistics_logcondens <- naive_monte_carlo_integrate_repeated(hfun_logcondens, 
                                                                generate_heteroskedastic_gaussian_samples_for_monte_carlo, K_samps, num_repeats_mc)

            exp_results_collect <- list("n"=n, "d"=d, "split_r"=all_split_r_vals[rind], "irep"=irep,
            "ip_mat"=ip_mat, "hell_err_mean_logcondens"=hellinger_error_estimate_statistics_logcondens$mean_val)

            all_results[[iexp]] <- exp_results_collect
            print("Done experiment")
            print(iexp)
            flush.console()
            Sys.sleep(0.2)
        }
    }

    # Save
    saveRDS(all_results, file=savefilename)

    return(all_results)
}





