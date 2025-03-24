### Import libraries

library(cubature)
library(MASS)
library(pracma)
library(mvtnorm)
library(LogConcDEAD)
library(logcondens) 
library(mclust)
library(ica)

### Basic utilities

fix_signs_fun <- function(my_mat) {
    # Make the signs of the first column positive
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

### Generate or package data

## The functions below generate synthetic data to be worked with. These should be used
## according to the desired application and data. These are based on matrices
## whose rows represent each data point.

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

get_uniform_data <- function(d, n, scalings_vec, true_mean_vec, W) {

    # W <- fix_signs_fun(randortho(d, type="orthonormal"))
    
    init_data <- matrix(0, nrow=n, ncol=d)

    # covariance_Z <- matrix(0, ncol = d, nrow = d)

    for (dind in 1:d){
        init_data[,dind] <- true_mean_vec[dind] + runif(n, min=-scalings_vec[dind], max=scalings_vec[dind])
        # covariance_Z[dind, dind] <- d+2-dind
    }
    # covariance_X <- t(W) %*% covariance_Z %*% W
    centered_init_data <- t(t(init_data) - true_mean_vec)
    
    pre_data <- centered_init_data %*% W
    pre_data <- t(t(pre_data) + true_mean_vec)

    my_output <- list("n"=n, "d"=d, "scalings_Z"=scalings_vec, "true_mean_vec"=true_mean_vec, "W"=W, "pre_data"=pre_data)
    return(my_output)

}

### Functions for contructing estimators

## Randomizes and splits the observations into two sets, with fraction r 
## (0 <= r <= 1) of them going into the first set.
## Takes as input the output of the data-generating functions

randomize_center_and_split <- function(SimData, r) {
    if (r == -1) {
        print("NOTE: r is -1. Will not split samples!")
        num_samples_1 <- SimData$n
        num_samples_2 <- SimData$n
        randomized_centered_data <- center_data_fun(SimData$pre_data[sample(SimData$n),])
        data_1 <- randomized_centered_data$centered_data
        data_2 <- randomized_centered_data$centered_data
    } else {
        num_samples_1 <- floor(r * SimData$n)
        num_samples_2 <- SimData$n - num_samples_1
        randomized_centered_data <- center_data_fun(SimData$pre_data[sample(SimData$n),])
        data_1 <- randomized_centered_data$centered_data[1:num_samples_1,]
        data_2 <- randomized_centered_data$centered_data[(num_samples_1+1):SimData$n,]
    }
    
    return(list("data_cov_estimation"=data_1, "data_marginal_estimation"=data_2, "mean_vec"=randomized_centered_data$mean_vec,
              "n_cov_estimation"=num_samples_1, "n_marginal_estimation"=num_samples_2))
}

## Generates the density estimate using the algorithm above, 
## based on the two sets of data from randomize_center_and_split

generate_estimator_with_logcondens <- function(SimData, r=0.9, use_ICA=FALSE, plotting=FALSE) {
    
    SplitData <- randomize_center_and_split(SimData, r)
    
    if (use_ICA) {
        my_output <- icafast(SplitData$data_cov_estimation, nc=d, center=FALSE, fun='kur')
        qrmat <- qr(t(my_output$W))
        W_hat <- t(qr.Q(qrmat))
        print("ICA done!")
    } else {
        cov_compute <- cov.wt(x=SplitData$data_cov_estimation, center=FALSE, method="ML")
        # weights_vec_renormalized <- weights_vec/(sum(weights_vec)*SplitData$n_cov_estimation)
        # weighted_cov_compute <- cov.wt(x=SplitData$data_cov_estimation, wt=weights_vec, center=center, method="ML")
        # weighted_emp_cov_mat <- cov_compute$cov
        # weighted_emp_mean <- cov_compute$center
        W_hat <- t(princomp(covmat=cov_compute$cov, fix_sign=TRUE)$loadings)
        print("PCA done!")
    }

    unmixed_obs <- SplitData$data_marginal_estimation %*% t(W_hat)
    
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

fit_marginals_given_W <- function(SimData, W_given, plotting=FALSE) {
    
    # CenteredData <- center_data_fun(SimData$pre_data[sample(SimData$n),])
    unmixed_obs <- SimData$pre_data %*% t(W_given)
    
    marginals <- list()
    sample_sorting_indices <- list()
    
    for (i in 1:SimData$d) {
        print("Marginal: ")
        print(i)
        sorted_unmixed_obs = sort(unmixed_obs[,i], index.return=TRUE)
        out1 <- logConDens(sorted_unmixed_obs$x, smoothed = FALSE)
        if (plotting){
            plot(out1)
        }
        marginals[[i]] <- out1
        sample_sorting_indices[[i]] <- sorted_unmixed_obs$ix
    }
    return(list("W_given"=W_given, "marginals"=marginals, "sample_sorting_indices"=sample_sorting_indices))
}




generate_estimator_with_logconcdead <- function(SimData, r=0.9, use_ICA=FALSE, plotting=FALSE) {
    
    SplitData <- randomize_center_and_split(SimData, r)
    
    if (use_ICA) {
        my_output <- icafast(SplitData$data_cov_estimation, nc=d, center=FALSE, fun='kur')
        qrmat <- qr(t(my_output$W))
        W_hat <- t(qr.Q(qrmat))
        print("ICA done!")
    } else {
        cov_compute <- cov.wt(x=SplitData$data_cov_estimation, center=FALSE, method="ML")
        # weights_vec_renormalized <- weights_vec/(sum(weights_vec)*SplitData$n_cov_estimation)
        # weighted_cov_compute <- cov.wt(x=SplitData$data_cov_estimation, wt=weights_vec, center=center, method="ML")
        # weighted_emp_cov_mat <- cov_compute$cov
        # weighted_emp_mean <- cov_compute$center
        W_hat <- t(princomp(covmat=cov_compute$cov, fix_sign=TRUE)$loadings)
        print("PCA done!")
    }

    unmixed_obs <- SplitData$data_marginal_estimation %*% t(W_hat)
    
    
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
    print("NA check:")
    print(any(is.na(unmixed_obs)))

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
    # Meant for logcondens 
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
    print("Not implemented yet!")
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


EM_with_lcic <- function(SimData, theta_mat_init, num_clusters, resample_factors=c(1,1), num_iter=10, theps=0) {
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
            try(cluster_densities[[kind]] <- generate_weighted_estimator_with_logcondens(SimData, theta_weights, resample_factors=resample_factors), silent=FALSE)
            
            # success_flag <- 1
            # redo_counter <- 1
            # while ((success_flag == 1) & (redo_counter < 10)) {

            #     success_flag <- tryCatch(
            #         {
            #             cluster_densities[[kind]] <<- generate_weighted_estimator_with_logcondens(SimData, theta_weights, resample_factors=resample_factors)
            #             # success return flag
            #             0
            #         },
            #         error = function(e) {
            #             print('Redo due to NA error')
            #             return(1)
            #         },
            #         finally = {
            #             print("Finally: Retry if needed")
            #         }

            #     )
            #     flush.console()
            #     Sys.sleep(0.2)
            #     redo_counter <- redo_counter + 1
            #     print(redo_counter)
            # }
            
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
        
        theta_mat <- pmax(theta_mat, theps)
        print("Theta min before = ")
        print(min(theta_mat))

        theta_mat <- theta_mat / rowSums(theta_mat)
        print("Theta min after = ")
        print(min(theta_mat))
        # print("Row sums:")
        # print(rowSums(theta_mat))
        # print("Theta mat")
        # print(theta_mat)

        # Collect densities and evaluate likelihood
        my_evaluations <- evaluate_mixture_density_vectorized(cluster_densities, pi_vec, SimData$pre_data)
        all_likelihoods[itind] <- (1/SimData$n) * sum(log(my_evaluations))
        print("Log-likelihood: ")
        print(all_likelihoods[itind])
        
        flush.console()
        Sys.sleep(0.2)
    }
    return(list("cluster_densities"=cluster_densities, "pi_vec"=pi_vec, "theta_mat"=theta_mat, "all_likelihoods"=all_likelihoods))
}

EM_with_fulldim <- function(SimData, num_clusters, max_iter=10, verbose=-1) {

    EM_results <- EMmixlcd(SimData$pre_data, k=num_clusters, max.iter=max_iter, verbose=verbose)
    
    n <- NROW(SimData$pre_data)

    pi_vec = EM_results$props
    evaluated_densities <- exp(EM_results$logf)
    
    theta_mat <- matrix(0, nrow=n, ncol=num_clusters)

    for (kind in 1:num_clusters) {
        theta_mat[,kind] <- pi_vec[kind] * evaluated_densities[,kind]
    }

    theta_mat <- theta_mat / rowSums(theta_mat)
    theta_mat[is.na(theta_mat)] <- 0

    return(list("cluster_densities"=NA, "pi_vec"=pi_vec, "theta_mat"=theta_mat, "all_likelihoods"=EM_results$lcdloglik))
}



heteroskedastic_gaussian_pdf_vectorized <- function(X_mat, SimData){
    # X_mat must be organized such that each row is a point in R^d
    my_mean <- SimData$true_mean_vec
    return(mvtnorm::dmvnorm(X_mat, my_mean, SimData$covariance_X))
}

axis_aligned_heteroskedastic_gamma_pdf_vectorized <- function(X_mat, SimData) {
    n_eval <- NROW(X_mat)
    d <- NCOL(X_mat)
    density_data <- rep(1,n_eval)
    for (dind in 1:d){
        density_data <- density_data * dgamma(X_mat[,dind], shape=SimData$covariance_Z[dind,dind]) 
    }
    return(density_data)
}

uniform_pdf_vectorized <- function(X_mat, scalings_vec, true_mean_vec, W) {

    n_eval <- NROW(X_mat)
    d <- NCOL(X_mat)

    X_mat_centered = t(t(X_mat) - true_mean_vec) # subtract the mean computed using the training data!

    Z_mat = X_mat_centered %*% t(W)
    
    result = rep(1.0, n_eval)
    for (dind in 1:d){
        result <- result * dunif(Z_mat[,dind], min=-scalings_vec[dind], max=scalings_vec[dind])
    }
    return(result)
}

### Functions for evaluating Monte Carlo integrals for Hellinger error computations

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

visualize_independent_directions_from_W <- function(W_true, W_estimate, xlim=c(-1,1), ylim=c(-1,1)){
    
    plot(NULL, xlim=xlim, ylim=ylim, xlab=expression(x[1]), ylab=expression(x[2]), cex.lab=1.5, cex.axis=1.5)
    arrows(rep(0,d), rep(0,d), W_true[,1], W_true[,2], col="blue")
    arrows(rep(0,d), rep(0,d), W_estimate[,1], W_estimate[,2], col="red")
    return(W_estimate %*% t(W_true))
}

direction_inner_products <- function(SimData, my_estimator){
    W = SimData$W
    W_hat = my_estimator$W_hat
    return(W_hat %*% t(W))
}


map_to_data_space <- function(scores_mat, pca_res) {
    return(t((pca_res$loadings %*% t(scores_mat)) + pca_res$center))
}



################
# TESTS


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


# Yuan and Samworth

#  See the function "fit_marginals_given_W" near line 200


compute_directional_derivative <- function(Y_tst, SimData, W_current, estimator_current) {

    # Looks correct

    n <- SimData$n
    d <- SimData$d

    K_inds = matrix(0, nrow=n, ncol=d) # K_ij
    slopes_b = matrix(0, nrow=n, ncol=d)

    WY <- W_current %*% Y_tst
    grad_g_val <- 0

    for (dind in 1:d) {
        density <- estimator_current$marginals[[dind]]
        sample_sorting_indices <- estimator_current$sample_sorting_indices[[dind]]
        sorted_samples <- SimData$pre_data[sample_sorting_indices, ]
        mask <- (density$IsKnot == 1)
        cumsum_knots <- cumsum(density$IsKnot)
        slopes_knotwise <- diff(density$phi[mask])/diff(density$x[mask])
        
        num_knots <- sum(mask)
        knot_indices <- which(mask)
        c_vec <- WY[dind,]
        transformed_samples <- sorted_samples %*% c_vec
        selection_indicators <- transformed_samples[mask] 

        K_inds[,dind] = cumsum_knots # handles most points correctly
        K_inds[knot_indices, dind] <- ifelse(selection_indicators<0, c(1:num_knots)-1, c(1:num_knots)) # corrects the knot points

        K_inds[1, dind] <- 1
        K_inds[n, dind] <- num_knots - 1 # corrects the first and the last knot points

        slopes_b[,dind] <- slopes_knotwise[K_inds[,dind]]

        grad_g_val <- grad_g_val + sum(transformed_samples*slopes_b[,dind])/n
    }
    return(grad_g_val)
}

compute_normalized_log_likelihood_g <- function(SimData, W_candidate, estimator_current) {

    # Looks correct, but check in a notebook! Test against alternative comparisons?
    # Compare with the likelihood value computed by the logConDens function?

    g_val <- 0

    Z_mat_candidate <- SimData$pre_data %*% t(W_candidate)

    for (dind in 1:SimData$d) {
        
        g_val <- g_val + (1/SimData$n)*sum(evaluateLogConDens(Z_mat_candidate[,dind], estimator_current$marginals[[dind]], which=1)[,2])
    }
    return(g_val)
}

update_W_given_marginals <- function(SimData, W_current, estimator_current, alp=0.3, gmm=0.5) {

    # Find steepest descent along the basis of the tangent space
    Y_tst <- matrix(0, nrow=SimData$d, ncol=SimData$d)
    grad_g_vals <- matrix(0, nrow=SimData$d, ncol=SimData$d)

    for (s in 1:(SimData$d - 1)){
        for (r in (s+1):SimData$d) {
            Y_tst[r,s] <- 1/sqrt(2)
            Y_tst[s,r] <- -1/sqrt(2)
            grad_g_vals[r,s] <- compute_directional_derivative(Y_tst, SimData, W_current, estimator_current)
            grad_g_vals[s,r] <- compute_directional_derivative(-Y_tst, SimData, W_current, estimator_current)
            Y_tst[r,s] <- 0
            Y_tst[s,r] <- 0
        }
    }
    grad_g_max <- max(grad_g_vals)
    max_inds <- which(grad_g_vals == grad_g_max, arr.ind = TRUE)
    r_max <- max_inds[1,1]
    s_max <- max_inds[1,2]

    # Construct Y_max
    Y_max <- matrix(0, nrow=SimData$d, ncol=SimData$d)
    Y_max[r_max,s_max] <- 1/sqrt(2)
    Y_max[s_max, r_max] <- -1/sqrt(2)

    # Line search
    ep <- 1.0
    W_candidate_update <- W_current %*% expm(ep*Y_max)
    g_val_current <- compute_normalized_log_likelihood_g(SimData, W_current, estimator_current)
    g_val_update <- compute_normalized_log_likelihood_g(SimData, W_candidate_update, estimator_current)
    counter <- 0
    while (g_val_update <= g_val_current + alp*ep*grad_g_max) {
        ep <- gmm*ep
        W_candidate_update <- W_current %*% expm(ep*Y_max)
        g_val_update <- compute_normalized_log_likelihood_g(SimData, W_candidate_update, estimator_current)
        counter <- counter + 1
        if (counter > 20) break
    }  
    return(W_candidate_update)
}


yuan_samworth_alternating <- function(SimData, W_init, alp=0.3, gmm=0.5, max_iter=20, converged_threshold=1e-7) {

    W_current <- W_init
    # visualize_independent_directions_from_W(SimData$W, W_current)
    estimator_current <- fit_marginals_given_W(SimData, W_current, plotting=FALSE)

    llh_all <- c()
    llh_all[1] <- compute_normalized_log_likelihood_g(SimData, W_current, estimator_current)
    cat("Initial likelihood =  ", llh_all[1])

    CONVERGED_FLAG <- FALSE

    for (itind in 1:max_iter) {


        W_current <- update_W_given_marginals(SimData, W_current, estimator_current, alp=alp, gmm=gmm)

        estimator_current <- fit_marginals_given_W(SimData, W_current, plotting=FALSE)

        llh_all[itind+1] <- compute_normalized_log_likelihood_g(SimData, W_current, estimator_current)
        cat("Current likelihood =  ", llh_all[itind+1])

        if ((llh_all[itind+1]-llh_all[itind])/abs(llh_all[itind]) < converged_threshold) {
            CONVERGED_FLAG <- TRUE
            cat("Converged by iteration ", itind)
            break
        }

        # cat("Iter: ")
        # print(itind)
        # print(llh_all[itind])
        # visualize_independent_directions_from_W(SimData$W, W_current)
    }

    return(list("estimator"=estimator_current, "llh_all"=llh_all, "converged_flag"=CONVERGED_FLAG))
}

yuan_samworth_alternating_with_repeats <- function(SimData, num_repeats=10, alp=0.3, gmm=0.5, max_iter=20, converged_threshold=1e-7) {

    all_algo_outputs <- list()
    all_final_log_likelihoods <- rep(0, num_repeats)

    for (repind in 1:num_repeats) {
        W_init <- fix_signs_fun(randortho(d, type="orthonormal"))
        all_algo_outputs[[repind]] <- yuan_samworth_alternating(SimData, W_init, alp=0.3, gmm=0.5, max_iter=20, converged_threshold=1e-7)
        all_final_log_likelihoods[repind] <- tail(all_algo_outputs[[repind]]$llh_all, n=1)
    }
    best_llh <- max(all_final_log_likelihoods)
    best_repeat <- which(all_final_log_likelihoods == best_llh)

    return(all_algo_outputs[[best_repeat]])
}


evaluate_yuan_samworth_estimator_vectorized <- function(X_mat, my_estimator){
    # X_mat must be organized such that each row is a point in R^d
    d <- NCOL(X_mat)
    n_eval <- NROW(X_mat)
    
    # X_mat_centered = t(t(X_mat) - my_estimator$mean_vec) # subtract the mean computed using the training data!

    Z_mat_hat = X_mat %*% t(my_estimator$W_given)
    
    result = rep(1.0, n_eval)
    for (dind in 1:d){
        result = result * evaluateLogConDens(Z_mat_hat[,dind], my_estimator$marginals[[dind]], which=2)[,3]
    }
    return(result)
}



