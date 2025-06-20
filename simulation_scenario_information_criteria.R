# The following function is the main part of the entire simulation study, which simulates a scenario for 
# specific copula, specific value of Kendall's tau and specific sample size

simulate_cop_name_tau_sample_size <- function(cop_name, cop_TYPES, true_tau, n, num_samples, d, log_file){
  
  # Preallocate memory for results
  cop_results <- data.frame(
    cop_name = character(num_samples),
    sample_size = integer(num_samples),
    iteration = integer(num_samples),
    
    AIC_clayton = numeric(num_samples),
    AIC_gumbel = numeric(num_samples),
    AIC_joe = numeric(num_samples ),
    AIC_frank = numeric(num_samples),
    AIC_gaussian = numeric(num_samples),
    
    xv1_clayton = numeric(num_samples),
    xv1_gumbel = numeric(num_samples),
    xv1_joe = numeric(num_samples ),
    xv1_frank = numeric(num_samples),
    xv1_gaussian = numeric(num_samples),
    
    xv_CIC_clayton = numeric(num_samples),
    xv_CIC_gumbel = numeric(num_samples),
    xv_CIC_joe = numeric(num_samples ),
    xv_CIC_frank = numeric(num_samples),
    xv_CIC_gaussian = numeric(num_samples),
    
    xv_nv_clayton = numeric(num_samples),
    xv_nv_gumbel = numeric(num_samples),
    xv_nv_joe = numeric(num_samples ),
    xv_nv_frank = numeric(num_samples),
    xv_nv_gaussian = numeric(num_samples),
    
    stringsAsFactors = FALSE
  )
  
  # All copulas are chosen by specifying different values of Kendall's tau,
  # via the function copula_choice_by_tau
  
  cop <- copula_choice_by_tau(cop_name, true_tau, d)
  true_param <- iTau(copula_funct_choice(cop_name), tau = true_tau)
  
  
  for (i in 1:num_samples) {
    
    # The following is a way to define unique seed values for each simulation of our algorithm
    seed_value <- which(cop_TYPES == cop_name)* 1e7 + as.integer(abs(true_tau) * 100) * 1e4 + n * 10 + i
    
    current_AIC_dict <- setNames(rep(0, 5), c("clayton", "gumbel", "joe", "frank", "gaussian"))
    
    current_xv1_dict <- setNames(rep(0, 5), c("clayton", "gumbel", "joe", "frank", "gaussian"))
    
    current_xv_CIC_dict <- setNames(rep(0, 5), c("clayton", "gumbel", "joe", "frank", "gaussian"))
    
    current_xv_nv_dict <- setNames(rep(0, 5), c("clayton", "gumbel", "joe", "frank", "gaussian"))
    
    set.seed(seed_value)
    
    data <- rCopula(n, cop)
    
    pseudo_obs <- pobs(data)
    
    
    for (cop_name_FIT in cop_TYPES) {
      
      param0 <- initial_param(cop_name_FIT, pseudo_obs)
      
      domain_param_FIT <- copula_param_domain(cop_name_FIT)
      
      # gradient of log pseudo-likelihood
      phi_copula <- phi_copula_LIST[[cop_name_FIT]]
      
      # As we have determined the reasonable domain for our parameters, the main optimization method is selected
      # as 'L-BFGS-B' (this method requires the function to have finite values on the searched interval)

      param_estim <- try(
        suppressWarnings(
          optim(
            par = param0,
            fn = function(param) neg_loglik_fun(param, pseudo_obs, cop_name_FIT),
            gr = function(param) sum(-phi_copula(pseudo_obs[,1], pseudo_obs[,2], param)),
            method = "L-BFGS-B",
            lower = domain_param_FIT[1],
            upper = domain_param_FIT[2]
          )$par
        ),
        silent = TRUE
      )
      
      # If the 'L-BFGS-B' method fails, we use the 'Brent' method, which is more robust for one-dimensional
      # optimization

      if (inherits(param_estim, "try-error")) {
        param_estim <- suppressWarnings(
          fitCopula(copula_funct_choice(cop_name_FIT), pseudo_obs, method = "mpl",
                    optim.method = "Brent", start = param0, lower = domain_param_FIT[1], upper = domain_param_FIT[2])@estimate
        )
      }
      
      ### Akaike Information Criterion
      current_AIC_dict[[cop_name_FIT]] <- 2*loglikCopula(param_estim, pseudo_obs, copula_funct_choice(cop_name_FIT)) - 2
      
      ### Leave one out cross-validation information criterion xv1
      
      param0_xv <- param_estim
      xv1 <- 0
      
      for (k in 1:n){
        pseudo_obs_ <- pseudo_obs[-k,]
        
        param_estim_ <- try(
          suppressWarnings(
            optim(
              par = param0_xv,
              fn = function(param) neg_loglik_fun(param, pseudo_obs_, cop_name_FIT),
              gr = function(param) sum(-phi_copula(pseudo_obs_[,1], pseudo_obs_[,2], param)),
              method = "L-BFGS-B",
              lower = domain_param_FIT[1],
              upper = domain_param_FIT[2]
            )$par
          ),
          silent = TRUE
        )
        
        if (inherits(param_estim_, "try-error")) {
          param_estim_ <- suppressWarnings(
            fitCopula(copula_funct_choice(cop_name_FIT), pseudo_obs_, method = "mpl",
                      optim.method = "Brent", start = param0_xv, lower = domain_param_FIT[1], upper = domain_param_FIT[2])@estimate
          )
        }
        
        xv1 <- xv1 + loglikCopula(param_estim_, pseudo_obs_, copula_funct_choice(cop_name_FIT))
      }
      
      xv1 <- xv1/n
      
      current_xv1_dict[[cop_name_FIT]] <- xv1
      
      ### Asymptotic equivalent version of leave one out cross validation information criterion (xv-CIC)
      
      # DEFINE FUNCTIONS
      log_copula <- log_copula_LIST[[cop_name_FIT]]
      d_log_copula <- d_log_copula_LIST[[cop_name_FIT]]
      #phi_copula <- phi_copula_LIST[[cop_name_FIT]]
      dparam_phi_copula <- dparam_phi_copula_LIST[[cop_name_FIT]]
      du_phi_copula <- du_phi_copula_LIST[[cop_name_FIT]]
      dv_phi_copula <- dv_phi_copula_LIST[[cop_name_FIT]]
      z_function_copula <- z_function_copula_LIST[[cop_name_FIT]]
      
      # Compute biases
      
      J_estim <- (-sum(dparam_phi_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim))/n)
      
      bias1_xv_CIC <- sum((phi_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim))^2)/(n*J_estim)
      
      du_phi_copula_vec <- du_phi_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim)
      dv_phi_copula_vec <- dv_phi_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim)
      
      bias2_xv_CIC <- sum(phi_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim)*
                            z_function_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim, pseudo_obs, du_phi_copula_vec, dv_phi_copula_vec))/(n*J_estim)
      
      bias3_xv_CIC <- sum(t(d_log_copula(pseudo_obs[,1], pseudo_obs[,2], param_estim))*cbind(1-pseudo_obs[,1], 1-pseudo_obs[,2]))/n
      
      
      current_xv_CIC_dict[[cop_name_FIT]] <- 2*loglikCopula(param_estim, pseudo_obs, copula_funct_choice(cop_name_FIT)) - 2*(bias1_xv_CIC + bias2_xv_CIC + bias3_xv_CIC)
      
      ### Leave n_{v} out cross validation information criterion
      
      nc <- round((n)^(0.9))
      nv <- n - nc
      b <- round(0.8*n)
      
      set.seed(seed_value)
      all_indices_nc <- replicate(b, sample(1:n, nc), simplify = FALSE)
      
      xv_nv <- 0
      
      for (j in 1:b){
        
        indices_nc <- all_indices_nc[[j]]
        pseudo_obs_nc <- pseudo_obs[indices_nc, ]
        pseudo_obs_nv <- pseudo_obs[-indices_nc, ]
        
        param_estim_nc <- try(
          suppressWarnings(
            optim(
              par = param0_xv,
              fn = function(param) neg_loglik_fun(param, pseudo_obs_nc, cop_name_FIT),
              gr = function(param) sum(-phi_copula(pseudo_obs_nc[,1], pseudo_obs_nc[,2], param)),
              method = "L-BFGS-B",
              lower = domain_param_FIT[1],
              upper = domain_param_FIT[2]
            )$par
          ),
          silent = TRUE
        )

        if (inherits(param_estim_nc, "try-error")) {
          param_estim_nc <- suppressWarnings(
            fitCopula(copula_funct_choice(cop_name_FIT), pseudo_obs_nc, method = "mpl",
                      optim.method = "Brent", start = param0_xv, lower = domain_param_FIT[1], upper = domain_param_FIT[2])@estimate
          )
        }
        
        xv_nv <- xv_nv + loglikCopula(param_estim_nc, pseudo_obs_nv, copula_funct_choice(cop_name_FIT))
      }
      
      xv_nv <- xv_nv/(nv*b)
      
      current_xv_nv_dict[[cop_name_FIT]] <- xv_nv
      
    }
    
    # Write progress to the log file
    log_data <- paste(cop_name, n, i, true_param, true_tau, 
                      current_AIC_dict[["clayton"]], current_AIC_dict[["gumbel"]], current_AIC_dict[["joe"]], current_AIC_dict[["frank"]],  current_AIC_dict[["gaussian"]], 
                      current_xv1_dict[["clayton"]], current_xv1_dict[["gumbel"]], current_xv1_dict[["joe"]], current_xv1_dict[["frank"]],  current_xv1_dict[["gaussian"]], 
                      current_xv_CIC_dict[["clayton"]], current_xv_CIC_dict[["gumbel"]], current_xv_CIC_dict[["joe"]], current_xv_CIC_dict[["frank"]], current_xv_CIC_dict[["gaussian"]], 
                      current_xv_nv_dict[["clayton"]], current_xv_nv_dict[["gumbel"]], current_xv_nv_dict[["joe"]], current_xv_nv_dict[["frank"]], current_xv_nv_dict[["gaussian"]], 
                      sep = ",")
    write(log_data, file = log_file, append = TRUE)
    
    # Add results
    cop_results[i, ] <- list(
      cop_name = cop_name
      ,sample_size = n
      ,iteration = i
      
      ,AIC_clayton = current_AIC_dict[["clayton"]]
      ,AIC_gumbel = current_AIC_dict[["gumbel"]]
      ,AIC_joe = current_AIC_dict[["joe"]]
      ,AIC_frank = current_AIC_dict[["frank"]]
      ,AIC_gaussian = current_AIC_dict[["gaussian"]]
      
      ,xv1_clayton = current_xv1_dict[["clayton"]]
      ,xv1_gumbel = current_xv1_dict[["gumbel"]]
      ,xv1_joe = current_xv1_dict[["joe"]]
      ,xv1_frank = current_xv1_dict[["frank"]]
      ,xv1_gaussian = current_xv1_dict[["gaussian"]]
      
      ,xv_CIC_clayton = current_xv_CIC_dict[["clayton"]]
      ,xv_CIC_gumbel = current_xv_CIC_dict[["gumbel"]]
      ,xv_CIC_joe = current_xv_CIC_dict[["joe"]]
      ,xv_CIC_frank = current_xv_CIC_dict[["frank"]]
      ,xv_CIC_gaussian = current_xv_CIC_dict[["gaussian"]]
      
      ,xv_nv_clayton = current_xv_nv_dict[["clayton"]]
      ,xv_nv_gumbel = current_xv_nv_dict[["gumbel"]]
      ,xv_nv_joe = current_xv_nv_dict[["joe"]]
      ,xv_nv_frank = current_xv_nv_dict[["frank"]]
      ,xv_nv_gaussian = current_xv_nv_dict[["gaussian"]]
    )
    
  }
  return(cop_results)
}

