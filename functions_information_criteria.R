# The following function constructs a dataframe of all possible copulas, values of Kendall's tau and sample sizes
construct_param_GRID <- function(cop_TYPES, sample_size) {
  # Initialize empty dataframes with correct column names
  param_GRID_indep <- data.frame(cop_name = character(), 
                                 tau = numeric(), 
                                 sample_size = integer())
  param_GRID_amh <- data.frame(cop_name = character(), 
                               tau = numeric(), 
                               sample_size = integer())
  param_GRID_other <- data.frame(cop_name = character(), 
                                 tau = numeric(), 
                                 sample_size = integer())
  
  for (cop_type in cop_TYPES) {
    if (cop_type %in% c("clayton", "gumbel", "joe", "frank", "gaussian")) {
      param_GRID_other <- rbind(param_GRID_other, 
                                expand.grid(cop_name = cop_type, 
                                            tau = c(0.25, 0.5, 0.75), 
                                            #tau = c(0.05, 0.1),
                                            sample_size = sample_size))
    }
  }
  
  # Combine and order the dataframes
  param_GRID_ <- rbind(param_GRID_other, param_GRID_amh)
  param_GRID_ <- param_GRID_[order(as.character(param_GRID_$cop_name), 
                                   param_GRID_$tau, 
                                   param_GRID_$sample_size), ]
  
  # Add independence copula at the top
  param_GRID <- rbind(param_GRID_indep, param_GRID_)
  
  return(param_GRID)
}



# The following function returns a copula object based on the specification of Kendall's tau (tau) and dimension (d)
copula_choice_by_tau <- function(cop_name, tau, d) {
  if (cop_name == "clayton") {
    cop <- claytonCopula(param=iTau(claytonCopula(), tau=tau), dim=d)} 
  else if (cop_name == "gumbel") {
    cop <- gumbelCopula(param=iTau(gumbelCopula(), tau=tau), dim=d)} 
  else if (cop_name == "gaussian") {
    cop <- normalCopula(param=iTau(normalCopula(), tau=rep(tau, d*(d-1)/2)), dim=d, dispstr="un")} 
  else if (cop_name == "joe") {
    cop <- joeCopula(param=iTau(joeCopula(), tau=tau), dim=d)}
  else if (cop_name == "frank") {
    cop <- frankCopula(param=iTau(frankCopula(), tau=tau), dim=d)}
  return(cop)
}



# The following function returns a copula object based on the specification of copula parameter (param) and dimension (d)
copula_choice_by_param <- function(cop_name, param, d) {
  if (cop_name == "clayton") {
    cop <- claytonCopula(param=param, dim=d)} 
  else if (cop_name == "gumbel") {
    cop <- gumbelCopula(param=param, dim=d)} 
  else if (cop_name == "gaussian") {
    cop <- normalCopula(param=param, dim=d, dispstr="un")} 
  else if (cop_name == "joe") {
    cop <- joeCopula(param=param, dim=d)}
  else if (cop_name == "frank") {
    cop <- frankCopula(param=param, dim=d)}
  return(cop)
}



# The following function returns a copula method WITHOUT specifying the copula parameter (param) or dimension (d)
copula_funct_choice <- function(cop_name){
  if (cop_name == "clayton"){
    return(claytonCopula())}
  else if (cop_name == "gumbel") {
    return(gumbelCopula())} 
  else if (cop_name == "gaussian") {
    return(normalCopula())}
  else if (cop_name == "joe") {
    return(joeCopula())}
  else if (cop_name == "frank") {
    return(frankCopula())}
}



# Our parameter domain contains parameter values for tau = 0.95 and tau = -0.95 (where it makes sense)
copula_param_domain <- function(cop_name){
  if (cop_name == "clayton"){
    upper_domain <- iTau(copula_funct_choice("clayton"), tau = 0.95)
    lower_domain <- iTau(copula_funct_choice("clayton"), tau = -0.95)
    return(c(lower_domain, upper_domain))}
  else if (cop_name == "gumbel") {
    upper_domain <- iTau(copula_funct_choice("gumbel"), tau = 0.95)
    lower_domain <- 1
    return(c(lower_domain, upper_domain))} 
  else if (cop_name == "gaussian") {
    upper_domain <- iTau(copula_funct_choice("gaussian"), tau = 0.95)
    lower_domain <- iTau(copula_funct_choice("gaussian"), tau = -0.95)
    return(c(lower_domain, upper_domain))}
  else if (cop_name == "joe") {
    upper_domain <- iTau(copula_funct_choice("joe"), tau = 0.95)
    lower_domain <- 1
    return(c(lower_domain, upper_domain))}
  else if (cop_name == "frank") {
    upper_domain <- iTau(copula_funct_choice("frank"), tau = 0.95)
    lower_domain <- iTau(copula_funct_choice("frank"), tau = -0.95)
    return(c(lower_domain, upper_domain))}
}



# Domain of Kendall's tau
copula_tau_domain <- function(cop_name){
  if (cop_name == "clayton"){
    return(c(-0.95, 0.95))}
  else if (cop_name == "gumbel") {
    return(c(0, 0.95))} 
  else if (cop_name == "gaussian") {
    return(c(-0.95, 0.95))}
  else if (cop_name == "t") {
    return(c(-0.95, 0.95))}
  else if (cop_name == "joe") {
    return(c(0, 0.95))}
  else if (cop_name == "frank") {
    return(c(-0.95, 0.95))}
}



# The following function compute initial estimate of our parameter
initial_param <- function(cop_name_fit, pseudo_obs){
  
  # Initial parameter estimate using the inverse tau method:
  
  if (is.na(suppressWarnings(cor(pseudo_obs[,1], pseudo_obs[,2], method = "kendall")))){ # This situation can occur, especially with small sample sizes ...
    # Jittering with a magnitude of 1e-6 should be small enough relative to pseudo_obs
    tau0 <- cor(pseudo_obs[,1] + runif(nrow(pseudo_obs), -1e-6, 1e-6), pseudo_obs[,2] + runif(nrow(pseudo_obs), -1e-6, 1e-6), method = "kendall")
  } else{
    tau0 <- cor(pseudo_obs[,1], pseudo_obs[,2], method = "kendall")
  }
  
  # Let determine the domain of Kendall's tau
  domain_tau <- copula_tau_domain(cop_name_fit)
  
  # Check if the initial estimate, tau0, lies within a reasonable domain
  if (tau0 < domain_tau[1]){
    tau0 <- domain_tau[1]
  }else if (domain_tau[2] < tau0){
    tau0 <- domain_tau[2]
  }
  
  # From tau0, one obtains the initial estimate of the parameter
  param0 <- iTau(copula_funct_choice(cop_name_fit), tau = tau0)
  
  # Let determine the domain of our parameter
  domain_param <- copula_param_domain(cop_name_fit)
  
  # One more check to ensure that param0 also lies within a reasonable parameter domain. This step may be 
  # unnecessary, as we have already checked if tau0 lies within a reasonable domain (for Kendall's tau), 
  # but is kept for safety
  
  if (param0 < domain_param[1]){
    param0 <- domain_param[1]
  }else if (domain_param[2] < param0){
    param0 <- domain_param[2]
  }
  
  return(param0)
}



neg_loglik_fun <- function(param, pseudo_obs, name) {
  return(-loglikCopula(param, pseudo_obs, copula_funct_choice(name)))
}

