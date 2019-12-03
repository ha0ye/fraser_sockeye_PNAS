compute_ricker_curve <- function(recruits, spawners, num_iter = 30000, num_burnin = 20000)
{
    jags.data <- list(recruits = recruits, spawners = spawners)
    jags.params <- c("alpha", "beta")
    
    # get estimates for params
    min_S_index <- which(jags.data$spawners == min(jags.data$spawners))
    alpha_hat <- log(jags.data$recruits[min_S_index] / jags.data$spawners[min_S_index])
    max_R_index <- which(jags.data$recruits == max(jags.data$recruits))
    beta_hat <- jags.data$spawners[max_R_index]
    
    # use param estimates to set initial values for chains
    jags.inits <- list(list(alpha = alpha_hat / 2.71828, beta = beta_hat / 2.71828, 
                            .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                       list(alpha = alpha_hat / 2.71828, beta = beta_hat * 2.71828, 
                            .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                       list(alpha = alpha_hat * 2.71828, beta = beta_hat / 2.71828, 
                            .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                       list(alpha = alpha_hat * 2.71828, beta = beta_hat * 2.71828, 
                            .RNG.seed = 1234, .RNG.name = "base::Super-Duper"))
    
    # run jags
    my_model <- jags.model(file = "ricker_model.txt", data = jags.data, inits = jags.inits, 
                           n.chains = length(jags.inits))
    update(my_model, n.iter = num_burnin)
    jags_output <- jags.samples(my_model, jags.params, n.iter = num_iter)
    posterior_alpha <- as.vector(jags_output[["alpha"]][1,,])
    posterior_beta <- as.vector(jags_output[["beta"]][1,,])
    
    return(data.frame(alpha = posterior_alpha, beta = posterior_beta))
}

extract_results_for_best_models <- function()
{
    simplex_univar_stats <- readRDS("stats_simple_EDM.RDS")
    simplex_univar_stats$data <- "univariate"
    simplex_univar_stats$columns <- "eff"
    simplex_univar_stats$model <- "EDM"
    simplex_univar_stats$E <- 1
    
    ricker_univar_stats <- readRDS("stats_standard_ricker.RDS")
    ricker_univar_stats$data <- "univariate"
    ricker_univar_stats$columns <- "eff"
    ricker_univar_stats$model <- "Ricker"
    ricker_univar_stats$E <- 1
    ricker_univar_stats$stk <- rownames(ricker_univar_stats)
    
    # load data    
    temp <- readRDS("stats_multivariate_EDM.RDS")
    
    simplex_multivar_stats <- do.call(rbind, temp)
    simplex_multivar_stats$model <- "EDM"
    simplex_multivar_stats$stk <- factor(simplex_multivar_stats$stk)
    
    ricker_multivar_stats <- readRDS("stats_extended_ricker.RDS")
    ricker_multivar_stats$model <- "Ricker"
    ricker_multivar_stats$E <- sapply(strsplit(ricker_multivar_stats$columns, ", "), length)
    
    # get results and save
    ricker_univar_preds <- readRDS("results_standard_ricker.RDS")
    simplex_univar_preds <- readRDS("results_simple_EDM.RDS")
    ricker_multivar_preds <- readRDS("results_extended_ricker.RDS")
    simplex_multivar_preds <- readRDS("results_multivariate_EDM.RDS")
    
    preds <- list()
    stats <- list()
    for(stk in names(ricker_univar_preds))
    {
        ricker_stats <- ricker_multivar_stats[ricker_multivar_stats$stk == stk,]
        ricker_multivar_model <- ricker_stats$columns[which.max(ricker_stats$rho)]
        ricker_stats$data <- "multivariate"
        
        simplex_stats <- simplex_multivar_stats[simplex_multivar_stats$stk == stk,]
        simplex_multivar_model <- simplex_stats$columns[which.max(simplex_stats$rho)]
        simplex_stats$data <- "multivariate"
        
        preds[[stk]] <- data.frame(year = ricker_univar_preds[[stk]]$year, 
                                   obs = ricker_univar_preds[[stk]]$obs, 
                                   ricker_univar_pred = ricker_univar_preds[[stk]][,"pred"], 
                                   simplex_univar_pred = simplex_univar_preds[[stk]][,"pred"], 
                                   ricker_multivar_pred = ricker_multivar_preds[[stk]][,ricker_multivar_model], 
                                   simplex_multivar_pred = simplex_multivar_preds[[stk]][,simplex_multivar_model])
        stats[[stk]] <- rbind(ricker_univar_stats[ricker_univar_stats$stk == stk,], 
                              simplex_univar_stats[simplex_univar_stats$stk == stk,], 
                              ricker_stats[which.max(ricker_stats$rho),],
                              simplex_stats[which.max(simplex_stats$rho),])
        
    }
    saveRDS(preds, file = "preds_combined.RDS")
    saveRDS(stats, file = "stats_combined.RDS")
    
    return()
}

compute_seymour_ricker_params <- function()
{
    compute_ricker_curve_params <- function(stock_df)
    {
        years <- stock_df$yr
        recruits <- stock_df$rec
        spawners <- stock_df$eff
        valid <- is.finite(recruits) & is.finite(spawners)
        years <- years[valid]
        recruits <- recruits[valid]
        spawners <- spawners[valid]
        
        first_half <- 1:floor(NROW(recruits)/2)
        second_half <- (floor(NROW(recruits)/2)+1):NROW(recruits)
        
        # compute params
        first_half_params <- compute_ricker_curve(recruits[first_half], spawners[first_half])
        second_half_params <- compute_ricker_curve(recruits[second_half], spawners[second_half])
        params <- compute_ricker_curve(recruits, spawners)
        
        return(list(first_half_params = first_half_params, 
                    second_half_params = second_half_params, 
                    params = params))
    }
    
    compute_ricker_curve <- function(recruits, spawners, num_iter = 30000, num_burnin = 20000)
    {
        jags.data <- list(recruits = recruits, spawners = spawners)
        jags.params <- c("alpha", "beta")
        
        # get estimates for params
        min_S_index <- which(jags.data$spawners == min(jags.data$spawners))
        alpha_hat <- log(jags.data$recruits[min_S_index] / jags.data$spawners[min_S_index])
        max_R_index <- which(jags.data$recruits == max(jags.data$recruits))
        beta_hat <- jags.data$spawners[max_R_index]
        
        # use param estimates to set initial values for chains
        jags.inits <- list(list(alpha = alpha_hat / 2.71828, beta = beta_hat / 2.71828, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                           list(alpha = alpha_hat / 2.71828, beta = beta_hat * 2.71828, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                           list(alpha = alpha_hat * 2.71828, beta = beta_hat / 2.71828, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                           list(alpha = alpha_hat * 2.71828, beta = beta_hat * 2.71828, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"))
        
        # run jags
        my_model <- jags.model(file = "ricker_model.txt", data = jags.data, inits = jags.inits, 
                               n.chains = length(jags.inits))
        update(my_model, n.iter = num_burnin)
        jags_output <- jags.samples(my_model, jags.params, n.iter = num_iter)
        posterior_alpha <- as.vector(jags_output[["alpha"]][1,,])
        posterior_beta <- as.vector(jags_output[["beta"]][1,,])
        
        return(data.frame(alpha = posterior_alpha, beta = posterior_beta))
    }
    
    load("block_data.Rdata")
    
    ricker_posteriors <- compute_ricker_curve_params(block_data[["Seymour"]])
    saveRDS(ricker_posteriors, file = "params_ricker_seymour.RDS")
    return()
}

compute_seymour_ricker_env_params <- function()
{
    fit_ricker_params <- function(df, 
                                  num_iter = 30000, num_burnin = 20000, 
                                  model_file = "ricker_model_env.txt", 
                                  quiet = FALSE)
    {        
        # fit model
        jags.params <- c("alpha", "beta", "g")
        jags.data <- df
        
        # get estimates for params
        min_S_index <- which.min(jags.data$spawners)
        alpha_hat <- log(jags.data$recruits[min_S_index] / jags.data$spawners[min_S_index])
        max_R_index <- which.max(jags.data$recruits)
        beta_hat <- 1/jags.data$spawners[max_R_index]
        
        # use param estimates to set initial values for chains
        jags.inits <- list(list(alpha = alpha_hat / 2.71828, beta = beta_hat / 2.71828, g = 0, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                           list(alpha = alpha_hat / 2.71828, beta = beta_hat * 2.71828, g = 0, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"), 
                           list(alpha = alpha_hat * 2.71828, beta = beta_hat / 2.71828, g = 0, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"),
                           list(alpha = alpha_hat * 2.71828, beta = beta_hat * 2.71828, g = 0, 
                                .RNG.seed = 1234, .RNG.name = "base::Super-Duper"))
        
        # run jags
        my_model <- jags.model(file = model_file, data = jags.data, inits = jags.inits, 
                               n.chains = length(jags.inits), quiet = quiet)
        update(my_model, n.iter = num_burnin)
        jags_output <- jags.samples(my_model, jags.params, n.iter = num_iter)
        return(data.frame(alpha = as.vector(jags_output[["alpha"]][1,,]), 
                          beta = as.vector(jags_output[["beta"]][1,,]), 
                          g = as.vector(jags_output[["g"]][1,,])))
    }
    
    load("block_data.Rdata")
    
    stk_name <- "Seymour"
    env_var <- "PT_apr"
    
    # setup data
    stock_df <- block_data[[stk_name]]
    valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
    recruits <- stock_df$rec45[valid]
    spawners <- stock_df$eff[valid]
    years <- stock_df$yr[valid]
    env <- stock_df[valid,env_var]
    
    params <- fit_ricker_params(data.frame(recruits, spawners, env))
    saveRDS(params, file = "params_ricker_env_seymour.RDS")
    return()
}

compute_chilko_smolt_forecasts <- function()
{
    cat("forecast with smolts for Chilko... ")
    start_time <- proc.time()
    
    load("block_data.Rdata")
    env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                   "ET_apr", "ET_may", "ET_jun", 
                   "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                   "PDO_win")
    
    stock_df <- block_data[["Chilko"]]
    
    temp <- normalize_by_cycle_line(stock_df$juv)
    stock_df$juv_n <- temp$ts
    
    # set up recruits and spawners
    valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
    years <- stock_df$yr[valid]
    returns <- stock_df$ret[valid]
    spawners <- stock_df$eff_n[valid]
    smolts <- stock_df$juv_n[valid]
    recruits_4 <- stock_df$rec4_n[valid]
    mu_4 <- stock_df$rec4_mu[valid]
    sigma_4 <- stock_df$rec4_sigma[valid]
    recruits_5 <- stock_df$rec5_n[valid]
    mu_5 <- stock_df$rec5_mu[valid]
    sigma_5 <- stock_df$rec5_sigma[valid]
    env <- normalize(stock_df[,env_names])
    
    # make block
    block <- data.frame(years = years, eff = spawners, juv = smolts, 
                        rec4 = recruits_4, rec5 = recruits_5)
    block <- cbind(block, env[valid, ])
    
    if(length(returns) < 2) # check for enough data
        return(data.frame(year = NaN, obs = NaN, pred = NaN))
    
    columns <- list()
    for(E in 1:2)
    {
        columns <- c(columns, combn(env_names, E, simplify = FALSE))
    }
    columns <- lapply(columns, function(embedding) c("eff", "juv", embedding))
    columns[[length(columns)+1]] <- c("eff", "juv")
    rec4_preds <- do.call(cbind, block_lnlp_4(block, target_column = 3, columns = columns))
    rec5_preds <- do.call(cbind, block_lnlp_4(block, target_column = 4, columns = columns))
    rec4_preds <- rec4_preds*sigma_4 + mu_4
    rec5_preds <- rec5_preds*sigma_5 + mu_5
    forecasts <- data.frame(rec4_preds + rbind(NA, rec5_preds[1:(NROW(block)-1),]))
    names(forecasts) <- lapply(columns, function(v) paste(v, sep = "", collapse = ", "))
    output <- cbind(year = years, obs = returns, forecasts)
    saveRDS(output, file = "results_chilko_smolts.RDS")
    
    stats <- do.call(rbind, lapply(3:NCOL(output), function(j) {
        compute_stats(output[,2], output[,j])
    }))
    stats$columns <- names(output)[3:NCOL(output)]
    saveRDS(stats, file = "stats_chilko_smolts.RDS")
    
    elapsed_time <- proc.time() - start_time
    cat("(", elapsed_time[3], " sec.)\n", sep = "")
    return()
}

compute_ccm <- function()
{
    load("block_data.Rdata")
    env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                   "ET_apr", "ET_may", "ET_jun", 
                   "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                   "PDO_win")
    
    ccm_table <- do.call(rbind, lapply(block_data, function(stock_df) {
        valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
        block <- stock_df[valid,]
        
        ccm_rhos <- do.call(cbind, lapply(env_names, function(env_var) {
            output <- block_lnlp(block, tp = 0, target_column = env_var, 
                                 columns = c("rec45_n", "eff_n"), silent = TRUE)
            return(output$rho)
        }))
        colnames(ccm_rhos) <- env_names
        ccm_rhos <- cbind(N = sum(valid), ccm_rhos)
        return(ccm_rhos)
    }))
    rownames(ccm_table) <- names(block_data)
    saveRDS(ccm_table, file = "results_ccm.RDS")
    return()
}