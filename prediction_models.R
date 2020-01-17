simple_EDM <- function()
{
    forecast <- function(stock_df)
    {  
        make_forecasts <- function(block, mu_4, sigma_4, mu_5, sigma_5)
        {
            rec4 <- block_lnlp_4(block, target_column = 2, columns = 1)
            rec5 <- block_lnlp_4(block, target_column = 3, columns = 1)
            
            rec4 <- rec4*sigma_4 + mu_4
            rec5 <- rec5*sigma_5 + mu_5
            return(rec4 + c(NA, rec5[1:(NROW(block)-1)]))
        }
        
        # set up recruits and spawners
        valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
        returns <- stock_df$ret[valid]
        years <- stock_df$yr[valid]
        spawners <- stock_df$eff_n[valid]
        recruits_4 <- stock_df$rec4_n[valid]
        mu_4 <- stock_df$rec4_mu[valid]
        sigma_4 <- stock_df$rec4_sigma[valid]
        recruits_5 <- stock_df$rec5_n[valid]
        mu_5 <- stock_df$rec5_mu[valid]
        sigma_5 <- stock_df$rec5_sigma[valid]
        
        # make block
        block <- data.frame(years = years, eff = spawners, 
                            rec4 = recruits_4, rec5 = recruits_5)
        
        if(length(returns) < 2) # check for enough data
            return(data.frame(year = NaN, obs = NaN, pred = NaN))
        
        forecasts <- make_forecasts(block, mu_4, sigma_4, mu_5, sigma_5)
        return(data.frame(year = years, obs = returns, pred = forecasts))
    }
    
    load("block_data.Rdata")
    
    # make forecasts for each stock
    results <- lapply(names(block_data), 
                      function(stk_name) {
                          cat("forecasting for ", stk_name, "... ", sep = "")
                          start_time <- proc.time()
                          output <- forecast(block_data[[stk_name]])
                          elapsed_time <- proc.time() - start_time
                          cat("(", elapsed_time[3], " sec.)\n", sep = "")
                          return(output)
                      })
    names(results) <- names(block_data)
    saveRDS(results, file = "results_simple_EDM.RDS")
    
    # compute stats
    stats <- do.call(rbind, lapply(results, function(stock_results) {
        compute_stats(stock_results$obs, stock_results$pred)
    }))
    stats$stk <- names(block_data)
    saveRDS(stats, file = "stats_simple_EDM.RDS")
    return()
}

multivariate_EDM <- function()
{
    forecast <- function(stock_df)
    {
        load("block_data.Rdata")
        env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                       "ET_apr", "ET_may", "ET_jun", 
                       "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                       "PDO_win")
        
        # set up recruits and spawners
        valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
        years <- stock_df$yr[valid]
        returns <- stock_df$ret[valid]
        spawners <- stock_df$eff_n[valid]
        recruits_4 <- stock_df$rec4_n[valid]
        mu_4 <- stock_df$rec4_mu[valid]
        sigma_4 <- stock_df$rec4_sigma[valid]
        recruits_5 <- stock_df$rec5_n[valid]
        mu_5 <- stock_df$rec5_mu[valid]
        sigma_5 <- stock_df$rec5_sigma[valid]
        env <- normalize(stock_df[,env_names])
        
        # make block
        block <- data.frame(years = years, eff = spawners, 
                            rec4 = recruits_4, rec5 = recruits_5)
        block <- cbind(block, env[valid, ])
        
        if(length(returns) < 2) # check for enough data
            return(data.frame(year = NaN, obs = NaN, pred = NaN))
        
        columns <- list()
        for(E in 1:2)
        {
            columns <- c(columns, combn(env_names, E, simplify = FALSE))
        }
        columns <- lapply(columns, function(embedding) c("eff", embedding))
        columns <- c(columns, "eff")
        rec4_preds <- do.call(cbind, block_lnlp_4(block, target_column = 2, columns = columns))
        rec5_preds <- do.call(cbind, block_lnlp_4(block, target_column = 3, columns = columns))
        rec4_preds <- rec4_preds*sigma_4 + mu_4
        rec5_preds <- rec5_preds*sigma_5 + mu_5
        forecasts <- data.frame(rec4_preds + rbind(NA, rec5_preds[1:NROW(block)-1,]))
        names(forecasts) <- lapply(columns, function(v) paste(v, sep = "", collapse = ", "))
        output <- cbind(year = years, obs = returns, forecasts)
        
        return(output)
    }
    
    load("block_data.Rdata")
    
    # make forecasts for each stock
    results <- lapply(names(block_data), 
                      function(stk_name) {
                          cat("forecasting for ", stk_name, "... ", sep = "")
                          start_time <- proc.time()
                          output <- forecast(block_data[[stk_name]])
                          elapsed_time <- proc.time() - start_time
                          cat("(", elapsed_time[3], " sec.)\n", sep = "")
                          return(output)
                      })
    names(results) <- names(block_data)
    saveRDS(results, file = "results_multivariate_EDM.RDS")
    
    # compute stats
    stats <- lapply(names(block_data), function(stk_name) {
        output <- results[[stk_name]]
        stats <- do.call(rbind, lapply(3:NCOL(output), function(j) {
            compute_stats(output[,2], output[,j])
        }))
        stats$columns <- names(output)[3:NCOL(output)]
        stats$stk <- stk_name
        return(stats)        
    })
    
    stats <- lapply(stats, function(stats_df) {
        stats_df$E <- sapply(strsplit(stats_df$columns, ", "), length)
        with_eff_only <- subset(stats_df, E == 1)
        with_one_env_var <- subset(stats_df, E == 2)
        if(max(with_one_env_var$rho) <= with_eff_only$rho)
            return(subset(stats_df, E <= 2))
        best_env_var <- strsplit(with_one_env_var$columns[which.max(with_one_env_var$rho)], 
                                 ", ")[[1]][2]
        with_two_env_var <- subset(stats_df, E == 3)
        idx <- grep(best_env_var, with_two_env_var$columns)
        return(rbind(with_eff_only, with_one_env_var, with_two_env_var[idx,]))
    })
    
    saveRDS(stats, file = "stats_multivariate_EDM.RDS")
    return()
}


write_model_files <- function()
{
    write(
        "model
        {
        for (i in 1:length(recruits))
        {
        y[i] <- alpha - beta * spawners[i] + log(spawners[i])
        recruits[i] ~ dlnorm(y[i], tau_r)
        }
        
        # priors
        alpha ~ dnorm(0, 1e-6)
        beta ~ dnorm(0, 1e-6)
        tau_r ~ dgamma(0.001, 0.001)
        sigma <- pow(tau_r, -2)
        }", file = "ricker_model.txt")

    write(
        "model
        {
        for (i in 1:length(recruits))
        {
        y[i] <- alpha - beta * spawners[i] + log(spawners[i]) + g * env[i]
        recruits[i] ~ dlnorm(y[i], tau_r)
        }
        
        # priors
        g ~ dnorm(0, 1e-6)
        alpha ~ dnorm(0, 1e-6)
        beta ~ dnorm(0, 1e-6)
        tau_r ~ dgamma(0.001, 0.001)
        sigma <- pow(tau_r, -2)
        }", file = "ricker_model_env.txt")

    return()
}

standard_ricker <- function()
{
    fit_ricker_to_stock <- function(stock_df)
    {
        # set up recruits and spawners
        valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
        recruits <- stock_df$rec45[valid]
        spawners <- stock_df$eff[valid]
        years <- stock_df$yr[valid]
        returns <- stock_df$ret[valid]
        p4 <- mean(stock_df$rec4 / stock_df$rec45, na.rm = TRUE)
        
        if(length(recruits) < 2) # check for enough data
            return(data.frame(year = NaN, obs = NaN, pred = NaN))
        
        # make ricker model predictions
        ricker_pred <- ricker_ret_model_4(data.frame(years, spawners, recruits), p4 = p4)
        
        return(data.frame(year = years, obs = returns, pred = ricker_pred))
    }
    
    load("block_data.Rdata")
    results_ricker <- lapply(block_data, fit_ricker_to_stock)
    saveRDS(results_ricker, file = "results_standard_ricker.RDS")
    
    stats_ricker <- do.call(rbind, lapply(results_ricker, function(results) {
        compute_stats(results$obs, results$pred)
    }))
    saveRDS(stats_ricker, file = "stats_standard_ricker.RDS")
    return()
}

extended_ricker <- function()
{
    fit_ricker_to_stock <- function(stock_df)
    {
        # set up recruits and spawners
        valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
        recruits <- stock_df$rec45[valid]
        spawners <- stock_df$eff[valid]
        years <- stock_df$yr[valid]
        returns <- stock_df$ret[valid]
        p4 <- mean(stock_df$rec4 / stock_df$rec45, na.rm = TRUE)
        
        if(length(recruits) < 2) # check for enough data
            return(data.frame(year = NaN, obs = NaN, pred = NaN))
        
        # make ricker model predictions
        forecasts <- lapply(env_names, function(env_var) {
            ricker_ret_model_env_4(data.frame(years, spawners, recruits, 
                                              env = stock_df[valid, env_var]), p4 = p4)
        })
        columns <- lapply(env_names, function(env_var) c("eff", env_var))
        forecasts <- data.frame(do.call(cbind, forecasts))
        names(forecasts) <- lapply(columns, function(v) paste(v, sep = "", collapse = ", "))
        
        return(cbind(year = years, obs = returns, forecasts))
    }
    
    env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                   "ET_apr", "ET_may", "ET_jun", 
                   "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                   "PDO_win")
    
    load("block_data.Rdata")
    results <- lapply(block_data, fit_ricker_to_stock)
    saveRDS(results, file = "results_extended_ricker.RDS")
    
    stats <- lapply(results, function(stock_results) {
        temp <- do.call(rbind, lapply(3:NCOL(stock_results), function(j) {
            compute_stats(stock_results[,2], stock_results[,j])
        }))
        temp$columns <- names(stock_results)[3:NCOL(stock_results)]
        return(temp)
    })
    for(stk_name in names(results))
    {
        stats[[stk_name]]$stk <- stk_name
    }
    
    stats <- do.call(rbind, stats)
    saveRDS(stats, file = "stats_extended_ricker.RDS")
    return()
}

block_lnlp_4 <- function(block, target_column, columns, norm = FALSE)
{
    if(norm)
    {
        block[,columns] <- normalize(block[,columns])
    }
    
    lib_segments <- matrix(NA, nrow = 4, ncol = 2)
    segment_size <- NROW(block)/4
    start_index <- 1
    for(i in 1:4)
    {
        lib_segments[i,1] <- floor(start_index)
        end_index <- start_index - 1 + segment_size
        lib_segments[i,2] <- floor(end_index)
        start_index <- end_index+1
    }
    
    if(is.list(columns))
    {
        preds <- lapply(1:length(columns), function(x) {rep.int(NA, times = NROW(block))})
        for(i in 1:4)
        {
            pred_index <- lib_segments[i,1]:lib_segments[i,2]
            
            temp <- block_lnlp(block, lib = lib_segments[-i,], pred = lib_segments[i,], 
                               target_column = target_column, tp = 0, 
                               first_column_time = TRUE, 
                               columns = columns, stats_only = FALSE)
            
            for(j in 1:length(columns))
                preds[[j]][pred_index] <- temp[[j, "model_output"]]$pred[seq_along(pred_index)]
        }
    }
    else
    {
        preds <- rep.int(NA, times = NROW(block))
        for(i in 1:4)
        {
            pred_index <- lib_segments[i,1]:lib_segments[i,2]
            
            temp <- block_lnlp(block, lib = lib_segments[-i,], pred = lib_segments[i,], 
                               target_column = target_column, tp = 0, 
                               first_column_time = TRUE, 
                               columns = columns, stats_only = FALSE)
            
            preds[pred_index] <- temp[[1, "model_output"]]$pred[seq_along(pred_index)]
        }
    }
    return(preds)
}

block_lnlp_4_v <- function(block, target_column, columns)
{
    lib_segments <- matrix(NA, nrow = 4, ncol = 2)
    segment_size <- NROW(block)/4
    start_index <- 1
    for(i in 1:4)
    {
        lib_segments[i,1] <- floor(start_index)
        end_index <- start_index - 1 + segment_size
        lib_segments[i,2] <- floor(end_index)
        start_index <- end_index+1
    }
    
    pred <- rep.int(NA, times = NROW(block))
    pred_var <- rep.int(NA, times = NROW(block))    
    for(i in 1:4)
    {
        pred_index <- lib_segments[i,1]:lib_segments[i,2]
        
        temp <- block_lnlp(block, lib = lib_segments[-i,], pred = lib_segments[i,], 
                           target_column = target_column, tp = 0, 
                           first_column_time = TRUE, 
                           columns = columns, stats_only = FALSE)
        
        pred[pred_index] <- temp[[1, "model_output"]]$pred[seq_along(pred_index)]
        pred_var[pred_index] <- temp[[1, "model_output"]]$pred_var[seq_along(pred_index)]
    }
    return(cbind(pred, pred_var))
}

ricker_ret_model_4 <- function(df, p4 = 0.95, num_iter = 30000, num_burnin = 20000, 
                               model_file = "ricker_model.txt")
{
    # setup data
    pred_recruits <- vector("list", NROW(df))
    jags.params <- c("alpha", "beta")
    
    lib_segments <- matrix(NA, nrow = 4, ncol = 2)
    segment_size <- NROW(df)/4
    start_index <- 1
    for(i in 1:4)
    {
        lib_segments[i,1] <- floor(start_index)
        end_index <- start_index - 1 + segment_size
        lib_segments[i,2] <- floor(end_index)
        start_index <- end_index+1
    }
    for(i in 1:4)
    {
        # setup lib and pred
        lib <- lib_segments[-i,]
        pred <- lib_segments[i,]
        lib_index <- do.call(c, lapply(1:NROW(lib), function(x) {lib[x,1]:lib[x,2]}))
        pred_index <- pred[1]:pred[2]
        
        jags.data <- list(recruits = df$recruits[lib_index], spawners = df$spawners[lib_index])
        
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
        my_model <- jags.model(file = model_file, data = jags.data, inits = jags.inits, 
                               n.chains = length(jags.inits))
        update(my_model, n.iter = num_burnin)
        jags_output <- jags.samples(my_model, jags.params, n.iter = num_iter)
        alpha <- as.vector(jags_output[["alpha"]][1,,])
        beta <- as.vector(jags_output[["beta"]][1,,])
        
        # make prediction
        for(k in pred_index)
        {
            pred_recruits[[k]] <- df$spawners[k] * exp(alpha - beta * df$spawners[k])
        }
    }
    
    r_pred <- rep.int(NaN, NROW(df))
    for(k in 2:NROW(df))
    {
        r_pred[k] <- median(pred_recruits[[k]] * p4 + 
                                pred_recruits[[k-1]] * (1-p4))
    }
    return(r_pred)
}

ricker_ret_model_env_4 <- function(df, p4 = 0.95, num_iter = 30000, num_burnin = 20000, 
                                   model_file = "ricker_model_env.txt")
{
    # setup data
    pred_recruits <- vector("list", NROW(df))
    jags.params <- c("alpha", "beta", "g")
    
    lib_segments <- matrix(NA, nrow = 4, ncol = 2)
    segment_size <- NROW(df)/4
    start_index <- 1
    for(i in 1:4)
    {
        lib_segments[i,1] <- floor(start_index)
        end_index <- start_index - 1 + segment_size
        lib_segments[i,2] <- floor(end_index)
        start_index <- end_index+1
    }
    
    for(i in 1:4)
    {
        # setup lib and pred
        lib <- lib_segments[-i,]
        pred <- lib_segments[i,]
        lib_index <- do.call(c, lapply(1:NROW(lib), function(x) {lib[x,1]:lib[x,2]}))
        pred_index <- pred[1]:pred[2]
        
        jags.data <- list(recruits = df$recruits[lib_index], 
                          spawners = df$spawners[lib_index], 
                          env = df$env[lib_index])
        
        # get estimates for params
        min_S_index <- which(jags.data$spawners == min(jags.data$spawners))
        alpha_hat <- log(jags.data$recruits[min_S_index] / jags.data$spawners[min_S_index])
        max_R_index <- which(jags.data$recruits == max(jags.data$recruits))
        beta_hat <- jags.data$spawners[max_R_index]
        
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
                               n.chains = length(jags.inits))
        update(my_model, n.iter = num_burnin)
        jags_output <- jags.samples(my_model, jags.params, n.iter = num_iter)
        alpha <- as.vector(jags_output[["alpha"]][1,,])
        beta <- as.vector(jags_output[["beta"]][1,,])
        g <- as.vector(jags_output[["g"]][1,,])
        
        # make prediction
        for(k in pred_index)
        {
            pred_recruits[[k]] <- df$spawners[k] * exp(alpha - 
                                                           beta * df$spawners[k] + 
                                                           g * df$env[k])
        }
    }
    r_pred <- rep.int(NaN, NROW(df))
    for(k in 2:NROW(df))
    {
        r_pred[k] <- median(pred_recruits[[k]] * p4 + 
                                pred_recruits[[k-1]] * (1-p4))
    }
    return(r_pred)
}