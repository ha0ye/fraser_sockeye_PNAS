preprocess_data <- function()
{
    preprocess_stock <- function(stock_df)
    {
        n <- NROW(stock_df)
        stock_df$rec45 <- stock_df$rec4 + stock_df$rec5
        stock_df$ret <- stock_df$rec4 + c(NA, stock_df$rec5[1:(n-1)]) # age-4 and age-5 fish (aligned to rec4)
        
        temp <- normalize_by_cycle_line(stock_df$rec45)
        stock_df$rec45_n <- temp$ts
        stock_df$rec45_mu <- temp$mu
        stock_df$rec45_sigma <- temp$sigma
        
        temp <- normalize_by_cycle_line(stock_df$rec4)
        stock_df$rec4_n <- temp$ts
        stock_df$rec4_mu <- temp$mu
        stock_df$rec4_sigma <- temp$sigma
        
        temp <- normalize_by_cycle_line(stock_df$rec5)
        stock_df$rec5_n <- temp$ts
        stock_df$rec5_mu <- temp$mu
        stock_df$rec5_sigma <- temp$sigma
        
        temp <- normalize_by_cycle_line(stock_df$eff)
        stock_df$eff_n <- temp$ts
        stock_df$eff_mu <- temp$mu
        stock_df$eff_sigma <- temp$sigma
        
        return(stock_df)
    }
    
    make_block <- function(stock_df, env_data)
    {
        discharge_names <- c("D_max", "D_apr", "D_may", "D_jun")
        temp_names <- c("ET_apr", "ET_may", "ET_jun", "PT_apr", "PT_may", "PT_jun", "PT_jul")
        pdo_names <- "PDO_win"
        discharge <- normalize(env_data[, discharge_names])
        temperature <- normalize(env_data[, temp_names])
        pdo <- normalize(env_data[, pdo_names])
        
        # line up environmental data
        # lag temperature and river discharge 2 years
        desired_years <- stock_df$yr + 2
        index_in_env_data <- match(desired_years, env_data$year)
        index_in_stock_df <- 1:length(desired_years)
        
        discharge_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = NCOL(discharge)))
        discharge_cols[index_in_stock_df,] <- discharge[index_in_env_data, ]
        stock_df[, discharge_names] <- discharge_cols
        
        temp_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = NCOL(temperature)))
        temp_cols[index_in_stock_df,] <- temperature[index_in_env_data, ]
        stock_df[, temp_names] <- temp_cols
        
        # lag PDO by 1 year (winter before smolt outmigration)
        desired_years <- stock_df$yr + 1
        index_in_env_data <- match(desired_years, env_data$year)
        pdo_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = 1))
        pdo_cols[index_in_stock_df,] <- pdo[index_in_env_data]
        stock_df[, pdo_names] <- pdo_cols
        
        return(stock_df)
    }
    
    data <- read.csv("sockeye_data.csv")
    
    # filter stocks we don't want
    stock_data <- split(data, data$stk)
    stock_data <- lapply(stock_data, preprocess_stock)
    
    # add env data
    env_data <- read.csv("env_data.csv")
    block_data <- lapply(stock_data, function(stock_df) { make_block(stock_df, env_data)})
    
    # save and return
    save(block_data, file = "block_data.Rdata")
    return()
}

compute_nonlinearity_aggregated <- function()
{
    load("block_data.Rdata")
    ret <- lapply(block_data, function(x) {
        temp <- x$ret
        temp <- temp[is.finite(temp)]
        return((temp - mean(temp)) / sd(temp))
    })
    x <- c()
    lib <- matrix(NA, nrow = 9, ncol = 2)
    last <- 0
    for(i in 1:9)
    {
        x <- c(x, ret[[i]])
        lib[i,] <- c(last+1, last + length(ret[[i]]))
        last <- lib[i,2]
    }
    simplex_output <- simplex(x, lib = lib, pred = lib, E = 1:6, exclusion_radius = 0, silent = TRUE)
    E <- simplex_output$E[which.max(simplex_output$rho)]
    smap_output <- s_map(x, lib = lib, pred = lib, E = E, exclusion_radius = 0, silent = TRUE)
    theta <- smap_output$theta[which.max(smap_output$rho)]
    
    save(simplex_output, E, smap_output, theta, file = "results_nonlinear_aggregated.Rdata")
    return()
}

test_nonlinearity_aggregated <- function(num_shuffles = 500)
{
    get_smap_stats <- function(x, lib, E = NULL)
    {
        if(is.null(E))
        {
            # compute E using simplex on recruits time series
            simplex_output <- simplex(x, E = 1:8, silent = TRUE)
            best_rho_E <- simplex_output$E[which.max(simplex_output$rho)]
            best_mae_E <- simplex_output$E[which.min(simplex_output$mae)]
            E <- min(best_rho_E, best_mae_E)
        }
        
        # compute theta using s-map and E 
        smap_output <- s_map(x, lib = lib, pred = lib, E = E, silent = TRUE)
        
        best_rho <- max(smap_output$rho)
        best_mae <- min(smap_output$mae)
        return(data.frame(delta_mae = best_mae - smap_output$mae[smap_output$theta == 0]))
    }
    
    load("block_data.Rdata")
    ret <- lapply(block_data, function(x) {
        temp <- x$ret
        temp <- temp[is.finite(temp)]
        return((temp - mean(temp)) / sd(temp))
    })
    x <- c()
    lib <- matrix(NA, nrow = 9, ncol = 2)
    last <- 0
    for(i in 1:9)
    {
        x <- c(x, ret[[i]])
        lib[i,] <- c(last+1, last + length(ret[[i]]))
        last <- lib[i,2]
    }
    E <- 4
    
    cat("calculating for actual data... ", sep = "")
    start_time <- proc.time()
    actual <- get_smap_stats(x, lib, E)
    delta_mae <- actual$delta_mae
    elapsed_time <- proc.time() - start_time
    cat("(", elapsed_time[3], " sec.)\n", sep = "")
    
    # null distribution
    cat("calculating for random shuffles... ", sep = "")
    start_time <- proc.time()
    null_dist <- do.call(rbind, lapply(1:num_shuffles, function(i) {
        x_shuffle <- c()
        for(i in 1:9)
        {
            n <- length(ret[[i]])
            x_shuffle <- c(x_shuffle, ret[[i]][sample(n, n)])
        }
        return(get_smap_stats(x_shuffle, lib, E))
    }))
    
    delta_mae_p = (sum(null_dist$delta_mae < delta_mae)+1) / num_shuffles
    elapsed_time <- proc.time() - start_time
    cat("(", elapsed_time[3], " sec.)\n", sep = "")
    
    save(delta_mae = delta_mae, delta_mae_p = delta_mae_p, 
         file = "test_nonlinear_aggregated.Rdata")
    return()
}

compute_nonlinearity_stock <- function()
{
    get_smap_stats <- function(x, E = NULL)
    {
        if(is.null(E))
        {
            # compute E using simplex on recruits time series
            simplex_output <- simplex(x, E = 1:8, silent = TRUE)
            best_rho_E <- simplex_output$E[which.max(simplex_output$rho)]
            best_mae_E <- simplex_output$E[which.min(simplex_output$mae)]
            E <- min(best_rho_E, best_mae_E)
        }
        
        # compute theta using s-map and E 
        smap_output <- s_map(x, E = E, silent = TRUE)
        
        best_rho <- max(smap_output$rho)
        best_mae <- min(smap_output$mae)
        return(data.frame(delta_mae = best_mae - smap_output$mae[smap_output$theta == 0]))
    }
    
    nonlinearity_for_stock <- function(stock_df, num_shuffles = 500, max_E = 8)
    {
        x <- stock_df$ret
        x <- x[is.finite(x)]
        n <- length(x)
        
        cat("calculating for actual data for ", as.character(stock_df$stk[1]), "... ", sep = "")
        start_time <- proc.time()
        simplex_output <- simplex(x, E = 1:8, silent = TRUE)
        best_rho_E <- simplex_output$E[which.max(simplex_output$rho)]
        best_mae_E <- simplex_output$E[which.min(simplex_output$mae)]
        E <- min(best_rho_E, best_mae_E)
        
        # compute theta using s-map and E 
        smap_output <- s_map(x, E = E, silent = TRUE)
        
        best_rho <- max(smap_output$rho)
        best_mae <- min(smap_output$mae)
        theta <- smap_output$theta[which.min(smap_output$mae)]
        delta_mae <- best_mae - smap_output$mae[smap_output$theta == 0]
        elapsed_time <- proc.time() - start_time
        cat("(", elapsed_time[3], " sec.)\n", sep = "")
        
        cat("calculating for random shuffles for ", as.character(stock_df$stk[1]), "... ", sep = "")
        start_time <- proc.time()
        null_dist <- do.call(rbind, lapply(1:num_shuffles, function(i) {
            x_shuffle <- x[sample(n, n)]
            return(get_smap_stats(x_shuffle, E))
        }))
        
        delta_mae_p = (sum(null_dist$delta_mae < delta_mae)+1) / num_shuffles
        elapsed_time <- proc.time() - start_time
        cat("(", elapsed_time[3], " sec.)\n", sep = "")
        
        return(list(simplex_output = simplex_output, 
                    smap_output = smap_output, 
                    E = E, 
                    theta = theta, 
                    delta_mae = delta_mae, 
                    delta_mae_p = delta_mae_p))
    }
    
    load("block_data.Rdata")
    nonlinearity_results <- lapply(block_data, nonlinearity_for_stock)
    saveRDS(nonlinearity_results, file = "results_nonlinearity_stock.RDS")
    return()
}