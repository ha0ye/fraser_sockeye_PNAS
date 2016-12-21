print_env_comparison_table <- function()
{
    stats <- readRDS("stats_combined.RDS")
    stats <- do.call(rbind, stats)
    stats <- subset(stats, data == "multivariate")
    stats_table <- data.frame(stock = stats$stk, model = stats$model, predictors = stats$columns, 
                              "num. predictions" = stats$N, rho = stats$rho, MAE = stats$mae)
    
    my_table <- xtable(stats_table, digits = 3)
    print(my_table, type = "html", file = "tables/Table_1.html")
    
    return()
}

print_nonlinearity_table <- function()
{
    nonlinear_results <- readRDS("results_nonlinearity_stock.RDS")
    temp_table <- do.call(rbind, lapply(nonlinear_results, function(res) {
        return(data.frame(E = res$E, theta = res$theta, 
                          delta_mae = res$delta_mae, p_value = res$delta_mae_p))
    }))
    temp_table$stock <- rownames(temp_table)
    temp_table <- temp_table[order(temp_table$stock), 
                             c("stock", "E", "theta", "delta_mae", "p_value")]
    temp_table$"significantly nonlinear?" <- temp_table$p_value <= 0.05
    
    my_table <- xtable(temp_table, digits = 3)
    print(my_table, type = "html", file = "tables/Table_S1.html", include.rownames = FALSE)
    
    return()    
}

print_comparison_table <- function()
{
    compute_p_values <- function(x1, x2, y)
    {
        index <- is.finite(x1) & is.finite(x2) & is.finite(y)
        x1 <- x1[index]
        x2 <- x2[index]
        y <- y[index]
        err1 <- abs(y - x1)
        err2 <- abs(y - x2)
        mae_ttest <- t.test(err1, err2, paired = TRUE, alternative = "less")
        mae_df <- mae_ttest$parameter
        mae_statistic <- mae_ttest$statistic
        mae_p <- mae_ttest$p.value
        rho_ttest <- rho_comp(x1, x2, y)
        rho_df <- rho_ttest$df
        rho_statistic <- rho_ttest$statistic
        rho_p <- rho_ttest$p.value
        return(data.frame(mae_df, mae_statistic, mae_p, rho_df, rho_statistic, rho_p))
    }
    
    preds <- readRDS("preds_combined.RDS")
    
    # normalize by mean obs value
    preds_n <- lapply(names(preds), function(stk_name) {
        df <- preds[[stk_name]]
        sigma <- sd(df$obs, na.rm = TRUE)
        mu <- mean(df$obs, na.rm = TRUE)
        df$obs <- (df$obs - sigma) / mu
        df$simplex_univar_pred <- (df$simplex_univar_pred - sigma) / mu
        df$ricker_univar_pred <- (df$ricker_univar_pred - sigma) / mu
        df$simplex_multivar_pred <- (df$simplex_multivar_pred - sigma) / mu
        df$ricker_multivar_pred <- (df$ricker_multivar_pred - sigma) / mu
        df$stk <- stk_name
        return(df)
    })
    preds_n <- do.call(rbind, preds_n)
    preds_n$stk <- factor(preds_n$stk)
    
    compare_from <- list(preds_n$simplex_univar_pred, 
                         preds_n$simplex_multivar_pred, 
                         preds_n$ricker_multivar_pred, 
                         preds_n$simplex_multivar_pred)
    compare_to <- list(preds_n$ricker_univar_pred, 
                       preds_n$ricker_multivar_pred, 
                       preds_n$ricker_univar_pred, 
                       preds_n$simplex_univar_pred)
    comparison_names <- list("simple EDM vs. Ricker", 
                             "multivariate EDM vs. extended Ricker", 
                             "extended Ricker vs. Ricker", 
                             "multivariate EDM vs. simple EDM")
    
    temp_table <- do.call(rbind, lapply(1:4, function(i) {
        temp <- compute_p_values(compare_from[[i]], compare_to[[i]], preds_n$obs)
        return(data.frame(comparison = comparison_names[[i]], 
                          performance_measure = c("rho", "MAE"), 
                          test_type = "t",
                          test_statistic = c(temp$rho_statistic, temp$mae_statistic), 
                          df = c(temp$rho_df, temp$mae_df), 
                          p_value = c(temp$rho_p, temp$mae_p)))
    }))
    my_table <- xtable(temp_table, digits = 3)
    print(my_table, type = "html", file = "tables/Table_S2.html", include.rownames = FALSE)
    
    return()
}

print_ccm_table <- function()
{
    ccm_table <- data.frame(readRDS("results_ccm.RDS"))
    ccm_table <- cbind("N" = ccm_table$N,
                       "95% p" = tanh(qnorm(0.95, sd = 1/sqrt(ccm_table$N - 3))), 
                       ccm_table[,2:NCOL(ccm_table)])
    my_table <- xtable(ccm_table, digits = 3)
    print(my_table, type = "html", file = "tables/Table_S3.html")
    return()
}

print_EDM_env_models <- function()
{
    stats <- readRDS("stats_multivariate_EDM.RDS")
    temp_table <- do.call(rbind, stats)
    temp_table <- data.frame(stock = temp_table$stk, 
                             predictors = temp_table$columns, 
                             num_predictions = temp_table$N, 
                             rho = temp_table$rho, 
                             MAE = temp_table$mae)
    my_table <- xtable(temp_table, digits = 3)
    print(my_table, type = "html", file = "tables/Table_S4.html", include.rownames = FALSE)
    
    return()
}