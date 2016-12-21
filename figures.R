make_ricker_curve_plot <- function(posterior, stock_df, point_size = 4, title = "")
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
    
    pred <- function(x) {sapply(x, function(s) 
        median(s * exp(posterior$params$alpha - posterior$params$beta * s)))}
    pred_1 <- function(x) {sapply(x, function(s) 
        median(s * exp(posterior$first_half_params$alpha - posterior$first_half_params$beta * s)))}
    pred_2 <- function(x) {sapply(x, function(s) 
        median(s * exp(posterior$second_half_params$alpha - posterior$second_half_params$beta * s)))}
    # make plot
    #output_file <- paste("figures/ricker_curve_", stock_df$stk[1], ".pdf", sep = "")
    #pdf(file = output_file, width = 6, height = 6)
    
    label_1 <- paste(years[head(first_half, 1)], " - ", years[tail(first_half, 1)], sep = "")
    label_2 <- paste(years[head(second_half, 1)], " - ", years[tail(second_half, 1)], sep = "")
    label_3 <- paste(years[1], " - ", years[length(years)], sep = "")
    my_labels <- c(label_1, label_2, label_3, label_1, label_2)
    index <- rep.int(NA, times = length(years))
    index[first_half] <- "1"
    index[second_half] <- "2"
    df <- data.frame(year = years, spawners = spawners, recruits = recruits, half = index)
    
    my_plot <- ggplot(data = df, aes(spawners, recruits, color = half, shape = half, linetype = half)) + 
        geom_point(size = point_size) + 
        stat_function(fun = pred_1, aes(color = "1", shape = "1", linetype = "1")) + 
        stat_function(fun = pred_2, aes(color = "2", shape = "2", linetype = "2")) + 
        stat_function(fun = pred, aes(color = "all", shape = "all", linetype = "all")) + 
        scale_color_manual(values = c("royalblue", "red3", "gray"), 
                           labels = c(label_1, label_2, label_3)) + 
        scale_shape_manual(values = c(24, 25, NA), 
                           labels = c(label_1, label_2, label_3)) + 
        scale_linetype_manual(values = c(1, 1, 1), 
                              labels = c(label_1, label_2, label_3)) + 
        # theme_classic() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              legend.background = element_rect(color = "black"), 
              legend.key = element_blank(), 
              legend.title = element_blank(), 
              legend.position = c(0, 1), 
              legend.justification = c(0, 1), 
              legend.margin = unit(0.0, "cm"), 
              legend.key.size = unit(1.05, "lines"), 
              legend.key.height = unit(1.05, "lines"), 
              panel.background = element_rect(color = "black", fill = NA)) +
        ggtitle(title) +
        xlab("Spawners (millions of fish)") + 
        ylab("Recruits (millions of fish)")
    
    return(my_plot)
}

plot_seymour_ricker_halves <- function()
{    
    load("block_data.Rdata")
    stock_df <- block_data[["Seymour"]]
    
    years <- stock_df$yr
    recruits <- stock_df$rec
    spawners <- stock_df$eff
    valid <- is.finite(recruits) & is.finite(spawners)
    years <- years[valid]
    recruits <- recruits[valid]
    spawners <- spawners[valid]
    
    first_half <- 1:floor(NROW(recruits)/2)
    second_half <- (floor(NROW(recruits)/2)+1):NROW(recruits)
    
    posteriors <- readRDS("params_ricker_seymour.RDS")
    pred <- function(x) {sapply(x, function(s) 
        median(s * exp(posteriors$params$alpha - posteriors$params$beta * s)))}
    pred_1 <- function(x) {sapply(x, function(s) 
        median(s * exp(posteriors$first_half_params$alpha - posteriors$first_half_params$beta * s)))}
    pred_2 <- function(x) {sapply(x, function(s) 
        median(s * exp(posteriors$second_half_params$alpha - posteriors$second_half_params$beta * s)))}
    
    label_1 <- paste(years[head(first_half, 1)], " - ", years[tail(first_half, 1)], sep = "")
    label_2 <- paste(years[head(second_half, 1)], " - ", years[tail(second_half, 1)], sep = "")
    label_3 <- paste(years[1], " - ", years[length(years)], sep = "")
    index <- rep.int(NA, times = length(years))
    index[first_half] <- "1"
    index[second_half] <- "2"
    df <- data.frame(year = years, spawners = spawners, recruits = recruits, half = index)
    
    my_plot <- ggplot(data = df, aes(spawners, recruits, color = half, shape = half, linetype = half)) + 
        geom_point(size = 4) + 
        stat_function(fun = pred_1, aes(color = "1", shape = "1", linetype = "1")) + 
        stat_function(fun = pred_2, aes(color = "2", shape = "2", linetype = "2")) + 
        stat_function(fun = pred, aes(color = "all", shape = "all", linetype = "all")) + 
        scale_color_manual(values = c("royalblue", "red3", "gray"), 
                           labels = c(label_1, label_2, label_3)) + 
        scale_shape_manual(values = c(24, 25, NA), 
                           labels = c(label_1, label_2, label_3)) + 
        scale_linetype_manual(values = c(1, 1, 1), 
                              labels = c(label_1, label_2, label_3)) + 
        # theme_classic() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              legend.background = element_rect(color = "black"), 
              legend.key = element_blank(), 
              legend.title = element_blank(), 
              legend.position = c(0, 1), 
              legend.justification = c(0, 1), 
              legend.margin = unit(0.0, "cm"), 
              legend.key.size = unit(1.05, "lines"), 
              legend.key.height = unit(1.05, "lines"), 
              panel.background = element_rect(color = "black", fill = NA)) +
        xlab("Spawners (millions of fish)") + 
        ylab("Recruits (millions of fish)")
    
    print(my_plot)
    return()
}


plot_rho_comparison <- function()
{
    stats <- readRDS("stats_combined.RDS")
    stats <- do.call(rbind, stats)
    rhos <- t(acast(stats, stk ~ model * data, value.var = "rho"))
    rhos <- rhos[c(4, 3, 2, 1),]
    
    par(lwd = 2, mar = c(8, 4, 1, 10), xpd = TRUE)
    barplot(rhos, beside = TRUE, density = c(20, -1), space = c(0, 0.5), 
            col = c("dodgerblue", "dodgerblue", "orange1", "orange1"), 
            las = 2, xlab = "", 
            ylab = expression("Forecast Accuracy (" ~ rho ~ ")"))
    legend(42, 0.5, c("Ricker", "extended Ricker", "simple EDM", "multivariate EDM"), 
           density = c(20, -1), cex = 1.2, 
           fill = c("dodgerblue", "dodgerblue", "orange1", "orange1"), 
           col = c("dodgerblue", "dodgerblue", "orange1", "orange1"))
    
    return()
}

plot_nonlinearity <- function()
{
    load("results_nonlinear_aggregated.Rdata")
    
    par(mfrow = c(1, 2), mar = c(3, 3, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
    max_rho <- max(simplex_output$rho, na.rm = TRUE)
    min_rho <- min(simplex_output$rho, na.rm = TRUE)
    y_limits <- c(max(0, 1.1*min_rho - 0.1*max_rho), min(1.0, 1.1*max_rho - 0.1*min_rho))
    plot(simplex_output$E, pmax(0, simplex_output$rho), type = "l", 
         lwd = 1.5, ylim = y_limits, 
         xlab = "E", ylab = expression(rho))
    max_rho <- max(smap_output$rho, na.rm = TRUE)
    min_rho <- min(smap_output$rho, na.rm = TRUE)
    y_limits <- c(max(0, 1.1*min_rho - 0.1*max_rho), min(1.0, 1.1*max_rho - 0.1*min_rho))
    plot(smap_output$theta, 1 - smap_output$mae, type = "l", 
         lwd = 1.5, xlim = c(0, 4),  ylim = c(0.44, 0.49), 
         xlab = expression(theta), ylab = "1 - MAE")
    text(max(smap_output$theta), y_limits[2], 
         paste("E = ", E, sep = ""), adj = c(1, 1))
    return()
}

plot_mae_comparison <- function()
{
    stats <- readRDS("stats_combined.RDS")
    stats <- do.call(rbind, stats)
    maes <- t(acast(stats, stk ~ model * data, value.var = "mae"))
    maes <- maes[c(4, 3, 2, 1),]
    
    par(lwd = 2, mar = c(8, 4, 1, 10), xpd = TRUE)
    barplot(maes, beside = TRUE, density = c(20, -1), space = c(0, 0.5), 
            col = c("dodgerblue", "dodgerblue", "orange1", "orange1"), 
            las = 2, xlab = "", 
            ylab = "Forecast Error (MAE)")
    legend(42, 0.5, c("Ricker", "extended Ricker", "simple EDM", "multivariate EDM"), 
           density = c(20, -1), cex = 1.2, 
           fill = c("dodgerblue", "dodgerblue", "orange1", "orange1"), 
           col = c("dodgerblue", "dodgerblue", "orange1", "orange1"))
    
    return()
    return()
}

plot_chilko_smolt_model <- function()
{
    univar_stats <- readRDS("stats_univariate_EDM.RDS")
    univar_stats <- univar_stats[univar_stats$stk == "Chilko",]
    univar_stats$columns <- "eff"
    univar_stats$E <- 0
    univar_stats$data <- "spawners"
    
    multivar_stats <- do.call(rbind, readRDS("stats_multivariate_EDM.RDS"))
    multivar_stats <- multivar_stats[multivar_stats$stk == "Chilko",]
    multivar_stats$E <- sapply(strsplit(multivar_stats$columns, ", "), length) - 1
    multivar_stats <- multivar_stats[which.max(multivar_stats$rho),]
    multivar_stats$data <- "+ environment"
    
    new_stats <- readRDS("stats_chilko_smolts.RDS")
    new_stats$data <- "+ environment \n& smolts"
    new_stats$E <- sapply(strsplit(new_stats$columns, ", "), length) - 1
    new_stats <- new_stats[which.max(new_stats$rho),]
    new_stats$stk <- "Chilko"
    
    stats <- rbind(univar_stats, multivar_stats, new_stats)
    
    stats <- within(stats, data <- factor(data, levels = c("spawners", "+ environment", 
                                                           "+ environment \n& smolts")))
    
    rho_plot <- ggplot(data = stats, aes(data, rho)) + 
        geom_bar(stat = "identity", position = "dodge", width = 1, show_guide = FALSE) +
        ylab("Forecast Accuracy (rho)") + coord_cartesian(ylim = c(0, 0.5))
    mae_plot <- ggplot(data = stats, aes(data, mae)) + 
        geom_bar(stat = "identity", position = "dodge", width = 1, show_guide = FALSE) +
        ylab("Forecast Error (MAE)") + coord_cartesian(ylim = c(0.5, 0.9))
    
    plots <- lapply(list(rho_plot, mae_plot), function(my_plot) {
        return(my_plot + geom_bar(stat = "identity", position = "dodge", color = "black", fill = "orange1", width = 1, show_guide = FALSE) + 
                   xlab("Data Usage") + 
                   theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.text = element_text(color = "black"), 
                         axis.text.x = element_text(angle = 90, hjust = 0.5), 
                         legend.background = element_rect(color = "black"), 
                         legend.key = element_blank(), 
                         legend.title = element_blank(), 
                         legend.position = c(0, 1), 
                         legend.justification = c(0, 1), 
                         legend.margin = unit(0.0, "cm"), 
                         legend.key.size = unit(1.05, "lines"), 
                         legend.key.height = unit(1.05, "lines"), 
                         panel.background = element_rect(color = "black", fill = NA)))
    })
    
    plots[["nrow"]] <- 1
    
    do.call(grid.arrange, plots)
    return()
}

plot_late_shuswap_CI <- function(file = NULL, width = 6, height = 4.5)
{
    stk_name <- "Late Shuswap"
    columns <- c("eff", "D_may", "PT_jul")
    
    load("block_data.Rdata")
    env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                   "ET_apr", "ET_may", "ET_jun", 
                   "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                   "PDO_win")
    
    stock_df <- block_data[[stk_name]]
    
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
    
    rec4_preds <- data.frame(block_lnlp_4_v(block, target_column = 2, columns = columns))
    rec5_preds <- data.frame(block_lnlp_4_v(block, target_column = 3, columns = columns))
    rec4_preds$pred <- rec4_preds$pred*sigma_4 + mu_4
    rec5_preds$pred <- rec5_preds$pred*sigma_5 + mu_5
    rec4_preds$pred_var <- rec4_preds$pred_var * sigma_4 * sigma_4
    rec5_preds$pred_var <- rec5_preds$pred_var * sigma_5 * sigma_5
    
    rets <- data.frame(pred = rec4_preds$pred + c(NA, rec5_preds$pred[1:NROW(block)-1]), 
                       pred_var = rec4_preds$pred_var + c(NA, rec5_preds$pred_var[1:NROW(block)-1]))
    rets$pred_std_err <- sqrt(rets$pred_var)
    
    if(!is.null(file))
    {
        pdf(file = file, width = width, height = height)
    }
    par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
    years <- years + 4
    plot(years, returns, type = "l", 
         #         ylim = c(0, 18), 
         xlab = "Year", ylab = "Returns (millions of fish)")
    points(years, rets$pred, col = "blue", pch = 1)
    for(i in seq_along(returns))
    {
        lines(c(years[i],years[i]), rets$pred[i] + c(rets$pred_std_err[i], -rets$pred_std_err[i]), 
              col = "blue")
    }
    legend(x = "topright", legend = c("Observed", "Predicted (+/- 1 SE)"), 
           col = c("black", "blue"), lwd = c(1,NA), pch = c(NA, 1), inset = 0.02)
    if(!is.null(file))
    {
        dev.off()
    }
    
    return()
}

plot_seymour_env_surface <- function(plot_ricker = FALSE)
{
    ricker_func <- function(S, E)
    {
        return(median(S * exp(params$alpha - params$beta * S + params$g * E)))
    }
    
    meshsurf3d <- function(x, y, z, color = "black")
    {
        for(i in seq_len(grid_size))
        {
            lines3d(x[i], y, z[i,], color = color)
            lines3d(x, y[i], z[,i], color = color)
        }
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
    returns <- stock_df$ret[valid]
    env <- stock_df[valid, env_var]
    
    params <- readRDS("params_ricker_env_seymour.RDS")
    
    # make predictions for surface plots
    grid_size <- 16
    x <- seq(from = min(spawners), to = max(spawners), length.out = grid_size)
    y <- seq(from = min(env), to = max(env), length.out = grid_size)
    x_mat <- matrix(rep.int(x, times = grid_size), nrow = grid_size)
    y_mat <- matrix(rep.int(y, times = grid_size), nrow = grid_size, byrow = TRUE)
    
    # plot
    plot3d(spawners, env, recruits, size = 5, 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    axes3d(edges=c("x+-", "y--", "z++"))
    mtext3d("Spawners", "x+-", line = 1.5)
    mtext3d("Temperature", "y--", line = 1.5)
    mtext3d("Recruitment", "z++", line = 1.5)
    box3d()
    sapply(1:length(recruits), function(i) {rgl.lines(rep(spawners[i], 2), 
                                                      rep(env[i], 2), 
                                                      c(min(recruits), recruits[i]), 
                                                      color = "gray", size = 1)
    })
    if(plot_ricker)
    {
        z_ricker <- mapply(ricker_func, x_mat, y_mat)
        meshsurf3d(x, y, matrix(z_ricker, nrow = grid_size), color = "gray20")
    }
    else
    {
        spawners_n <- (spawners - mean(spawners)) / sd(spawners)
        env_n <- (env - mean(env)) / sd(env)
        x_n <- (x - mean(spawners)) / sd(spawners)
        y_n <- (y - mean(env)) / sd(env)
        block <- rbind(data.frame(recruits = recruits, 
                                  spawners = spawners_n, 
                                  env = env_n), 
                       data.frame(recruits = 1, 
                                  spawners = rep(x_n, each = grid_size), 
                                  env = rep(y_n, grid_size)))
        
        out <- block_lnlp(block, lib = c(1, length(recruits)), 
                          pred = c(length(recruits)+1, NROW(block)), 
                          tp = 0, columns = c(2, 3), stats_only = FALSE)
        z_simplex <- out[[1]]$model_output[(length(recruits)+1):NROW(block), "pred"]
        z_simplex <- matrix(z_simplex, nrow = grid_size, byrow = TRUE)
        meshsurf3d(x, y, z_simplex, color = "gray20")
    }
    p <- matrix(c(-0.8, -0.6, 0, 0, 
                  0.2, -0.3, 0.9, 0, 
                  -0.6, 0.7, 0.3, 0, 
                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(userMatrix = p)
    return()
}



plot_total_returns <- function()
{
    df <- read.csv("sockeye_ret_data.csv")
    df$cycle <- factor(df$yr %% 4)
    
    my_plot <- ggplot(data = df, aes(yr, ret, fill = cycle)) + 
        geom_bar(stat = "identity", color = "black") + 
        scale_fill_manual(values = c("white", "white", "black", "white")) + 
        scale_x_continuous(breaks = seq(1950, 2010, by = 10)) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none", 
              axis.text = element_text(color = "black"), 
              panel.background = element_rect(color = "black", fill = NA)) + 
        xlab("Year") + ylab("Returns (millions of fish)")
    
    return()
}