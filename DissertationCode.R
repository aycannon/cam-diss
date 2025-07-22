### Dissertation Code

# premable
source("preamble.R")

# data
spot_df <- read_xlsx("data/data2.xlsx", sheet = "spot")
fwd_df <- read_xlsx("data/data2.xlsx", sheet = "forward")
join <- spot_df %>%
    inner_join(fwd_df, by = "DATES")

df <- join %>%
    mutate(
     Date = as.Date(DATES),
     .keep = "unused",
     .before = "CHF"
    ) %>%
    pivot_longer(-Date,names_to = "Instrument",
                 values_to = "value") %>%
    mutate(
        Type = ifelse(grepl("1M", Instrument), "fwd", "spot")
    )

currencies <- df %>%
    filter(Type == "spot") %>%
    distinct(Instrument) %>%
    as.list()

monthly <- df %>%
    mutate(YM = floor_date(Date, "month")) %>%
    group_by(YM, Instrument) %>%
    summarise(mid = mean(value, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(
        names_from = Instrument,
        values_from = mid
    ) %>%
    mutate(
        EUR1M = 1/EUR1M,
        GBP1M = 1/GBP1M,
        .keep = "unused"
    ) %>%
    mutate(
        across(-YM, log)
    )
daily <- join %>%
    mutate(
        Date = as.Date(DATES),
        EUR1M = 1/EUR1M,
        GBP1M = 1/GBP1M,
        .keep = "unused",
        .before = "CHF"
    ) %>%
    mutate(
        across(-Date, log)
    )



# lag the spot to create excess and premium

lag <- monthly %>%
    transmute(
        YM = YM %m+% months(1), 
        across(-YM, ~ ., .names = "{.col}_t1"))
lag_daily <- daily %>%
    transmute(
        Date = Date %m+% months(1),
        across(-Date, ~ ., .names = "{.col}_t1")
    )

monthly_df <- monthly %>%
    inner_join(lag, by = "YM")

daily_df <- daily %>%
    inner_join(lag_daily, by = "Date")

spot_cols <- names(monthly)[str_detect(names(monthly), "^[A-Z]{3}$")]
forward_cols <- names(monthly)[str_detect(names(monthly), "1M$")]

forward_map <- c(
    "GBP" = "GBP1M",
    "EUR" = "EUR1M",
    "JPY" = "JPY1M",
    "CAD" = "CAD1M",
    "KRW" = "KWN1M",  
    "VND" = "VDN1M",
    "SGD" = "SGD1M",
    "CHF" = "CHF1M",
    "THB" = "THB1M"
)

monthly_df1 <- monthly_df

# creating the excess and prem variables

for (cur in spot_cols) {
    forward_name <- forward_map[[cur]]
    if (!is.null(forward_name)) {
        monthly_df1 <- monthly_df1 %>%
            mutate(
                !!paste0(cur, "_excess") := .data[[cur]] - .data[[paste0(cur, "_t1")]],
                !!paste0(cur, "_prem") := .data[[paste0(cur, "_t1")]] - .data[[paste0(forward_name, "_t1")]]
            )
    } else {
        warning("No forward mapping found for: ", cur)
    }
}

daily2 <- daily_df

for (cur in spot_cols) {
    forward_name <- forward_map[[cur]]
    if (!is.null(forward_name)) {
        daily2 <- daily2 %>%
            mutate(
                !!paste0(cur, "_excess") := .data[[cur]] - .data[[paste0(cur, "_t1")]],
                !!paste0(cur, "_prem") := .data[[paste0(cur, "_t1")]] - .data[[paste0(forward_name, "_t1")]]
            )
    } else {
        warning("No forward mapping found for: ", cur)
    }
}

# loop through the variables for the OLS Forward Premium Bias Model
lm_data <- list()

for (cur in spot_cols) {
    excess_var <- paste0(cur, "_excess")
    prem_var   <- paste0(cur, "_prem")
    
    df <- monthly_df1 %>%
        dplyr::select(all_of(c(excess_var, prem_var))) %>%
        filter(!is.na(.data[[excess_var]]), !is.na(.data[[prem_var]]))
    
    
    model <- lm(as.formula(paste0(excess_var, " ~ ", prem_var)), data = df)
    
    summary_stats <- tibble(
        Currency   = cur,
        Intercept  = coef(model)[1],
        SE.Intercept = summary(model)$coefficients[1, 2],
        Slope      = coef(model)[2],
        SE.Slope   = summary(model)$coefficients[2, 2],
        R2         = summary(model)$r.squared,
        N          = nobs(model)
    )
    
    fitted_df <- tibble(
        Currency = cur,
        Fitted   = fitted(model),
        Residual = resid(model)
    )
    
    lm_data[[cur]] <- list(summary = summary_stats, residuals = fitted_df)
    print(linearHypothesis(model, paste0(prem_var, " = 1")))
}

lmresults_df <- bind_rows(list(lm_data$CAD$summary, lm_data$CHF$summary,
                               lm_data$EUR$summary, lm_data$GBP$summary,
                               lm_data$JPY$summary, lm_data$KRW$summary,
                               lm_data$SGD$summary, lm_data$THB$summary,
                               lm_data$VND$summary))

print(lmresults_df)

# on daily data
lm_daily <- list()

for (cur in spot_cols) {
    excess_var <- paste0(cur, "_excess")
    prem_var   <- paste0(cur, "_prem")
    
    df <- daily2 %>%
        dplyr::select(all_of(c(excess_var, prem_var))) %>%
        filter(!is.na(.data[[excess_var]]), !is.na(.data[[prem_var]]))
    
    
    model <- lm(as.formula(paste0(excess_var, " ~ ", prem_var)), data = df)
    ht <- linearHypothesis(model, paste0(prem_var, " = 1"))
    
    summary_stats <- tibble(
        Currency   = cur,
        Intercept  = coef(model)[1],
        SE.Intercept = summary(model)$coefficients[1, 2],
        Slope      = coef(model)[2],
        SE.Slope   = summary(model)$coefficients[2, 2],
        R2         = summary(model)$r.squared,
        N          = nobs(model),
        p_val_HT   = ht$`Pr(>F)`[2]
    )
    
    fitted_df <- tibble(
        Currency = cur,
        Fitted   = fitted(model),
        Residual = resid(model)
    )
    
    
    
    lm_daily[[cur]] <- list(summary = summary_stats, residuals = fitted_df)
    cat("-------------", cur, "----------- \n")
    print(linearHypothesis(model, paste0(prem_var, " = 1")))
}

lm_daily_df <- bind_rows(list(lm_daily$CAD$summary, lm_daily$CHF$summary,
                              lm_daily$EUR$summary, lm_daily$GBP$summary,
                              lm_daily$JPY$summary, lm_daily$KRW$summary,
                              lm_daily$SGD$summary, lm_daily$THB$summary,
                              lm_daily$VND$summary))

lm_daily_df <- lm_daily_df %>%
    mutate(
        Sig = case_when(
            p_val_HT < 0.01 ~ "***",  # significant at 1%
            p_val_HT < 0.05 ~ "**",   # significant at 5%
            p_val_HT < 0.10 ~ "*",    # significant at 10%
            TRUE            ~ ""      # not significant
        )
    )

print(lm_daily_df)

ggplot(lm_daily_df[1:nrow(lm_daily_df)-1,], aes(x = Currency, y = Slope, color = Currency)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_errorbar(aes(ymin = Slope - 1.96 * SE.Slope,
                      ymax = Slope + 1.96 * SE.Slope),
                  width = 0.25, size = 1) +
    labs(x = "Currency", 
         y = "Slope Coefficient") +
    theme_minimal() +
    theme(legend.position = "none")

ggsave( file = "betaCoefficients.png",
        path = "plots/")


###------------------------------ GMM Alpha Estimation---------------------------

plan(multisession, workers = parallel::detectCores() - 1)

roll_alpha_gmm_parallel <- function(data, currencies, window = 60, p = 2, Inst = 0) {
    n <- nrow(data)
    date_seq <- data$YM[(window + 1):n]
    inv_logit <- function(x) 1 / (1 + exp(-x))
    
    
    # Inner function for one currency
    estimate_alpha_for_currency <- function(cur) {
        e_col <- paste0(cur, "_excess")
        f_col <- paste0(cur, "_prem")
        
        alpha_vec <- rep(NA_real_, length(date_seq))
        alpha_se_vec <- rep(NA_real_, length(date_seq))
        
        for (i in seq_along(date_seq)) {
            idx_start <- i
            idx_end <- i + window - 1
            sub_data <- data[idx_start:idx_end, ]
            
            dat_cur <- sub_data %>%
                transmute(
                    forecast_error = .data[[e_col]] - .data[[f_col]],
                    lag_error = lag(.data[[e_col]] - .data[[f_col]]),
                    lag2_error = lag(lag(.data[[e_col]] - .data[[f_col]])),
                    ALE = abs(lag(.data[[e_col]] - .data[[f_col]])),
                    const = 1
                ) %>% na.omit()
            
            Inst_local <- Inst
            
            if (nrow(dat_cur) > 10) {
                gmm_moments <- function(theta, data, Inst) {
                    alpha <- theta[1]
                    e <- data$forecast_error
                    z <- switch(
                        as.character(Inst_local),
                        "0" = matrix(1, nrow = length(e), ncol = 1),
                        "1" = as.matrix(data[, c("const", "lag_error")]),
                        "2" = as.matrix(data[, c("const", "lag_error", "lag2_error")]),
                        "3" = as.matrix(data[, c("const", "ALE")]),
                        stop("Invalid instrument set. Choose Inst = 0, 1, 2, or 3.")
                    )
                    lambda <- (e < 0) - alpha
                    g <- z * lambda * abs(e)^(p - 1)
                    return(g)
                }
                
                res <- tryCatch({
                    gmm(g = gmm_moments, x = dat_cur, t0 = c(0),  # t0 is on logit scale
                        type = "cue", method = "BFGS")
                }, error = function(e) NULL)
                
                if (!is.null(res)) {
                    alpha_logit     <- coef(res)[1]
                    alpha_hat       <- inv_logit(alpha_logit)
                    
                    # Delta method: SE(logit) × derivative of inverse logit
                    vcov_mat <- tryCatch(vcov(res), error = function(e) NA)
                    
                    if (!anyNA(vcov_mat)) {
                        se_logit     <- sqrt(vcov_mat[1, 1])
                        se_alpha     <- se_logit * alpha_hat * (1 - alpha_hat)
                    } else {
                        se_alpha     <- NA
                    }
                    
                    alpha_vec[i]    <- alpha_hat
                    alpha_se_vec[i] <- se_alpha
                } else {
                    alpha_vec[i]    <- NA
                    alpha_se_vec[i] <- NA
                }
                
            }
        }
        
        return(list(alpha = alpha_vec, se = alpha_se_vec))
    }
    
    # Run in parallel
    alpha_list       <- future_map(currencies, estimate_alpha_for_currency, .progress = TRUE)
    names(alpha_list) <- currencies
    
    alpha_matrix     <- do.call(cbind, lapply(alpha_list, `[[`, "alpha"))
    alpha_se_matrix  <- do.call(cbind, lapply(alpha_list, `[[`, "se"))
    rownames(alpha_matrix)    <- as.character(date_seq)
    rownames(alpha_se_matrix) <- as.character(date_seq)
    
    return(list(alpha = alpha_matrix, se = alpha_se_matrix))
}

alpha_results <- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 0)
alpha_results_i1<- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 1)
alpha_results_i2<- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 2)
alpha_results_i3<- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 3)

## we want to remove Vietnamese Dong as there are so many NAs

alpha_matrix <- alpha_results$alpha[,1:ncol(alpha_results$alpha)-1]
alpha_matrix_i1 <- alpha_results_i1$alpha[,1:ncol(alpha_results_i1$alpha)-1]
alpha_matrix_i2 <- alpha_results_i2$alpha[,1:ncol(alpha_results_i2$alpha)-1]
alpha_matrix_i3 <- alpha_results_i3$alpha[,1:ncol(alpha_results_i3$alpha)-1]

alpha_se <- alpha_results$se[,1:ncol(alpha_results$se)-1]
alpha_se_i1 <- alpha_results_i1$se[,1:ncol(alpha_results_i1$se)-1]
alpha_se_i2 <- alpha_results_i2$se[,1:ncol(alpha_results_i2$se)-1]
alpha_se_i3 <- alpha_results_i3$se[,1:ncol(alpha_results_i3$se)-1]

## Plotting the alphas

recessions <- data.frame(
    start = as.Date(c("2007-12-01", "2020-02-01")),
    end   = as.Date(c("2009-06-01", "2022-04-01")),
    label = c("GFC", "COVID")
)

key_dates <- data.frame(
    start = as.Date(c("2014-06-01")),
    end = as.Date(c("2016-02-01")),
    label = c("Oil Crash")
)


smooth_alpha <- as.data.frame(alpha_matrix) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date))

alpha_ewma <- smooth_alpha %>%
    mutate(across(-Date, ~ EMA(., n = 10)))
alpha_sma <- smooth_alpha %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))

alpha_vals <- alpha_sma %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    pull(Alpha)

alpha_min <- min(alpha_vals, na.rm = TRUE)
alpha_max <- max(alpha_vals, na.rm = TRUE)

smooth_alpha_i1 <- as.data.frame(alpha_matrix_i1) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date))
alpha_ewma_i1 <- smooth_alpha_i1 %>%
    mutate(across(-Date, ~ EMA(., n = 10)))
alpha_sma_i1 <- smooth_alpha_i1 %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))

smooth_alpha_i2 <- as.data.frame(alpha_matrix_i2) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date))
alpha_ewma_i2 <- smooth_alpha_i2 %>%
    mutate(across(-Date, ~ EMA(., n = 10)))
alpha_sma_i2 <- smooth_alpha_i2 %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))

smooth_alpha_i3 <- as.data.frame(alpha_matrix_i3) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date))
alpha_ewma_i3 <- smooth_alpha_i3 %>%
    mutate(across(-Date, ~ EMA(., n = 10)))
alpha_sma_i3 <- smooth_alpha_i3 %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))


ggplot(smooth_alpha %>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())

    ggsave("alpha time series.png",
           path = "plots/")

ggplot(alpha_ewma%>% pivot_longer(-Date, names_to = "Currency",
                                  values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha EWMA time series.png",
           path = "plots/")
    
    ggplot(alpha_sma %>% pivot_longer(-Date, names_to = "Currency",
                                            values_to = "Alpha"),
           aes(x = Date, y = Alpha, color = Currency)) +
        geom_line() +
        labs(x = "Date",
             y = "Alpha") +
        theme_minimal() +
        theme(legend.position = "bottom", 
              legend.title = element_blank())
    ggsave("alpha SMA time series.png",
           path = "plots/")

ggplot(alpha_sma%>% pivot_longer(-Date, names_to = "Currency",
                                 values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_rect(data = recessions, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0.57, ymax = 0.69),
              fill = "darkgrey", alpha = 0.3) +
    geom_rect(data = key_dates, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0.57, ymax = 0.69),
              fill = "salmon", alpha = 0.3) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    annotate("text", x = as.Date("2008-01-11"), y = .62, label = "Great Recession",
             angle=90, vjust=-.5, size=4, alpha=0.6) +
    annotate("text", x = as.Date("2020-03-11"), y = .66, label = "Pandemic Announced", 
             angle=90, vjust=-.5, size=4, alpha=0.6) +
    annotate("text", x = as.Date("2014-07-11"), y = .665, label = "Oil Crash",
             angle=90, vjust=-.5, size=4, alpha=0.6) +
    theme(legend.position = "none")

    ggsave("alpha SMA w Recession time series.png",
           path = "plots/")
    
### smoothing the standard errors for plotting:
    
alpha_se_sma <- as.data.frame(alpha_se) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(Date) %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))

plotting_alpha_se_df <- alpha_sma %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    left_join(alpha_se_sma %>% pivot_longer(-Date, names_to = "Currency", values_to = "SE"),
              by = c("Date", "Currency")) %>%
    mutate(
        Lower = Alpha - 1.96 * SE,
        Upper = Alpha + 1.96 * SE
    )

ggplot(plotting_alpha_se_df, aes(x = Date, y = Alpha, color = Currency, fill = Currency)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1, color = NA) +
    geom_line() +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    labs(
         x = "Date", y = expression(alpha)) +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave("alpha with SE.png",
       path = "plots/",
       width = 10, height = 6)

plotting_alpha_se_df %>%
    summarise(
        prop_0.5 = mean(Lower < 0.5, na.rm = TRUE),
    )


# plotting alpha i1

ggplot(smooth_alpha_i1 %>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 1 time series.png",
           path = "plots/")

ggplot(alpha_ewma_i1%>% pivot_longer(-Date, names_to = "Currency",
                                  values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 1 EWMA time series.png",
           path = "plots/")

ggplot(alpha_sma_i1%>% pivot_longer(-Date, names_to = "Currency",
                                 values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 1 SMA time series.png",
           path = "plots/")

# plottiong alpha i2

ggplot(smooth_alpha_i2 %>% pivot_longer(-Date, names_to = "Currency",
                                        values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 2 time series.png",
           path = "plots/")

ggplot(alpha_ewma_i2%>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 2 EWMA time series.png",
           path = "plots/")

ggplot(alpha_sma_i2%>% pivot_longer(-Date, names_to = "Currency",
                                    values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 2 SMA time series.png",
           path = "plots/")

# plotting i3

ggplot(smooth_alpha_i3 %>% pivot_longer(-Date, names_to = "Currency",
                                        values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 3 time series.png",
           path = "plots/")

ggplot(alpha_ewma_i3%>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 3 EWMA time series.png",
           path = "plots/")

ggplot(alpha_sma_i3%>% pivot_longer(-Date, names_to = "Currency",
                                    values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha Inst 3 SMA time series.png",
           path = "plots/")

## --- Portfolio Allocation: No instruments ------------------------------------ 

### -------- Creating Weights ----
alpha_long <- as.data.frame(alpha_matrix) %>%
    rownames_to_column("Date") %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
           Alpha = as.numeric(Alpha),
           Deviation = abs(Alpha - 0.5),                    
           Direction = ifelse(Alpha < 0.5, "Long", "Short"),
           Delta = Alpha - dplyr::lag(Alpha)
    )

alpha_weights <- alpha_long %>%
    group_by(Date) %>%
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n / 2) ~ +1 / floor(n / 2),  # Long lowest alphas
            rank >  floor(n / 2) ~ -1 / ceiling(n / 2), # Short highest alphas
            TRUE                 ~ 0
        )
    ) %>%
    ungroup()

wide_smooth_ewma <- alpha_ewma %>%  
    arrange(Date) %>%
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ EMA(.x, n = 10),
            .names = "{.col}_EWMA"
        )
    ) %>%
    group_by() %>% # if you ever have multiple symbols, do group_by(symbol)
    mutate(
        across(
            .cols = ends_with("_EWMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

wide_smooth_sma <- alpha_sma %>%     # or whatever your SMA table is
    arrange(Date) %>%
    # 1a) compute a 10‐month EWMA of *each* currency's alpha
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ SMA(.x, n = 10),
            .names = "{.col}_SMA"
        )
    ) %>%
    group_by() %>% 
    mutate(
        across(
            .cols = ends_with("_SMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

alpha_ewma_long <- wide_smooth_ewma %>%
    pivot_longer(
        cols      = ends_with("_EWMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_EWMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n / 2) ~ +1 / floor(n / 2),  # Long lowest alphas
            rank >  floor(n / 2) ~ -1 / ceiling(n / 2), # Short highest alphas
            TRUE                 ~ 0
        )
    ) %>%
    ungroup()

alpha_sma_long <- wide_smooth_sma %>%
    pivot_longer(
        cols      = ends_with("_SMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_SMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n / 2) ~ +1 / floor(n / 2),  # Long lowest alphas
            rank >  floor(n / 2) ~ -1 / ceiling(n / 2), # Short highest alphas
            TRUE                 ~ 0
        )
    ) %>%
    ungroup()



## try to align the start dates so that we can compare the strategies fairly
start_dates <- list(
    NS = min(alpha_weights$Date, na.rm = TRUE),
    EWMA = alpha_ewma_long %>% filter(!is.na(Alpha)) %>% summarise(min(Date)) %>% pull(),
    SMA  = alpha_sma_long  %>% filter(!is.na(Alpha)) %>% summarise(min(Date)) %>% pull(),
    EWMA_d = alpha_ewma_long %>% filter(!is.na(w_d)) %>% summarise(min(Date)) %>% pull(),
    SMA_d  = alpha_sma_long  %>% filter(!is.na(w_d)) %>% summarise(min(Date)) %>% pull()
)

common_start <- max(unlist(start_dates))
alpha_weights_filtered <- alpha_weights %>% filter(Date >= common_start)

### -------- Calculating Returns of simple strats ----

returns_long <- monthly_df1 %>%
    mutate(Date = YM) %>%
    dplyr::select(Date, ends_with("_excess")) %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "ExcessReturn") %>%
    mutate(Currency = sub("_excess", "", Currency))

alpha_portfolios <- alpha_weights_filtered %>%
    left_join(returns_long, by = c("Date", "Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )


## plotting simple strat - non-smoothed
ggplot(alpha_portfolios, aes(x = Date, y = Cume_a)) +
    geom_line() +
    labs(title = "Cumulative Returns of Alpha-Based Portfolio",
         x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom") 


### smoothed data

a_ewma_portfolios <- alpha_ewma_long %>% filter(Date >= common_start) %>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    ) 

a_sma_portfolios <- alpha_sma_long %>% filter(Date >= common_start)%>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE)
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    ) 

## plotting all strats

non_smoothed_long <- alpha_portfolios %>%
    mutate(Source = "Non-smoothed") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

ewma_long <- a_ewma_portfolios %>%
    mutate(Source = "EWMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

sma_long <- a_sma_portfolios %>%
    mutate(Source = "SMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)


# Combine all into one data frame
combined_cume <- bind_rows(non_smoothed_long, ewma_long, sma_long) %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )


# Plot
ggplot(combined_cume, aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Source)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Date", 
         y = "Cumulative Return") +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(2, "cm"))

ggsave(file = "allRegionsAllStrats.png",
           path = "plots/")

#### ------- Geographical Strategies ----------------------------

west <- c("CAD", "EUR", "GBP", "CHF")
east <- c("JPY", "SGD", "KRW", "THB")
g7 <- c("CAD", "EUR", "JPY", "GBP")
all <- c("CAD", "EUR", "JPY", "GBP", "CHF", "SGD", "KRW", "THB")

compute_strategy <- function(wide_data, returns_long, cols, suffix = NULL) {
    suffix_pattern <- if (!is.null(suffix)) paste0("_", suffix) else ""
    
    alpha_long <- wide_data %>% filter(Date >= common_start) %>%
        pivot_longer(
            cols = any_of(cols),
            names_to = "Currency",
            values_to = "Alpha"
        ) %>%
        group_by(Currency) %>%
        arrange(Date) %>%
        mutate(
            Delta = Alpha - lag(Alpha)
        ) %>%
        ungroup() %>%
        mutate(
            dev   = abs(Alpha - 0.5),
            dir   = if_else(Alpha >  0.5, -1, +1),
            raw_a = dir * dev
        ) %>%
        group_by(Date) %>%
        mutate(
            w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
            n_pos = sum(Delta >  0, na.rm = TRUE),
            n_neg = sum(Delta <  0, na.rm = TRUE),
            w_d   = case_when(
                Delta >  0 ~ -1/n_pos,
                Delta <  0 ~ +1/n_neg,
                TRUE       ~ 0
            )
        ) %>% 
        ungroup() %>%
        group_by(Date) %>%
        arrange(Date, Alpha) %>%
        mutate(
            rank = row_number(),
            n = n(),
            w_r = case_when(
                rank <= floor(n/2) ~ +1/floor(n/2),
                rank > floor(n/2) ~ -1/ceiling(n/2),
                T ~ 0
            )
            ) %>%
        ungroup()
    
    alpha_long %>%
        left_join(returns_long, by = c("Date", "Currency")) %>%
        group_by(Date) %>%
        summarise(
            Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
            Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
            Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            Cume_a = cumprod(1 + Ret_a) - 1,
            Cume_d = cumprod(1 + Ret_d) - 1,
            Cume_r = cumprod(1 + Ret_r) - 1,
            Date   = as.Date(Date)
        )
}
compute_strategy2 <- function(wide_data, returns_long, cols, suffix = NULL) {
    suffix_pattern <- if (!is.null(suffix)) paste0("_", suffix) else ""
    
    alpha_long <- wide_data %>%
        filter(Date >= common_start) %>%
        pivot_longer(
            cols = any_of(cols),
            names_to = "Currency",
            values_to = "Alpha"
        ) %>%
        group_by(Currency) %>%
        arrange(Date) %>%
        mutate(
            Delta = Alpha - lag(Alpha)
        ) %>%
        ungroup() %>%
        mutate(
            dev   = abs(Alpha - 0.5),
            dir   = if_else(Alpha > 0.5, -1, +1),
            raw_a = dir * dev
        ) %>%
        group_by(Date) %>%
        mutate(
            w_a = raw_a / sum(abs(raw_a), na.rm = TRUE),
            n_pos = sum(Delta > 0, na.rm = TRUE),
            n_neg = sum(Delta < 0, na.rm = TRUE),
            w_d = case_when(
                Delta > 0 ~ -1 / n_pos,
                Delta < 0 ~ +1 / n_neg,
                TRUE ~ 0
            )
        ) %>%
        ungroup() %>%
        group_by(Date) %>%
        arrange(Date, Alpha) %>%
        mutate(
            rank = row_number(),
            n = n(),
            w_r = case_when(
                rank <= floor(n / 2) ~ +1 / floor(n / 2),
                rank >  floor(n / 2) ~ -1 / ceiling(n / 2),
                TRUE ~ 0
            )
        ) %>%
        ungroup()
    
    # Join returns
    full_df <- alpha_long %>%
        left_join(returns_long, by = c("Date", "Currency"))
    
    # Portfolio returns
    portfolio_returns <- full_df %>%
        group_by(Date) %>%
        summarise(
            Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
            Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
            Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            Cume_a = cumprod(1 + Ret_a) - 1,
            Cume_d = cumprod(1 + Ret_d) - 1,
            Cume_r = cumprod(1 + Ret_r) - 1,
            Date = as.Date(Date)
        )
    
    return(list(
        weights = alpha_long %>% select(Date, Currency, Alpha, w_a, w_d, w_r),
        returns = portfolio_returns
    ))
}



### respecifying the non-smoothed alpha data

wide_non_smooth <- alpha_weights %>%
    dplyr::select(Date, Currency, Alpha) %>%
    pivot_wider(names_from = Currency, values_from = Alpha)

WEST_EWMA <- compute_strategy(wide_smooth_ewma, 
                              returns_long, 
                              suffix = "_EWMA", 
                              cols = west)
WEST_SMA <- compute_strategy(wide_smooth_sma, 
                              returns_long, 
                              suffix = "_SMA", 
                              cols = west)
WEST_NS <- compute_strategy(wide_non_smooth, 
                             returns_long,
                             cols = west)

ASIA_EWMA <- compute_strategy(wide_smooth_ewma, 
                              returns_long, 
                              suffix = "_EWMA", 
                              cols = east)
ASIA_SMA <- compute_strategy(wide_smooth_sma, 
                             returns_long, 
                             suffix = "_SMA", 
                             cols = east)
ASIA_NS <- compute_strategy(wide_non_smooth, 
                            returns_long,
                            cols = east)

ASIA_NS.w <- compute_strategy2(wide_non_smooth, 
                            returns_long,
                            cols = east)
G7_EWMA <- compute_strategy(wide_smooth_ewma, 
                              returns_long, 
                              suffix = "_EWMA", 
                              cols = g7)
G7_SMA <- compute_strategy(wide_smooth_sma, 
                             returns_long, 
                             suffix = "_SMA", 
                             cols = g7)
G7_NS <- compute_strategy(wide_non_smooth, 
                            returns_long,
                            cols = g7)
G7_NS.w <- compute_strategy2(wide_non_smooth, 
                          returns_long,
                          cols = g7)
ALL_EWMA <- compute_strategy(wide_smooth_ewma, 
                            returns_long, 
                            suffix = "_EWMA", 
                            cols = all)
ALL_SMA <- compute_strategy(wide_smooth_sma, 
                           returns_long, 
                           suffix = "_SMA", 
                           cols = all)
ALL_NS <- compute_strategy(wide_non_smooth, 
                          returns_long,
                          cols = all)
 
### combining returns

strat_dfs <- list(
    G7_NS = G7_NS,
    G7_SMA = G7_SMA,
    G7_EWMA = G7_EWMA,
    ASIA_NS = ASIA_NS,
    ASIA_SMA = ASIA_SMA,
    ASIA_EWMA = ASIA_EWMA,
    WEST_NS = WEST_NS,
    WEST_SMA = WEST_SMA,
    WEST_EWMA = WEST_EWMA,
    ALL_NS = ALL_NS,
    ALL_SMA = ALL_SMA,
    ALL_EWMA = ALL_EWMA
)

combined_geo_strats <- bind_rows(strat_dfs, .id = "Source") %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    separate(Source, into = c("Region", "Smoothing"), sep = "_") %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )

combined_geo_strats %>%
    filter(Region == "G7") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "G7I0.png",
       path = "plots/")

combined_geo_strats %>%
    filter(Region == "ASIA") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "ASIAI0.png",
       path = "plots/")

combined_geo_strats %>%
    filter(Region == "WEST") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "WESTI0.png",
       path = "plots/")

# using peerperformance package

wide_perf <- bind_rows(strat_dfs, .id = "Source") %>%
    pivot_longer(cols = c(Ret_a, Ret_d, Ret_r), names_to = "Type", values_to = "Return") %>%
    mutate(Strategy = paste(Source, Type, sep = "_")) %>%
    select(Date, Strategy, Return) %>%
    pivot_wider(names_from = Strategy, values_from = Return) %>%
    arrange(Date)

ret_xts <- xts(wide_perf[,-1], order.by = wide_perf$Date)

alphaScreen <- alphaScreening(X = ret_xts, control = c("nCore" = parallel::detectCores() - 1))
sharpeScreen <- sharpeScreening(X = ret_xts, control = c("nCore" = parallel::detectCores() - 1))
msharpeScreen <- msharpeScreening(X = ret_xts, na.neg = F, control = c("nCore" = parallel::detectCores() - 1))

Sharpe <- sharpe(wide_perf[,-1])
mSharpe <- msharpe(wide_perf[,-1])

## plotting weights

alpha_ewma_long %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs( y = "Weight", x = "Date") +
    theme(legend.position = "bottom")

ggsave(file = "WeightsEWMAI0.png",
       path = "plots/")

G7_NS.w$weights %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weights") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weights, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")

ggsave(file = "WeightsG7I0.png",
       path = "plots/")

ASIA_NS.w$weights %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weights") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weights, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")

ggsave(file = "WeightsASIAI0.png",
       path = "plots/")

alpha_sma_long %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")

ggsave(file = "WeightsSMAI0.png",
       path = "plots/")

alpha_weights %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")

ggsave(file = "WeightsNSI0.png",
       path = "plots/")

## --- Portfolio Allocation: Instrument 1 (constant + lagged Error) ----

### -------- Creating Weights ----
alpha_long_i1 <- as.data.frame(alpha_matrix_i1) %>%
    rownames_to_column("Date") %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
           Alpha = as.numeric(Alpha),
           Deviation = abs(Alpha - 0.5),                    
           Direction = ifelse(Alpha < 0.5, "Long", "Short"),
           Delta = Alpha - dplyr::lag(Alpha)
    )

alpha_weights_i1 <- alpha_long_i1 %>%
    group_by(Date) %>%
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()

wide_smooth_ewma_i1 <- alpha_ewma_i1 %>%  
    arrange(Date) %>%
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ EMA(.x, n = 10),
            .names = "{.col}_EWMA"
        )
    ) %>%
    group_by() %>% # if you ever have multiple symbols, do group_by(symbol)
    mutate(
        across(
            .cols = ends_with("_EWMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

wide_smooth_sma_i1 <- alpha_sma_i1 %>%     # or whatever your SMA table is
    arrange(Date) %>%
    # 1a) compute a 10‐month EWMA of *each* currency's alpha
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ SMA(.x, n = 10),
            .names = "{.col}_SMA"
        )
    ) %>%
    group_by() %>% 
    mutate(
        across(
            .cols = ends_with("_SMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

alpha_ewma_long_i1 <- wide_smooth_ewma_i1 %>%
    pivot_longer(
        cols      = ends_with("_EWMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_EWMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()

alpha_sma_long_i1 <- wide_smooth_sma_i1 %>%
    pivot_longer(
        cols      = ends_with("_SMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_SMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()


alpha_weights_filtered_i1 <- alpha_weights_i1 %>% filter(Date >= common_start)

### -------- Calculating Returns of simple strats ----

alpha_portfolios_i1 <- alpha_weights_filtered_i1 %>%
    left_join(returns_long, by = c("Date", "Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )


## plotting simple strat - non-smoothed
ggplot(alpha_portfolios_i1, aes(x = Date, y = Cume_a)) +
    geom_line() +
    labs(title = "Cumulative Returns of Alpha-Based Portfolio",
         x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom") 


### smoothed data

a_ewma_portfolios_i1 <- alpha_ewma_long_i1 %>% filter(Date >= common_start) %>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    ) 

a_sma_portfolios_i1 <- alpha_sma_long_i1 %>% filter(Date >= common_start)%>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )

## plotting all strats

non_smoothed_long_i1 <- alpha_portfolios_i1 %>%
    mutate(Source = "Non-smoothed") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

ewma_long_i1 <- a_ewma_portfolios_i1 %>%
    mutate(Source = "EWMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

sma_long_i1 <- a_sma_portfolios_i1 %>%
    mutate(Source = "SMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)


# Combine all into one data frame
combined_cume_i1 <- bind_rows(non_smoothed_long_i1, ewma_long_i1, sma_long_i1) %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )


# Plot
ggplot(combined_cume_i1, aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Source)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Date", 
         y = "Cumulative Return") +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(2, "cm"))
    ggsave("allRegionsAllStratsi1.png",
           path = "plots/")

#### ------- Geographical Strategies ----------------------------

wide_non_smooth_i1 <- alpha_weights_i1 %>%
    dplyr::select(Date, Currency, Alpha) %>%
    pivot_wider(names_from = Currency, values_from = Alpha)

WEST_EWMA_i1 <- compute_strategy(wide_smooth_ewma_i1, 
                                 returns_long, 
                                 suffix = "_EWMA", 
                                 cols = west)
WEST_SMA_i1 <- compute_strategy(wide_smooth_sma_i1, 
                                returns_long, 
                                suffix = "_SMA", 
                                cols = west)
WEST_NS_i1 <- compute_strategy(wide_non_smooth_i1, 
                               returns_long,
                               cols = west)

ASIA_EWMA_i1 <- compute_strategy(wide_smooth_ewma_i1, 
                                 returns_long, 
                                 suffix = "_EWMA", 
                                 cols = east)
ASIA_SMA_i1 <- compute_strategy(wide_smooth_sma_i1, 
                                returns_long, 
                                suffix = "_SMA", 
                                cols = east)
ASIA_NS_i1 <- compute_strategy(wide_non_smooth_i1, 
                               returns_long,
                               cols = east)
G7_EWMA_i1 <- compute_strategy(wide_smooth_ewma_i1, 
                               returns_long, 
                               suffix = "_EWMA", 
                               cols = g7)
G7_SMA_i1 <- compute_strategy(wide_smooth_sma_i1, 
                              returns_long, 
                              suffix = "_SMA", 
                              cols = g7)
G7_NS_i1 <- compute_strategy(wide_non_smooth_i1, 
                             returns_long,
                             cols = g7)

### combining returns

strat_dfs_i1 <- list(
    G7_NS_i1 = G7_NS_i1,
    G7_SMA_i1 = G7_SMA_i1,
    G7_EWMA_i1 = G7_EWMA_i1,
    ASIA_NS_i1 = ASIA_NS_i1,
    ASIA_SMA_i1 = ASIA_SMA_i1,
    ASIA_EWMA_i1 = ASIA_EWMA_i1,
    WEST_NS_i1 = WEST_NS_i1,
    WEST_SMA_i1 = WEST_SMA_i1,
    WEST_EWMA_i1 = WEST_EWMA_i1
)

combined_geo_strats_i1 <- bind_rows(strat_dfs_i1, .id = "Source") %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    separate(Source, into = c("Region", "Smoothing"), sep = "_", extra = "merge") %>%
    mutate(
        Smoothing = case_when(
            Smoothing == "NS_i1" ~ "Non-Smoothed",
            Smoothing == "EWMA_i1" ~ "EWMA",
            Smoothing == "SMA_i1" ~ "SMA",
            TRUE ~ Smoothing
        )
    ) %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )

combined_geo_strats_i1 %>%
    filter(Region == "G7") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom") 
ggsave(file = "G7I1.png",
       path = "plots/")

combined_geo_strats_i1 %>%
    filter(Region == "ASIA") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "ASIAI1.png",
       path = "plots/")

combined_geo_strats_i1 %>%
    filter(Region == "WEST") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "WESTI1.png",
       path = "plots/")

## performance metrics

wide_perf_i1 <- bind_rows(strat_dfs_i1, .id = "Source") %>%
    pivot_longer(cols = c(Ret_a, Ret_d, Ret_r), names_to = "Type", values_to = "Return") %>%
    mutate(Strategy = paste(Source, Type, sep = "_")) %>%
    select(Date, Strategy, Return) %>%
    pivot_wider(names_from = Strategy, values_from = Return) %>%
    arrange(Date)

ret_xts_i1 <- xts(wide_perf_i1[,-1], order.by = wide_perf_i1$Date)

alphaScreen_i1 <- alphaScreening(X = ret_xts_i1, control = c("nCore" = parallel::detectCores() - 1))
sharpeScreen_i1 <- sharpeScreening(X = ret_xts_i1, control = c("nCore" = parallel::detectCores() - 1))
msharpeScreen_i1 <- msharpeScreening(X = ret_xts_i1, control = c("nCore" = parallel::detectCores() - 1))

Sharpe_i1 <- sharpe(wide_perf_i1[,-1])
mSharpe_i1 <- msharpe(wide_perf_i1[,-1])


## plotting weights
alpha_sma_long_i1 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs( y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsSMAI1.png",
       path = "plots/")

alpha_ewma_long_i1 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsEWMAI1.png",
       path = "plots/")

alpha_weights_i1 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsNSI1.png",
       path = "plots/")

## --- Portfolio Allocation: Instrument 2 (constant + 2 lagged Errors) ---------
### -------- Creating Weights ----
alpha_long_i2 <- as.data.frame(alpha_matrix_i2) %>%
    rownames_to_column("Date") %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
           Alpha = as.numeric(Alpha),
           Deviation = abs(Alpha - 0.5),                    
           Direction = ifelse(Alpha < 0.5, "Long", "Short"),
           Delta = Alpha - dplyr::lag(Alpha)
    )

alpha_weights_i2 <- alpha_long_i2 %>%
    group_by(Date) %>%
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()

wide_smooth_ewma_i2 <- alpha_ewma_i2 %>%  
    arrange(Date) %>%
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ EMA(.x, n = 10),
            .names = "{.col}_EWMA"
        )
    ) %>%
    group_by() %>% # if you ever have multiple symbols, do group_by(symbol)
    mutate(
        across(
            .cols = ends_with("_EWMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

wide_smooth_sma_i2 <- alpha_sma_i2 %>%     # or whatever your SMA table is
    arrange(Date) %>%
    # 1a) compute a 10‐month EWMA of *each* currency's alpha
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ SMA(.x, n = 10),
            .names = "{.col}_SMA"
        )
    ) %>%
    group_by() %>% 
    mutate(
        across(
            .cols = ends_with("_SMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

alpha_ewma_long_i2 <- wide_smooth_ewma_i2 %>%
    pivot_longer(
        cols      = ends_with("_EWMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_EWMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()

alpha_sma_long_i2 <- wide_smooth_sma_i2 %>%
    pivot_longer(
        cols      = ends_with("_SMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_SMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()


alpha_weights_filtered_i2 <- alpha_weights_i2 %>% filter(Date >= common_start)

### -------- Calculating Returns of simple strats ----

alpha_portfolios_i2 <- alpha_weights_filtered_i2 %>%
    left_join(returns_long, by = c("Date", "Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )


## plotting simple strat - non-smoothed
ggplot(alpha_portfolios_i2, aes(x = Date, y = Cume_a)) +
    geom_line() +
    labs(title = "Cumulative Returns of Alpha-Based Portfolio",
         x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_colour_viridis_d(option = "D") 


### smoothed data

a_ewma_portfolios_i2 <- alpha_ewma_long_i2 %>% filter(Date >= common_start) %>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    ) 

a_sma_portfolios_i2 <- alpha_sma_long_i2 %>% filter(Date >= common_start)%>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )

## plotting all strats

non_smoothed_long_i2 <- alpha_portfolios_i2 %>%
    mutate(Source = "Non-smoothed") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

ewma_long_i2 <- a_ewma_portfolios_i2 %>%
    mutate(Source = "EWMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

sma_long_i2 <- a_sma_portfolios_i2 %>%
    mutate(Source = "SMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)


# Combine all into one data frame
combined_cume_i2 <- bind_rows(non_smoothed_long_i2, ewma_long_i2, sma_long_i2) %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )


# Plot
ggplot(combined_cume_i2, aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Source)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Date", 
         y = "Cumulative Return") +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(2, "cm"))
    ggsave("allRegionsAllStratsi2.png",
           path = "plots/")

#### ------- Geographical Strategies ----------------------------

wide_non_smooth_i2 <- alpha_weights_i2 %>%
    dplyr::select(Date, Currency, Alpha) %>%
    pivot_wider(names_from = Currency, values_from = Alpha)

WEST_EWMA_i2 <- compute_strategy(wide_smooth_ewma_i2, 
                                 returns_long, 
                                 suffix = "_EWMA", 
                                 cols = west)
WEST_SMA_i2 <- compute_strategy(wide_smooth_sma_i2, 
                                returns_long, 
                                suffix = "_SMA", 
                                cols = west)
WEST_NS_i2 <- compute_strategy(wide_non_smooth_i2, 
                               returns_long,
                               cols = west)

ASIA_EWMA_i2 <- compute_strategy(wide_smooth_ewma_i2, 
                                 returns_long, 
                                 suffix = "_EWMA", 
                                 cols = east)
ASIA_SMA_i2 <- compute_strategy(wide_smooth_sma_i2, 
                                returns_long, 
                                suffix = "_SMA", 
                                cols = east)
ASIA_NS_i2 <- compute_strategy(wide_non_smooth_i2, 
                               returns_long,
                               cols = east)
G7_EWMA_i2 <- compute_strategy(wide_smooth_ewma_i2, 
                               returns_long, 
                               suffix = "_EWMA", 
                               cols = g7)
G7_SMA_i2 <- compute_strategy(wide_smooth_sma_i2, 
                              returns_long, 
                              suffix = "_SMA", 
                              cols = g7)
G7_NS_i2 <- compute_strategy(wide_non_smooth_i2, 
                             returns_long,
                             cols = g7)

### combining returns

strat_dfs_i2 <- list(
    G7_NS_i2 = G7_NS_i2,
    G7_SMA_i2 = G7_SMA_i2,
    G7_EWMA_i2 = G7_EWMA_i2,
    ASIA_NS_i2 = ASIA_NS_i2,
    ASIA_SMA_i2 = ASIA_SMA_i2,
    ASIA_EWMA_i2 = ASIA_EWMA_i2,
    WEST_NS_i2 = WEST_NS_i2,
    WEST_SMA_i2 = WEST_SMA_i2,
    WEST_EWMA_i2 = WEST_EWMA_i2
)

combined_geo_strats_i2 <- bind_rows(strat_dfs_i2, .id = "Source") %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    separate(Source, into = c("Region", "Smoothing"), sep = "_", extra = "merge") %>%
    mutate(
        Smoothing = case_when(
            Smoothing == "NS_i2" ~ "Non-Smoothed",
            Smoothing == "EWMA_i2" ~ "EWMA",
            Smoothing == "SMA_i2" ~ "SMA",
            TRUE ~ Smoothing
        )
    ) %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha Based",
            Strategy == "Cume_d" ~ "Delta-Based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )

combined_geo_strats_i2 %>%
    filter(Region == "G7") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "G7I2.png",
       path = "plots/")

combined_geo_strats_i2 %>%
    filter(Region == "ASIA") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "ASIAI2.png",
       path = "plots/")

combined_geo_strats_i2 %>%
    filter(Region == "WEST") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(file = "WESTI2.png",
       path = "plots/")


## performance metrics

wide_perf_i2 <- bind_rows(strat_dfs_i2, .id = "Source") %>%
    pivot_longer(cols = c(Ret_a, Ret_d, Ret_r), names_to = "Type", values_to = "Return") %>%
    mutate(Strategy = paste(Source, Type, sep = "_")) %>%
    select(Date, Strategy, Return) %>%
    pivot_wider(names_from = Strategy, values_from = Return) %>%
    arrange(Date)

ret_xts_i2 <- xts(wide_perf_i2[,-1], order.by = wide_perf_i2$Date)

alphaScreen_i2 <- alphaScreening(X = ret_xts_i2, control = c("nCore" = parallel::detectCores() - 1))
sharpeScreen_i2 <- sharpeScreening(X = ret_xts_i2, control = c("nCore" = parallel::detectCores() - 1))
msharpeScreen_i2 <- msharpeScreening(X = ret_xts_i2, control = c("nCore" = parallel::detectCores() - 1))

sharpe(wide_perf_i2[,-1])
msharpe(wide_perf_i2[,-1])


## plotting weights

alpha_ewma_long_i2 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsEWMAI2.png",
       path = "plots/")

alpha_sma_long_i2 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsSMAI2.png",
       path = "plots/")

alpha_weights_i2 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsNSI2.png",
       path = "plots/")


## ---- Portfolio Allocation: Instrument 3 (constant + Absolute Lagged Errors) ----

### -------- Creating Weights ----
alpha_long_i3 <- as.data.frame(alpha_matrix_i3) %>%
    rownames_to_column("Date") %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
           Alpha = as.numeric(Alpha),
           Deviation = abs(Alpha - 0.5),                    
           Direction = ifelse(Alpha < 0.5, "Long", "Short"),
           Delta = Alpha - dplyr::lag(Alpha)
    )

alpha_weights_i3 <- alpha_long_i3 %>%
    group_by(Date) %>%
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()

wide_smooth_ewma_i3 <- alpha_ewma_i3 %>%  
    arrange(Date) %>%
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ EMA(.x, n = 10),
            .names = "{.col}_EWMA"
        )
    ) %>%
    group_by() %>% # if you ever have multiple symbols, do group_by(symbol)
    mutate(
        across(
            .cols = ends_with("_EWMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

wide_smooth_sma_i3 <- alpha_sma_i3 %>%     # or whatever your SMA table is
    arrange(Date) %>%
    # 1a) compute a 10‐month EWMA of *each* currency's alpha
    mutate(
        across(
            .cols = -Date,
            .fns  = ~ SMA(.x, n = 10),
            .names = "{.col}_SMA"
        )
    ) %>%
    group_by() %>% 
    mutate(
        across(
            .cols = ends_with("_SMA"),
            .fns  = ~ .x - lag(.x),
            .names = "{.col}.d"
        )
    ) %>%
    ungroup()

alpha_ewma_long_i3 <- wide_smooth_ewma_i3 %>%
    pivot_longer(
        cols      = ends_with("_EWMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_EWMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()

alpha_sma_long_i3 <- wide_smooth_sma_i3 %>%
    pivot_longer(
        cols      = ends_with("_SMA"),
        names_to  = "Currency",
        values_to = "Alpha"
    ) %>%
    mutate(
        Currency = sub("_SMA$", "", Currency)
    ) %>% 
    group_by(Currency) %>% 
    arrange(Date) %>% 
    mutate(
        Delta = Alpha - lag(Alpha)
    ) %>% 
    ungroup() %>% 
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a   = raw_a / sum(abs(raw_a), na.rm = TRUE),
        n_pos = sum(Delta >  0, na.rm = TRUE),
        n_neg = sum(Delta <  0, na.rm = TRUE),
        w_d   = case_when(
            Delta >  0 ~  -1/n_pos,
            Delta <  0 ~  +1/n_neg,
            TRUE       ~   0
        )
    ) %>%
    ungroup() %>%
    group_by(Date) %>%
    arrange(Date, Alpha) %>%
    mutate(
        rank = row_number(),
        n = n(),
        w_r = case_when(
            rank <= floor(n/2) ~ +1/floor(n/2),
            rank > floor(n/2) ~ -1/ceiling(n/2),
            T ~ 0
        )
    ) %>%
    ungroup()


alpha_weights_filtered_i3 <- alpha_weights_i3 %>% filter(Date >= common_start)

### -------- Calculating Returns of simple strats ----

alpha_portfolios_i3 <- alpha_weights_filtered_i3 %>%
    left_join(returns_long, by = c("Date", "Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )


## plotting simple strat - non-smoothed
ggplot(alpha_portfolios_i3, aes(x = Date, y = Cume_a)) +
    geom_line() +
    labs(title = "Cumulative Returns of Alpha-Based Portfolio",
         x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_colour_viridis_d(option = "D") 


### smoothed data

a_ewma_portfolios_i3 <- alpha_ewma_long_i3 %>% filter(Date >= common_start) %>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    ) 

a_sma_portfolios_i3 <- alpha_sma_long_i3 %>% filter(Date >= common_start)%>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
        Ret_r = sum(w_r * ExcessReturn, na.rm = TRUE),
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_r = cumprod(1 + Ret_r) - 1,
        Date = as.Date(Date)
    )

## plotting all strats

non_smoothed_long_i3 <- alpha_portfolios_i3 %>%
    mutate(Source = "Non-smoothed") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

ewma_long_i3 <- a_ewma_portfolios_i3 %>%
    mutate(Source = "EWMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)

sma_long_i3 <- a_sma_portfolios_i3 %>%
    mutate(Source = "SMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Cume_r, Source)


# Combine all into one data frame
combined_cume_i3 <- bind_rows(non_smoothed_long_i3, ewma_long_i3, sma_long_i3) %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )


# Plot
ggplot(combined_cume_i3, aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Source)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Date", 
         y = "Cumulative Return") +
    #scale_color_viridis_d(option = "E") +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(2, "cm"))
    ggsave("allRegionsAllStratsi3.png",
           path = "plots/")

#### ------- Geographical Strategies ----------------------------

wide_non_smooth_i3 <- alpha_weights_i3 %>%
    dplyr::select(Date, Currency, Alpha) %>%
    pivot_wider(names_from = Currency, values_from = Alpha)

WEST_EWMA_i3 <- compute_strategy(wide_smooth_ewma_i3, 
                                 returns_long, 
                                 suffix = "_EWMA", 
                                 cols = west)
WEST_SMA_i3 <- compute_strategy(wide_smooth_sma_i3, 
                                returns_long, 
                                suffix = "_SMA", 
                                cols = west)
WEST_NS_i3 <- compute_strategy(wide_non_smooth_i3, 
                               returns_long,
                               cols = west)

ASIA_EWMA_i3 <- compute_strategy(wide_smooth_ewma_i3, 
                                 returns_long, 
                                 suffix = "_EWMA", 
                                 cols = east)
ASIA_SMA_i3 <- compute_strategy(wide_smooth_sma_i3, 
                                returns_long, 
                                suffix = "_SMA", 
                                cols = east)
ASIA_NS_i3 <- compute_strategy(wide_non_smooth_i3, 
                               returns_long,
                               cols = east)
G7_EWMA_i3 <- compute_strategy(wide_smooth_ewma_i3, 
                               returns_long, 
                               suffix = "_EWMA", 
                               cols = g7)
G7_SMA_i3 <- compute_strategy(wide_smooth_sma_i3, 
                              returns_long, 
                              suffix = "_SMA", 
                              cols = g7)
G7_NS_i3 <- compute_strategy(wide_non_smooth_i3, 
                             returns_long,
                             cols = g7)

### combining returns

strat_dfs_i3 <- list(
    G7_NS_i3 = G7_NS_i3,
    G7_SMA_i3 = G7_SMA_i3,
    G7_EWMA_i3 = G7_EWMA_i3,
    ASIA_NS_i3 = ASIA_NS_i3,
    ASIA_SMA_i3 = ASIA_SMA_i3,
    ASIA_EWMA_i3 = ASIA_EWMA_i3,
    WEST_NS_i3 = WEST_NS_i3,
    WEST_SMA_i3 = WEST_SMA_i3,
    WEST_EWMA_i3 = WEST_EWMA_i3
)

combined_geo_strats_i3 <- bind_rows(strat_dfs_i3, .id = "Source") %>%
    pivot_longer(cols = c(Cume_a, Cume_d, Cume_r), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    separate(Source, into = c("Region", "Smoothing"), sep = "_", extra = "merge") %>%
    mutate(
        Smoothing = case_when(
            Smoothing == "NS_i3" ~ "Non-Smoothed",
            Smoothing == "EWMA_i3" ~ "EWMA",
            Smoothing == "SMA_i3" ~ "SMA",
            TRUE ~ Smoothing
        )
    ) %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            Strategy == "Cume_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    )

combined_geo_strats_i3 %>%
    filter(Region == "G7") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(file = "G7I3.png",
       path = "plots/")

combined_geo_strats_i3 %>%
    filter(Region == "ASIA") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(file = "ASIAI3.png",
       path = "plots/")

combined_geo_strats_i3 %>%
    filter(Region == "WEST") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(file = "WESTI3.png",
       path = "plots/")


## performance metrics

wide_perf_i3 <- bind_rows(strat_dfs_i3, .id = "Source") %>%
    pivot_longer(cols = c(Ret_a, Ret_d, Ret_r), names_to = "Type", values_to = "Return") %>%
    mutate(Strategy = paste(Source, Type, sep = "_")) %>%
    select(Date, Strategy, Return) %>%
    pivot_wider(names_from = Strategy, values_from = Return) %>%
    arrange(Date)

ret_xts_i3 <- xts(wide_perf_i3[,-1], order.by = wide_perf_i3$Date)

alphaScreen_i3 <- alphaScreening(X = ret_xts_i3, control = c("nCore" = parallel::detectCores() - 1))
sharpeScreen_i3 <- sharpeScreening(X = ret_xts_i3, control = c("nCore" = parallel::detectCores() - 1))
msharpeScreen_i3 <- msharpeScreening(X = ret_xts_i3, control = c("nCore" = parallel::detectCores() - 1))

sharpe(wide_perf_i3[,-1])
msharpe(wide_perf_i3[,-1])


## plotting weights

alpha_ewma_long_i3 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsEWMAI3.png",
       path = "plots/")

alpha_sma_long_i3 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsSMAI3.png",
       path = "plots/")

alpha_weights_i3 %>%
    pivot_longer(cols = c(w_a, w_d, w_r), names_to = "Strategy", values_to = "Weight") %>%
    mutate(
        Strategy = case_when(
            Strategy == "w_a" ~ "Alpha-based",
            Strategy == "w_d" ~ "Delta-based",
            Strategy == "w_r" ~ "Ranked Alphas",
            TRUE ~ Strategy
        )
    ) %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    theme_minimal() +
    labs(y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
ggsave(file = "WeightsNSI3.png",
       path = "plots/")



##### ----- Expectile Regression Strategies -----

## first we want to join the alpha and the monthly data
monthly_long <- monthly_df1 %>%
    dplyr::select(YM, ends_with("_prem"), ends_with("_excess")) %>%
    pivot_longer(-YM, names_to = c("Currency", ".value"), names_sep = "_") %>%
    rename(Date = YM)

joined <- monthly_long %>%
    inner_join(alpha_long %>% dplyr::select(Date, Currency, Alpha), by = c("Date", "Currency")) %>%
    arrange(Currency, Date)

### Expectile regression:

### --- Expectile Functions ----
expectile_loss <- function(beta, X, y, tau) {
    u <- y - X %*% beta
    weights <- ifelse(u >= 0, tau, 1 - tau)
    sum(weights * u^2)
}
expectile_optim <- function(X, y, tau) {
    X <- as.matrix(cbind(1, X))
    if (any(is.na(X)) || any(is.na(y))) {
        message("→ Skipping: NA in input")
        return(rep(NA_real_, ncol(X)))
    }
    
    start <- tryCatch(coef(lm.fit(X, y)), error = function(e) rep(NA_real_, ncol(X)))
    if (any(is.na(start))) {
        message("→ Skipping: OLS failed")
        return(rep(NA_real_, ncol(X)))
    }
    
    fit <- tryCatch(
        optim(par = start,
              fn = expectile_loss,
              X = X, y = y, tau = tau,
              method = "BFGS",
              control = list(fnscale = 1)),
        error = function(e) NA
    )
    
    if (is.list(fit) && !is.null(fit$par) && all(is.finite(fit$par))) {
        message("✓ Optim succeeded. β₀ = ", round(fit$par[1], 6),
                ", β₁ = ", round(fit$par[2], 6), ", tau = ", round(tau, 3))
        return(as.numeric(fit$par))
    } else {
        message("→ Optimization failed or returned NA at tau = ", round(tau, 3))
        return(rep(NA_real_, ncol(X)))
    }
    return(fit$par)
}
run_expectile_forecast <- function(monthly_df, alpha_df, window = 60) {
    # Reshape prem and excess into long format
    data_long <- monthly_df %>%
        select(YM, ends_with("_prem"), ends_with("_excess")) %>%
        pivot_longer(-YM, names_to = c("Currency", ".value"), names_sep = "_") %>%
        rename(Date = YM)
    
    # Join alpha source
    df <- left_join(data_long, alpha_df, by = c("Date", "Currency")) %>%
        arrange(Currency, Date)
    
    all_results <- list()
    currencies <- unique(df$Currency)
    
    for (curr in currencies) {
        message("Running expectile model for ", curr)
        d <- df %>% filter(Currency == curr) %>% arrange(Date)
        
        n <- nrow(d)
        d$Beta_0 <- NA_real_
        d$Beta_1 <- NA_real_
        d$ExpectileForecast <- NA_real_
        
        for (i in (window + 1):(n - 1)) {
            train <- d[(i - window):(i - 1), ]
            test <- d[i + 1, ]
            
            tau <- train[nrow(train), ]$Alpha  # lagged alpha for asymmetry
            
            if (anyNA(train$prem) || anyNA(train$excess) || is.na(tau)) next
            
            X_train <- matrix(train$prem, ncol = 1)
            y_train <- train$excess
            beta <- as.numeric(expectile_optim(X_train, y_train, tau))
            if (length(beta) != 2 || any(is.na(beta))) next
            
            d$Beta_0[i + 1] <- beta[1]
            d$Beta_1[i + 1] <- beta[2]
            
            x_test <- test$prem
            if (!is.na(x_test)) {
                d$ExpectileForecast[i + 1] <- c(1, x_test) %*% beta
            }
        }
        
        all_results[[curr]] <- d
    }
    
    return(bind_rows(all_results))
}

### Expectile Prelim Results ----

expectile_results <- run_expectile_forecast(monthly_df1, alpha_long, window = 60)

expectile_results <- expectile_results %>%
    mutate(
        ForecastError = excess - ExpectileForecast,
        ForecastLoss = ifelse(
            ForecastError < 0,
            2 * (1 - Alpha) * ForecastError^2,
            2 * Alpha * ForecastError^2
        ),
        NaiveError = excess - prem,
        NaiveLoss = ifelse(
            NaiveError <0,
            2 * (1- Alpha) * NaiveError^2,
            2 * Alpha * NaiveError^2
        )
    )

expectile_results %>%
    summarise(
        MeanLoss_Expectile = mean(ForecastLoss, na.rm = TRUE),
        MeanLoss_Naive = mean(NaiveLoss, na.rm = TRUE),
        RMSE_Expectile = sqrt(mean(ForecastError^2, na.rm = TRUE)),
        RMSE_Naive = sqrt(mean(NaiveError^2, na.rm = TRUE))
    )
expectile_results %>%    
    # Diebold-Marino Test under asymmetric loss
    mutate(loss_diff = ForecastLoss - NaiveLoss) %>%
    summarise(
        DM_stat = mean(loss_diff, na.rm = TRUE) / (sd(loss_diff, na.rm = TRUE) / sqrt(n())),
        p_value = 2 * (1 - pnorm(abs(DM_stat)))
    )

dm_by_currency <- expectile_results %>%
    filter(!is.na(ForecastLoss), !is.na(NaiveLoss)) %>%
    group_by(Currency) %>%
    summarise(
        DM_stat = mean(ForecastLoss - NaiveLoss, na.rm = TRUE) /
            (sd(ForecastLoss - NaiveLoss, na.rm = TRUE) / sqrt(n())),
        p_value = 2 * (1 - pnorm(abs(DM_stat))),
        .groups = "drop"
    ) %>%
    arrange(p_value)

rolling_dm <- function(df, window = 60) {
    n <- nrow(df)
    dm_stats <- rep(NA_real_, n - window + 1)
    p_values <- rep(NA_real_, n - window + 1)
    dates <- df$Date[window:n]
    
    for (i in seq_len(n - window + 1)) {
        sub_df <- df[i:(i + window - 1), ]
        d <- sub_df$ForecastLoss - sub_df$NaiveLoss
        d <- d[!is.na(d)]
        if (length(d) > 10) {
            dm_stat <- mean(d) / (sd(d) / sqrt(length(d)))
            p_val <- 2 * (1 - pnorm(abs(dm_stat)))
            
            dm_stats[i] <- dm_stat
            p_values[i] <- p_val
        }
    }
    
    return(tibble(
        Date = dates,
        DM_stat = dm_stats,
        p_value = p_values
    ))
}

# Function to extract runs of significance
get_sig_blocks <- function(df) {
    v       <- df$Significant 
    rle_out <- rle(v)
    lengths <- rle_out$lengths
    values  <- rle_out$values
    ends    <- cumsum(lengths)
    starts  <- c(1, head(ends + 1, -1))
    
    sig_runs <- tibble(
        Currency = df$Currency[1],
        start = df$Date[starts],
        end   = df$Date[ends],
        sig   = values
    ) %>% filter(sig)
    
    return(sig_runs)
}


dm_rolling_by_currency <- expectile_results %>%
    filter(!is.na(ForecastLoss), !is.na(NaiveLoss)) %>%
    group_by(Currency) %>%
    group_modify(~ rolling_dm(.x, window = 60)) %>%
    ungroup() %>%
    mutate(Significant = abs(DM_stat) > qnorm(0.975))

sig_blocks <- dm_rolling_by_currency %>%
    group_by(Currency) %>%
    group_modify(get_sig_blocks(.x)) %>%
    ungroup()

ggplot(dm_rolling_by_currency, aes(x = Date, y = DM_stat, color = Currency)) +
    geom_line() +
    geom_hline(yintercept = qnorm(0.975), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -qnorm(0.975), linetype = "dashed", color = "grey40") +
    labs(y = "Rolling Diebold–Mariano Statistic", x = "Date") +
    facet_wrap(~ Currency, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none")
ggsave("RollingDMtestCurrency.png",
       path = "plots/",
       width = 10, height = 6)

dm_rolling_aggregate <- expectile_results %>%
    filter(!is.na(ForecastLoss), !is.na(NaiveLoss)) %>%
    arrange(Date) %>%
    group_by(Date) %>%
    summarise(
        ForecastLoss = mean(ForecastLoss, na.rm = TRUE),
        NaiveLoss = mean(NaiveLoss, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    rolling_dm(window = 60)

ggplot(dm_rolling_aggregate, aes(x = Date, y = DM_stat)) +
    geom_line() +
    geom_hline(yintercept = qnorm(0.975), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -qnorm(0.975), linetype = "dashed", color = "grey40") +
    labs(y = "Rolling Diebold–Mariano Statistic", x = "Date") +
    theme_minimal()
ggsave("RollingDMtestAggregate.png",
       path = "plots/",
       width = 10, height = 6)

# Create a logical column to flag significance
dm_agg_df <- dm_rolling_aggregate %>%
    mutate(significant = DM_stat > 1.96)

# Identify contiguous significant regions
rle_sig <- rle(dm_agg_df$significant)
lengths <- rle_sig$lengths
values <- rle_sig$values
ends <- cumsum(lengths)
starts <- c(1, head(ends + 1, -1))

rects <- tibble(start = starts, end = ends, is_sig = values) %>%
    filter(is_sig) %>%
    mutate(
        xmin = dm_agg_df$Date[start],
        xmax = dm_agg_df$Date[end]
    )

# Plot
ggplot(dm_agg_df, aes(x = Date, y = DM_stat)) +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
    geom_line(color = "black") +
    geom_hline(yintercept = qnorm(0.975), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -qnorm(0.975), linetype = "dashed", color = "grey40") +
    theme_minimal() +
    labs(y = "Rolling DM Statistic")
ggsave("RollingDMtestAggregate2.png",
       path = "plots/",
       width = 10, height = 6)






ggplot(expectile_results, aes(x = ExpectileForecast, y = excess)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(
        title = "Actual vs Forecasted Excess Returns",
        x = "Forecasted Excess Return",
        y = "Actual Excess Return"
    ) +
    theme_minimal()


ggplot(expectile_results, aes(x = ExpectileForecast, y = excess)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    facet_wrap(~ Currency, scales = "free") +
    labs(
        title = "Forecast vs Actual Excess Returns by Currency",
        x = "Forecasted",
        y = "Actual"
    ) +
    theme_minimal()


ggplot(expectile_results, aes(x = Date)) +
    geom_line(aes(y = excess, color = "Actual")) +
    geom_line(aes(y = ExpectileForecast, color = "Forecast")) +
    facet_wrap(~ Currency, scales = "free_y") +
    labs(
        title = "Forecasted vs Actual Excess Returns Over Time",
        y = "Excess Return",
        color = "Series"
    ) +
    theme_minimal()

expectile_results %>%
    filter(Currency == "CAD") %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = excess, color = "Actual")) +
    geom_line(aes(y = ExpectileForecast, color = "Forecast")) +
    labs(
        title = "CAD: Forecasted vs Actual Excess Returns",
        y = "Excess Return",
        color = "Series"
    ) +
    theme_minimal()

mz_table <- expectile_results %>%
    filter(!is.na(ExpectileForecast), !is.na(excess)) %>%
    group_by(Currency) %>%
    group_modify(~ {
        model <- lm(excess ~ ExpectileForecast, data = .x)
        tidy_res <- tidy(model)
        glance_res <- glance(model)
        linear_test <- car::linearHypothesis(model, c("ExpectileForecast = 1", "(Intercept) = 0"))
        
        tibble(
            Intercept     = tidy_res$estimate[1],
            SE_Intercept  = tidy_res$std.error[1],
            Beta          = tidy_res$estimate[2],
            SE_Beta       = tidy_res$std.error[2],
            t_Beta        = tidy_res$statistic[2],
            p_Beta        = tidy_res$p.value[2],
            F_optimality  = linear_test$F[2],
            p_optimality  = linear_test$`Pr(>F)`[2],
            R2            = glance_res$r.squared,
            N             = glance_res$nobs
        )
    }) %>%
    ungroup()
print(mz_table)


###  ---- Comparing Expectile vs Naive (prem) -----

strategy_df <- expectile_results %>%
    filter(Currency != "VND") %>%
    filter(!is.na(prem), !is.na(ExpectileForecast), !is.na(excess)) %>%
    mutate(
        Position_Naive = sign(prem),
        Position_Expectile = sign(ExpectileForecast),
        Return_Naive = Position_Naive * excess,
        Return_Expectile = Position_Expectile * excess,
        weight_forecast = ifelse(!is.na(ExpectileForecast),
                                 ExpectileForecast, 0),
        weight_forecast = weight_forecast / sum(abs(weight_forecast), na.rm = TRUE),
        Return_Weight = weight_forecast * excess
    )

portfolio_returns_exp <- strategy_df %>%
    group_by(Date) %>%
    summarise(
        PortReturn_Naive = mean(Return_Naive, na.rm = TRUE),
        PortReturn_Expectile = mean(Return_Expectile, na.rm = TRUE),
        PortReturn_Weight = sum(Return_Weight, na.rm = TRUE),
        Active_Naive = sum(!is.na(Return_Naive)),
        Active_Expectile = sum(!is.na(Return_Expectile))
    ) %>%
    filter(!is.na(PortReturn_Naive), !is.na(PortReturn_Expectile))

portfolio_returns_exp <- portfolio_returns_exp %>%
    mutate(
        Cum_Naive = cumsum(PortReturn_Naive),
        Cum_Expectile = cumsum(PortReturn_Expectile),
        Cum_Weight = cumsum(PortReturn_Weight),
    )

# Plot
ggplot(portfolio_returns_exp, aes(x = Date)) +
    geom_line(aes(y = Cum_Naive, color = "Naive")) +
    geom_line(aes(y = Cum_Expectile, color = "Expectile")) +
    labs(title = "Cumulative Portfolio Returns",
         y = "Cumulative Return",
         color = "Strategy") +
    theme_minimal() +
    theme(legend.position = "bottom")







### --- Compute Expectile Strategies ----

compute_strategy_expectile <- function(forecast_df, return_df, suffix = NULL, cols = NULL) {
    suffix_pattern <- if (!is.null(suffix)) paste0("_", suffix) else ""
    
    # Clean and align
    strategy_df <- forecast_df %>%
        filter(!is.na(ExpectileForecast), !is.na(excess)) %>%
        filter(Currency != "VND") %>%
        mutate(
            Position = sign(ExpectileForecast),
            ForecastError = excess - ExpectileForecast,
            ForecastLoss = ifelse(
                ForecastError < 0,
                2 * (1 - Alpha) * ForecastError^2,
                2 * Alpha * ForecastError^2
            ),
            PositionNaive = sign(prem)
        ) %>%
        left_join(return_df, by = c("Date", "Currency")) %>%
        mutate(
            Return = Position * ExcessReturn,
            ReturnNaive = PositionNaive * ExcessReturn
        )
    
    if (!is.null(cols)) {
        strategy_df <- strategy_df %>% filter(Currency %in% cols)
    }
    
    # Portfolio returns (equal weight across currencies)
    port_returns <- strategy_df %>%
        group_by(Date) %>%
        summarise(
            PortReturn = mean(Return, na.rm = TRUE),
            PortLoss = mean(ForecastLoss, na.rm = TRUE),
            n_assets = sum(!is.na(Return)),
            PortReturn_naive = mean(ReturnNaive, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            Cumulative = cumprod(1 + PortReturn) - 1,
            Cum_Naive = cumprod(1 + PortReturn_naive) - 1
        )
    
    return(list(
        weights = strategy_df %>%
            select(Date, Currency, Alpha, ExpectileForecast, Position, ForecastLoss),
        returns = port_returns
    ))
}


#### Implementing the expectile strats function ----

forecast_df_exp <- run_expectile_forecast(monthly_df1, alpha_long, window = 60)
forecast_exp_ewma <- run_expectile_forecast(monthly_df1, alpha_ewma_long, window = 60)
forecast_exp_sma <- run_expectile_forecast(monthly_df1, alpha_sma_long, window = 60)

all_exp <- compute_strategy_expectile(forecast_df_exp, returns_long)
all_exp_ewma <- compute_strategy_expectile(forecast_exp_ewma, returns_long)
all_exp_sma <- compute_strategy_expectile(forecast_exp_sma, returns_long)

g7_exp <- compute_strategy_expectile(forecast_df_exp, returns_long, cols = g7)
g7_exp_ewma <- compute_strategy_expectile(forecast_exp_ewma, returns_long, cols = g7)
g7_exp_sma <- compute_strategy_expectile(forecast_exp_sma, returns_long, cols = g7)

asia_exp <- compute_strategy_expectile(forecast_df_exp, returns_long, cols = east)
asia_exp_ewma <- compute_strategy_expectile(forecast_exp_ewma, returns_long, cols = east)
asia_exp_sma <- compute_strategy_expectile(forecast_exp_sma, returns_long, cols = east)

west_exp <- compute_strategy_expectile(forecast_df_exp, returns_long, cols = west)
west_exp_ewma <- compute_strategy_expectile(forecast_exp_ewma, returns_long, cols = west)
west_exp_sma <- compute_strategy_expectile(forecast_exp_sma, returns_long, cols = west)

### Combine all strategies into one data frame
combined_exp_strats <- bind_rows(
    all_exp$returns %>% mutate(Source = "All",
                               Smoothing = "NS"),
    all_exp_ewma$returns %>% mutate(Source = "All",
                                    Smoothing = "EWMA"),
    all_exp_sma$returns %>% mutate(Source = "All",
                                   Smoothing = "SMA"),
    g7_exp$returns %>% mutate(Source = "G7",
                               Smoothing = "NS"),
    g7_exp_ewma$returns %>% mutate(Source = "G7",
                                   Smoothing = "EWMA"),
    g7_exp_sma$returns %>% mutate(Source = "G7",
                                   Smoothing = "SMA"),
    asia_exp$returns %>% mutate(Source = "Asia",
                                Smoothing = "NS"),
    asia_exp_ewma$returns %>% mutate(Source = "Asia",
                                     Smoothing = "EWMA"),
    asia_exp_sma$returns %>% mutate(Source = "Asia",
                                     Smoothing = "SMA"),
    west_exp$returns %>% mutate(Source = "West",
                                Smoothing = "NS"),
    west_exp_ewma$returns %>% mutate(Source = "West",
                                    Smoothing = "EWMA"),
    west_exp_sma$returns %>% mutate(Source = "West",
                                     Smoothing = "SMA")
)
# Reshape for plotting
combined_exp_strats_long <- combined_exp_strats %>%
    pivot_longer(cols = c(PortReturn, PortLoss, Cumulative, Cum_Naive), names_to = "Metric", values_to = "Value") %>%
    mutate(
        Strategy = case_when(
            Metric == "PortReturn" ~ "Portfolio Return",
            Metric == "PortLoss" ~ "Portfolio Loss",
            Metric == "Cumulative" ~ "Expectile Forecast",
            Metric == "Cum_Naive" ~ "Naive Forecast",
            TRUE ~ Metric
        )
    ) %>%
    mutate(
        Smoothing = factor(Smoothing, levels = c("NS", "EWMA", "SMA")),
        Source = factor(Source, levels = c("All", "G7", "Asia", "West"))
    )

# Plotting the cumulative returns for all strategies
ggplot(combined_exp_strats_long %>% filter(
    Strategy %in% c("Expectile Forecast", "Naive Forecast"),
    Smoothing == "NS"), 
       aes(x = Date, y = Value, color = Source, linetype = Strategy)) +
    geom_line() +
    labs(x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom")

ggplot(combined_exp_strats_long %>% filter(
    Metric %in% c("Cumulative Return", "Naive Cumulative Return"),
    Smoothing == "SMA"), 
    aes(x = Date, y = Value, color = Source, linetype = Metric)) +
    geom_line() +
    labs(x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom")

ggplot(combined_exp_strats_long %>% filter(
    Metric %in% c("Cumulative Return", "Naive Cumulative Return"),
    Smoothing == "EWMA"), 
    aes(x = Date, y = Value, color = Source, linetype = Metric)) +
    geom_line() +
    labs(x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom")



ggplot(
    combined_exp_strats_long %>%
        filter(Strategy %in% c("Expectile Forecast", "Naive Forecast")),
    aes(x = Date, y = Value, color = Smoothing, linetype = Strategy)
) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon", alpha = 0.5) +
    geom_line() +
    facet_wrap(~ Source, scales = "fixed", nrow = 2) +
    labs(
        x = "Date", y = "Cumulative Return"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom")
ggsave("ExpectileStrats.png",
       path = "plots/")


combined_geo_strats %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Smoothing, linetype = Strategy)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon", alpha = 0.5) +
    geom_line() +
    facet_wrap(~ Region, scales = "free_y", nrow = 2) +
    labs(
        x = "Date", y = "Cumulative Return"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom")
ggsave("CombinedGeoStrats.png",
       path = "plots/")




##### --- Comparing Expectile with all others ----


# rearrange expectile results:

exp_returns <- combined_exp_strats_long %>%
    filter(Metric == "PortReturn") %>%
    mutate(
        Strategy = paste(Source, Smoothing, "Expectile", sep = "_")
    ) %>%
    dplyr::select(Date, Strategy, Return = Value)

naive_returns <- combined_exp_strats_long %>%
    dplyr::select(Date, Source, Smoothing, PortReturn_naive) %>%
    distinct(Date, Source, Smoothing, PortReturn_naive) %>%
    rename(Value = PortReturn_naive) %>%
    mutate(
        Strategy = paste(Source, Smoothing, "Naive", sep = "_")
    ) %>%
    dplyr::select(Date, Strategy, Return = Value)


geo_returns <- combined_geo_strats %>%
    rename(Alpha = Ret_a,
           Delta = Ret_d,
           Ranked = Ret_r) %>%
    pivot_longer(cols = c(Alpha, Delta, Ranked), names_to = "strat", values_to = "Return") %>%
    mutate(Strategy = paste(Region, Smoothing, strat, sep = "_")) %>%
    dplyr::select(Date, Strategy, Return)

all_portfolios <- bind_rows(
    exp_returns,
    naive_returns,
    geo_returns
) %>%
    arrange(Date) %>%
    group_by(Date, Strategy) %>%
    summarise(Return = mean(Return, na.rm = TRUE), .groups = "drop")

all_portfolios_wide <- all_portfolios %>%
    pivot_wider(names_from = Strategy, values_from = Return) 

min_common_date <- all_portfolios %>%
    group_by(Strategy) %>%
    summarise(min_date = min(Date)) %>%
    summarise(max(min_date)) %>%
    pull()

all_portfolios_wide <- all_portfolios_wide %>%
    filter(Date >= min_common_date) %>%
    arrange(Date)

all_portfolios_wide %>% dplyr::select(sort(names(.)))
all_portfolios_wide <- all_portfolios_wide[, c("Date", 
                                               "All_NS_Naive", "All_NS_Expectile", "ALL_NS_Alpha", "ALL_NS_Delta", "ALL_NS_Ranked",
                                               "All_SMA_Naive", "All_SMA_Expectile", "ALL_SMA_Alpha", "ALL_SMA_Delta", "ALL_SMA_Ranked",
                                               "All_EWMA_Naive", "All_EWMA_Expectile", "ALL_EWMA_Alpha", "ALL_EWMA_Delta", "ALL_EWMA_Ranked",
                                               "G7_NS_Naive", "G7_NS_Expectile", "G7_NS_Alpha", "G7_NS_Delta", "G7_NS_Ranked",
                                               "G7_SMA_Naive", "G7_SMA_Expectile", "G7_SMA_Alpha", "G7_SMA_Delta", "G7_SMA_Ranked",
                                               "G7_EWMA_Naive", "G7_EWMA_Expectile", "G7_EWMA_Alpha", "G7_EWMA_Delta", "G7_EWMA_Ranked",
                                               "Asia_NS_Naive", "Asia_NS_Expectile", "ASIA_NS_Alpha", "ASIA_NS_Delta", "ASIA_NS_Ranked",
                                               "Asia_SMA_Naive", "Asia_SMA_Expectile", "ASIA_SMA_Alpha", "ASIA_SMA_Delta", "ASIA_SMA_Ranked",
                                               "Asia_EWMA_Naive", "Asia_EWMA_Expectile", "ASIA_EWMA_Alpha", "ASIA_EWMA_Delta", "ASIA_EWMA_Ranked",
                                               "West_NS_Naive", "West_NS_Expectile", "WEST_NS_Alpha", "WEST_NS_Delta", "WEST_NS_Ranked",
                                               "West_SMA_Naive", "West_SMA_Expectile", "WEST_SMA_Alpha", "WEST_SMA_Delta", "WEST_SMA_Ranked",
                                               "West_EWMA_Naive", "West_EWMA_Expectile", "WEST_EWMA_Alpha", "WEST_EWMA_Delta", "WEST_EWMA_Ranked")]

all_portfolios_xts <- xts(all_portfolios_wide[,-1], order.by = all_portfolios_wide$Date)

all_alphaScreen <- alphaScreening(all_portfolios_xts, control = c("nCore" = parallel::detectCores() - 1))
all_sharpeScreen <- sharpeScreening(all_portfolios_xts, control = c("nCore" = parallel::detectCores() - 1))
all_msharpeScreen <- msharpeScreening(all_portfolios_xts, na.neg = F,control = c("nCore" = parallel::detectCores() - 1))





