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
                
                alpha_vec[i] <- if (!is.null(res)) inv_logit(coef(res)[1]) else NA
            }
        }
        
        return(alpha_vec)
    }
    
    # Run in parallel
    alpha_list <- future_map(currencies, estimate_alpha_for_currency, .progress = TRUE)
    names(alpha_list) <- currencies
    
    alpha_matrix <- do.call(cbind, alpha_list)
    rownames(alpha_matrix) <- as.character(date_seq)
    return(alpha_matrix)
}

alpha_matrix <- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 0)
alpha_matrix_i1<- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 1)
alpha_matrix_i2<- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 2)
alpha_matrix_i3<- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p, Inst = 3)

## we want to remove Vietnamese Dong as there are so many NAs

alpha_matrix <- alpha_matrix[,1:ncol(alpha_matrix)-1]
alpha_matrix_i1 <- alpha_matrix_i1[,1:ncol(alpha_matrix_i1)-1]
alpha_matrix_i2 <- alpha_matrix_i2[,1:ncol(alpha_matrix_i2)-1]
alpha_matrix_i3 <- alpha_matrix_i3[,1:ncol(alpha_matrix_i3)-1]

## Plotting the alphas

smooth_alpha <- as.data.frame(alpha_matrix) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date))

alpha_ewma <- smooth_alpha %>%
    mutate(across(-Date, ~ EMA(., n = 10)))
alpha_sma <- smooth_alpha %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))

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

ggplot(alpha_sma%>% pivot_longer(-Date, names_to = "Currency",
                                 values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank())
    ggsave("alpha SMA time series.png",
           path = "plots/")

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
