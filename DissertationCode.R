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

ggplot(lm_daily_df[1:nrow(lm_daily_df)-1,], aes(x = Currency, y = Slope)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_errorbar(aes(ymin = Slope - 1.96 * SE.Slope,
                      ymax = Slope + 1.96 * SE.Slope),
                  width = 0.1) +
    labs(x = "Currency", 
         y = "Slope Coefficient") +
    theme_minimal()


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
            
            
            if (nrow(dat_cur) > 10) {
                gmm_moments <- function(theta, data, Inst) {
                    alpha <- theta[1]
                    e <- data$forecast_error
                    z <- as.matrix(data[, c("const", "ALE")])
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
alpha_sma_i3 <- smooth_alpha_i1 %>%
    mutate(across(-Date, ~ rollmean(., k = 10, fill = NA, align = "right")))


ggplot(smooth_alpha %>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha time series.png")

ggplot(alpha_ewma%>% pivot_longer(-Date, names_to = "Currency",
                                  values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D")  +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha EWMA time series.png")

ggplot(alpha_sma%>% pivot_longer(-Date, names_to = "Currency",
                                 values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha SMA time series.png")

# plotting alpha i1

ggplot(smooth_alpha_i1 %>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 1 time series.png")

ggplot(alpha_ewma_i1%>% pivot_longer(-Date, names_to = "Currency",
                                  values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D")  +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 1 EWMA time series.png")

ggplot(alpha_sma_i1%>% pivot_longer(-Date, names_to = "Currency",
                                 values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 1 SMA time series.png")

# plottiong alpha i2

ggplot(smooth_alpha_i2 %>% pivot_longer(-Date, names_to = "Currency",
                                        values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 2 time series.png")

ggplot(alpha_ewma_i2%>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D")  +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 2 EWMA time series.png")

ggplot(alpha_sma_i2%>% pivot_longer(-Date, names_to = "Currency",
                                    values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 2 SMA time series.png")

# plotting i3

ggplot(smooth_alpha_i3 %>% pivot_longer(-Date, names_to = "Currency",
                                        values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 3 time series.png")

ggplot(alpha_ewma_i3%>% pivot_longer(-Date, names_to = "Currency",
                                     values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D")  +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 3 EWMA time series.png")

ggplot(alpha_sma_i3%>% pivot_longer(-Date, names_to = "Currency",
                                    values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank()) +
    ggsave("alpha Inst 3 SMA time series.png")

## Simple Portfolio Allocation:

alpha_long <- as.data.frame(alpha_matrix) %>%
    rownames_to_column("Date") %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "Alpha") %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
           Alpha = as.numeric(Alpha),
           Deviation = abs(Alpha - 0.5),                    
           Direction = ifelse(Alpha < 0.5, "Long", "Short"),
           Delta = Alpha - dplyr::lag(Alpha)
    )

# creating the weights
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
    ungroup()

alpha_weights <- alpha_weights %>%
    mutate(
        dev   = abs(Alpha - 0.5),
        dir   = if_else(Alpha >  0.5, -1, +1),
        raw_a = dir * dev
    ) %>%
    group_by(Date) %>%
    mutate(
        w_a = raw_a / sum(abs(raw_a), na.rm = TRUE)  # alpha-based weights
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

## creating the weights for smoothed simple portfolios
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
        .groups = "drop"
    ) %>%
    mutate(
        Cume_d = cumprod(1 + Ret_d) - 1,
        Cume_a = cumprod(1 + Ret_a) - 1,
        Date = as.Date(Date)
    )


## plotting simple strat - non-smoothed
ggplot(alpha_portfolios, aes(x = Date, y = Cume_a)) +
    geom_line() +
    labs(title = "Cumulative Returns of Alpha-Based Portfolio",
         x = "Date",
         y = "Cumulative Return") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_colour_viridis_d(option = "D") 


## simple on smoothed data


a_ewma_portfolios <- alpha_ewma_long %>% filter(Date >= common_start) %>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE)
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1
    ) %>%
    mutate(Date = as.Date(Date))

a_sma_portfolios <- alpha_sma_long %>% filter(Date >= common_start)%>%
    left_join(returns_long, by = c("Date","Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE)
    ) %>%
    mutate(
        Cume_a = cumprod(1 + Ret_a) - 1,
        Cume_d = cumprod(1 + Ret_d) - 1
    )  %>%
    mutate(Date = as.Date(Date))

## plotting all strats

non_smoothed_long <- alpha_portfolios %>%
    mutate(Source = "Non-smoothed") %>%
    dplyr::select(Date, Cume_a, Cume_d, Source)

ewma_long <- a_ewma_portfolios %>%
    mutate(Source = "EWMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Source)

sma_long <- a_sma_portfolios %>%
    mutate(Source = "SMA") %>%
    dplyr::select(Date, Cume_a, Cume_d, Source)

# Combine all into one data frame
combined_cume <- bind_rows(non_smoothed_long, ewma_long, sma_long) %>%
    pivot_longer(cols = c(Cume_a, Cume_d), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    mutate(
        Strategy = case_when(
            Strategy == "Cume_a" ~ "Alpha-based",
            Strategy == "Cume_d" ~ "Delta-based",
            TRUE ~ Strategy
        )
    )


# Plot
ggplot(combined_cume, aes(x = Date, y = CumulativeReturn, color = Source, linetype = Strategy)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Date", 
         y = "Cumulative Return") +
    #scale_color_viridis_d(option = "E") +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(2, "cm"))



#### ------------------ Geographical Strategies ----------------------------

west <- c("CAD", "EUR", "GBP", "CHF")
east <- c("JPY", "SGD", "KRW", "THB")
g7 <- c("CAD", "EUR", "JPY", "GBP")

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
        ungroup()
    
    alpha_long %>%
        left_join(returns_long, by = c("Date", "Currency")) %>%
        group_by(Date) %>%
        summarise(
            Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),
            Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            Cume_a = cumprod(1 + Ret_a) - 1,
            Cume_d = cumprod(1 + Ret_d) - 1,
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
    WEST_EWMA = WEST_EWMA
)

combined_geo_strats <- bind_rows(strat_dfs, .id = "Source") %>%
    pivot_longer(cols = c(Cume_a, Cume_d), names_to = "Strategy", values_to = "CumulativeReturn") %>%
    separate(Source, into = c("Region", "Smoothing"), sep = "_")

combined_geo_strats %>%
    filter(Region == "G7") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")

combined_geo_strats %>%
    filter(Region == "ASIA") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")

combined_geo_strats %>%
    filter(Region == "WEST") %>%
    ggplot(aes(x = Date, y = CumulativeReturn, color = Strategy, linetype = Smoothing)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom")


## return profiles - volatility, drawdown, Sharpe, etc
# using peerperformance package

wide_perf <- bind_rows(strat_dfs, .id = "Source") %>%
    pivot_longer(cols = c(Ret_a, Ret_d), names_to = "Type", values_to = "Return") %>%
    mutate(Strategy = paste(Source, Type, sep = "_")) %>%
    select(Date, Strategy, Return) %>%
    pivot_wider(names_from = Strategy, values_from = Return) %>%
    arrange(Date)

ret_xts <- xts(wide_perf[,-1], order.by = wide_perf$Date)

alphaScreening(X = ret_xts)

sharpe(wide_perf[,-1])
msharpe(wide_perf[,-1])
alphaTesting(wide_perf[,"G7_NS_Ret_d"], wide_perf[,"ASIA_NS_Ret_d"])
alphaScreening()

## plotting weights

alpha_ewma_long %>%
    pivot_longer(cols = c(w_a, w_d), names_to = "Strategy", values_to = "Weight") %>%
    ggplot(aes(x = Date, y = Weight, fill = Currency)) +
    geom_area(position = "stack") +
    facet_wrap(~ Strategy, ncol = 1) +
    scale_fill_viridis_d(option = "D") +
    theme_minimal() +
    labs(title = "Portfolio Weights Over Time", y = "Weight", x = "Date") +
    theme(legend.position = "bottom")
