### Dissertation Code

# premable
source("preamble.R")

# data
spot_df <- read_xlsx("data/data2.xlsx", sheet = "spot")
fwd_df <- read_xlsx("data/data2.xlsx", sheet = "forward")

df <- spot_df %>%
    inner_join(fwd_df, by = "DATES") %>%
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

# lag the spot to create excess and premium

lag <- monthly %>%
    transmute(
        YM = YM %m+% months(1), 
        across(-YM, ~ ., .names = "{.col}_t1"))

monthly_df <- monthly %>%
    inner_join(lag, by = "YM")

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


###------------------------------ GMM Alpha Estimation---------------------------

plan(multisession, workers = parallel::detectCores() - 1)

roll_alpha_gmm_parallel <- function(data, currencies, window = 60, p = 2) {
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
                    const = 1
                ) %>% na.omit()
            
            
            if (nrow(dat_cur) > 10) {
                gmm_moments <- function(theta, data) {
                    alpha <- theta[1]
                    e <- data$forecast_error
                    z <- matrix(1, nrow = length(e), ncol = 1) #z <- as.matrix(data[, c("const", "lag_error")])
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

alpha_matrix <- roll_alpha_gmm_parallel(monthly_df1, spot_cols, window = 60, p = p)

## we want to remove Vietnamese Dong as there are so many NAs

alpha_matrix <- alpha_matrix[,1:ncol(alpha_matrix)-1]
## Plotting the alphas

smooth_alpha <- as.data.frame(alpha_matrix) %>%
    rownames_to_column("Date") %>%
    mutate(Date = as.Date(Date))

alpha_ewma <- smooth_alpha %>%
    mutate(across(-Date, ~ EMA(., n = 10)))
alpha_sma <- smooth_alpha %>%
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
          legend.title = element_blank())

ggplot(alpha_ewma%>% pivot_longer(-Date, names_to = "Currency",
                                  values_to = "Alpha"),
       aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D")  +
    theme(legend.position = "bottom", 
          legend.title = element_blank())

ggplot(alpha_sma%>% pivot_longer(-Date, names_to = "Currency",
                                 values_to = "Alpha"), aes(x = Date, y = Alpha, color = Currency)) +
    geom_line() +
    labs(x = "Date",
         y = "Alpha") +
    theme_minimal() +
    scale_colour_viridis_d("Currency", option = "D") +
    theme(legend.position = "bottom", 
          legend.title = element_blank())

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

returns_long <- monthly_df1 %>%
    mutate(Date = YM) %>%
    dplyr::select(Date, ends_with("_excess")) %>%
    pivot_longer(-Date, names_to = "Currency", values_to = "ExcessReturn") %>%
    mutate(Currency = sub("_excess", "", Currency))

alpha_portfolios <- alpha_weights %>%
    left_join(returns_long, by = c("Date", "Currency")) %>%
    group_by(Date) %>%
    summarise(
        Ret_d = sum(w_d * ExcessReturn, na.rm = TRUE),   # delta-based return
        Ret_a = sum(w_a * ExcessReturn, na.rm = TRUE),    # alpha-based return
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

a_ewma_portfolios <- alpha_ewma_long %>%
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

a_sma_portfolios <- alpha_sma_long %>%
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
    
    alpha_long <- wide_data %>%
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

## Regressions for CAPM? check alpha
## return profiles - volatility, drawdown, Sharpe, etc


