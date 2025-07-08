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

lm_data_monthly <- list()

for (cur in currencies) {
    excess <- paste0(cur, "_excess")
    prem   <- paste0(cur, "_prem")
    
    if (any(!c(excess, prem) %in% names(monthly_df1))) {
        message("Skipping ", cur, ": missing columns.")
        next
    }
    
    df <- monthly_df1 %>%
        dplyr::select(all_of(c(excess, prem))) %>%
        filter(!is.na(.data[[excess]]), !is.na(.data[[prem]]))
    
    f <- as.formula(paste0(excess, " ~ ", prem))
    m <- lm(f, data = df)
    sum <- summary(m)
    se  <- sum$coefficients[, "Std. Error"]
    
    nw_vcov <- NeweyWest(m, lag = 14, prewhite = FALSE)
    nw_se   <- coeftest(m, vcov. = nw_vcov)
    
    stats <- tibble(
        Currency         = cur,
        Intercept        = coef(m)[1],
        SE_Intercept     = se[1],
        SE_NW_Intercept  = nw_se[1, 2],
        Slope            = coef(m)[2],
        SE_Slope         = se[2],
        SE_NW_Slope      = nw_se[2, 2],
        R2               = sum$r.squared,
        AIC              = AIC(m),
        BIC              = BIC(m),
        N                = nobs(m)
    )
    
    resid_df <- list(
        Currency = cur,
        fitted   = fitted(m),
        resid    = resid(m)
    )
    
    lm_data_monthly[[cur]] <- list(
        model     = m,
        summary   = stats,
        residuals = resid_df
    )
}







