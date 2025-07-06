```{r estimate ARIMA for residuals}
library(forecast)

for (cur in currencies) {
    excess_var <- paste0(cur, "_excess")
    prem_var   <- paste0(cur, "_prem")
    
    df <- wide_data %>%
        dplyr::select(all_of(c(excess_var, prem_var))) %>%
        dplyr::filter(!is.na(.data[[excess_var]]), !is.na(.data[[prem_var]]))
    
    model <- lm(as.formula(paste0(excess_var, " ~ ", prem_var)), data = df)
    res <- resid(model)
    
    cat("-----", cur, "-----\n")
    auto_model <- auto.arima(res, max.p = 0, max.q = 30, seasonal = FALSE, ic = "bic")
    print(auto_model)
}

```

```{r Newey-West Standard Errors - everything insig}

lmnw_data <- list()
for (cur in currencies) {
    excess_var <- paste0(cur, "_excess")
    prem_var   <- paste0(cur, "_prem")
    
    df <- wide_data %>%
        dplyr::select(all_of(c(excess_var, prem_var))) %>%
        filter(!is.na(.data[[excess_var]]), !is.na(.data[[prem_var]]))
    
    model <- lm(as.formula(paste0(excess_var, " ~ ", prem_var)), data = df)
    nw_vcov <- NeweyWest(model, lag = 14, prewhite = FALSE)
    nw_se <- coeftest(model, vcov. = nw_vcov)
    
    summary_stats <- tibble(
        Currency   = cur,
        Intercept  = coef(model)[1],
        SE.Intercept = nw_se[1, "Std. Error"],
        Slope      = coef(model)[2],
        SE.Slope   = nw_se[2, "Std. Error"],
        R2         = summary(model)$r.squared,
        N          = nobs(model)
    )
    
    fitted_df <- tibble(
        Currency = cur,
        Fitted   = fitted(model),
        Residual = resid(model)
    )
    
    lmnw_data[[cur]] <- list(summary = summary_stats, residuals = fitted_df)
}
lmnwresults_df <- bind_rows(list(lmnw_data$GBP$summary, lmnw_data$EUR$summary, lmnw_data$CAD$summary, lmnw_data$JPY$summary))
print(lmnwresults_df)

```