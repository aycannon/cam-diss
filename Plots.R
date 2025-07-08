quad_quad_loss <- function(eps, alpha){
    ifelse(
        eps >= 0,
        alpha        *  eps^2,     # positive side, weight = alpha
        (1 - alpha)  *  (-eps)^2   # negative side, weight = 1–alpha
    )
}

# build a data.frame of epsilons and α’s
loss_df <- expand.grid(
    eps   = seq(-1, 1, length.out = 501),
    alpha = c(0.1, 0.3, 0.5, 0.7, 0.9)
)

# compute the loss
loss_df <- loss_df %>%
    mutate(L = quad_quad_loss(eps, alpha))

# plot
ggplot(loss_df, aes(x = eps, y = L, colour = factor(alpha))) +
    geom_line(alpha = 0.8) +
    scale_colour_viridis_d(
        name   = "alpha",
        option = "D"
    ) +
    labs(
        x     = expression(epsilon),
        y     = expression(L(epsilon))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

### Needs main code to be ran first

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








#### animated alphas
library(gganimate)
eps_seq <- seq(-1,1,length.out = 500)

animate_alphas <- alpha_ewma_long %>%
    filter(!is.na(Alpha)) %>%
    dplyr::select(Date, Currency, Alpha) %>%
    crossing(eps = eps_seq) %>%
    mutate(Loss = quad_quad_loss(eps, Alpha))

anim <- ggplot(animate_alphas, aes(x = eps, y = Loss, color = Currency)) +
    geom_line() +
    labs(
        title = "Dynamic Asymmetric Loss Functions",
        subtitle = "Date: {frame_time}",
        x = "Error",
        y = "L(Error)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    transition_time(Date) +
    ease_aes("linear") 

animate(anim, nframes = 200, fps = 10, renderer = gifski_renderer("loss_animation.gif"))
