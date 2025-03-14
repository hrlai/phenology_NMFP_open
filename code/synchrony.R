require(bpnreg)

climate_dat <- 
    mod$data %>% 
    data_grid(rain_wet = seq_range(rain_wet, n = 20),
              Tmin = 0) %>% 
    mutate(scenario = "rain_wet") %>% 
    bind_rows(
        mod$data %>% 
            data_grid(rain_wet = 0,
                      Tmin = seq_range(Tmin, n = 20)) %>% 
            mutate(scenario = "Tmin")
    )


# Calculate synchrony over years -----------------------------------------

synchrony_trend <- 
    coefs_species %>% 
    expand_grid(
        climate_dat %>% 
            rename(X_Tmin = Tmin,
                   X_rain_wet = rain_wet)) %>% 
    mutate(
        alpha = 
            Intercept +
            rain_wet*X_rain_wet + Tmin*X_Tmin,
        beta1 = 
            C1 +
            `C1:rain_wet`*X_rain_wet + `C1:Tmin`*X_Tmin,
        beta2 = 
            C2 +
            `C2:rain_wet`*X_rain_wet + `C2:Tmin`*X_Tmin,
        gamma1 = 
            S1 +
            `S1:rain_wet`*X_rain_wet + `S1:Tmin`*X_Tmin,
        gamma2 = 
            S2 +
            `S2:rain_wet`*X_rain_wet + `S2:Tmin`*X_Tmin
    ) %>% 
    # convert to amplitude and phase 
    mutate(
        d1 = atan(abs(gamma1 / beta1)),
        d2 = atan(abs(gamma2 / beta2)),
        # correction following Bingham et al. (1982) Chronobiologia
        # but I multiplied by -1
        # because the acrophase has negative sign? See Fig. 4a 
        d1 = ifelse(beta1 >= 0 & gamma1 > 0, d1,
                    ifelse(beta1 < 0  & gamma1 >= 0, pi - d1,
                           ifelse(beta1 <= 0 & gamma1 < 0, pi + d1,
                                  2 * pi - d1))),
        d2 = ifelse(beta2 >= 0 & gamma2 > 0, d2,
                    ifelse(beta2 < 0  & gamma2 >= 0, pi - d2,
                           ifelse(beta2 <= 0 & gamma2 < 0, pi + d2,
                                  2 * pi - d2)))
    ) %>% 
    # convert phase to circular object
    mutate_at(vars(starts_with("d")), ~ circular(., units = "radians", modulo = "2pi")) %>% 
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    left_join(A_dominant_newleaf) %>% 
    mutate(d = ifelse(Dominant == 1, d1, d2)) %>% 
    group_by(phenology, scenario, X_rain_wet, X_Tmin, .chain, .iteration, .draw) %>% 
    summarise(sd_d = sd.circular(d)) %>% 
    ungroup %>% 
    mutate(sync = 1 / sd_d) %>% 
    group_by(phenology, scenario, X_rain_wet, X_Tmin) %>% 
    median_qi(sync, .width = 0.89) %>% 
    ungroup %>% 
    mutate(phenology = fct_relevel(phenology,
                                   "noleaf", "newleaf", "flowers", "maturefruit"),
           phenology = fct_recode(phenology,
                                  `Leaf shedding` = "noleaf",
                                  `Leaf flush` = "newleaf", 
                                  `Flowering` = "flowers", 
                                  `Fruiting` = "maturefruit")) %>% 
    mutate(rain_wet = X_rain_wet * sd(dat$rain_wet) + mean(dat$rain_wet),
           Tmin = X_Tmin * sd(dat$Tmin) + mean(dat$Tmin),
           .keep = "unused")






# Plot --------------------------------------------------------------------

p_sync_rain <- 
    synchrony_trend %>% 
    filter(scenario == "rain_wet") %>% 
    ggplot() +
    facet_wrap(~ phenology, nrow = 1) +
    geom_ribbon(aes(rain_wet, ymin = .lower, ymax = .upper),
                fill = "grey") +
    geom_line(aes(rain_wet, sync)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 1), expand = c(0, 0)) +
    labs(x = expression(paste(R[wet], " (mm)")), 
         y = "Synchrony") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          plot.margin = margin(0, 1, 1, 1),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 11))
p_sync_temp <- 
    synchrony_trend %>% 
    filter(scenario == "Tmin") %>% 
    ggplot() +
    facet_wrap(~ phenology, nrow = 1) +
    geom_ribbon(aes(Tmin, ymin = .lower, ymax = .upper),
                fill = "grey") +
    geom_line(aes(Tmin, sync)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 1), expand = c(0, 0)) +
    labs(x = expression(paste(T[min], " (°C)")), 
         y = "Synchrony") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          plot.margin = margin(0, 1, 1, 1),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 11))
