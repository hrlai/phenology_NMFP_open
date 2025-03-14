fourier_cf <- 
    coefs_species %>% 
    expand_grid(
        mod$data %>% 
            distinct(year, Tmin, rain_wet) %>% 
            rename(X_Tmin = Tmin,
                   X_rain_wet = rain_wet)) %>% 
    left_join(coefs_year) %>% 
    mutate(
        # full submodel
        alpha = 
            Intercept +
            rain_wet*X_rain_wet + Tmin*X_Tmin + 
            Intercept_year,
        beta1 = 
            C1 +
            `C1:rain_wet`*X_rain_wet + `C1:Tmin`*X_Tmin + 
            C1_year,
        beta2 =
            C2 +
            `C2:rain_wet`*X_rain_wet + `C2:Tmin`*X_Tmin + 
            C2_year,
        gamma1 = 
            S1 +
            `S1:rain_wet`*X_rain_wet + `S1:Tmin`*X_Tmin + 
            S1_year,
        gamma2 = 
            S2 +
            `S2:rain_wet`*X_rain_wet + `S2:Tmin`*X_Tmin + 
            S2_year,
        # climate components only
        alpha_clim = 
            Intercept +
            rain_wet*X_rain_wet + Tmin*X_Tmin,
        beta1_clim = 
            C1 +
            `C1:rain_wet`*X_rain_wet + `C1:Tmin`*X_Tmin,
        beta2_clim = 
            C2 +
            `C2:rain_wet`*X_rain_wet + `C2:Tmin`*X_Tmin,
        gamma1_clim = 
            S1 +
            `S1:rain_wet`*X_rain_wet + `S1:Tmin`*X_Tmin,
        gamma2_clim = 
            S2 +
            `S2:rain_wet`*X_rain_wet + `S2:Tmin`*X_Tmin,
        .keep = "unused"
    ) %>% 
    # convert to amplitude and phase 
    mutate(
        A1 = sqrt(beta1^2 + gamma1^2),
        A2 = sqrt(beta2^2 + gamma2^2),
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
                                  2 * pi - d2))) ,
        A1_clim = sqrt(beta1_clim^2 + gamma1_clim^2),
        A2_clim = sqrt(beta2_clim^2 + gamma2_clim^2),
        d1_clim = atan(abs(gamma1_clim / beta1_clim)),
        d2_clim = atan(abs(gamma2_clim / beta2_clim)),
        d1_clim = ifelse(beta1_clim >= 0 & gamma1_clim > 0, d1_clim,
                         ifelse(beta1_clim < 0  & gamma1_clim >= 0, pi - d1_clim,
                                ifelse(beta1_clim <= 0 & gamma1_clim < 0, pi + d1_clim,
                                       2 * pi - d1_clim))),
        d2_clim = ifelse(beta2_clim >= 0 & gamma2_clim > 0, d2_clim,
                         ifelse(beta2_clim < 0  & gamma2_clim >= 0, pi - d2_clim,
                                ifelse(beta2_clim <= 0 & gamma2_clim < 0, pi + d2_clim,
                                       2 * pi - d2_clim)))
    ) %>% 
    # convert phase to circular object
    mutate_at(vars(starts_with("d")), ~ circular(., units = "radians", modulo = "2pi")) %>%
    group_by(phenology, Taxa, year) %>% 
    summarise_at(vars(alpha, alpha_clim, A1:d2_clim), median)

fourier_cf_plot <- 
    fourier_cf %>% 
    ungroup %>% 
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    left_join(A_dominant_newleaf) %>% 
    mutate(A = ifelse(Dominant == 1, A1, A2),
           d = ifelse(Dominant == 1, d1, d2),
           .keep = "unused") %>% 
    select(phenology, Name, year, alpha, A, d) %>% 
    pivot_longer(cols = c(alpha, A, d),
                 names_to = "component",
                 values_to = "median") %>% 
    mutate(phenology = fct_relevel(phenology,
                                   "noleaf", "newleaf", "flowers", "maturefruit"),
           phenology = fct_recode(phenology,
                                  `Leaf shedding` = "noleaf",
                                  `Leaf flush` = "newleaf", 
                                  `Flowering` = "flowers", 
                                  `Fruiting` = "maturefruit")) 

p_ts_midline <- 
    fourier_cf_plot %>% 
    filter(component == "alpha") %>% 
    ggplot() +
    facet_wrap(~ phenology, 
               strip.position = "left",
               ncol = 1) +
    geom_jitter(aes(year, median),
                pch = 16, alpha = 0.2, size = 1, 
                height = 0, width = 0.2) +
    scale_x_continuous(breaks = c(2004, 2010, 2014, 2020)) +
    scale_y_continuous(breaks = c(-4, 0, 4)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "in"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black", angle = 90, hjust = 0.5),
          axis.title = element_blank(),
          legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          plot.margin = margin(20, 1, 1, 1),
          plot.title.position = "plot")

p_ts_amplitude <- 
    fourier_cf_plot %>% 
    filter(component == "A") %>% 
    ggplot() +
    facet_wrap(~ phenology, 
               strip.position = "left",
               ncol = 1) +
    geom_jitter(aes(year, median),
                pch = 16, alpha = 0.2, size = 1, 
                height = 0, width = 0.2) +
    scale_x_continuous(breaks = c(2004, 2010, 2014, 2020)) +
    scale_y_sqrt(breaks = c(2, 10, 20)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "in"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black", angle = 90, hjust = 0.5),
          axis.title = element_blank(),
          legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          plot.margin = margin(20, 1, 1, 1),
          plot.title.position = "plot")

p_ts_phase <- 
    fourier_cf_plot %>% 
    filter(component == "d") %>% 
    ggplot() +
    facet_wrap(~ phenology, 
               strip.position = "left",
               ncol = 1) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = 2*pi/12*5, ymax = 2*pi/12*9,
             fill = "blue", alpha = 0.2) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = 2*pi/12*11, ymax = 2*pi/12*12,
             fill = "red", alpha = 0.2) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = 2*pi/12*0, ymax = 2*pi/12*2,
             fill = "red", alpha = 0.2) +
    geom_jitter(aes(year, median),
                pch = 16, alpha = 0.2, size = 1,
                height = 0, width = 0.2) +
    scale_x_continuous(breaks = c(2004, 2010, 2014, 2020)) +
    scale_y_continuous(breaks = seq(0, 2*pi, length.out = 5)[-1],
                       labels = seq(0, 12, length.out = 5)[-1]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "in"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black", angle = 90, hjust = 0.5),
          axis.title = element_blank(),
          legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          plot.margin = margin(20, 1, 1, 1),
          plot.title.position = "plot")
