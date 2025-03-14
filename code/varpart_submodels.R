
# species coefs (fixed + random)
coefs_species <-
    mod %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(.variable = str_remove(.variable, "b_")) %>%
    separate(.variable, c("phenology", ".variable"), sep = "_", extra = "merge") %>%
    filter(!str_detect(.variable, "Intercept")) %>%
    right_join(
        mod %>%
            gather_draws(`r_Taxa__.*`[Taxa, Param], regex = TRUE) %>%
            mutate(.variable = str_remove(.variable, "r_Taxa__")) %>%
            rename(ranef = .value,
                   phenology = .variable,
                   .variable = Param)
    ) %>%
    mutate(.value = ifelse(is.na(.value) & .variable == "Intercept", 0, .value)) %>%
    mutate(coef = .value + ranef,
           .keep = "unused") %>%
    pivot_wider(names_from = .variable,
                values_from = coef)

# year ranefs
coefs_year <- 
    mod %>% 
    gather_draws(`r_year_.*`[year, Param], regex = TRUE) %>% 
    mutate(.variable = str_remove(.variable, "r_year__")) %>%
    rename(phenology = .variable) %>% 
    pivot_wider(names_from = Param,
                values_from = .value,
                names_glue = "{Param}_year")

# calculate equations 6-8
varpart_submodels <- 
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
            rain_wet*X_rain_wet + Tmin*X_Tmin + 
            Intercept_year,
        beta1 = 
            `C1:rain_wet`*X_rain_wet + `C1:Tmin`*X_Tmin + 
            C1_year,
        beta2 =
            `C2:rain_wet`*X_rain_wet + `C2:Tmin`*X_Tmin + 
            C2_year,
        gamma1 = 
            `S1:rain_wet`*X_rain_wet + `S1:Tmin`*X_Tmin + 
            S1_year,
        gamma2 = 
            `S2:rain_wet`*X_rain_wet + `S2:Tmin`*X_Tmin + 
            S2_year,
        # climate components only
        alpha_clim = 
            rain_wet*X_rain_wet + Tmin*X_Tmin,
        beta1_clim = 
            `C1:rain_wet`*X_rain_wet + `C1:Tmin`*X_Tmin,
        beta2_clim = 
            `C2:rain_wet`*X_rain_wet + `C2:Tmin`*X_Tmin,
        gamma1_clim = 
            `S1:rain_wet`*X_rain_wet + `S1:Tmin`*X_Tmin,
        gamma2_clim = 
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
    group_by(.chain, .iteration, .draw, phenology) %>% 
    summarise(r2_alpha = cor(alpha, alpha_clim)^2,
              r2_beta1 = cor(beta1, beta1_clim)^2,
              r2_beta2 = cor(beta2, beta2_clim)^2,
              r2_gamma1 = cor(gamma1, gamma1_clim)^2,
              r2_gamma2 = cor(gamma2, gamma2_clim)^2,
              r2_A1 = cor(A1, A1_clim)^2,
              r2_A2 = cor(A2, A2_clim)^2,
              r2_d1 = cor.circular(d1, d1_clim)^2,
              r2_d2 = cor.circular(d2, d2_clim)^2) %>% 
    pivot_longer(cols = starts_with("r2_"),
                 names_to = "Param",
                 names_prefix = "r2_",
                 values_to = "r2") %>% 
    group_by(phenology, Param) %>% 
    median_qi(r2, .width = 0.9)

# plot
p_varpart <- 
    varpart_submodels %>% 
    filter(Param %in% c("alpha", "A1", "A2", "d1", "d2")) %>% 
    mutate(phenology = fct_relevel(phenology,
                                   "noleaf", "newleaf", "flowers", "maturefruit"),
           phenology = fct_recode(phenology,
                                  `Leaf shedding` = "noleaf",
                                  `Leaf flush` = "newleaf", 
                                  `Flowering` = "flowers", 
                                  `Fruiting` = "maturefruit"),
           Param = fct_relevel(Param,
                               "alpha", "A1", "A2", "d1", "d2"),
           Param = fct_recode(Param,
                              `Midline` = "alpha",
                              `Amplitude\n(period 12)` = "A1", 
                              `Amplitude\n(period 6)` = "A2", 
                              `Phase\n(period 12)` = "d1", 
                              `Phase\n(period 6)` = "d2")) %>% 
    ggplot() +
    facet_wrap(~ Param, nrow = 1) +
    geom_errorbar(aes(phenology, ymin = .lower, ymax = .upper),
                  width = 0) +
    geom_point(aes(phenology, r2)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(y = "Proportion of variance\nexplained by climate") +
    coord_cartesian(clip = "off") +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, size = 11, colour = "black", hjust = 1, vjust = 0.5),
          strip.text.x = element_text(size = 11),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "in"))
