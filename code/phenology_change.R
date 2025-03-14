require(bpnreg)


# Phenology change with rainfall ------------------------------------------

# Change from min to max climate values
rain_cf <- 
    mod$data %>% 
    data_grid(rain_wet = seq_range(rain_wet, n = 2),
              Tmin = 0)
temp_cf <- 
    mod$data %>% 
    data_grid(rain_wet = 0,
              Tmin = seq_range(Tmin, n = 2))

# calculate species responses
newleaf_rain_cf <- 
    mod %>% 
    spread_draws(r_Taxa__newleaf[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__newleaf = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_newleaf_C1, 
                         b_newleaf_S1, 
                         `b_newleaf_C1:rain_wet`, `b_newleaf_C1:Tmin`,
                         `b_newleaf_S1:rain_wet`, `b_newleaf_S1:Tmin`,
                         b_newleaf_C2, 
                         b_newleaf_S2, 
                         `b_newleaf_C2:rain_wet`, `b_newleaf_C2:Tmin`,
                         `b_newleaf_S2:rain_wet`, `b_newleaf_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_newleaf_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__newleaf) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(rain_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, rain_wet, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(rain_wet = cut(rain_wet, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = rain_wet,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))

noleaf_rain_cf <- 
    mod %>% 
    spread_draws(r_Taxa__noleaf[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__noleaf = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_noleaf_C1, 
                         b_noleaf_S1, 
                         `b_noleaf_C1:rain_wet`, `b_noleaf_C1:Tmin`,
                         `b_noleaf_S1:rain_wet`, `b_noleaf_S1:Tmin`,
                         b_noleaf_C2, 
                         b_noleaf_S2, 
                         `b_noleaf_C2:rain_wet`, `b_noleaf_C2:Tmin`,
                         `b_noleaf_S2:rain_wet`, `b_noleaf_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_noleaf_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__noleaf) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(rain_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, rain_wet, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(rain_wet = cut(rain_wet, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = rain_wet,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))

flowers_rain_cf <- 
    mod %>% 
    spread_draws(r_Taxa__flowers[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__flowers = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_flowers_C1, 
                         b_flowers_S1, 
                         `b_flowers_C1:rain_wet`, `b_flowers_C1:Tmin`,
                         `b_flowers_S1:rain_wet`, `b_flowers_S1:Tmin`,
                         b_flowers_C2, 
                         b_flowers_S2, 
                         `b_flowers_C2:rain_wet`, `b_flowers_C2:Tmin`,
                         `b_flowers_S2:rain_wet`, `b_flowers_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_flowers_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__flowers) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(rain_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, rain_wet, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(rain_wet = cut(rain_wet, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = rain_wet,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))

maturefruit_rain_cf <- 
    mod %>% 
    spread_draws(r_Taxa__maturefruit[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__maturefruit = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_maturefruit_C1, 
                         b_maturefruit_S1, 
                         `b_maturefruit_C1:rain_wet`, `b_maturefruit_C1:Tmin`,
                         `b_maturefruit_S1:rain_wet`, `b_maturefruit_S1:Tmin`,
                         b_maturefruit_C2, 
                         b_maturefruit_S2, 
                         `b_maturefruit_C2:rain_wet`, `b_maturefruit_C2:Tmin`,
                         `b_maturefruit_S2:rain_wet`, `b_maturefruit_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_maturefruit_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__maturefruit) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(rain_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, rain_wet, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(rain_wet = cut(rain_wet, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = rain_wet,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))




# Phenology change with temperature ---------------------------------------

newleaf_temp_cf <- 
    mod %>% 
    spread_draws(r_Taxa__newleaf[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__newleaf = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_newleaf_C1, 
                         b_newleaf_S1, 
                         `b_newleaf_C1:rain_wet`, `b_newleaf_C1:Tmin`,
                         `b_newleaf_S1:rain_wet`, `b_newleaf_S1:Tmin`,
                         b_newleaf_C2, 
                         b_newleaf_S2, 
                         `b_newleaf_C2:rain_wet`, `b_newleaf_C2:Tmin`,
                         `b_newleaf_S2:rain_wet`, `b_newleaf_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_newleaf_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__newleaf) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(temp_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, Tmin, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(Tmin = cut(Tmin, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = Tmin,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))

noleaf_temp_cf <- 
    mod %>% 
    spread_draws(r_Taxa__noleaf[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__noleaf = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_noleaf_C1, 
                         b_noleaf_S1, 
                         `b_noleaf_C1:rain_wet`, `b_noleaf_C1:Tmin`,
                         `b_noleaf_S1:rain_wet`, `b_noleaf_S1:Tmin`,
                         b_noleaf_C2, 
                         b_noleaf_S2, 
                         `b_noleaf_C2:rain_wet`, `b_noleaf_C2:Tmin`,
                         `b_noleaf_S2:rain_wet`, `b_noleaf_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_noleaf_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__noleaf) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(temp_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, Tmin, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(Tmin = cut(Tmin, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = Tmin,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))

flowers_temp_cf <- 
    mod %>% 
    spread_draws(r_Taxa__flowers[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__flowers = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_flowers_C1, 
                         b_flowers_S1, 
                         `b_flowers_C1:rain_wet`, `b_flowers_C1:Tmin`,
                         `b_flowers_S1:rain_wet`, `b_flowers_S1:Tmin`,
                         b_flowers_C2, 
                         b_flowers_S2, 
                         `b_flowers_C2:rain_wet`, `b_flowers_C2:Tmin`,
                         `b_flowers_S2:rain_wet`, `b_flowers_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_flowers_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__flowers) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(temp_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, Tmin, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(Tmin = cut(Tmin, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = Tmin,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))

maturefruit_temp_cf <- 
    mod %>% 
    spread_draws(r_Taxa__maturefruit[Taxa,Param]) %>% 
    filter(Param %in% c("C1", "S1", 
                        "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                        "C2", "S2", 
                        "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
    # add average species with zero random effects
    bind_rows(
        data.frame(
            Taxa = "Community", 
            .chain = rep(1:4, each = 1000),
            .iteration = rep(1:1000, 4),
            .draw = 1:4000
        ) %>% 
            expand_grid(Param = c("C1", "S1", 
                                  "C1:rain_wet", "C1:Tmin", "S1:rain_wet", "S1:Tmin",
                                  "C2", "S2", 
                                  "C2:rain_wet", "C2:Tmin", "S2:rain_wet", "S2:Tmin")) %>% 
            mutate(r_Taxa__maturefruit = 0)
    ) %>% 
    # join population-level effects
    left_join(
        mod %>% 
            gather_draws(b_maturefruit_C1, 
                         b_maturefruit_S1, 
                         `b_maturefruit_C1:rain_wet`, `b_maturefruit_C1:Tmin`,
                         `b_maturefruit_S1:rain_wet`, `b_maturefruit_S1:Tmin`,
                         b_maturefruit_C2, 
                         b_maturefruit_S2, 
                         `b_maturefruit_C2:rain_wet`, `b_maturefruit_C2:Tmin`,
                         `b_maturefruit_S2:rain_wet`, `b_maturefruit_S2:Tmin`) %>% 
            mutate(Param = str_remove(.variable, "b_maturefruit_"))
    ) %>% 
    mutate(Estimate = .value + r_Taxa__maturefruit) %>% 
    ungroup %>% 
    select(Taxa, Param, .chain, .iteration, .draw, Estimate) %>% 
    pivot_wider(names_from = Param,
                values_from = Estimate) %>% 
    # filter(.draw %in% 1:1000) %>% 
    expand_grid(temp_cf) %>% 
    mutate(C1_cf = C1 + `C1:rain_wet`*rain_wet + `C1:Tmin`*Tmin,
           S1_cf = S1 + `S1:rain_wet`*rain_wet + `S1:Tmin`*Tmin,
           C2_cf = C2 + `C2:rain_wet`*rain_wet + `C2:Tmin`*Tmin,
           S2_cf = S2 + `S2:rain_wet`*rain_wet + `S2:Tmin`*Tmin,
           # convert to amplitude and phase 
           A1 = sqrt(C1_cf^2 + S1_cf^2),
           A2 = sqrt(C2_cf^2 + S2_cf^2),
           d1 = atan(abs(S1_cf / C1_cf)),
           d2 = atan(abs(S2_cf / C2_cf)),
           # correction following Bingham et al. (1982) Chronobiologia
           # but I multiplied by -1
           # because the acrophase has negative sign? See Fig. 4a 
           d1 = ifelse(C1_cf >= 0 & S1_cf > 0, d1,
                       ifelse(C1_cf < 0  & S1_cf >= 0, pi - d1,
                              ifelse(C1_cf <= 0 & S1_cf < 0, pi + d1,
                                     2 * pi - d1))),
           d2 = ifelse(C2_cf >= 0 & S2_cf > 0, d2,
                       ifelse(C2_cf < 0  & S2_cf >= 0, pi - d2,
                              ifelse(C2_cf <= 0 & S2_cf < 0, pi + d2,
                                     2 * pi - d2)))) %>% 
    # calculate change
    select(Taxa, .chain, .iteration, .draw, Tmin, A1, A2, d1, d2) %>% 
    # select the more important climate effect (previous vs. current)
    mutate(Taxa = str_replace(Taxa, "\\.", " ")) %>% 
    mutate(Tmin = cut(Tmin, 2, c("low", "high"))) %>% 
    pivot_wider(names_from = Tmin,
                values_from = c(A1, A2, d1, d2)) %>% 
    mutate(A1_change = A1_high - A1_low,
           A2_change = A2_high - A2_low,
           d1_change = d1_high - d1_low,
           d2_change = d2_high - d2_low) %>% 
    mutate(d1_change_circ = circular(d1_change, units = "radians", modulo = "asis"),
           d2_change_circ = circular(d2_change, units = "radians", modulo = "asis")) %>% 
    group_by(Taxa) %>% 
    summarise(A1_change.lower = quantile(A1_change, 0.055),
              A1_change.upper = quantile(A1_change, 0.945),
              A1_change = median(A1_change),
              A2_change.lower = quantile(A2_change, 0.055),
              A2_change.upper = quantile(A2_change, 0.945),
              A2_change = median(A2_change),
              d1_change.lower = quantile(d1_change, 0.055),
              d1_change.upper = quantile(d1_change, 0.945),
              d1_change = median(d1_change),
              d2_change.lower = quantile(d2_change, 0.055),
              d2_change.upper = quantile(d2_change, 0.945),
              d2_change = median(d2_change),
              # d1_change_circ.lower = quantile(d1_change_circ, 0.055),
              # d1_change_circ.upper = quantile(d1_change_circ, 0.945),
              d1_change_circ.lower = hpd_est(d1_change_circ)[1],
              d1_change_circ.upper = hpd_est(d1_change_circ)[2],
              d1_change_circ = median(d1_change_circ),
              d2_change_circ.lower = hpd_est(d2_change_circ)[1],
              d2_change_circ.upper = hpd_est(d2_change_circ)[2],
              d2_change_circ = median(d2_change_circ)) %>%
    # median_qi(A1_change, A2_change, d1_change, d2_change, d1_change_circ, d2_change_circ,
    #           .width = 0.89) %>%
    mutate_at(vars(d1_change_circ, d1_change_circ.lower, d1_change_circ.upper), ~(. * 12 / (2 * pi * 2))) %>%
    mutate_at(vars(d2_change_circ, d2_change_circ.lower, d2_change_circ.upper), ~(. * 12 / (2 * pi * 4)))






# Combine -----------------------------------------------------------------

newleaf_cf <- 
    bind_rows(newleaf_rain_cf %>% mutate(Climate = "rain_wet"),
              newleaf_temp_cf %>% mutate(Climate = "Tmin")) 
noleaf_cf <- 
    bind_rows(noleaf_rain_cf %>% mutate(Climate = "rain_wet"),
              noleaf_temp_cf %>% mutate(Climate = "Tmin")) 
flowers_cf <- 
    bind_rows(flowers_rain_cf %>% mutate(Climate = "rain_wet"),
              flowers_temp_cf %>% mutate(Climate = "Tmin"))
maturefruit_cf <- 
    bind_rows(maturefruit_rain_cf %>% mutate(Climate = "rain_wet"),
              maturefruit_temp_cf %>% mutate(Climate = "Tmin"))





# Amplitude change --------------------------------------------------------

A_newleaf_change <- 
    newleaf_cf %>% 
    left_join(A_dominant_newleaf) %>% 
    mutate(A_change = ifelse(Dominant == 1, A1_change, A2_change),
           A_change.lower = ifelse(Dominant == 1, A1_change.lower, A2_change.lower),
           A_change.upper = ifelse(Dominant == 1, A1_change.upper, A2_change.upper)) %>% 
    select(Taxa, Name, starts_with("A_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("A_change")) 
A_noleaf_change <- 
    noleaf_cf %>% 
    left_join(A_dominant_noleaf) %>% 
    mutate(A_change = ifelse(Dominant == 1, A1_change, A2_change),
           A_change.lower = ifelse(Dominant == 1, A1_change.lower, A2_change.lower),
           A_change.upper = ifelse(Dominant == 1, A1_change.upper, A2_change.upper)) %>% 
    select(Taxa, Name, starts_with("A_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("A_change")) 
A_flowers_change <- 
    flowers_cf %>% 
    left_join(A_dominant_flowers) %>% 
    mutate(A_change = ifelse(Dominant == 1, A1_change, A2_change),
           A_change.lower = ifelse(Dominant == 1, A1_change.lower, A2_change.lower),
           A_change.upper = ifelse(Dominant == 1, A1_change.upper, A2_change.upper)) %>% 
    select(Taxa, Name, starts_with("A_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("A_change")) 
A_maturefruit_change <- 
    maturefruit_cf %>% 
    left_join(A_dominant_maturefruit) %>% 
    mutate(A_change = ifelse(Dominant == 1, A1_change, A2_change),
           A_change.lower = ifelse(Dominant == 1, A1_change.lower, A2_change.lower),
           A_change.upper = ifelse(Dominant == 1, A1_change.upper, A2_change.upper)) %>% 
    select(Taxa, Name, starts_with("A_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("A_change")) 


p_A_noleaf_change <- 
    ggplot(A_noleaf_change) +
    geom_errorbar(aes(x = -A_change_rain_wet,
                      ymin = A_change.lower_Tmin,
                      ymax = A_change.upper_Tmin),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = A_change_Tmin,
                       xmin = -A_change.lower_rain_wet,
                       xmax = -A_change.upper_rain_wet),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-A_change_rain_wet, A_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    # geom_errorbar(data = A_noleaf_change %>% filter(Taxa == "Community"),
    #               aes(x = A_change_rain_wet,
    #                   ymin = A_change.lower_Tmin,
    #                   ymax = A_change.upper_Tmin),
    #               width = 0, size = 1,
    #               colour = "red") +
    # geom_errorbarh(data = A_noleaf_change %>% filter(Taxa == "Community"),
    #                aes(y = A_change_Tmin,
    #                    xmin = A_change.lower_rain_wet,
    #                    xmax = A_change.upper_rain_wet),
    #                height = 0, size = 1,
    #                colour = "red") +
    # geom_point(data = A_noleaf_change %>% filter(Taxa == "Community"),
    #            aes(A_change_rain_wet, A_change_Tmin),
    #            fill = "red",
    #            size = 3,
    #            pch = 23) +
    coord_cartesian(ylim = c(-15, 15), xlim = c(-15, 15), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on amplitude")),
         y = expression(paste("Warming ", T[min], " effect on amplitude")),
         title = "(a) Leaf shedding") +
    theme_classic() +
    theme(legend.position = "none")

p_A_newleaf_change <- 
    ggplot(A_newleaf_change) +
    geom_errorbar(aes(x = -A_change_rain_wet,
                      ymin = A_change.lower_Tmin,
                      ymax = A_change.upper_Tmin),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = A_change_Tmin,
                       xmin = -A_change.lower_rain_wet,
                       xmax = -A_change.upper_rain_wet),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-A_change_rain_wet, A_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-15, 15), xlim = c(-15, 15), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on amplitude")),
         y = expression(paste("Warming ", T[min], " effect on amplitude")),
         title = "(b) Leaf flush") +
    theme_classic() +
    theme(legend.position = "none")

p_A_flowers_change <- 
    ggplot(A_flowers_change) +
    geom_errorbar(aes(x = -A_change_rain_wet,
                      ymin = A_change.lower_Tmin,
                      ymax = A_change.upper_Tmin),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = A_change_Tmin,
                       xmin = -A_change.lower_rain_wet,
                       xmax = -A_change.upper_rain_wet),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-A_change_rain_wet, A_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) + 
    coord_cartesian(ylim = c(-15, 15), xlim = c(-15, 15), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on amplitude")),
         y = expression(paste("Warming ", T[min], " effect on amplitude")),
         title = "(c) Flowering") +
    theme_classic() +
    theme(legend.position = "none")

p_A_maturefruit_change <- 
    ggplot(A_maturefruit_change) +
    geom_errorbar(aes(x = -A_change_rain_wet,
                      ymin = A_change.lower_Tmin,
                      ymax = A_change.upper_Tmin),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = A_change_Tmin,
                       xmin = -A_change.lower_rain_wet,
                       xmax = -A_change.upper_rain_wet),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-A_change_rain_wet, A_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-15, 15), xlim = c(-15, 15), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on amplitude")),
         y = expression(paste("Warming ", T[min], " effect on amplitude")),
         title = "(d) Fruiting") +
    theme_classic() +
    theme(legend.position = "none")




# Phase change --------------------------------------------------------

d_noleaf_change <- 
    noleaf_cf %>% 
    left_join(A_dominant_noleaf) %>% 
    mutate(d_change = ifelse(Dominant == 1, d1_change_circ, d2_change_circ),
           d_change.lower = ifelse(Dominant == 1, d1_change_circ.lower, d2_change_circ.lower),
           d_change.upper = ifelse(Dominant == 1, d1_change_circ.upper, d2_change_circ.upper)) %>% 
    select(Taxa, Name, starts_with("d_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("d_change")) %>% 
    mutate(colour = ifelse(Dominant == 2, "white", "black")) %>% 
    # for illustration purpose
    mutate(lci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.lower_rain_wet, d_change.upper_rain_wet),
           uci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.upper_rain_wet, d_change.lower_rain_wet),
           lci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, 6),
           uci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, -6),
           lci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.lower_Tmin, d_change.upper_Tmin),
           uci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.upper_Tmin, d_change.lower_Tmin),
           lci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, 6),
           uci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, -6)) 
d_newleaf_change <- 
    newleaf_cf %>% 
    left_join(A_dominant_newleaf) %>% 
    mutate(d_change = ifelse(Dominant == 1, d1_change_circ, d2_change_circ),
           d_change.lower = ifelse(Dominant == 1, d1_change_circ.lower, d2_change_circ.lower),
           d_change.upper = ifelse(Dominant == 1, d1_change_circ.upper, d2_change_circ.upper)) %>% 
    select(Taxa, Name, starts_with("d_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("d_change")) %>% 
    mutate(colour = ifelse(Dominant == 2, "white", "black")) %>% 
    # for illustration purpose
    mutate(lci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.lower_rain_wet, d_change.upper_rain_wet),
           uci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.upper_rain_wet, d_change.lower_rain_wet),
           lci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, 6),
           uci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, -6),
           lci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.lower_Tmin, d_change.upper_Tmin),
           uci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.upper_Tmin, d_change.lower_Tmin),
           lci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, 6),
           uci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, -6))  
d_flowers_change <- 
    flowers_cf %>% 
    left_join(A_dominant_flowers) %>% 
    mutate(d_change = ifelse(Dominant == 1, d1_change_circ, d2_change_circ),
           d_change.lower = ifelse(Dominant == 1, d1_change_circ.lower, d2_change_circ.lower),
           d_change.upper = ifelse(Dominant == 1, d1_change_circ.upper, d2_change_circ.upper)) %>% 
    select(Taxa, Name, starts_with("d_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("d_change")) %>% 
    mutate(colour = ifelse(Dominant == 2, "white", "black")) %>% 
    # for illustration purpose
    mutate(lci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.lower_rain_wet, d_change.upper_rain_wet),
           uci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.upper_rain_wet, d_change.lower_rain_wet),
           lci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, 6),
           uci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, -6),
           lci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.lower_Tmin, d_change.upper_Tmin),
           uci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.upper_Tmin, d_change.lower_Tmin),
           lci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, 6),
           uci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, -6))  
d_maturefruit_change <- 
    maturefruit_cf %>% 
    left_join(A_dominant_maturefruit) %>% 
    mutate(d_change = ifelse(Dominant == 1, d1_change_circ, d2_change_circ),
           d_change.lower = ifelse(Dominant == 1, d1_change_circ.lower, d2_change_circ.lower),
           d_change.upper = ifelse(Dominant == 1, d1_change_circ.upper, d2_change_circ.upper)) %>% 
    select(Taxa, Name, starts_with("d_change"), Climate, Dominant) %>% 
    pivot_wider(names_from = Climate,
                values_from = starts_with("d_change")) %>% 
    mutate(colour = ifelse(Dominant == 2, "white", "black")) %>% 
    # for illustration purpose
    mutate(lci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.lower_rain_wet, d_change.upper_rain_wet),
           uci_rain_wet = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet,
                      d_change.upper_rain_wet, d_change.lower_rain_wet),
           lci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, 6),
           uci_rain_wet_plot = 
               ifelse(d_change.lower_rain_wet <= d_change_rain_wet & d_change.upper_rain_wet >= d_change_rain_wet, 
                      d_change_rain_wet, -6),
           lci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.lower_Tmin, d_change.upper_Tmin),
           uci_Tmin = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin,
                      d_change.upper_Tmin, d_change.lower_Tmin),
           lci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, 6),
           uci_Tmin_plot = 
               ifelse(d_change.lower_Tmin <= d_change_Tmin & d_change.upper_Tmin >= d_change_Tmin, 
                      d_change_Tmin, -6))  


p_d_noleaf_change <- 
    ggplot(d_noleaf_change) +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = lci_Tmin, yend = lci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = uci_Tmin, yend = uci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -lci_rain_wet, xend = -lci_rain_wet_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -uci_rain_wet, xend = -uci_rain_wet_plot),
                 colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-d_change_rain_wet, d_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-6, 6), xlim = c(-6, 6), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on phase")),
         y = expression(paste("Warming ", T[min], " effect on phase")),
         title = "(a) Leaf shedding") +
    theme_classic() +
    theme(legend.position = "none")

p_d_newleaf_change <- 
    ggplot(d_newleaf_change) +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = lci_Tmin, yend = lci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = uci_Tmin, yend = uci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -lci_rain_wet, xend = -lci_rain_wet_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -uci_rain_wet, xend = -uci_rain_wet_plot),
                 colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-d_change_rain_wet, d_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-6, 6), xlim = c(-6, 6), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on phase")),
         y = expression(paste("Warming ", T[min], " effect on phase")),
         title = "(b) Leaf flush") +
    theme_classic() +
    theme(legend.position = "none")

p_d_flowers_change <- 
    ggplot(d_flowers_change) +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = lci_Tmin, yend = lci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = uci_Tmin, yend = uci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -lci_rain_wet, xend = -lci_rain_wet_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -uci_rain_wet, xend = -uci_rain_wet_plot),
                 colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-d_change_rain_wet, d_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-6, 6), xlim = c(-6, 6), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on phase")),
         y = expression(paste("Warming ", T[min], " effect on phase")),
         title = "(c) Flowering") +
    theme_classic() +
    theme(legend.position = "none")

p_d_maturefruit_change <- 
    ggplot(d_maturefruit_change) +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = lci_Tmin, yend = lci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(x = -d_change_rain_wet, xend = -d_change_rain_wet, y = uci_Tmin, yend = uci_Tmin_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -lci_rain_wet, xend = -lci_rain_wet_plot),
                 colour = "grey") +
    geom_segment(aes(y = d_change_Tmin, yend = d_change_Tmin, x = -uci_rain_wet, xend = -uci_rain_wet_plot),
                 colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-d_change_rain_wet, d_change_Tmin,
                   fill = as.factor(Dominant)),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-6, 6), xlim = c(-6, 6), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on phase")),
         y = expression(paste("Warming ", T[min], " effect on phase")),
         title = "(d) Fruiting") +
    theme_classic() +
    theme(legend.position = "none")
