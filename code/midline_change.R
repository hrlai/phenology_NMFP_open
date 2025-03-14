midline_change_spp <- 
    A_delta_taxa$Taxa[, , c("newleaf_Tmin", "newleaf_rain_wet",
                            "noleaf_Tmin", "noleaf_rain_wet",
                            "flowers_Tmin", "flowers_rain_wet",
                            "maturefruit_Tmin", "maturefruit_rain_wet")] %>% 
    reshape2::melt() %>% 
    group_by(Var2, Var3) %>% 
    median_qi(value, .width = 0.89) %>% 
    separate(Var3, c("Phenology", "Climate"), sep = "_", extra = "merge") %>% 
    select(Taxa = Var2, Phenology, Climate, value, .lower, .upper) %>% 
    pivot_wider(names_from = Climate,
                values_from = c(value, .lower, .upper),
                names_glue = "{Climate}_{.value}")

midline_change <- 
    fixef(mod, probs = c(0.055, 0.945)) %>% 
    as.data.frame() %>% 
    rownames_to_column("Param") %>% 
    filter(Param %in% c("newleaf_Tmin", "newleaf_rain_wet",
                        "noleaf_Tmin", "noleaf_rain_wet",
                        "flowers_Tmin", "flowers_rain_wet",
                        "maturefruit_Tmin", "maturefruit_rain_wet")) %>% 
    separate(Param, c("Phenology", "Climate"), sep = "_", extra = "merge") %>% 
    pivot_wider(names_from = Climate,
                values_from = c(Estimate, Est.Error, Q5.5, Q94.5),
                names_glue = "{Climate}_{.value}") 

# Plot --------------------------------------------------------------------

p_midline_noleaf <- 
    midline_change_spp %>% 
    filter(Phenology == "noleaf") %>% 
    left_join(A_delta_noleaf %>% select(Taxa, A_diff)) %>% 
    ggplot() +
    geom_errorbar(aes(x = -rain_wet_value,
                      ymin = Tmin_.lower,
                      ymax = Tmin_.upper),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = Tmin_value,
                       xmin = -rain_wet_.lower,
                       xmax = -rain_wet_.upper),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-rain_wet_value, Tmin_value,
                   fill = A_diff < 0),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-1.25, 1.4), xlim = c(-1.25, 2), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on midline")),
         y = expression(paste("Warming ", T[min], " effect on midline")),
         title = "(a) Leaf shedding") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title.position = "plot")

p_midline_newleaf <- 
    midline_change_spp %>% 
    filter(Phenology == "newleaf") %>% 
    left_join(A_delta_newleaf %>% select(Taxa, A_diff)) %>% 
    ggplot() +
    geom_errorbar(aes(x = -rain_wet_value,
                      ymin = Tmin_.lower,
                      ymax = Tmin_.upper),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = Tmin_value,
                       xmin = -rain_wet_.lower,
                       xmax = -rain_wet_.upper),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-rain_wet_value, Tmin_value,
                   fill = A_diff < 0),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-1.25, 1.4), xlim = c(-1.25, 2), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on midline")),
         y = expression(paste("Warming ", T[min], " effect on midline")),
         title = "(b) Leaf flush") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title.position = "plot")

p_midline_flowers <- 
    midline_change_spp %>% 
    filter(Phenology == "flowers") %>% 
    left_join(A_delta_flowers %>% select(Taxa, A_diff)) %>% 
    ggplot() +
    geom_errorbar(aes(x = -rain_wet_value,
                      ymin = Tmin_.lower,
                      ymax = Tmin_.upper),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = Tmin_value,
                       xmin = -rain_wet_.lower,
                       xmax = -rain_wet_.upper),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-rain_wet_value, Tmin_value,
                   fill = A_diff < 0),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-1.25, 1.4), xlim = c(-1.25, 2), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on midline")),
         y = expression(paste("Warming ", T[min], " effect on midline")),
         title = "(c) Flowering") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title.position = "plot")

p_midline_maturefruit <- 
    midline_change_spp %>% 
    filter(Phenology == "maturefruit") %>% 
    left_join(A_delta_maturefruit %>% select(Taxa, A_diff)) %>% 
    ggplot() +
    geom_errorbar(aes(x = -rain_wet_value,
                      ymin = Tmin_.lower,
                      ymax = Tmin_.upper),
                  width = 0,
                  colour = "grey") +
    geom_errorbarh(aes(y = Tmin_value,
                       xmin = -rain_wet_.lower,
                       xmax = -rain_wet_.upper),
                   height = 0,
                   colour = "grey") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
    geom_point(aes(-rain_wet_value, Tmin_value,
                   fill = A_diff < 0),
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    coord_cartesian(ylim = c(-1.25, 1.4), xlim = c(-1.25, 2), expand = FALSE) +
    labs(x = expression(paste("Drying ", R[wet], " effect on midline")),
         y = expression(paste("Warming ", T[min], " effect on midline")),
         title = "(d) Fruiting") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title.position = "plot")
