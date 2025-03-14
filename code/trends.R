# Counterfactual comparison between average year and drought year -------------

# climate full series
rain <- 
    read_csv("data/Gembu Rain 5x3.csv") %>% 
    # correct 2019 rainfall according to Silvio
    mutate(Jun = ifelse(year == 2019, 309.9, Jun),
           Aug = ifelse(year == 2019, 253.4, Aug)) %>%
    mutate_at(vars(Jan:Dec), ~ifelse(is.na(.), 0, .)) %>% 
    mutate(Wet = (Jun + Jul + Aug + Sep)/4) %>% 
    select(year, rain_wet = Wet)
Tmin <- 
    read_csv("data/ERA5 (1976_2023) MinT 5x3.csv") %>%
    rowwise() %>% 
    mutate(Tmin = min(c(dry, intermediate, wet))) %>% 
    ungroup %>% 
    select(year, Tmin)
climate_full <- 
    left_join(rain, Tmin) %>% 
    # remove 1995 (no data) and 2023 (incomplete data)
    filter(!year %in% c(1995, 2023))

ann_df <- 
    climate_full %>% 
    filter(year == "2014") %>% 
    mutate(year_lab = paste0("Year ", as.character(year))) %>%
    bind_rows(
        climate_full %>% 
            summarise_all(mean, na.rm = TRUE) %>% 
            mutate(year = NA,
                   year_lab = "Average year")
    )

# new data
dat_eg_spp <- 
    ann_df %>% 
    mutate(rain_wet = (rain_wet - mean(dat$rain_wet)) / sd(dat$rain_wet),
           Tmin = (Tmin - mean(dat$Tmin)) / sd(dat$Tmin)) %>% 
    expand_grid(Taxa = unique(mod$data$Taxa),
                month = seq(0, 12, length.out = 50),
                uid = NA,
                transect = NA)
cs_new <- make_c_s_mat(time_steps = dat_eg_spp$month, P = period)
cs_new <- do.call(cbind, lapply(cs_new, \(x) x[, 1:n_cs]))
colnames(cs_new) <- paste0(rep(c("C", "S"), each = n_cs), 
                           rep(1:n_cs, n_cs))
dat_eg_spp <- bind_cols(dat_eg_spp, cs_new)
pred_eg_spp <- fitted(mod, 
                      dpar = "mu",
                      scale = "linear",
                      newdata = dat_eg_spp,
                      allow_new_levels = TRUE,
                      summary = TRUE,
                      robust = TRUE,
                      ndraws = 1000)
dat_eg_spp <- 
    bind_cols(dat_eg_spp, pred_eg_spp[, "Estimate", ]) %>% 
    pivot_longer(cols = c(newleaf, noleaf, flowers, maturefruit),
                 names_to = "phenology",
                 values_to = "value") %>% 
    mutate(phenology = fct_relevel(phenology,
                                   "noleaf", "newleaf", "flowers", "maturefruit"),
           phenology = fct_recode(phenology,
                                  `Leaf shedding` = "noleaf",
                                  `Leaf flush` = "newleaf", 
                                  `Flowering` = "flowers", 
                                  `Fruiting` = "maturefruit")) 

# panel background
p_ts_bg <- data.frame(
    year_lab = unique(dat_eg_spp$year_lab),
    fill = c("springgreen", "white"),
    xmin = 0, xmax = 12,
    ymin = -Inf, ymax = Inf
)

# plot
p_ts_cf_year <- 
    ggplot(dat_eg_spp) +
    facet_grid(phenology ~ year_lab,
               switch = "y") +
    geom_rect(data = p_ts_bg, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = fill)) +
    geom_line(aes(month, value, 
                  group = Taxa), 
              alpha = 0.2) +
    scale_x_continuous(expand = c(0, 0), 
                       breaks = 1:12) +
    scale_fill_identity() +
    labs(x = "Month", 
         y = expression(paste("Phenology, ", eta))) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.y = element_text(size = 9),
          strip.text.x = element_text(size = 11, hjust = 0),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9),
          axis.text.x = element_text(colour = "black", size = 8),
          axis.ticks.y = element_blank(),
          panel.spacing.y = unit(0, "in"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())



# Fitted values over the study period -------------------------------------

dat_new <- 
    dat_brms %>% 
    distinct(date, rain_wet, Tmin, Taxa) %>% 
    # add average species (treat it as community mean)
    bind_rows(
        dat_brms %>% 
            distinct(date, rain_wet, Tmin) %>% 
            mutate(Taxa = NA)
    ) %>% 
    mutate(year = year(date),
           month = month(date),
           uid = NA,
           transect = NA)
cs_new <- make_c_s_mat(time_steps = dat_new$month, P = period)
cs_new <- do.call(cbind, lapply(cs_new, \(x) x[, 1:n_cs]))
colnames(cs_new) <- paste0(rep(c("C", "S"), each = n_cs), 
                           rep(1:n_cs, n_cs))
dat_new <- bind_cols(dat_new, cs_new)
ts_pred <- fitted(mod, 
                  dpar = "mu",
                  scale = "linear",
                  newdata = dat_new,
                  allow_new_levels = TRUE,
                  summary = TRUE,
                  robust = TRUE,
                  ndraws = 1000)
dat_new <- 
    bind_cols(dat_new, ts_pred[, "Estimate", ]) %>% 
    pivot_longer(cols = c(newleaf, noleaf, flowers, maturefruit),
                 names_to = "phenology",
                 values_to = "value") %>% 
    mutate(phenology = fct_relevel(phenology,
                                   "noleaf", "newleaf", "flowers", "maturefruit"),
           phenology = fct_recode(phenology,
                                  `Leaf shedding` = "noleaf",
                                  `Leaf flush` = "newleaf", 
                                  `Flowering` = "flowers", 
                                  `Fruiting` = "maturefruit")) 

# Plot
p_ts <- 
    ggplot() +
    facet_wrap(~ phenology, 
               strip.position = "left",
               ncol = 1) +
    geom_vline(xintercept = 
                   ym(paste(min(dat_new$year):(max(dat_new$year)+1), 1, sep = "-")),
               colour = "lightgrey", 
               linewidth = 0.3) +
    geom_line(data = dat_new %>% filter(!is.na(Taxa)),
              aes(date, value, group = Taxa), 
              colour = "black",
              alpha = 0.1) +
    scale_x_date(expand = c(0, 0),
                 date_labels = "%Y",
                 breaks = ym(paste(c(2004, 2010, 2014, 2020), 7, sep = "-")),
                 labels = ym(paste(c(2004, 2010, 2014, 2020), 1, sep = "-"))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "in"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          plot.margin = margin(20, 1, 1, 10),
          plot.title.position = "plot")





# Climate trends ----------------------------------------------------------

slope_pval <- 
    data.frame(climate = c("paste(R[wet], ' (mm)')", "paste(T[min], ' (°C)')"),
               text = c(paste("Slope = -0.65\np = 0.19"),
                        paste("Slope = 0.017\np < 0.001")))

p_clim <- 
    climate_full %>% 
    select(year, rain_wet, Tmin) %>% 
    pivot_longer(cols = -year,
                 names_to = "climate",
                 values_to = "value") %>% 
    mutate(climate = fct_recode(climate, 
                                `paste(R[wet], ' (mm)')` = "rain_wet",
                                `paste(T[min], ' (°C)')` = "Tmin")) %>% 
    ggplot() +
    facet_wrap(~ climate, 
               strip.position = "left",
               scales = "free_y", 
               labeller = "label_parsed",
               ncol = 1) +
    annotate("rect", 
             xmin = 2004, xmax = 2020, ymin = -Inf, ymax = Inf,
             fill = "snow2") +
    geom_text(data = slope_pval,
              aes(min(climate_full$year), Inf, label = text),
              size = 3, hjust = 0, vjust = 1.1) +
    geom_smooth(aes(year, value),
                method = "lm",
                formula = y ~ x,
                fill = "grey",
                colour = "black",
                alpha = 1) +
    geom_point(aes(year, value),
               pch = 21, fill = "white") +
    scale_y_continuous(n.breaks = 4) +
    scale_x_continuous(expand = c(0, 0.5),
                       breaks = c(1976, 1990, 2004, 2014, 2020)) +
    scale_fill_identity() +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "in"),
          axis.text = element_text(colour = "black"),
          axis.title = element_blank(),
          legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          plot.margin = margin(20, 1, 1, 1),
          plot.title.position = "plot")






# Trends in periodic components -------------------------------------------

source("code/trends_component.R")
