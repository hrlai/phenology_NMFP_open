A_delta_taxa <- 
    coef(mod,
         groups = "Taxa",
         summary = FALSE)

# calculate difference between the amplitudes of different fourier components
# to find the dominant component (large amplitude)

A_delta_newleaf <- 
    A_delta_taxa$Taxa[, , c("newleaf_C1", "newleaf_S1", "newleaf_C2", "newleaf_S2")] %>% 
    reshape2::melt() %>% 
    pivot_wider(names_from = Var3,
                values_from = value) %>% 
    mutate(A1 = sqrt(newleaf_C1^2 + newleaf_S1^2),
           A2 = sqrt(newleaf_C2^2 + newleaf_S2^2)) %>% 
    rename(Taxa = Var2) %>% 
    mutate(A_diff = A1 - A2) %>% 
    group_by(Taxa) %>% 
    median_qi(A_diff, .width = 0.89) %>% 
    left_join(spp_names)

A_delta_noleaf <- 
    A_delta_taxa$Taxa[, , c("noleaf_C1", "noleaf_S1", "noleaf_C2", "noleaf_S2")] %>% 
    reshape2::melt() %>% 
    pivot_wider(names_from = Var3,
                values_from = value) %>% 
    mutate(A1 = sqrt(noleaf_C1^2 + noleaf_S1^2),
           A2 = sqrt(noleaf_C2^2 + noleaf_S2^2)) %>% 
    rename(Taxa = Var2) %>% 
    mutate(A_diff = A1 - A2) %>% 
    group_by(Taxa) %>% 
    median_qi(A_diff, .width = 0.89) %>% 
    left_join(spp_names)

A_delta_flowers <- 
    A_delta_taxa$Taxa[, , c("flowers_C1", "flowers_S1", "flowers_C2", "flowers_S2")] %>% 
    reshape2::melt() %>% 
    pivot_wider(names_from = Var3,
                values_from = value) %>% 
    mutate(A1 = sqrt(flowers_C1^2 + flowers_S1^2),
           A2 = sqrt(flowers_C2^2 + flowers_S2^2)) %>% 
    rename(Taxa = Var2) %>% 
    mutate(A_diff = A1 - A2) %>% 
    group_by(Taxa) %>% 
    median_qi(A_diff, .width = 0.89) %>% 
    left_join(spp_names)

A_delta_maturefruit <- 
    A_delta_taxa$Taxa[, , c("maturefruit_C1", "maturefruit_S1", "maturefruit_C2", "maturefruit_S2")] %>% 
    reshape2::melt() %>% 
    pivot_wider(names_from = Var3,
                values_from = value) %>% 
    mutate(A1 = sqrt(maturefruit_C1^2 + maturefruit_S1^2),
           A2 = sqrt(maturefruit_C2^2 + maturefruit_S2^2)) %>% 
    rename(Taxa = Var2) %>% 
    mutate(A_diff = A1 - A2) %>% 
    group_by(Taxa) %>% 
    median_qi(A_diff, .width = 0.89) %>% 
    left_join(spp_names)

# determine species' dominant fourier component by phenology
A_dominant_newleaf <- 
    A_delta_newleaf %>% 
    mutate(Dominant = ifelse(A_diff > 0, 1, 2)) %>% 
    select(Taxa, Name, Dominant) %>% 
    bind_rows(
        data.frame(Taxa = "Community", Name = "Community", Dominant = 1)
    )
A_dominant_noleaf <- 
    A_delta_noleaf %>% 
    mutate(Dominant = ifelse(A_diff > 0, 1, 2)) %>% 
    select(Taxa, Name, Dominant) %>% 
    bind_rows(
        data.frame(Taxa = "Community", Name = "Community", Dominant = 1)
    )
A_dominant_flowers <- 
    A_delta_flowers %>% 
    mutate(Dominant = ifelse(A_diff > 0, 1, 2)) %>% 
    select(Taxa, Name, Dominant) %>% 
    bind_rows(
        data.frame(Taxa = "Community", Name = "Community", Dominant = 1)
    )
A_dominant_maturefruit <- 
    A_delta_maturefruit %>% 
    mutate(Dominant = ifelse(A_diff > 0, 1, 2)) %>% 
    select(Taxa, Name, Dominant) %>% 
    bind_rows(
        data.frame(Taxa = "Community", Name = "Community", Dominant = 1)
    )




# Plot --------------------------------------------------------------------
p_A_dominant_newleaf <- 
    A_delta_newleaf %>% 
    mutate(Taxa = fct_reorder(Name, A_diff)) %>% 
    ggplot() + 
    geom_vline(xintercept = 0, colour = "grey") +
    geom_errorbarh(aes(y = Taxa, xmin = .lower, xmax = .upper),
                   height = 0) +
    geom_point(aes(A_diff, Taxa, fill = A_diff < 0), 
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face = 3, size = 5),
          plot.margin = margin(20, 1, 1, 1),
          legend.position = "none")
p_A_dominant_noleaf <- 
    A_delta_noleaf %>% 
    mutate(Taxa = fct_reorder(Name, A_diff)) %>% 
    ggplot() + 
    geom_vline(xintercept = 0, colour = "grey") +
    geom_errorbarh(aes(y = Taxa, xmin = .lower, xmax = .upper),
                   height = 0) +
    geom_point(aes(A_diff, Taxa, fill = A_diff < 0), 
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face = 3, size = 5),
          plot.margin = margin(20, 1, 1, 1),
          legend.position = "none")
p_A_dominant_flowers <- 
    A_delta_flowers %>% 
    mutate(Taxa = fct_reorder(Name, A_diff)) %>% 
    ggplot() + 
    geom_vline(xintercept = 0, colour = "grey") +
    geom_errorbarh(aes(y = Taxa, xmin = .lower, xmax = .upper),
                   height = 0) +
    geom_point(aes(A_diff, Taxa, fill = A_diff < 0), 
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face = 3, size = 5),
          plot.margin = margin(20, 1, 1, 1),
          legend.position = "none")
p_A_dominant_maturefruit <- 
    A_delta_maturefruit %>% 
    mutate(Taxa = fct_reorder(Name, A_diff)) %>% 
    ggplot() + 
    geom_vline(xintercept = 0, colour = "grey") +
    geom_errorbarh(aes(y = Taxa, xmin = .lower, xmax = .upper),
                   height = 0) +
    geom_point(aes(A_diff, Taxa, fill = A_diff < 0), 
               pch = 21) +
    scale_fill_grey(start = 0.5, end = 1) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(face = 3, size = 5),
          plot.margin = margin(20, 1, 1, 1),
          legend.position = "none")
