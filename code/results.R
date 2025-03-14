library(tidybayes)
library(modelr)
library(ggpubr)
library(scales)
library(circular)





# Amplitude and phase shift -----------------------------------------------

source("code/dominant_amplitude.R")





# Midline, amplitude and phase change with climate ------------------------

source("code/midline_change.R")
source("code/phenology_change.R")

png("fig/phenology_change.png", width = 9, height = 12, units = "in", res = 300)
ggarrange(
    p_midline_noleaf,
    p_A_noleaf_change + ggtitle(""),
    p_d_noleaf_change + ggtitle(""),
    p_midline_newleaf,
    p_A_newleaf_change + ggtitle(""),
    p_d_newleaf_change + ggtitle(""),
    p_midline_flowers,
    p_A_flowers_change + ggtitle(""),
    p_d_flowers_change + ggtitle(""),
    p_midline_maturefruit,
    p_A_maturefruit_change + ggtitle(""),
    p_d_maturefruit_change + ggtitle(""),
    nrow = 4,
    ncol = 3
)
dev.off()





# Variance partitioning ---------------------------------------------------

source("code/varpart_submodels.R")

pdf("fig/varpart_submodels.pdf", width = 7, height = 3.5)
p_varpart
dev.off()





# Synchrony ---------------------------------------------------------------

source("code/synchrony.R")

png("fig/synchrony.png", width = 7.5, height = 4, units = "in", res = 300)
ggarrange(p_sync_rain, p_sync_temp,
          nrow = 2)
dev.off()




# Temporal trends ---------------------------------------------------------

source("code/trends.R")

p_clim_ts <- 
    ggarrange(p_clim, p_ts, 
              ncol = 2,
              widths = c(0.8, 1),
              labels = c("(a) Climate", "(b) Phenology"),
              font.label = list(size = 12),
              hjust = 0)
p_ts_component <- 
    ggarrange(p_ts_midline, p_ts_amplitude, p_ts_phase,
              nrow = 1,
              labels = c("(c) Midline (mean intensity)",
                         "(d) Amplitude",
                         "(e) Phase (timing)"),
              font.label = list(size = 12),
              hjust = 0)

png("fig/trends.png", width = 7.5, height = 8, units = "in", res = 300)
ggarrange(p_clim_ts, p_ts_component,
          heights = c(0.8, 1),
          ncol = 1)
dev.off()
