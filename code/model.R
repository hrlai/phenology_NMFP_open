
library(brms)


# Data --------------------------------------------------------------------

dat <- read.csv("data/data_open.csv")




# Model -------------------------------------------------------------------

mod <- brm(
    mvbind(new_leaf, no_leaf, flowers, mature_fruit) ~ 
        1 + (C1 + S1 + C2 + S2) * (rain_wet + Tmin) + 
        (1 + (C1 + S1 + C2 + S2) * (rain_wet + Tmin) || Taxa) +
        (1 + C1 + S1 + C2 + S2 || year) +
        (1 || transect) +
        (1 || uid),
    data = dat,
    family = cumulative(link = "logit"),
    control = list(adapt_delta = 0.95),
    cores = 4
)
