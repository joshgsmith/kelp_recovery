
library(dplyr)

set.seed(123)  # set random seed for reproducibility

n <- 1200  # number of urchins
diam_mean <- 40  # mean diameter
diam_min <- 25  # minimum diameter
diam_max <- 60  # maximum diameter

# generate random diameters
diam <- rnorm(n, mean = diam_mean, sd = 4)
# adjust diameters to be within desired range
diam[diam < diam_min] <- diam_min
diam[diam > diam_max] <- diam_max

# create data frame
df <- data.frame(test_diameter = diam) %>%
      mutate(biomass_g = -22.45 + 12.23*exp(0.0394*diam),
             biomass_kg = biomass_g /1000 )


total_required_pounds <- sum(df$biomass_kg)*2.20462


plot(df$biomass_g, df$test_diameter)




