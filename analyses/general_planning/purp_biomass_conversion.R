

rm(list = ls())

librarian::shelf(tidyverse)

#load mlpa size fq data
size_fq <- read.csv("/Users/jossmith/Desktop/size_year.csv") #aggregated mlpa data

################################################################################
#determine no. urchins per 20 meter swath
# average urchin density is 11 per m2 -- (reported in Smith et al. 2023 ProcB)

swath_den <- 11*2*10 #two by 10 meter swath clearing

#determine size frequency
 DatU.raw <- size_fq %>%
   filter(year == 2021)%>%
   dplyr::group_by(size) %>%
   dplyr::summarize(Freq = sum(count))

colnames(DatU.raw)[1] = "UrchSize"
FrequencyTable <- DatU.raw
DatU <- vcdExtra::expand.dft(FrequencyTable, freq="Freq")


# Calculate the relative proportions of each size class
size_class_proportions <- table(DatU$UrchSize) / length(DatU$UrchSize)

# Estimate the number of individuals for each size class in a sample size of 220

size_class_counts_in_sample <- round(size_class_proportions * swath_den)

#expand
size_class_long <- vcdExtra::expand.dft(size_class_counts_in_sample, freq="Freq")


# Convert size from centimeters to millimeters 
size_class_long$Var1_mm <- size_class_long$Var1 * 10

################################################################################
# Define the biomass conversion function
biomass_conversion <- function(diameter_mm) {
  biomass_g <- -22.45 + 12.23 * exp(0.0394 * diameter_mm)
  return(biomass_g)
}

# Apply the biomass conversion function to the size column in size_class_long
size_class_long$Biomass <- biomass_conversion(size_class_long$Var1_mm)

# Calculate the total biomass for all individuals in the sample
total_biomass_in_sample <- sum(size_class_long$Biomass) / 1000 # convert to kg


################################################################################
# plot size fq

plot <- ggplot(size_class_long, aes(x = Var1)) +
  geom_histogram(binwidth = 1, fill = "purple", color = "black") +
  labs(
    x = "Size (cm)",
    y = "Frequency per 20m² swath"
  ) +
  theme_bw()

plot


################################################################################
# plot biomass relationship

# Create a sequence of size 
size_values_cm <- seq(0, 10, by = 0.1)  

# Calculate the corresponding biomass values using the biomass conversion function
biomass_values <- biomass_conversion(size_values_cm)

# Combine size and biomass values into a dataframe
biomass_df <- data.frame(Size_cm = size_values_cm, Biomass = biomass_values)

# Create a scatterplot with a line plot overlay
plot <- ggplot(biomass_df, aes(x = Size_cm, y = Biomass)) +
  geom_point(color = "blue") +
  geom_line(color = "red") +
  labs(
    title = "Biomass Conversion Function",
    x = "Size (cm)",
    y = "Biomass (g)"
  ) +
  theme_minimal()

plot









biomass_g <- -22.45 + 12.23 * exp(0.0394 * 40)





