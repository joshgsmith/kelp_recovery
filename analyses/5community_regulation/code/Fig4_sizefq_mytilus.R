

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages

librarian::shelf(tidyverse, sf, raster, terra, janitor)

basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","5community_regulation","figures")

#read mussel size fq. 
mus_size_orig <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mussel_size_fq.xlsx"),sheet = 1)


################################################################################
#Step 1 - filter to focal study area
mus_build1 <- mus_size_orig %>% 
  filter(marine_site_name %in% c("Point Lobos","Hopkins","Stillwater","Point Pinos"))%>%
  rename(year = marine_common_year)%>%
  mutate(ssw_period = ifelse(year <=2013, "before","after")) 

DatU <- vcdExtra::expand.dft(mus_build1, freq="total")


# Calculate the weighted mean size bin for each year using row count as frequency
mean_size_by_year <- DatU %>%
  group_by(year) %>%
  summarize(mean_size = mean(size_bin))
            
            # Create a histogram-like plot with bars touching and vertical bars for the mean
     ggplot(DatU, aes(x = size_bin)) +
      geom_bar(binwidth = 2, fill = "gray40", color = "gray40") +
    geom_vline(data = mean_size_by_year, aes(xintercept = mean_size), linetype = "dotted", color = "red", size = 1) +
       facet_wrap(~ year, ncol = 1) +
          labs(x = "Size Bin", y = "Frequency") +
              ggtitle("Size Frequency Distribution by Year") +
              theme_bw()
            
