#Joshua G. Smith 
#jossmith@mbayaq.org


# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(rerddap, lubridate, tidyverse, beepr, data.table, purrr)

# Directories
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
figdir <- here::here("analyses","4patch_drivers","Figures")

#read SST
sst_dat <- readRDS("2003_2022_daily_SST.Rds")

#determine monthly baseline mean 2003-2012
################################################################################

# Group your data by year_month and calculate the mean analysed_sst for each group
monthly_baseline_sst <- sst_dat %>% 
  filter(year < 2013)%>%
  group_by(month) %>% 
  summarize(monthly_baseline_sst = mean(analysed_sst),
            sd_sst = sd(analysed_sst),
            one_sd_above = monthly_baseline_sst + sd_sst,
            one_sd_below = monthly_baseline_sst - sd_sst,
            two_sd_above = monthly_baseline_sst +  2*sd_sst,
            two_sd_below = monthly_baseline_sst -  2*sd_sst)


ggplot(monthly_baseline_sst, aes(x = month, y = monthly_baseline_sst)) +
  geom_point() +
  geom_line()+
  geom_line(aes(y = one_sd_above),color="orange", linetype="dashed")+
  geom_line(aes(y = one_sd_below),color="orange", linetype="dashed")+
  geom_line(aes(y = two_sd_below), color = "indianred", linetype="dashed")+
  geom_line(aes(y = two_sd_above), color = "indianred", linetype="dashed")+
  labs(x = "Month", y = "Baseline SST (2003-2012)") +
  scale_x_discrete(limits = month.name) +
  theme_classic()


