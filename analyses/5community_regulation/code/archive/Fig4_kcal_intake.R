

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages

librarian::shelf(tidyverse, sf, raster, terra, janitor)

basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","5community_regulation","figures")

#read sofa output
mass_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 4, skip=1) %>% clean_names()
dietcomp_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 5, skip = 1) %>% clean_names()
kcal_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 6, skip=1) %>% clean_names()

################################################################################

kcal_sum <- kcal_class %>%
  filter(period <=2019)%>%
  rename(year = period) %>%
  mutate(period = ifelse(year < 2013, "Pre-SSW", "Post-SSW")) %>%
  group_by(period) %>%
  summarize(across(.cols = 2:25, .fns = mean, na.rm = TRUE)) %>%
  pivot_longer(cols = 2:25, names_to = "prey", values_to = "kcal") %>%
  filter(!prey %in% c("lobster", "fish", "fishegg", "small_kelp_invert")) %>%
  mutate(predator = "Sea otter",
         period = factor(period, levels = c("Pre-SSW","Post-SSW")),
         prey_kcal = paste(prey,round(kcal,2))) %>%
  filter(!(kcal<0.01))

hist(kcal_sum$kcal)

################################################################################
#plot

# Theme
my_theme <-  theme(axis.text=element_text(size=8,color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=10,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=8,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=8, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())



# Calculate percent change and format prey
kcal_sum_perc <- kcal_sum %>%
  dplyr::select(period, prey, kcal) %>%
  pivot_wider(names_from = period, values_from = kcal) %>%
  filter(!(prey %in% c('sand_crab'))) %>%
  mutate(
    prey = str_replace_all(prey, "_", " "),  # Remove underscores
    prey = tolower(prey),                    # Convert to lowercase
    prey = str_to_sentence(prey),             # Convert to sentence case
    perc_change = ((`Post-SSW` - `Pre-SSW`) / `Pre-SSW`) * 100,
    label = paste0("(", round(`Pre-SSW`, 2), ", ", round(`Post-SSW`, 2), ")"),
    #assign functional group
    functional_group = case_when(
      prey == "Mussel" ~ "Planktivore",
      prey == "Urchin" ~ "Grazer",
      prey == "Snail" ~ "Grazer",
      prey == "Abalone" ~ "Grazer",
      prey == "Clam" ~ "Planktivore",
      prey == "Octopus" ~ "Macroinvertivore",
      prey == "Squid" ~ "Macroinvertivore",
      prey == "Kelp crab" ~ "Grazer",
      prey == "Worm" ~ "Detritivore",
      prey == "Cancer crab" ~ "Macroinvertivore",
      prey == "Sand dollar" ~ "Planktivore",
      prey == "Crab other" ~ "Macroinvertivore",
      prey == "Chiton" ~ "Grazer",
      prey == "Star" ~ "Macroinvertivore"
    )
  )


# Create a dumbbell plot
g1 <- ggplot(kcal_sum_perc, aes(x = perc_change, y = reorder(prey, perc_change)
                                           )) +
  geom_point(aes(color = functional_group, fill = functional_group)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = prey,color = functional_group), linetype = "solid", size = 1) +
  geom_text(aes( x = ifelse(kcal_sum_perc$perc_change >= 0, -50, 50),
                label = label), size=2.5) +  
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Percent change", y = "") +
  scale_x_continuous(limits=c(-100,350))+
  theme_bw() + my_theme



ggsave(g1, filename = file.path(figdir, "Fig4_diet_change.png"), 
       width = 5, height = 6, units = "in", dpi = 600)






################################################################################
#old toy code

# Create the plot
g <- ggplot(kcal_sum, aes(x = predator, xend = prey, y = 0, yend = kcal, color = prey)) +
  geom_segment(arrow = arrow(type = "closed", length = unit(0.01, "npc")))+
  geom_text(aes(x = prey, y = kcal, label = prey_kcal), vjust = 1, size=3) +
  labs(x = "", y = "") +
  facet_wrap(~period) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_y_reverse() +
  theme_bw() +
  guides(color = FALSE)+  # Hide the color legend
  my_theme 

g


#ggsave(g, filename = file.path(figdir, "Fig4_kcal_intake.png"), 
 #      width = 7, height = 3.5, units = "in", dpi = 600, bg = "white")


library(ggplot2)


# Create the circular rosette plot
g <- ggplot(kcal_sum, aes(x = factor(1), y = kcal, fill = prey)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = prey_kcal), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar(theta = "y") +
  labs(x = "", y = "") +
  facet_wrap(~period) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_fill_manual(values = rainbow(length(unique(kcal_sum$prey)))) +
  guides(fill = FALSE) +
  my_theme

g

# Calculate the width for each slice based on kcal
kcal_sum <- kcal_sum %>%
  group_by(period) %>%
  mutate(width = kcal / sum(kcal, na.rm = TRUE))

# Create the circular rosette plot
g <- ggplot(kcal_sum, aes(x = 1, y = width, fill = prey)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = prey_kcal, x = 1.2), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar(theta = "y") +
  labs(x = "", y = "") +
  facet_wrap(~period) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_fill_manual(values = rainbow(length(unique(kcal_sum$prey)))) +
  guides(fill = FALSE) +
  my_theme

g





# Function to create pie slice
pie.slice <- function(center, radius, arc1, arc2, color, label = NULL, label.dist = 1) {
  stopifnot(length(center) == 2)
  stopifnot(length(radius) <= 2)
  if (length(radius) == 1) radius <- c(radius, radius)
  
  xb <- seq(arc1, arc2, length.out = (arc2 * 180 - arc1 * 180) * 1)  # fine = 1
  
  x <- center[1] + c(0, cospi(2 * xb), 0) * radius[1]
  y <- center[2] + c(0, sinpi(2 * xb), 0) * radius[2]
  
  polygon(x, y, col = color)
  
  mid <- (arc1 + arc2) / 2
  anchor <- c(center[1] + cospi(2 * mid) * radius[1], center[2] + sinpi(2 * mid) * radius[2])
  
  if (!is.null(label))
    text(anchor[1], anchor[2], label, adj = c(-0.5, 1.0) * label.dist, srt = mid * 360)
  
  anchor
}

# Create the circular rosette plot
pdf("pieslice.pdf", width = 12, height = 12, pointsize = 12 / 0.75)
plot(0, xlim = c(-1, 2), ylim = c(-1, 2), xlab = "", ylab = "", type = "n")

radius <- 0.04
sl.out <- seq(0.0, 1.0, 0.1)

my.center <- c(-0.25, 1)
for (i in 2:(length(sl.out) - 1)) {
  prey_data <- filter(kcal_sum, period == "Pre-SSW")
  kcal_value <- prey_data$kcal[i]
  color <- rainbow(length(sl.out) - 2)[i - 1]
  pie.slice(my.center, radius * i, sl.out[i - 1], sl.out[i], color, label = kcal_value, label.dist = 1)
}

my.center <- c(1.0, 0.0)
for (i in 2:length(sl.out)) {
  prey_data <- filter(kcal_sum, period == "Post-SSW")
  kcal_value <- prey_data$kcal[i - 1]
  color <- rainbow(length(sl.out) - 1)[i - 1]
  pie.slice(my.center, 1.5 * radius * (length(sl.out) - i + 1), sl.out[i - 1], sl.out[i],
            color, label = kcal_value, label.dist = 0.25)
}

dev.off()



###########



# Function to create pie slice
pie.slice <- function(center, radius, arc1, arc2, color, label = NULL, label.dist = 1) {
  stopifnot(length(center) == 2)
  stopifnot(length(radius) <= 2)
  if (length(radius) == 1) radius <- c(radius, radius)
  
  xb <- seq(arc1, arc2, length.out = (arc2 * 180 - arc1 * 180) * 1)  # fine = 1
  
  x <- center[1] + c(0, cospi(2 * xb), 0) * radius[1]
  y <- center[2] + c(0, sinpi(2 * xb), 0) * radius[2]
  
  polygon(x, y, col = color)
  
  mid <- (arc1 + arc2) / 2
  anchor <- c(center[1] + cospi(2 * mid) * radius[1], center[2] + sinpi(2 * mid) * radius[2])
  
  if (!is.null(label))
    text(anchor[1], anchor[2], label, adj = c(-0.5, 1.0) * label.dist, srt = mid * 360)
  
  anchor
}

# Create the circular rosette plot
plot(0, xlim = c(-1, 2), ylim = c(-1, 2), xlab = "", ylab = "", type = "n")

radius <- 0.04
sl.out <- seq(0.0, 1.0, 0.1)

my.center <- c(-0.25, 1)
for (i in 2:(length(sl.out) - 1)) {
  prey_data <- filter(kcal_sum, period == "Pre-SSW")
  kcal_value <- prey_data$kcal[i]
  color <- rainbow(length(sl.out) - 2)[i - 1]
  prey_name <- prey_data$prey[i - 1]
  pie.slice(my.center, radius * i, sl.out[i - 1], sl.out[i], color, label = prey_name, label.dist = 1)
}

my.center <- c(1.0, 0.0)
for (i in 2:length(sl.out)) {
  prey_data <- filter(kcal_sum, period == "Post-SSW")
  kcal_value <- prey_data$kcal[i - 1]
  color <- rainbow(length(sl.out) - 1)[i - 1]
  prey_name <- prey_data$prey[i - 1]
  pie.slice(my.center, 1.5 * radius * (length(sl.out) - i + 1), sl.out[i - 1], sl.out[i],
            color, label = prey_name, label.dist = 0.25)
}





#combined rosette



# Function to create pie slice
pie.slice <- function(center, radius, arc1, arc2, color, label = NULL, label.dist = 1) {
  stopifnot(length(center) == 2)
  stopifnot(length(radius) <= 2)
  if (length(radius) == 1) radius <- c(radius, radius)
  
  xb <- seq(arc1, arc2, length.out = (arc2 * 180 - arc1 * 180) * 1)  # fine = 1
  
  x <- center[1] + c(0, cospi(2 * xb), 0) * radius[1]
  y <- center[2] + c(0, sinpi(2 * xb), 0) * radius[2]
  
  polygon(x, y, col = color)
  
  mid <- (arc1 + arc2) / 2
  anchor <- c(center[1] + cospi(2 * mid) * radius[1], center[2] + sinpi(2 * mid) * radius[2])
  
  if (!is.null(label))
    text(anchor[1], anchor[2], label, adj = c(-0.5, 1.0) * label.dist, srt = mid * 360)
  
  anchor
}

# Create the circular rosette plot
plot(0, xlim = c(-1, 2), ylim = c(-1, 2), xlab = "", ylab = "", type = "n")

radius <- 0.04
sl.out <- seq(0.0, 1.0, 0.1)

my.center <- c(0, 0)  # Center point for the combined rosette

prey_data_pre <- filter(kcal_sum, period == "Pre-SSW")
prey_data_post <- filter(kcal_sum, period == "Post-SSW")

for (i in 2:(length(sl.out) - 1)) {
  # Pre-SSW prey data
  kcal_value_pre <- prey_data_pre$kcal[i]
  color_pre <- rainbow(length(sl.out) - 2)[i - 1]
  prey_name_pre <- prey_data_pre$prey[i - 1]
  
  # Post-SSW prey data
  kcal_value_post <- prey_data_post$kcal[i - 1]
  color_post <- rainbow(length(sl.out) - 1)[i - 1]
  prey_name_post <- prey_data_post$prey[i - 1]
  
  # Calculate the radius for Pre-SSW and Post-SSW leaves
  radius_pre <- radius * i
  radius_post <- 1.5 * radius * (length(sl.out) - i + 1)
  
  # Calculate the start and end angles for Pre-SSW and Post-SSW leaves
  start_angle_pre <- sl.out[i - 1]
  end_angle_pre <- sl.out[i]
  
  start_angle_post <- sl.out[i - 1] + 0.05  # Offset Post-SSW angle to avoid overlap
  end_angle_post <- sl.out[i] + 0.05  # Offset Post-SSW angle to avoid overlap
  
  # Create Pre-SSW leaf
  pie.slice(my.center, radius_pre, start_angle_pre, end_angle_pre, color_pre, label = prey_name_pre, label.dist = 1)
  
  # Create Post-SSW leaf
  pie.slice(my.center, radius_post, start_angle_post, end_angle_post, color_post, label = prey_name_post, label.dist = 0.25)
}















###cobined test 2 

# Function to create pie slice
# Function to create pie slice
pie.slice <- function(center, radius, arc1, arc2, color, label = NULL, label.dist = 1) {
  stopifnot(length(center) == 2)
  stopifnot(length(radius) <= 2)
  if (length(radius) == 1) radius <- c(radius, radius)
  
  xb <- seq(arc1, arc2, length.out = 100)
  
  x <- center[1] + c(0, cospi(2 * xb), 0) * radius[1]
  y <- center[2] + c(0, sinpi(2 * xb), 0) * radius[2]
  
  polygon(x, y, col = color)
  
  mid <- (arc1 + arc2) / 2
  anchor <- c(center[1] + cospi(2 * mid) * radius[1], center[2] + sinpi(2 * mid) * radius[2])
  
  if (!is.null(label))
    text(anchor[1], anchor[2], label, adj = c(-0.5, 1.0) * label.dist, srt = mid * 360)
  
  anchor
}

# Create the circular rosette plot
plot(0, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", type = "n")

radius <- 0.04

my.center <- c(0, 0)  # Center point for the combined rosette

# Filter Pre-SSW and Post-SSW data
prey_data_pre <- filter(kcal_sum, period == "Pre-SSW")
prey_data_post <- filter(kcal_sum, period == "Post-SSW")

num_prey <- nrow(prey_data_pre)

for (i in 1:num_prey) {
  # Pre-SSW prey data
  kcal_value_pre <- prey_data_pre$kcal[i]
  color_pre <- rainbow(num_prey)[i]
  prey_name_pre <- prey_data_pre$prey[i]
  
  # Post-SSW prey data
  kcal_value_post <- prey_data_post$kcal[i]
  color_post <- rainbow(num_prey)[i]
  prey_name_post <- prey_data_post$prey[i]
  
  # Calculate the radius for Pre-SSW and Post-SSW leaves
  radius_pre <- 1.2  # Fixed radius for Pre-SSW leaves
  radius_post <- 1.6  # Fixed radius for Post-SSW leaves
  
  # Calculate the angles for Pre-SSW and Post-SSW leaves
  angle_step <- 2 * pi / num_prey
  start_angle_pre <- (i - 1) * angle_step
  end_angle_pre <- i * angle_step
  
  start_angle_post <- start_angle_pre
  end_angle_post <- end_angle_pre
  
  # Create Pre-SSW leaf
  pie.slice(my.center, radius_pre, start_angle_pre, end_angle_pre, color_pre, label = prey_name_pre, label.dist = 1)
  
  # Create Post-SSW leaf
  pie.slice(my.center, radius_post, start_angle_post, end_angle_post, color_post, label = prey_name_post, label.dist = 0.25)
}
