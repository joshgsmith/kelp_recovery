#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce, reshape2)


################################################################################
#set directories and load data
figdir <- here::here("analyses","4patch_drivers","Figures")
basedir <- here::here("analyses","4patch_drivers","Output")

#load multivariate data
load(file.path(basedir, "multivariate_data.Rdata"))

#load standardized data
stan_dat <- read.csv(file.path(basedir, "kelp_stan_CC.csv")) 

#load raw dat
fish_raw <- read.csv(file.path(basedir, "kelp_fish_counts_CC.csv")) 

upc_raw <- read.csv(file.path(basedir, "kelp_upc_cov_CC.csv")) 

swath_raw <- read.csv(file.path(basedir, "kelp_swath_counts_CC.csv")) %>%
  #remove kelps -- we will handle them as their own group
  dplyr::select(-macrocystis_pyrifera,
                -pterygophora_californica,
                -eisenia_arborea,
                -stephanocystis_osmundacea,
                -laminaria_setchellii,
                -nereocystis_luetkeana,
                -alaria_marginata,
                -costaria_costata,
                -laminaria_farlowii,
                -pleurophycus_gardneri)

kelp_raw <- read.csv(file.path(basedir, "kelp_swath_counts_CC.csv")) %>% 
  #extract kelps as their own group              
          dplyr::select(1:11, macrocystis_pyrifera,
                              pterygophora_californica,
                              eisenia_arborea,
                              stephanocystis_osmundacea,
                              laminaria_setchellii,
                              nereocystis_luetkeana,
                              alaria_marginata,
                              costaria_costata,
                              laminaria_farlowii,
                              pleurophycus_gardneri)


################################################################################
#process data

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0)))

fish_sum <- fish_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(10:114, mean, na.rm = TRUE))

swath_sum <- swath_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:57, mean, na.rm = TRUE))

upc_sum <- upc_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:64, mean, na.rm = TRUE))

kelp_sum <- kelp_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:18, mean, na.rm = TRUE))


################################################################################
#Step 2 - determine optimal centroid clustering

#indexing
species_names <- colnames(stan_ord_dat)
# Find the index of the species of interest
species_index <- match("strongylocentrotus_purpuratus", species_names)
species_scores <- scores(stan_ord, display = "sites")[,] 


scrs<- as.data.frame(scores(stan_ord, display="site"))
scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year) #to facilitate computing centroids, add group var
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year, data = scrs, FUN = mean) #computes centroids by MHW and region

# determine optimal number of clusters using elbow method
wss <- (nrow(cent)-1)*sum(apply(cent[,2:3],2,var))
for(i in 1:13) wss[i] <- sum(kmeans(cent[,2:3],centers=i)$withinss)
plot(1:13, wss, type="b", xlab="Number of clusters", ylab="Within groups sum of squares")

# Create a data frame for plotting
df <- data.frame(NumClusters = 1:13, WSS = wss)

# Plot using ggplot2
elbow_plot <- ggplot(df, aes(x = NumClusters, y = WSS)) +
  geom_line(size=0.5) +
  geom_point(size=1) +
  labs(x = "Number of clusters", y = "Within groups sum of squares", size=2) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     axis.text = element_text(size = 6),
                     axis.title = element_text(size=8))


#ggsave(elbow_plot, filename=file.path(figdir, "FigSX_elbow_plot.png"), bg = "white",
 #      width=4, height=4, units="in", dpi=600) 


# choose number of clusters based on elbow in plot
k <- 2

# perform k-means clustering and add cluster label to data frame
set.seed(123)
cent$cluster <- kmeans(cent[,2:3], k)$cluster

################################################################################
#Step 3 - process diversity
#process swath data

swath_dat <- swath_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
swath_groups <- swath_sum %>% ungroup() %>% dplyr::select(1:3)

swath_richness <- data.frame(S.obs = apply(swath_dat[,1:49]>0, 1, sum)) #count each unique species once
swath_evenness <- diversity(swath_dat)/log(specnumber(swath_dat))
swath_shannon <- diversity(swath_dat, index="shannon")
swath_simpson <- diversity(swath_dat, index="simpson")
swath_abund <- rowSums(swath_dat[,1:49])

swath_alphadiv <- cbind(swath_groups, swath_richness, swath_shannon, swath_simpson, swath_evenness, swath_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))


#process upc data
upc_dat <- upc_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
upc_groups <- upc_sum %>% ungroup() %>% dplyr::select(1:3)

total_area_covered <- rowSums(upc_dat)

abundance <- upc_dat / total_area_covered * 166.67 #convert %cov to abundance, multiply by scalar, in this case 60m2 transect, 166.67 (which is 10,000 divided by 60)

upc_dat[,1:56] <- abundance

#convert perc cov to relative abundance

upc_richness <- data.frame(S.obs = apply(upc_dat[,1:56]>0, 1, sum))
upc_evenness <- diversity(upc_dat)/log(specnumber(upc_dat))
upc_shannon <- diversity(upc_dat, index="shannon")
upc_simpson <- diversity(upc_dat, index="simpson")
upc_abund <- rowSums(upc_dat[,1:56])

upc_alphadiv <- cbind(upc_groups, upc_richness, upc_shannon, upc_simpson, upc_evenness, upc_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))


#process fish data

fish_dat <- fish_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
fish_groups <- fish_sum %>% ungroup() %>% dplyr::select(1:3)

fish_richness <- data.frame(S.obs = apply(fish_dat[,1:105]>0, 1, sum))
fish_evenness <- diversity(fish_dat)/log(specnumber(fish_dat))
fish_shannon <- diversity(fish_dat, index="shannon")
fish_simpson <- diversity(fish_dat, index="simpson")
fish_abund <- rowSums(fish_dat[,1:105])

fish_alphadiv <- cbind(fish_groups, fish_richness, fish_shannon, fish_simpson, fish_evenness, fish_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")),
         fish_evenness = ifelse(fish_evenness == "NaN",NA,fish_evenness))


#process kelps

kelp_dat <- kelp_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
kelp_groups <- kelp_sum %>% ungroup() %>% dplyr::select(1:3)

kelp_richness <- data.frame(S.obs = apply(kelp_dat[,1:10]>0, 1, sum)) #count each unique species once
kelp_evenness <- diversity(kelp_dat)/log(specnumber(kelp_dat))
kelp_shannon <- diversity(kelp_dat, index="shannon")
kelp_simpson <- diversity(kelp_dat, index="simpson")
kelp_abund <- rowSums(kelp_dat[,1:10])

kelp_alphadiv <- cbind(kelp_groups, kelp_richness, kelp_shannon, kelp_simpson, kelp_evenness, kelp_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))

################################################################################
#Step 4 determine spp that explain changes over time

fish_join <- cbind(fish_alphadiv, fish_dat) %>% rename(richness=S.obs) %>%
              #calculate annual mean
              group_by(year)%>%
              summarize(across(3:112,mean))

fish_rich <- fish_join %>%
  select(year, 7:111) %>%
  mutate(across(2:106, ~ifelse(. > 0, 1, 0))) 

# Create a data frame for species changes (added or lost)
species_changes <- data.frame(year = numeric(0), species_added = character(0), species_lost = character(0))

# Sort the data frame by year
fish_rich <- fish_rich %>%
  arrange(year)

# Iterate through years starting from the second year (2008)
for (i in 2:nrow(fish_rich)) {
  current_year <- fish_rich$year[i]
  previous_year <- fish_rich$year[i - 1]
  
  current_species <- fish_rich[i, 3:106]  # Exclude the year column from comparison
  previous_species <- fish_rich[i - 1, 3:106]
  
  species_added <- colnames(current_species)[current_species > previous_species]
  species_lost <- colnames(previous_species)[previous_species > current_species]
  
  species_changes <- rbind(species_changes, data.frame(year = current_year, species_added = paste(species_added, collapse = ", "), species_lost = paste(species_lost, collapse = ", ")))
}

# Print or work with species_changes to see which species were added or lost each year relative to the previous year
print(species_changes)


# Reshape the data to long format
fish_rich_long <- fish_rich %>%
  pivot_longer(cols = -year, names_to = "species", values_to = "presence")

# Filter out species that were never observed (presence == 0) in any year
fish_rich_filtered <- fish_rich_long %>%
  group_by(species) %>%
  filter(any(presence == 1))

# Calculate richness (number of species observed) for each year based on the filtered dataset
richness_data <- fish_rich_filtered %>%
  group_by(year) %>%
  summarize(richness = sum(presence))

# Filter out species that were observed in ALL years
fish_rich_reduced <- fish_rich_filtered %>%
  group_by(species) %>%
  filter(sum(presence) != n_distinct(year))

# Create a heatmap
g <- ggplot(data = fish_rich_reduced, aes(x = year, y = species, fill = factor(presence))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "white", "1" = "green"), labels = c("0" = "Absent", "1" = "Present")) +
  labs(title = "Species Presence/Absence Over Years", x = "Year", y = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Create a colored line plot with a gradient color transition between years
ggplot(data = richness_data, aes(x = year, y = 0, color = richness)) +
  geom_line(size = 20, aes(group = 1)) +
  labs(title = "Richness Over Years", x = "Year", y = NULL) +
  scale_color_gradientn(
    colors = c("navyblue", "indianred"),
    values = scales::rescale(c(min(richness_data$year), max(richness_data$year))),
    guide = "none"
  ) +  # Adjust the gradient colors and values as needed
  theme_minimal()

################################################################################
#Step 5 - examine site cohesion and distance changes over time

stan_group_vars2 <- stan_group_vars %>% mutate(period = ifelse(year < 2013,"Before","After"))

#use betadisper to reduce vegdist to principal coords
disper_mat <- betadisper(stan_max_distmat, type="centroid",
                         group = stan_group_vars2$period)

plot(disper_mat)





shift_dist <- reshape2::melt(as.matrix(sqrt(dist(disper_mat$centroids[,disper_mat$eig>0]^2)-
                                        dist(disper_mat$centroids[,disper_mat$eig<0]^2))))%>%
            tibble::rownames_to_column("distance")










# get betadisper dataframes ####
# have written functions to grab the necessary data from the betadisper object

# functions ####
# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

dist_explore <- betadisper_distances(disper_mat)

betadisper_dat <- get_betadisper_data(disper_mat)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(dplyr::select(betadisper_dat$centroids, group, PCoA1, PCoA2), 
                          dplyr::select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group')) %>%
  mutate(distance = sqrt((PCoA1.x - PCoA1.y)^2 + (PCoA2.x - PCoA2.y)^2))



# do some transformations on the data
#betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig))




# add convex hull points ####
# this could be put in a function
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()



# Now the dataframes are all ready to be completely customisable in ggplot
# plot betadispersion plot
ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 4) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$eigenvector) +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull ) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = group), betadisper_lines) +
  theme_bw(base_size = 12, base_family = 'Helvetica') 




# plot distances from centroid
ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.2)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  scale_color_manual('', values = c('black', 'grey'), labels = c("Grazed", 'Ungrazed')) +
  scale_fill_manual('', values = c('black', 'grey'), labels = c("Grazed", 'Ungrazed')) +
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('') +
  scale_x_discrete(labels = c('Grazed', 'Ungrazed'))







####################################
#betadisper mod



betadisper_mod <-
  function(d, group, type = c("median","centroid"), bias.adjust=FALSE,
           sqrt.dist = FALSE, add = FALSE)
  {
    ## inline function for double centring. We used .C("dblcen", ...,
    ## PACKAGE = "stats") which does not dublicate its argument, but
    ## it was removed from R in r60360 | ripley | 2012-08-22 07:59:00
    ## UTC (Wed, 22 Aug 2012) "more conversion to .Call, clean up".
    dblcen <- function(x, na.rm = TRUE) {
      cnt <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2L, cnt, check.margin = FALSE)
      cnt <- rowMeans(x, na.rm = na.rm)
      sweep(x, 1L, cnt, check.margin = FALSE)
    }
    ## inline function for spatial medians
    spatialMed <- function(vectors, group, pos) {
      axes <- seq_len(NCOL(vectors))
      spMedPos <- ordimedian(vectors, group, choices = axes[pos])
      spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
      cbind(spMedPos, spMedNeg)
    }
    ## inline function for centroids
    centroidFUN <- function(vec, group) {
      cent <- apply(vec, 2,
                    function(x, group) tapply(x, INDEX = group, FUN = mean),
                    group = group)
      if(!is.matrix(cent)) { ## if only 1 group, cent is vector
        cent <- matrix(cent, nrow = 1,
                       dimnames = list(as.character(levels(group)),
                                       paste0("Dim", seq_len(NCOL(vec)))))
      }
      cent
    }
    ## inline function for distance computation
    Resids <- function(x, c) {
      if(is.matrix(c))
        d <- x - c
      else
        d <- sweep(x, 2, c)
      rowSums(d^2)
    }
    ## Tolerance for zero Eigenvalues
    TOL <- sqrt(.Machine$double.eps)
    ## uses code from stats:::cmdscale by R Core Development Team
    if(!inherits(d, "dist"))
      stop("distances 'd' must be a 'dist' object")
    ## Someone really tried to analyse correlation like object in range -1..+1
    if (any(d < -TOL, na.rm = TRUE))
      stop("dissimilarities 'd' must be non-negative")
    ## adjust to avoid negative eigenvalues (if they disturb you)
    if (sqrt.dist)
      d <- sqrt(d)
    if (is.logical(add) && isTRUE(add))
      add <- "lingoes"
    if (is.character(add)) {
      add <- match.arg(add, c("lingoes", "cailliez"))
      if (add == "lingoes") {
        ac <- addLingoes(as.matrix(d))
        d <- sqrt(d^2 + 2 * ac)
      }
      else if (add == "cailliez") {
        ac <- addCailliez(as.matrix(d))
        d <- d + ac
      }
    }
    if(missing(type))
      type <- "median"
    type <- match.arg(type)
    ## checks for groups - need to be a factor for later
    group <- if(!is.factor(group)) {
      as.factor(group)
    } else { ## if already a factor, drop empty levels
      droplevels(group, exclude = NA) # need exclude = NA under Rdevel r71113
    }
    n <- attr(d, "Size")
    x <- matrix(0, ncol = n, nrow = n)
    x[row(x) > col(x)] <- d^2
    ## site labels
    labs <- attr(d, "Labels")
    ## remove NAs in group
    if(any(gr.na <- is.na(group))) {
      group <- group[!gr.na]
      x <- x[!gr.na, !gr.na]
      ## update n otherwise C call crashes
      n <- n - sum(gr.na)
      ## update labels
      labs <- labs[!gr.na]
      message("missing observations due to 'group' removed")
    }
    ## remove NA's in d
    if(any(x.na <- apply(x, 1, function(x) any(is.na(x))))) {
      x <- x[!x.na, !x.na]
      group <- group[!x.na]
      ## update n otherwise C call crashes
      n <- n - sum(x.na)
      ## update labels
      labs <- labs[!x.na]
      message("missing observations due to 'd' removed")
    }
    x <- x + t(x)
    x <- dblcen(x)
    e <- eigen(-x/2, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    ## Remove zero eigenvalues
    eig <- eig[(want <- abs(eig) > max(TOL, TOL * eig[1L]))]
    ## scale Eigenvectors
    vectors <- vectors[, want, drop = FALSE] %*% diag(sqrt(abs(eig)),
                                                      nrow = length(eig))
    ## store which are the positive eigenvalues
    pos <- eig > 0
    ## group centroids in PCoA space
    centroids <-
      switch(type,
             centroid = centroidFUN(vectors, group),
             median = spatialMed(vectors, group, pos)
      )
    ## for each of the groups, calculate distance to centroid for
    ## observation in the group
    ## Uses in-line Resids function as we want LAD residuals for
    ## median method, and LSQ residuals for centroid method
    dist.pos <- Resids(vectors[, pos, drop=FALSE],
                       centroids[group, pos, drop=FALSE])
    dist.neg <- 0
    if(any(!pos))
      dist.neg <- Resids(vectors[, !pos, drop=FALSE],
                         centroids[group, !pos, drop=FALSE])
    
    
    # Calculate distances to centroids of all groups
    dist_to_other_centroids <- matrix(0, ncol = n, nrow = n)
    for (i in 1:length(levels(group))) {
      if (type == "centroid" && i == 1) {
        # If calculating centroids, only need to do it once
        dist_to_other_centroids <- Resids(vectors, centroids[group == levels(group)[i],, drop = FALSE])
      } else {
        # Calculate distances to centroids of other groups
        other_group_indices <- group != levels(group)[i]
        dist_to_other_centroids[other_group_indices] <- Resids(vectors[other_group_indices, , drop = FALSE], centroids[group == levels(group)[i],, drop = FALSE])
      }
    }
    
    # zij are the distances of each point to the other group centroids
    zij <- sqrt(dist_to_other_centroids)
    
    if (bias.adjust) {
      n.group <- as.vector(table(group))
      zij <- zij * sqrt(n.group[group] / (n.group[group] - 1))
    }
    
    
    ## zij are the distances of each point to its group centroid
    if (any(dist.neg > dist.pos)) {
      ## Negative squared distances give complex valued distances:
      ## take only the real part (which is zero). Github issue #306.
      warning("some squared distances are negative and changed to zero")
      zij <- Re(sqrt(as.complex(dist.pos - dist.neg)))
    } else {
      zij <- sqrt(dist.pos - dist.neg)
    }
    if (bias.adjust) {
      n.group <- as.vector(table(group))
      zij <- zij*sqrt(n.group[group]/(n.group[group]-1))
    }
    ## pre-compute group mean distance to centroid/median for `print` method
    grp.zij <- tapply(zij, group, "mean")
    ## add in correct labels
    if (any(want))
      colnames(vectors) <- names(eig) <-
      paste("PCoA", seq_along(eig), sep = "")
    if(is.matrix(centroids))
      colnames(centroids) <- names(eig)
    else
      names(centroids) <- names(eig)
    rownames(vectors) <- names(zij) <- labs
    retval <- list(eig = eig, vectors = vectors, distances = zij,
                   group = group, centroids = centroids,
                   group.distances = grp.zij, call = match.call())
    class(retval) <- "betadisper"
    attr(retval, "method") <- attr(d, "method")
    attr(retval, "type") <- type
    attr(retval, "bias.adjust") <- bias.adjust
    
    # Return the modified results
    retval$distances_other_centroids <- zij
    
    retval
  }




################################################################################
#Step 4 - plot


my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_text(size=8, color = "black"),
                   plot.title =element_text(size=7, face="bold", color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6, color = "black"),
                   legend.title = element_text(size = 7, color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold", color ="black"),
)


#plot NMDS
stan_trajectory <- ggplot(data = cent %>% mutate(basin = ifelse(year < 2012, "before",ifelse(year > 2014, "after",NA))),
                          aes(NMDS1, NMDS2)) +
  geom_segment(aes(x = NMDS1, y = NMDS2, xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               size = 0.5, color = "black",
               data = cent, size=4) +
  #geom_path(size = 1, lineend = "round") + #use geom_path to validate the segment
  stat_ellipse(aes(fill = factor(cluster)), type = "norm", level = 0.8, geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values = c("forestgreen", "purple")) +
  #labs(shape="Heatwave period",color='Site type')+
  ggrepel::geom_label_repel(aes(label = year),
                            box.padding   = 1, 
                            point.padding = 0.1,
                            segment.color = 'grey',
                            max.overlaps=Inf,
                            size=2) +
  ggtitle("Kelp forest community structure") +
  labs(color = 'K-means \ncluster', fill = "K-means \ncluster", tag = "A")+
  theme_bw() + my_theme

#stan_trajectory


ymax <- max(swath_alphadiv$swath_shannon, na.rm=TRUE)
shannon <- ggplot() +
  #kelp
  geom_point(aes(x = year, y = kelp_shannon, color = "Kelps", fill = "Kelps"), alpha = 0.1, size = 0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = kelp_shannon, color = "Kelps", fill = "Kelps"), data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y = swath_shannon, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), alpha = 0.1, size = 0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = swath_shannon, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #upc
  geom_point(aes(x = year + 0.1, y = upc_shannon, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), alpha = 0.1, size = 0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = upc_shannon, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year - 0.1, y = fish_shannon, color = "Fishes", fill = "Fishes"), alpha = 0.1, size = 0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = fish_shannon, color = "Fishes", fill = "Fishes"), data = fish_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa") +
  ylab("Shannon diversity") +
  xlab("Year") +
  ggtitle("Shannon diversity") +
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))


#shannon

ymax <- max(swath_alphadiv$swath_simpson, na.rm=TRUE)
simpson <- ggplot() +
  #kelp
  geom_point(aes(x = year, y = kelp_simpson, color = "Kelps", fill = "Kelps"), alpha = 0.1, size = 0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = kelp_simpson, color = "Kelps", fill = "Kelps"), data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y = swath_simpson, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), alpha = 0.1, size = 0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = swath_simpson, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #upc
  geom_point(aes(x = year + 0.1, y = upc_simpson, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), alpha = 0.1, size = 0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = upc_simpson, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year - 0.1, y = fish_simpson, color = "Fishes", fill = "Fishes"), alpha = 0.1, size = 0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = fish_simpson, color = "Fishes", fill = "Fishes"), data = fish_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa") +
  ylab("Simpson diversity") +
  xlab("Year") +
  ggtitle("Simpson diversity") +
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))
#simpson


richness <- ggplot()+
  #kelp
  geom_point(aes(x = year, y=S.obs, color = "Kelps", fill="Kelps"), alpha = 0.1, size=0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Kelps", fill="Kelps"),data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y=S.obs, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=S.obs, color = "Fishes", fill="Fishes"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Fishes", fill="Fishes"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=S.obs, color = "Sessile inverts and \nmacroalgae",fill="Sessile inverts and \nmacroalgae"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa")+
  ylab("Species richness (n)")+
  xlab("Year")+
  ggtitle("Species richness")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))

#richness



evenness <- ggplot()+
  #kelp
  geom_point(aes(x = year, y=kelp_evenness, color = "Kelps", fill="Kelps"), alpha = 0.1, size=0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=kelp_evenness, color = "Kelps", fill="Kelps"),data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y=swath_evenness, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=swath_evenness, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=fish_evenness, color = "Fishes", fill="Fishes"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=fish_evenness, color = "Fishes", fill="Fishes"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=upc_evenness, color = "Sessile inverts and \nmacroalgae",fill="Sessile inverts and \nmacroalgae"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=upc_evenness, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa")+
  ylab("Species evenness (n)")+
  xlab("Year")+
  ggtitle("Species evenness")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))

#evenness

# Combine ggplots with ggarrange()
combined_plot <- ggpubr::ggarrange(shannon, simpson, richness, evenness, ncol = 2, nrow=2, common.legend = TRUE, legend = "right") 

# Display the plot
combined_plot <- combined_plot + theme(plot.margin = margin(5, 0, 5, 20, "pt")) + 
  labs(tag = "B") + theme(plot.tag=element_text(size=8, color = "black", face = "plain")) 

#ggsave(combined_plot, filename=file.path(figdir, "Fig3_richness_diversity_evenness.png"), bg = "white",
 #      width=5.5, height=4.5, units="in", dpi=600) 


region_wide_plot <- ggpubr::ggarrange(stan_trajectory, combined_plot, ncol=1)
region_wide_plot 

ggsave(region_wide_plot, filename=file.path(figdir, "Fig3_regional_metrics_new6.png"), bg = "white",
      width=5.5, height=8, units="in", dpi=600) 





