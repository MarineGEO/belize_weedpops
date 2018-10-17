############################################################################
###                   BELIZE SQUIDPOP ANALYSIS                           ###
############################################################################

# Author: Jon Lefcheck
# Last updated: 17 Oct 2018
# Contact: LefcheckJ@si.edu

#---------------------------------------------------------------------------

# Load required libraries
library(cowplot)
library(effects)
library(ggplot2)
library(lme4)
library(lmtest)
library(lubridate)
library(piecewiseSEM)
library(readxl)
library(tidyverse)

setwd("C:/Users/lefcheckj/OneDrive - Smithsonian Institution/Documents/GitHub/belize_weedpops")

#---------------------------------------------------------------------------

# Read in data
wp <- read_excel("RITTER UPDATE 09_2018_CBC_weedpop-vid_analysis.xlsx", sheet = 2)

vidlengths <- read.csv("weedpopsVidLengths.csv")

# Merge weedpop and vidlengths data
wp <- left_join(wp, vidlengths)

# Code habitats
wp$habitat <- recode(wp$locality,
                     "Curlew Patch Reef" = "Patch Reef",
                     "House Reef" = "Reef",
                     "South Reef" = "Reef",
                     "Tobacco Reef" = "Reef",
                     "CBC Central 30" = "Reef",
                     "CBC Patch Reef" = "Patch Reef")

# Code sampling location
wp$location <- recode(wp$locality,
                      "Curlew Patch Reef" = "Curlew",
                      "House Reef" = "House",
                      "South Reef" = "South",
                      "Tobacco Reef" = "Tobacco",
                      "CBC Central 30" = "Carrie Bow",
                      "CBC Patch Reef" = "Carrie Bow")
                      
# # Convert measurementValue to seconds
# td1 <- as.numeric(seconds_to_period(na.omit(as.numeric(wp$measurementValue))))
# 
# # Convert BaitEatenTime to seconds
# td2 <- na.omit(as.numeric(ms(wp$BaitEatenTime)) - as.numeric(ms(wp$clockstart)))
# 
# wp$timestamp <- c(td1, td2)

# Code species
wp$species <- paste(wp$genus, wp$specificEpithet)

wp$species <- as.factor(wp$species)

wp$algalType <- as.factor(wp$algalType)

# Assign bites
wp$bite <- 1

# Fix typos
wp$algalType <- recode(wp$algalType, "acanthophora" = "Acanthophora")

wp$specificEpithet <- recode(wp$specificEpithet, "coerulus" = "coeruleus")

#---------------------------------------------------------------------------

# Analyze richness effects

# Summarize data by replicate
wp.summary <- wp %>% group_by(habitat, location, videoNumber, algalType) %>%
  
  summarize(
    richness = length(unique(species)),
    bites = sum(bite)
  )

# Model relationship between richness and bites based on algal type
richness.model <- glmer(bites ~ richness * algalType + (1|location/videoNumber) + (1|habitat), family = "poisson", data = wp.summary)

summary(richness.model)

# Get standardized coefficients
# coefs(model)

# Check assumptions
hist(resid(richness.model))

plot(richness.model)

# Get model R^2
rsquared(richness.model)

# Get partial effects
richness.model.effects <- allEffects(richness.model)

# Extract predictor scores for plotting
richness.predict.df <- cbind.data.frame(
  richness.model.effects$`richness:algalType`$x,
  fit = richness.model.effects$`richness:algalType`$fit,
  lower = richness.model.effects$`richness:algalType`$lower,
  upper = richness.model.effects$`richness:algalType`$upper
)

# Plot predicted scores
ggplot(richness.predict.df, aes(x = richness, y = fit, group = algalType, color = algalType, fill = algalType)) +
  # geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.2, col = NA) +
  geom_line(lwd = 1) +
  scale_fill_manual(values = c("firebrick3", "coral1", "chartreuse3", "burlywood4", "forestgreen"), name = "Type") +
  scale_color_manual(values = c("firebrick3", "coral1", "chartreuse3", "burlywood4", "forestgreen"), name = "Type") +
  labs(x = "Richness", y = "E(Bites)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#---------------------------------------------------------------------------

# Fit species-interaction model

# Summarize by species
wp.species.summary <- wp %>% group_by(habitat, location, videoNumber, algalType, species) %>%
  
  summarize(bites = sum(bite)) %>%
  
  spread(species, bites, fill = 0)

# Subset species with >= bites
remove <- names(which(colSums(wp.species.summary[, 5:ncol(wp.species.summary)]) < 5))

wp.species.summary <- wp.species.summary[, !colnames(wp.species.summary) %in% remove]

# Summarize total bites
wp.species.summary$total.bites <- rowSums(wp.species.summary[, 5:ncol(wp.species.summary)])

# Convert species to presence/absence
wp.species.summary[, 5:(ncol(wp.species.summary) - 1)][wp.species.summary[, 5:(ncol(wp.species.summary) - 1)] > 0] <- 1

# Fit model for each algalType
sig.species.df <- do.call(rbind, lapply(unique(wp.species.summary$algalType), function(i) {
  
  # Subset data
  dat <- subset(wp.species.summary, algalType == i)
  
  # Remove species with no bites
  remove <- names(which(colSums(dat[, 5:(ncol(dat) - 1)]) == 0))
  
  dat <- dat[, !colnames(dat) %in% remove]
  
  # Create model formula
  form <- paste("total.bites ~ ", paste0("`", colnames(dat[, 5:(ncol(dat) - 1)]), "`", collapse = " + "))
  
  # Fit model
  species.model <- glmer(formula(paste(form, "+ (1|location) + (1|habitat)")), family = "poisson", data = dat)

  # Extract significant species (P < 0.05)
  sig.species <- names(which(summary(species.model)$coefficients[, 4] < 0.05))
  
  sig.species <- sig.species[sig.species != "(Intercept)"]
  
  sig.species <- gsub("`", "", sig.species)
  
  data.frame(
    algalType = i,
    sig.species = sig.species
  )
  
} ) ) 

# Plot results
ggplot(sig.species.df, aes(x = sig.species, y = forcats::fct_rev(algalType), fill = algalType)) +
  geom_tile(col = "white") +
  scale_fill_manual(values = c("firebrick3", "coral1", "chartreuse3", "burlywood4", "forestgreen"), name = "Type") +
  # scale_fill_manual(values = c("forestgreen", "burlywood4", "chartreuse3", "coral1", "firebrick3"), name = "Type") +
  labs(x = "", y = "Algal type") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#---------------------------------------------------------------------------

# Generate bite accumulation curves

#Summarize data
wp.bites.cum <- wp %>%
  
  cbind.data.frame(., bites = 1) %>% 
  
  group_by(location, habitat) %>%
  
  arrange(measurementValue) %>%
  
  mutate(total.bites = cumsum(bites))

# Plot bite accumulation curves
ggplot(wp.bites.cum, aes(x = measurementValue/60, y = total.bites, group = location, col = habitat)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = c("goldenrod", "dodgerblue", "dodgerblue4", "forestgreen"), guide = F) +
  facet_grid( ~ habitat) + #, scales = "free_x") +
  labs(x = "Time (min)", y = "Cumulative bites") +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

# Plot time series of bites and species' identities

# Assign colors by genus
colors <- c(
  paste0("darkgoldenrod", 1:4),
  "purple",
  "cadetblue",
  "deeppink2",
  "aquamarine",
  "coral1",
  paste0("dodgerblue", 2:3),
  "darkslateblue",
  paste0("darkseagreen", 1:3),
  paste0("firebrick", 1:4),
  "darkolivegreen3")

names(colors) <- unique(wp$species)[order(unique(wp$species))]

# Plot results
# out.sig <- subset(out, sig == "Y")

plotlist <- lapply(unique(wp$habitat), function(i) {
  
  dat <- subset(wp, habitat == i)
  
  legend.plot <- ggplot(dat, aes(x = measurementValue, y = bite, fill = species)) +
    geom_tile() +
    scale_fill_manual(values = colors, name = "") 
  
  legend <- get_legend(legend.plot)
  
  plist <- lapply(unique(dat$location), function(j) {
    
    dat1 <- subset(dat, location == j)
    
    # dat1$sig <- ifelse(dat1$species %in% unlist(subset(out.sig, location == j & habitat == i)[, c("species1", "species2")]), "Y", "N")
    
    # dat1 <- subset(dat1, sig == "Y")
    
    dat1$new.value <- 1:nrow(dat1)
    
    p <- ggplot(dat1, aes(x = measurementValue/60, #new.value, 
                          y = bite - 0.5, group = species, fill = species)) + #, alpha = sig)) +
      geom_tile() + 
      # scale_alpha_manual(values = c(0, 1), guide = F) +
      scale_fill_manual(values = colors[names(colors) %in% dat$species], guide = F) +
      geom_hline(yintercept = 0, col = "white", lwd = 1.1) +
      labs(x = "Time", y = "", title = paste(j)) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    return(p)
    
  } )
  
  pcol <- cowplot::plot_grid(plotlist = plist, ncol = 1) 
  
  finalp <- plot_grid(pcol, legend, ncol = 2, rel_widths = c(1, 0.3))
  
  title <- ggdraw() + draw_label(paste(i), fontface = "bold", size = 20)
  
  plot_grid(title, finalp, ncol = 1, rel_heights = c(0.1, 1))
  
} )

pdf("Weedpops Patch Reef species time series.PDF", width = 11, height = 9)
bquiet <- print(plotlist[[1]])
dev.off()    

pdf("Weedpops Reef species time series.PDF", width = 10, height = 12)
bquiet <- print(plotlist[[2]])
dev.off()   

#---------------------------------------------------------------------------

# Univariate Tests of Granger Causality

#' cgrangertest: a function to apply Granger Causality to communities
#' 
#' @param species a vector of species names to compare
#' @param data a time-by-species data.frame
#' 
#' @return a data.frame of comparisons and associated P-values
#' 
cgrangertest <- function(species, data) {
  
  # Compute comparisons
  combos <- combn(species, 2)
  
  # Loop over combos
  res <- do.call(rbind, lapply(1:ncol(combos), function(i) {
    
    print(i)
    
    sp <- combos[, i]
  
    dat <- data[, sp]
  
    model <- try(grangertest(dat), silent = TRUE)
    
    if(all(class(model) == "try-error")) data.frame() else 
      
      data.frame(species1 = sp[1], species2 = sp[2], P.value = round(model$`Pr(>F)`[2], 4))
    
  } ) )
  
  # Adjust alpha
  alpha <- 0.05 / ncol(combos)
  
  # Add significance
  res$sig <- ifelse(res$P.value < alpha, "Y", "N")
  
  return(res)
  
}

# Apply to weedpop data
wp.long <- wp %>% group_by(location, habitat, algalType, species) %>%
  
  complete(measurementValue = min(measurementValue):max(measurementValue)) %>%
  
  spread(species, bite, fill = 0) # %>%
  
  # select(location, habitat, algalType, measurementValue, `Acanthurus bahianus`:`Thalassoma bifasciatum`)

species <- names(wp.long)[32:51] # [4:23]

out <- do.call(rbind, lapply(unique(wp.long$habitat), function(i) {
    
    dat <- subset(wp.long, habitat == i)
    
    do.call(rbind, lapply(unique(dat$location), function(j) {
      
      dat1 <- subset(dat, location == j)
      
      # lapply(unique(dat1$algalType), function(k) {
      #   
      #   dat2 <- subset(dat1, algalType == k)
        
      dat2 <- dat1
      
        # Remove columns with zeros
        data <- dat2[, 32:51][, !colSums(dat2[, 32:51]) == 0]
        
        print(paste(i, j))
        
        data.frame(location = j, habitat = i, cgrangertest(species = colnames(data), data)) # alga = k, 
        
      # } ) )
      
    } ) )

} ) )

lapply(unique(out$habitat), function(i) {
  
  dat <- subset(out, habitat == i)
  
  plist <- lapply(unique(dat$location), function(j) {
    
    dat2 <- subset(dat, location == j)
    
    dat2$species1 <- factor(dat2$species1, levels = levels(dat2$species1)[order(levels(dat2$species1))])
  
    dat2$species2 <- factor(dat2$species2, levels = levels(dat2$species2)[rev(order(levels(dat2$species2)))])
    
    # Plot results
    ggplot(dat2, aes(x = species1, y = species2, fill = sig)) +
      geom_raster() +
      scale_fill_manual(values = c("red", "grey80"), guide = F) +
      theme_bw(base_size = 12) +
      labs(x = "", y = "", title = j) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
  } )

  pfinal <- plot_grid(plotlist = plist, ncol = 1)
  
  title <- ggdraw() + draw_label(paste(i), fontface = "bold", size = 20)
  
  plot_grid(title, pfinal, ncol = 1, rel_heights = c(0.1, 1))
  
} )

