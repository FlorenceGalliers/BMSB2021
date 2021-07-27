## HEADER ####
## Florence Galliers
## Direct Comparison of Climate Variables
## Last Edited: 2021-07-22

## Install Packages ####
library(dplyr)
library(agricolae)
library(raster)
library(tidyverse)
library(car)
install.packages("multcompView")
library(multcompView)


## Table of Contents ####
## 1.0 Import Current Climate Data
## 2.0 Crop Climatic Data into Extents
## 3.0 Import Occurrences
## 4.0 Crop Occurrences into extents
## 5.0 Extract values 
## 6.0 Combine results
## 7.0 Box Plot
## 8.0 ANOVA 
##     BIO 02
##     BIO 03
##     BIO 08
##     BIO 09
##     BIO 13
##     BIO 14
##     BIO 15
##     BIO 18
##     BIO 19

## 1.0 Import Current Climate Data ####
## Select only variables used in modeling
world_preds <- stack("data/world_preds_red_2.5.grd")

## 2.0 Crop Climatic Data into Extents ####
## USA Region
usa_extent <- extent(c(-128.2961, -55.3049, 20.79696, 55.83368))
usa_preds <- raster::crop(world_preds, usa_extent)
# Re-stack so it changes to raster stack rather than raster brick
usa_preds <- stack(usa_preds)
## Europe Region
europe_extent <- extent(c(-12, 45, 34, 72))
europe_preds <- raster::crop(world_preds, europe_extent)
europe_preds <- stack(europe_preds)
## Native Region
native_extent <- extent(c(92.54915, 149.632, 1.643559, 47.89202))
nat_preds <- raster::crop(world_preds, native_extent)
nat_preds <- stack(nat_preds)

## 3.0 Import Occurrences ####
occur <- read.csv("data/01_world_occur.csv")
## Assign coordinates and CRS to occurrence data
coordinates(occur) <- ~lon+lat
projection(occur) <- CRS('+proj=longlat +datum=WGS84')

## 4.0 Crop Occurrences into extents ####
usa_occur <- raster::crop(occur, usa_extent)
europe_occur <- raster::crop(occur, europe_extent)
nat_occur <- raster::crop(occur, native_extent)

## 5.0 Extract values ####
## Extract values of each bioclimatic variable for each occurrence point
# Native
nat_values <- raster::extract(nat_preds, nat_occur, 
                              method = "bilinear",
                              df = TRUE)
# USA
usa_values <- raster::extract(usa_preds, usa_occur, 
                              method = "bilinear",
                              df = TRUE)
# Europe
europe_values <- raster::extract(europe_preds, europe_occur, 
                              method = "bilinear",
                              df = TRUE)

## 6.0 Combine results ####

# Make new column giving the region of the values and remove the ID column
# For native values
nat_values$region <- "native"
nat_values <- nat_values[,-1]
# For usa values
usa_values$region <- "usa"
usa_values <- usa_values[,-1]
# For europe values
europe_values$region <- "europe"
europe_values <- europe_values[,-1]
# Bind together the three regions data into one data frame
mid <- rbind(nat_values, usa_values) # first bind two together
all_values <- rbind(mid, europe_values) # then all three
# Check all looks okay
head(all_values)
# Make region a factor
all_values$region <- factor(all_values$region)
levels(all_values$region)
# Put factor levels in order - native, usa, europe
all_values$region <- ordered(all_values$region,
                         levels = c("native", "usa", "europe"))

## 7.0 Box Plot ####
# Pivot longer to get into tidy format to make box plot
all_val <- pivot_longer(all_values,
                        cols = 1:9,
                        names_to = "BIOCLIMATIC",
                        values_to = "value")
head(all_val)

# Make region a factor with levels in correct order
all_val$region <- factor(all_val$region, levels = c("native", "usa", "europe"))
# Make levels have capital letters
levels(all_val$region) <- c("Native", "USA", "Europe")

# Make bioclimatic labels a factor
all_val$BIOCLIMATIC <- factor(all_val$BIOCLIMATIC)
# Rename factor levels to full names
levels(all_val$BIOCLIMATIC) <- c("02: Mean Diurnal Range ",
                                 "03: Isothermality",
                                 "08: Mean Temp. of Wettest Quarter",
                                 "09: Mean Temp. of Driest Quarter",
                                 "13: Precip. of Wettest Month",
                                 "14: Precip. of Driest Month",
                                 "15: Precip. Seasonality",
                                 "18: Precip. of Warmest Quarter", 
                                 "19: Precip. of Coldest Quarter")

# Make the plot
p2 <- ggplot(all_val, 
             aes(x = region, 
                 y = value, 
                 fill = region)) + # fill by region
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0) +
  scale_fill_manual(values = c("grey40", "grey90", "grey65")) +
  # facet for each variable
  facet_wrap(~BIOCLIMATIC, 
             scale = "free") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.spacing = unit(1, "cm"),
    axis.title = element_blank(),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12)
    )
# view plot
p2
# Save plot
ggsave("output/fig2.png", last_plot(),
       width = 30, height = 30, units = "cm", dpi = 600)
  

## 8.0 ANOVA ####

# Make separate data frames for each of the bioclimatic variables
bio_02 <- all_values[, c(1, 10)]
bio_03 <- all_values[, c(2, 10)]
bio_08 <- all_values[, c(3, 10)]
bio_09 <- all_values[, c(4, 10)]
bio_13 <- all_values[, c(5, 10)]
bio_14 <- all_values[, c(6, 10)]
bio_15 <- all_values[, c(7, 10)]
bio_18 <- all_values[, c(8, 10)]
bio_19 <- all_values[, c(9, 10)]

## BIO 02 ####

bio02_stats <- group_by(bio_02, region) %>%
  summarise(
    mean = mean(bio02, na.rm = TRUE),
    SE = sd(bio02, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio02_model <- lm(bio02 ~ region, data = bio_02)
# Summary of the analysis
anova(bio02_model)
bio02_aov <- aov(bio02_model)
TukeyHSD(bio02_aov, ordered = TRUE)
HSD.test(bio02_aov, "region", console = TRUE)

# test for normality 
plot(bio02_model, which = 2)
# test for equal variance
plot(bio02_model, which = 3)

ggplot(data = bio02_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO02") +
  coord_flip() +
  theme_bw()

## BIO 03 ####
bio03_stats <- group_by(bio_03, region) %>%
  summarise(
    mean = mean(bio03, na.rm = TRUE),
    SE = sd(bio03, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio03_model <- lm(bio03 ~ region, data = bio_03)
# Summary of the analysis
anova(bio03_model)
bio03_aov <- aov(bio03_model)
TukeyHSD(bio03_aov, ordered = TRUE)
HSD.test(bio03_aov, "region", console = TRUE)

# test for normality 
plot(bio03_model, which = 2)
# test for equal variance
plot(bio03_model, which = 3)

ggplot(data = bio03_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO03") +
  coord_flip() +
  theme_bw()


## BIO 08 ####
bio08_stats <- group_by(bio_08, region) %>%
  summarise(
    mean = mean(bio08, na.rm = TRUE),
    SE = sd(bio08, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio08_model <- lm(bio08 ~ region, data = bio_08)
# Summary of the analysis
anova(bio08_model)
bio08_aov <- aov(bio08_model)
TukeyHSD(bio08_aov, ordered = TRUE)
HSD.test(bio08_aov, "region", console = TRUE)

# test for normality 
plot(bio08_model, which = 2)
# test for equal variance
plot(bio08_model, which = 3)

ggplot(data = bio08_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO08") +
  coord_flip() +
  theme_bw()


## BIO 09 ####
bio09_stats <- group_by(bio_09, region) %>%
  summarise(
    mean = mean(bio09, na.rm = TRUE),
    SE = sd(bio09, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio09_model <- lm(bio09 ~ region, data = bio_09)
# Summary of the analysis
anova(bio09_model)
bio09_aov <- aov(bio09_model)
TukeyHSD(bio09_aov, ordered = TRUE)
HSD.test(bio09_aov, "region", console = TRUE)

# test for normality 
plot(bio09_model, which = 2)
# test for equal variance
plot(bio09_model, which = 3)

ggplot(data = bio09_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO09") +
  coord_flip() +
  theme_bw()


## BIO 13 ####
bio13_stats <- group_by(bio_13, region) %>%
  summarise(
    mean = mean(bio13, na.rm = TRUE),
    SE = sd(bio13, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio13_model <- lm(bio13 ~ region, data = bio_13)
# Summary of the analysis
anova(bio13_model)
bio13_aov <- aov(bio13_model)
TukeyHSD(bio13_aov, ordered = TRUE)
HSD.test(bio13_aov, "region", console = TRUE)

# test for normality 
plot(bio13_model, which = 2)
# test for equal variance
plot(bio13_model, which = 3)

ggplot(data = bio13_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO13") +
  coord_flip() +
  theme_bw()

## BIO 14 ####
bio14_stats <- group_by(bio_14, region) %>%
  summarise(
    mean = mean(bio14, na.rm = TRUE),
    SE = sd(bio14, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio14_model <- lm(bio14 ~ region, data = bio_14)
# Summary of the analysis
anova(bio14_model)
bio14_aov <- aov(bio14_model)
TukeyHSD(bio14_aov, ordered = TRUE)
HSD.test(bio14_aov, "region", console = TRUE)

# test for normality 
plot(bio14_model, which = 2)
# test for equal variance
plot(bio14_model, which = 3)

ggplot(data = bio14_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO14") +
  coord_flip() +
  theme_bw()


## BIO 15 ####
bio15_stats <- group_by(bio_15, region) %>%
  summarise(
    mean = mean(bio15, na.rm = TRUE),
    SE = sd(bio15, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio15_model <- lm(bio15 ~ region, data = bio_15)
# Summary of the analysis
anova(bio15_model)
bio15_aov <- aov(bio15_model)
TukeyHSD(bio15_aov, ordered = TRUE)
HSD.test(bio15_aov, "region", console = TRUE)

# test for normality 
plot(bio15_model, which = 2)
# test for equal variance
plot(bio15_model, which = 3)

ggplot(data = bio15_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO15") +
  coord_flip() +
  theme_bw()

## BIO 18 ####
bio18_stats <- group_by(bio_18, region) %>%
  summarise(
    mean = mean(bio18, na.rm = TRUE),
    SE = sd(bio18, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio18_model <- lm(bio18 ~ region, data = bio_18)
# Summary of the analysis
anova(bio18_model)
bio18_aov <- aov(bio18_model)
TukeyHSD(bio18_aov, ordered = TRUE)
HSD.test(bio18_aov, "region", console = TRUE)

# test for normality 
plot(bio18_model, which = 2)
# test for equal variance
plot(bio18_model, which = 3)

ggplot(data = bio18_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO18") +
  coord_flip() +
  theme_bw()

## BIO 19 ####
bio19_stats <- group_by(bio_19, region) %>%
  summarise(
    mean = mean(bio19, na.rm = TRUE),
    SE = sd(bio19, na.rm = TRUE)/sqrt(n())
  )

# Compute the analysis of variance
bio19_model <- lm(bio19 ~ region, data = bio_19)
# Summary of the analysis
anova(bio19_model)
bio19_aov <- aov(bio19_model)
TukeyHSD(bio19_aov, ordered = TRUE)
HSD.test(bio19_aov, "region", console = TRUE)

# test for normality 
plot(bio19_model, which = 2)
# test for equal variance
plot(bio19_model, which = 3)

ggplot(data = bio19_stats, 
       aes(x = region, y = mean, ymin = mean - SE, ymax = mean + SE)) + 
  geom_point(colour = "blue", size = 1) + 
  geom_errorbar(width = 0.1, colour = "blue") + 
  xlab("Region") + ylab("BIO19") +
  coord_flip() +
  theme_bw()



