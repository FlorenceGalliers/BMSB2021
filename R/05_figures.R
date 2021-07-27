## HEADER ####
## Figures for Report
## Florence Galliers
## Last Edited: 2021-07-22

# Install Packages
library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyverse)

## Table of Contents ####
## 1.0 Fig.1 = Occurrences Plot
## 2.0 Fig.3 = Five Models
## 3.0 Fig.6 = Future Climate Plots
## 4.0 Fig.5 = Habitat-Climatic Plots
## 5.0 Table 3 = Evaluation Metrics
## 6.0 Fig.4 = Individual Algorithm Metrics

## crs for plots
land_crs <- crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80
                 +units=m +no_defs")


## 1.0 Fig.1 = Occurrences Plot ####
## Import Occurrences
points <- read.csv("data/01_world_occur.csv")
points <- points[, -c(1, 2, 5)]

## Get world map data from rnaturalearth package
worldmap <- ne_countries(scale = "medium", returnclass = "sf")

# Plot
occurrenceplot <- ggplot() + geom_sf(data = worldmap, 
                                     fill = "white", # background of countries colour
                                     size = 0.3, # linewidth
                                     colour = "grey40") +  # line colour
  coord_sf(xlim = c(-180, 180), ylim = c(-65, 90.68569), expand = FALSE) +
  geom_point(data = points, aes(x = lon, y = lat), size = 0.5, 
             shape = 16, col = "black") +
  # add USA extent rectangle
  geom_rect(aes(xmin=-128.2961, xmax=-55.3049, ymin=20.79696, ymax=55.83368), 
            inherit.aes = FALSE, 
            color = "blue", fill = "transparent") +
  # add europe extent rectangle
  geom_rect(aes(xmin=-12, xmax=45, ymin=34, ymax=72), 
            inherit.aes = FALSE, 
            color = "red", fill = "transparent") +
  # add native extent rectangle
  geom_rect(aes(xmin=92,54915, xmax=149.632, ymin=1.643559, ymax=47.89202), 
            inherit.aes = FALSE, 
            color = "seagreen", fill = "transparent") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8),
    plot.title = element_text(vjust = -7.5, hjust = 0.005, size = 11)
  ) +
  annotate("text", x = -118, y = 61, 
           label = "(a) USA",
           color = "blue") +
  annotate("text", x = 3, y = 78, 
           label = "(b) Europe",
           color = "red") +
  annotate("text", x = 105, y = 54, 
           label = "(c) Native",
           color = "seagreen") +
  ggtitle("(d) Global")

ggsave("output/fig1.png", plot = occurrenceplot,
       width = 21, height = 10, unit = "cm", dpi = 600)

## 2.0 Fig.3 = Five Models ####
# 2.0 Import Projected Models
# Native
native_proj <- stack("BMSB/proj_nat_onto_uk/proj_nat_onto_uk_BMSB_ensemble.grd")
n <- native_proj[[2]]
plot(n)
summary(n)
# Europe
europe_proj <- stack("BMSB/proj_europe_onto_uk/proj_europe_onto_uk_BMSB_ensemble.grd")
eu <- europe_proj[[2]]
plot(eu)
summary(eu)
# USA
usa_proj <- stack("BMSB/proj_usa_onto_uk/proj_usa_onto_uk_BMSB_ensemble.grd")
u <- usa_proj[[2]]
plot(u)
summary(u)
# Worldwide
world_proj <- stack("BMSB/proj_world_onto_uk/proj_world_onto_uk_BMSB_ensemble.grd")
w <- world_proj[[2]]
plot(w)
summary(w)
# Invaded
invaded_proj <- stack("BMSB/proj_invaded_onto_uk/proj_invaded_onto_uk_BMSB_ensemble.grd")
inv <- invaded_proj[[2]]
plot(inv)
summary(inv)
## Reproject projections into different CRS 
n <- projectRaster(from = n, crs = land_crs, method = "bilinear")
eu <- projectRaster(from = eu, crs = land_crs, method = "bilinear")
u <- projectRaster(from = u, crs = land_crs, method = "bilinear")
w <- projectRaster(from = w, crs = land_crs, method = "bilinear")
inv <- projectRaster(from = inv, crs = land_crs, method = "bilinear")

#fixInNamespace("draw.colorkey", "lattice")
my.at <- seq(0, 1000, by = 50)
labs <- seq(0, 1000, by = 100)

## Plot 

png(file="output/fig3.png", width=21, height=21, units = "cm", res=600)

fourmodels <- rasterVis::levelplot(
  stack(n, u, eu, w, inv), 
  at = my.at,
  col.regions = rev(terrain.colors(20)), 
  colorkey = list(space = "bottom",
                  labels = list(labels = seq(0, 1, by = 0.1), 
                                at = labs, cex = 0.8),
                  width = 1, height = 1),
  names.attr = c("(a) Native", "(b) USA", "(c) Europe", "(d) Global", "(e) Invaded"),
  layout = c(2, 3),
  scales = list(tck = c(0,0),
                draw = FALSE),
  xlab = list("Climatic Suitability", cex = 0.8, vjust = 7),
  ylab = NULL,
  par.settings = list(axis.line = list(col = "transparent"),
                      strip.background = list(col = "transparent"),
                      strip.border = list(col = "transparent"),
                      layout.heights = list(bottom.padding = 3)),
) 

print(fourmodels)

dev.off()


## 3.0 Fig.6 = Future Climate Plots ####
# Import Projections
future_crs <- crs(native_proj)

## Invaded Model
en_26_2060_inv <- stack("BMSB/proj_inv_uk_26_2060/proj_inv_uk_26_2060_BMSB_ensemble.grd")
en_26_2060_inv <- en_26_2060_inv[[2]]
crs(en_26_2060_inv) <- future_crs
en_26_2060_inv <- projectRaster(from = en_26_2060_inv, crs = land_crs, method = "bilinear")

en_26_2080_inv <- stack("BMSB/proj_inv_uk_26_2080/proj_inv_uk_26_2080_BMSB_ensemble.grd")
en_26_2080_inv <- en_26_2080_inv[[2]]
crs(en_26_2080_inv) <- future_crs
en_26_2080_inv <- projectRaster(from = en_26_2080_inv, crs = land_crs, method = "bilinear")

en_85_2060_inv <- stack("BMSB/proj_inv_uk_85_2060/proj_inv_uk_85_2060_BMSB_ensemble.grd")
en_85_2060_inv <- en_85_2060_inv[[2]]
crs(en_85_2060_inv) <- future_crs
en_85_2060_inv <- projectRaster(from = en_85_2060_inv, crs = land_crs, method = "bilinear")

en_85_2080_inv <- stack("BMSB/proj_inv_uk_85_2080/proj_inv_uk_85_2080_BMSB_ensemble.grd")
en_85_2080_inv <- en_85_2080_inv[[2]]
crs(en_85_2080_inv) <- future_crs
en_85_2080_inv <- projectRaster(from = en_85_2080_inv, crs = land_crs, method = "bilinear")

summary(en_26_2060_inv)
summary(en_26_2080_inv)
summary(en_85_2060_inv)
summary(en_85_2080_inv)


my.at <- seq(0, 1000, by = 50)
labs <- seq(0, 1000, by = 100)

png(file="output/fig6.png", width=21, height=21, units = "cm", res=600)

future <- rasterVis::levelplot(stack(en_26_2060_inv, en_26_2080_inv, 
                           en_85_2060_inv, en_85_2080_inv), 
                     at = my.at,
                     col.regions = rev(terrain.colors(20)), 
                     colorkey = list(space = "bottom",
                                     labels = list(labels = seq(0, 1, by = 0.1),
                                                   at = labs, cex = 0.8),
                                     width = 1, height = 1),
                     names.attr = c("RCP2.6, 2041-2060", 
                                    "RCP2.6, 2061-2080",
                                    "RCP8.5, 2041-2060", 
                                    "RCP8.5, 2061-2080"),
                     layout = c(2,2),
                     scales = list(tck = c(0,0),
                                   draw = FALSE),
                     xlab = list("Climatic Suitability", cex = 0.8, vjust = 7),
                     ylab = NULL,
                     par.settings = list(axis.line = list(col = "transparent"),
                                         strip.background = list(col = "transparent"),
                                         strip.border = list(col = "transparent"),
                                         layout.heights = list(bottom.padding = 3)),
)

print(future)

dev.off()


## 4.0 Fig.5 = Habitat-Climatic Plots ####

uk_clim <- raster("output/03_uk_clim.grd")
uk_hab <- raster("output/03_habitat_only_uk.grd")
hab_clim_uk <- raster("output/03_hab_clim_uk.grd")
summary(uk_clim)
summary(uk_hab)
summary(hab_clim_uk)


at_at <- seq(0, 1, by = 0.05)
hab_labs <- seq(0, 1, by = 0.1)

png(file="output/fig5.png", width=21, height=21, units = "cm", res=600)

habclim <- rasterVis::levelplot(stack(uk_clim, uk_hab, hab_clim_uk), 
                     at = at_at,
                     col.regions = rev(terrain.colors(20)), 
                     colorkey = list(space = "bottom",
                                     labels = list(labels = seq(0, 1, by = 0.1),
                                                   at = hab_labs, cex = 0.8),
                                     width = 1, height = 1),
                     names.attr = c("(a) Climatic Suitability", 
                                    "(b) Habitat Suitability",
                                    "(c) Combined Suitability"),
                     layout = c(2, 2),
                     scales = list(tck = c(0,0),
                                   draw = FALSE),
                     xlab = list("Suitability", cex = 0.8, vjust = 7),
                     ylab = NULL,
                     par.settings = list(axis.line = list(col = "transparent"),
                                         strip.background = list(col = "transparent"),
                                         strip.border = list(col = "transparent"),
                                         layout.heights = list(bottom.padding = 3)),
)

print(habclim)

dev.off()

## 5.0 Table 3 = Evaluation Metrics ####

## Native Model
native <- load("BMSB/BMSB.native_modelensemble.models.out")
native <- get(native)
native_metrics <- get_evaluations(native)
native_metrics <- native_metrics$BMSB_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData
## USA Model
usa <- load("BMSB/BMSB.usaensemble.models.out")
usa <- get(usa)
usa_metrics <- get_evaluations(usa)
usa_metrics <- usa_metrics$BMSB_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData
## Europe Model
europe <- load("BMSB/BMSB.europe_modelensemble.models.out")
europe <- get(europe)
europe_metrics <- get_evaluations(europe)
europe_metrics <- europe_metrics$BMSB_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData
## World Model
world <- load("BMSB/BMSB.world_modelensemble.models.out")
world <- get(world)
world_metrics <- get_evaluations(world)
world_metrics <- world_metrics$BMSB_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData

TSS_metrics <- c(native_metrics[1], usa_metrics[1], europe_metrics[1], world_metrics[1])
ROC_metrics <- c(native_metrics[5], usa_metrics[5], europe_metrics[5], world_metrics[5])


## 6.0 Fig.4 = Individual Algorithm Metrics ####

## 6.1 Algo Metrics Invaded model ####
world2 <- load("BMSB/BMSB.invaded_model.models.out")
world2 <- get(world2)
EM_algo <- BIOMOD_EnsembleModeling(modeling.output = world2,
                                           chosen.models = 'all', # all models
                                           em.by = 'algo', # a total consensus model
                                           eval.metric = c('TSS'), # use TSS for evaluation
                                           eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                           models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                           prob.mean.weight = TRUE,
                                           prob.mean.weight.decay = 'proportional')
world_algo_metrics <- get_evaluations(EM_algo) 
world_algo_metrics

world_algos <- c("GLM", "GBM", "GAM", "CTA", "ANN",
                 "SRE", "FDA", "MARS", "RF", "MAXENT")
world_algo_TSS <- c("0.836", "0.891", "0.897", "0.935", "0.886",
                    "0.698", "0.841", "0.859", "0.974", "0.882")
world_algo_ROC <- c("0.966", "0.982", "0.981", "0.982", "0.978",
                    "0.849", "0.969", "0.972", "1.000", "0.982")

world_TSS <- as.data.frame(cbind(world_algos, world_algo_TSS))
world_TSS$world_algo_TSS <- as.numeric(world_TSS$world_algo_TSS)

world_ROC <- as.data.frame(cbind(world_algos, world_algo_ROC))
world_ROC$world_algo_ROC <- as.numeric(world_ROC$world_algo_ROC)

world_TSS$metric <- "TSS"
world_ROC$metric <- "ROC"

world_TSS <- world_TSS %>%
  rename(
    score = world_algo_TSS,
    algo = world_algos
  )

world_ROC <- world_ROC %>%
  rename(
    score = world_algo_ROC,
    algo = world_algos
  )

all_world_metrics <- rbind(world_TSS, world_ROC)

# create dummy data frame for 0.7 threshold  of TSS
dummy2 <- data.frame(X = c("ROC", "TSS"),
                     Z = c(0.9, 0.7))

# Change the names of facet labels for plot
metric.labs <- c("(a) ROC", "(b) TSS")
names(metric.labs) <- c("ROC", "TSS")

# Plot
ggplot(all_world_metrics) +
  geom_bar(aes(score, reorder(algo, score)), 
           stat = "identity",
           fill = "white",
           color = "black") +
  xlab("Accuracy Score") +
  ylab("Algorithm") +
  facet_grid(~metric,
             labeller = labeller(metric = metric.labs)) +
  geom_vline(data = subset(all_world_metrics, metric == "TSS"), aes(xintercept = 0.71),
             color = "red", lty = 2, lwd = 0.8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
    )

ggsave("metrics-plot.png", plot = last_plot(), 
       width = 15, height = 10, units = "cm")

## 6.2 Algo Metrics Native Model ####
native <- load("BMSB/BMSB.native_model.models.out")
native <- get(native)
native_EM_algo <- BIOMOD_EnsembleModeling(modeling.output = native,
                                   chosen.models = 'all', # all models
                                   em.by = 'algo', # a total consensus model
                                   eval.metric = c('TSS'), # use TSS for evaluation
                                   eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                   models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                   prob.mean.weight = TRUE,
                                   prob.mean.weight.decay = 'proportional')
native_algo_metrics <- get_evaluations(native_EM_algo) 
native_algo_metrics

native_algos <- c("GLM", "GBM", "GAM", "CTA", "ANN",
                 "SRE", "FDA", "MARS", "RF", "MAXENT")
native_algo_TSS <- c("0.814", "0.881", "0.878", "0.901", "0.841",
                    "0.626", "0.770", "0.847", "0.993", "0.851")
native_algo_ROC <- c("0.948", "0.981", "0.978", "0.977", "0.962",
                    "0.813", "0.945", "0.968", "1.000", "0.974")

native_TSS <- as.data.frame(cbind(native_algos, native_algo_TSS))
native_TSS$native_algo_TSS <- as.numeric(native_TSS$native_algo_TSS)

native_ROC <- as.data.frame(cbind(native_algos, native_algo_ROC))
native_ROC$native_algo_ROC <- as.numeric(native_ROC$native_algo_ROC)

native_TSS$metric <- "TSS"
native_ROC$metric <- "ROC"

native_TSS <- native_TSS %>%
  rename(
    score = native_algo_TSS,
    algo = native_algos
  )

native_ROC <- native_ROC %>%
  rename(
    score = native_algo_ROC,
    algo = native_algos
  )

all_native_metrics <- rbind(native_TSS, native_ROC)

# create dummy data frame for 0.7 threshold  of TSS
dummy2 <- data.frame(X = c("ROC", "TSS"),
                     Z = c(0.9, 0.7))

# Change the names of facet labels for plot
metric.labs <- c("(a) ROC", "(b) TSS")
names(metric.labs) <- c("ROC", "TSS")

# Plot
ggplot(all_native_metrics) +
  geom_bar(aes(score, reorder(algo, score)), 
           stat = "identity",
           fill = "white",
           color = "black") +
  xlab("Accuracy Score") +
  ylab("Algorithm") +
  facet_grid(~metric,
             labeller = labeller(metric = metric.labs)) +
  geom_vline(data = subset(all_native_metrics, metric == "TSS"), aes(xintercept = 0.71),
             color = "red", lty = 2, lwd = 0.8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("native-metrics-plot.png", plot = last_plot(), 
       width = 15, height = 10, units = "cm")


## 6.3 Algo Metrics Europe Model ####
europe <- load("BMSB/BMSB.europe_model.models.out")
europe <- get(europe)
europe_EM_algo <- BIOMOD_EnsembleModeling(modeling.output = europe,
                                          chosen.models = 'all', # all models
                                          em.by = 'algo', # a total consensus model
                                          eval.metric = c('TSS'), # use TSS for evaluation
                                          eval.metric.quality.threshold = c(0.5), # remove models that score below 0.7
                                          models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                          prob.mean.weight = TRUE,
                                          prob.mean.weight.decay = 'proportional')
europe_algo_metrics <- get_evaluations(europe_EM_algo) 
europe_algo_metrics

europe_algos <- c("GLM", "GBM", "GAM", "CTA", "ANN",
                  "SRE", "FDA", "MARS", "RF", "MAXENT")
europe_algo_TSS <- c("0.702", "0.805", "0.773", "0.879", "0.794",
                     "0.638", "0.716", "0.747", "0.972", "0.731")
europe_algo_ROC <- c("0.910", "0.958", "0.945", "0.976", "0.953",
                     "0.828", "0.924", "0.935", "0.999", "0.941")

europe_TSS <- as.data.frame(cbind(europe_algos, europe_algo_TSS))
europe_TSS$europe_algo_TSS <- as.numeric(europe_TSS$europe_algo_TSS)

europe_ROC <- as.data.frame(cbind(europe_algos, europe_algo_ROC))
europe_ROC$europe_algo_ROC <- as.numeric(europe_ROC$europe_algo_ROC)

europe_TSS$metric <- "TSS"
europe_ROC$metric <- "ROC"

europe_TSS <- europe_TSS %>%
  rename(
    score = europe_algo_TSS,
    algo = europe_algos
  )

europe_ROC <- europe_ROC %>%
  rename(
    score = europe_algo_ROC,
    algo = europe_algos
  )

all_europe_metrics <- rbind(europe_TSS, europe_ROC)

# create dummy data frame for 0.7 threshold  of TSS
dummy2 <- data.frame(X = c("ROC", "TSS"),
                     Z = c(0.9, 0.7))

# Change the names of facet labels for plot
metric.labs <- c("(a) ROC", "(b) TSS")
names(metric.labs) <- c("ROC", "TSS")

# Plot
ggplot(all_europe_metrics) +
  geom_bar(aes(score, reorder(algo, score)), 
           stat = "identity",
           fill = "white",
           color = "black") +
  xlab("Accuracy Score") +
  ylab("Algorithm") +
  facet_grid(~metric,
             labeller = labeller(metric = metric.labs)) +
  geom_vline(data = subset(all_europe_metrics, metric == "TSS"), aes(xintercept = 0.71),
             color = "red", lty = 2, lwd = 0.8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("europe-metrics-plot.png", plot = last_plot(), 
       width = 15, height = 10, units = "cm")


## 6.4 Algo Metrics USA Model ####
usa <- load("BMSB/BMSB.usa.models.out")
usa <- get(usa)
usa_EM_algo <- BIOMOD_EnsembleModeling(modeling.output = usa,
                                          chosen.models = 'all', # all models
                                          em.by = 'algo', # a total consensus model
                                          eval.metric = c('TSS'), # use TSS for evaluation
                                          eval.metric.quality.threshold = c(0.5), # remove models that score below 0.7
                                          models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                          prob.mean.weight = TRUE,
                                          prob.mean.weight.decay = 'proportional')
usa_algo_metrics <- get_evaluations(usa_EM_algo) 
usa_algo_metrics

usa_algos <- c("GLM", "GBM", "GAM", "CTA", "ANN",
                  "SRE", "FDA", "MARS", "RF", "MAXENT")
usa_algo_TSS <- c("0.663", "0.772", "0.756", "0.844", "0.759",
                     "0.607", "0.708", "0.712", "0.948", "0.763")
usa_algo_ROC <- c("0.885", "0.942", "0.935", "0.958", "0.935",
                     "0.809", "0.913", "0.918", "0.998", "0.947")

usa_TSS <- as.data.frame(cbind(usa_algos, usa_algo_TSS))
usa_TSS$usa_algo_TSS <- as.numeric(usa_TSS$usa_algo_TSS)

usa_ROC <- as.data.frame(cbind(usa_algos, usa_algo_ROC))
usa_ROC$usa_algo_ROC <- as.numeric(usa_ROC$usa_algo_ROC)

usa_TSS$metric <- "TSS"
usa_ROC$metric <- "ROC"

usa_TSS <- usa_TSS %>%
  rename(
    score = usa_algo_TSS,
    algo = usa_algos
  )

usa_ROC <- usa_ROC %>%
  rename(
    score = usa_algo_ROC,
    algo = usa_algos
  )

all_usa_metrics <- rbind(usa_TSS, usa_ROC)

# create dummy data frame for 0.7 threshold  of TSS
dummy2 <- data.frame(X = c("ROC", "TSS"),
                     Z = c(0.9, 0.7))

# Change the names of facet labels for plot
metric.labs <- c("(a) ROC", "(b) TSS")
names(metric.labs) <- c("ROC", "TSS")

# Plot
ggplot(all_usa_metrics) +
  geom_bar(aes(score, reorder(algo, score)), 
           stat = "identity",
           fill = "white",
           color = "black") +
  xlab("Accuracy Score") +
  ylab("Algorithm") +
  facet_grid(~metric,
             labeller = labeller(metric = metric.labs)) +
  geom_vline(data = subset(all_usa_metrics, metric == "TSS"), aes(xintercept = 0.71),
             color = "red", lty = 2, lwd = 0.8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("usa-metrics-plot.png", plot = last_plot(), 
       width = 15, height = 10, units = "cm")


## 6.5 Algo Metrics World Model ####
world <- load("BMSB/BMSB.world_model.models.out")
world <- get(world)
world_EM_algo <- BIOMOD_EnsembleModeling(modeling.output = world,
                                       chosen.models = 'all', # all models
                                       em.by = 'algo', # a total consensus model
                                       eval.metric = c('TSS'), # use TSS for evaluation
                                       eval.metric.quality.threshold = c(0.5), # remove models that score below 0.7
                                       models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
world_algo_metrics <- get_evaluations(world_EM_algo) 
world_algo_metrics

world_algos <- c("GLM", "GBM", "GAM", "CTA", "ANN",
               "SRE", "FDA", "MARS", "RF", "MAXENT")
world_algo_TSS <- c("0.807", "0.887", "0.880", "0.929", "0.866",
                  "0.710", "0.814", "0.837", "0.973", "0.868")
world_algo_ROC <- c("0.958", "0.979", "0.977", "0.981", "0.975",
                  "0.863", "0.964", "0.968", "1.000", "0.978")

world_TSS <- as.data.frame(cbind(world_algos, world_algo_TSS))
world_TSS$world_algo_TSS <- as.numeric(world_TSS$world_algo_TSS)

world_ROC <- as.data.frame(cbind(world_algos, world_algo_ROC))
world_ROC$world_algo_ROC <- as.numeric(world_ROC$world_algo_ROC)

world_TSS$metric <- "TSS"
world_ROC$metric <- "ROC"

world_TSS <- world_TSS %>%
  rename(
    score = world_algo_TSS,
    algo = world_algos
  )

world_ROC <- world_ROC %>%
  rename(
    score = world_algo_ROC,
    algo = world_algos
  )

all_world_metrics <- rbind(world_TSS, world_ROC)

# create dummy data frame for 0.7 threshold  of TSS
dummy2 <- data.frame(X = c("ROC", "TSS"),
                     Z = c(0.9, 0.7))

# Change the names of facet labels for plot
metric.labs <- c("(a) ROC", "(b) TSS")
names(metric.labs) <- c("ROC", "TSS")

# Plot
ggplot(all_world_metrics) +
  geom_bar(aes(score, reorder(algo, score)), 
           stat = "identity",
           fill = "white",
           color = "black") +
  xlab("Accuracy Score") +
  ylab("Algorithm") +
  facet_grid(~metric,
             labeller = labeller(metric = metric.labs)) +
  geom_vline(data = subset(all_world_metrics, metric == "TSS"), aes(xintercept = 0.71),
             color = "red", lty = 2, lwd = 0.8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("world-metrics-plot.png", plot = last_plot(), 
       width = 15, height = 10, units = "cm")


