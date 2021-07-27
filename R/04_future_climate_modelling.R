## HEADER ####
## Florence Galliers
## Future Climate Predictions
## Last Edited 2021-06-09

## Install Packages ####
library(biomod2)
library(dismo)
library(openxlsx)
library(maptools)
library(usdm)
library(tictoc)
library(tidyverse)
library(rasterVis)

## Table of Contents ####
## 1.0 Variable Importance ####

## Invaded Model Variable Importance
invaded <- load("BMSB/BMSB.invaded_model.models.out")
invaded_model <- get(invaded)
invaded <- BIOMOD_LoadModels(get(invaded))
invaded

expl_data <- get_formal_data(invaded_model,'expl.var')

vi <- list()
for(mod in invaded){
  cat("\n> variables importance of ", mod)
  vi <- c(vi, variables_importance(model = get(mod), 
                                   data = expl_data,
                                   method = "full_rand",
                                   nb_rand = 2))
}
names(vi) <- invaded

vi
# write to csv file
write.csv(vi, "output/variables_imp.csv")

read.csv("output/variables_imp.csv")

uk_extent <- extent(c(-10.67, 1.87, 49.92, 59.48))

## 2.0 Load Models and Extents ####
# Get Invasive Model Object from files directory
world <- load("BMSB/BMSB.invaded_model.models.out")
world <- get(world)
# Get Ensemble Invasive Model Object from file directory
ensemble <- load("BMSB/BMSB.invaded_modelensemble.models.out")
ensemble <- get(ensemble)
# UK Extent
uk_extent <- extent(c(-10.67, 1.87, 49.92, 59.48))


# RCP26 2060 ####
path_26_2060 <- file.path("rcp26_2060")
files_26_2060 <- list.files(path_26_2060, full.names=TRUE)
preds_26_2060 <- stack(files_26_2060)
# Crop predictors to UK extent
preds_26_2060 <- raster::crop(preds_26_2060, uk_extent)
preds_26_2060 <- aggregate(preds_26_2060, fact = 10, fun = mean)
preds_26_2060 <- stack(preds_26_2060)
# Project using future climate data
proj_26_2060 <- BIOMOD_Projection(modeling.output = world,
                                   new.env = preds_26_2060,
                                   proj.name = 'inv_uk_26_2060',
                                   selected.models = "all",
                                   build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_26_2060 <- BIOMOD_EnsembleForecasting(projection.output = proj_26_2060,
                                               EM.output = ensemble,
                                               total.consensus = TRUE)

ensemble_26_2060 <- get_predictions(EF_26_2060)
ensemble_26_2060_w <- ensemble_26_2060[[2]] 
plot(ensemble_26_2060_w)

# RCP45 2060 ####
path_45_2060 <- file.path("rcp45_2060")
files_45_2060 <- list.files(path_45_2060, full.names=TRUE)
preds_45_2060 <- stack(files_45_2060)
# Crop predictors to UK extent
preds_45_2060 <- raster::crop(preds_45_2060, uk_extent)
preds_45_2060 <- aggregate(preds_45_2060, fact = 10, fun = mean)
preds_45_2060 <- stack(preds_45_2060)
# Project using future climate data
proj_45_2060 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_45_2060,
                                  proj.name = 'uk_45_2060',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_45_2060 <- BIOMOD_EnsembleForecasting(projection.output = proj_45_2060,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_45_2060 <- get_predictions(EF_45_2060)
ensemble_45_2060_w <- ensemble_45_2060[[2]] 
plot(ensemble_45_2060_w)

# RCP60 2060 ####
path_60_2060 <- file.path("rcp60_2060")
files_60_2060 <- list.files(path_60_2060, full.names=TRUE)
preds_60_2060 <- stack(files_60_2060)
# Crop predictors to UK extent
preds_60_2060 <- raster::crop(preds_60_2060, uk_extent)
preds_60_2060 <- aggregate(preds_60_2060, fact = 10, fun = mean)
preds_60_2060 <- stack(preds_60_2060)
# Project using future climate data
proj_60_2060 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_60_2060,
                                  proj.name = 'inv_uk_60_2060',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_60_2060 <- BIOMOD_EnsembleForecasting(projection.output = proj_60_2060,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_60_2060 <- get_predictions(EF_60_2060)
ensemble_60_2060_w <- ensemble_60_2060[[2]] 
plot(ensemble_60_2060_w)


# RCP85 2060 ####
path_85_2060 <- file.path("rcp85_2060")
files_85_2060 <- list.files(path_85_2060, full.names=TRUE)
preds_85_2060 <- stack(files_85_2060)
# Crop predictors to UK extent
preds_85_2060 <- raster::crop(preds_85_2060, uk_extent)
preds_85_2060 <- aggregate(preds_85_2060, fact = 10, fun = mean)
preds_85_2060 <- stack(preds_85_2060)
# Project using future climate data
proj_85_2060 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_85_2060,
                                  proj.name = 'inv_uk_85_2060',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_85_2060 <- BIOMOD_EnsembleForecasting(projection.output = proj_85_2060,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_85_2060 <- get_predictions(EF_85_2060)
ensemble_85_2060_w <- ensemble_85_2060[[2]] 
plot(ensemble_85_2060_w)

# RCP26 2080 ####
path_26_2080 <- file.path("rcp26_2080")
files_26_2080 <- list.files(path_26_2080, full.names=TRUE)
preds_26_2080 <- stack(files_26_2080)
# Crop predictors to UK extent
preds_26_2080 <- raster::crop(preds_26_2080, uk_extent)
preds_26_2080 <- aggregate(preds_26_2080, fact = 10, fun = mean)
preds_26_2080 <- stack(preds_26_2080)
# Project using future climate data
proj_26_2080 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_26_2080,
                                  proj.name = 'inv_uk_26_2080',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_26_2080 <- BIOMOD_EnsembleForecasting(projection.output = proj_26_2080,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_26_2080 <- get_predictions(EF_26_2080)
ensemble_26_2080_w <- ensemble_26_2080[[2]] 
plot(ensemble_26_2080_w)

# RCP45 2080 ####
path_45_2080 <- file.path("rcp45_2080")
files_45_2080 <- list.files(path_45_2080, full.names=TRUE)
preds_45_2080 <- stack(files_45_2080)
# Crop predictors to UK extent
preds_45_2080 <- raster::crop(preds_45_2080, uk_extent)
preds_45_2080 <- aggregate(preds_45_2080, fact = 10, fun = mean)
preds_45_2080 <- stack(preds_45_2080)
# Project using future climate data
proj_45_2080 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_45_2080,
                                  proj.name = 'uk_45_2080',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_45_2080 <- BIOMOD_EnsembleForecasting(projection.output = proj_45_2080,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_45_2080 <- get_predictions(EF_45_2080)
ensemble_45_2080_w <- ensemble_45_2080[[2]] 
plot(ensemble_45_2080_w)

# RCP60 2080 ####
path_60_2080 <- file.path("rcp60_2080")
files_60_2080 <- list.files(path_60_2080, full.names=TRUE)
preds_60_2080 <- stack(files_60_2080)
# Crop predictors to UK extent
preds_60_2080 <- raster::crop(preds_60_2080, uk_extent)
preds_60_2080 <- aggregate(preds_60_2080, fact = 10, fun = mean)
preds_60_2080 <- stack(preds_60_2080)
# Project using future climate data
proj_60_2080 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_60_2080,
                                  proj.name = 'inv_uk_60_2080',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_60_2080 <- BIOMOD_EnsembleForecasting(projection.output = proj_60_2080,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_60_2080 <- get_predictions(EF_60_2080)
ensemble_60_2080_w <- ensemble_60_2080[[2]] 
plot(ensemble_60_2080_w)


# RCP85 2080 ####
path_85_2080 <- file.path("rcp85_2080")
files_85_2080 <- list.files(path_85_2080, full.names=TRUE)
preds_85_2080 <- stack(files_85_2080)
# Crop predictors to UK extent
preds_85_2080 <- raster::crop(preds_85_2080, uk_extent)
preds_85_2080 <- aggregate(preds_85_2080, fact = 10, fun = mean)
preds_85_2080 <- stack(preds_85_2080)
# Project using future climate data
proj_85_2080 <- BIOMOD_Projection(modeling.output = world,
                                  new.env = preds_85_2080,
                                  proj.name = 'inv_uk_85_2080',
                                  selected.models = "all",
                                  build.clamping.mask = FALSE)

# Ensemble Forecasting
EF_85_2080 <- BIOMOD_EnsembleForecasting(projection.output = proj_85_2080,
                                         EM.output = ensemble,
                                         total.consensus = TRUE)

ensemble_85_2080 <- get_predictions(EF_85_2080)
ensemble_85_2080_w <- ensemble_85_2080[[2]] 
plot(ensemble_85_2080_w)



