## HEADER ####
## Florence Galliers
## Climatic and Habitat Suitability over Europe - BMSB
## Last Edited: 2021-07-22

## Install packages ####
library(climateStability)
library(biomod2)
library(dismo)
library(openxlsx)
library(maptools)
library(usdm)
library(tictoc)
library(tidyverse)
library(rasterVis)

## Table of Contents ####
## 1.0 EUROPE HABITAT SUITABILITY
## 1.1 Project Invaded Model over Europe 
## 1.2 Land Cover Data (CORINE)
## 1.3 NDVI Data
## 2.0 UK Habitat Suitability
## 2.1 Corine Land Cover
## 2.2 Binary Land Data
## 2.3 NDVI Data
## 2.4 Binary Land * NDVI = HABITAT
## 2.5 HABITAT * CLIMATIC = HAB-CLIM

## 1.0 EUROPE HABITAT SUITABILITY #####

## 1.1 Project Invaded Model over Europe ####
# Get Invaded Model Object from files directory
invaded <- load("BMSB/BMSB.invaded_model.models.out")
invaded <- get(invaded)
# Get Ensemble Model Object from file directory
ensemble <- load("BMSB/BMSB.invaded_modelensemble.models.out")
ensemble <- get(ensemble)

# Get European extent
europe_extent <- extent(c(-12, 45, 34, 72))

# Get Climate Variables for European extent
world_preds <- stack("data/world_preds_red_2.5.grd")
# Crop reduced variables to European Extent
europe_preds <- raster::crop(world_preds, europe_extent)
europe_preds <- stack(europe_preds)

# Project model over Europe
world_onto_europe <- BIOMOD_Projection(modeling.output = invaded, # invaded model object
                                   new.env = europe_preds, # european predictors
                                   proj.name = 'world_onto_europe',
                                   selected.models = "all",
                                   build.clamping.mask = FALSE)

# Ensemble Forecasting for final output
EF_world_onto_europe <- BIOMOD_EnsembleForecasting(projection.output = world_onto_europe,
                                               EM.output = ensemble,
                                               total.consensus = TRUE)

EF_europe <- get_predictions(EF_world_onto_europe)
ensemble_preds_w <- EF_europe[[2]] 
plot(ensemble_preds_w)

## 1.2 Land Cover Data (CORINE) ####
# Import land cover data file
corine <- raster("data/02_corine_land.tif")
land_crs <- crs(corine)
# Import climatic suitability predictions
europe_proj<- stack("BMSB/proj_world_onto_europe/proj_world_onto_europe_BMSB_ensemble.grd")
europe_proj <- europe_proj[[2]]
# reproject climatic suitability into crs that matches land cover
europe_climatic <- projectRaster(from = europe_proj, crs = land_crs, method = "bilinear")
# get extent of land cover
extent2 <- extent(corine)
# crop climatic to same extent
europe_clim <- raster::crop(europe_climatic, extent2)
# crop land to climatic extent
land <- raster::crop(corine, extent(europe_clim))
hist(land, breaks = 128)
# reclassify land cover classes into 1 and 0 for suitable and unsuitable
rcl <- matrix(c(1, 1, 2, 1, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 1,
                11, 0, 12, 1, 13, 1,14, 0, 15, 1,16, 1,17, 1,18, 1,19, 1,
                20, 1,21, 1,22, 1,23, 1,24, 0,25, 1, 26, 0, 27, 0, 28, 0, 29, 1, 
                30, 0, 31, 0, 32, 0, 33, 0, 34, 0, 35, 0, 36, 0,37, 0, 38, 0, 39, 0,
                40, 0, 41, 0,42, 0,43, 0,44, 0, 48, 0, 49, 0, 50, 0, 128, NA),
              ncol = 2,
              byrow = TRUE)
recl_land <- reclassify(land, 
                        rcl = rcl)
# put binary land into a raster file
writeRaster(recl_land, 
            filename = "output/03_binary_land_europe.grd", 
            options = "INTERLEAVE=BAND", 
            overwrite = TRUE)
# reimport raster file
binary_land <- raster("output/binary_land_europe.grd")

## 1.3 UK NDVI DATA ####
# April
april <- raster("data/ndvi/april.tiff")
# find 0 values (these = 0.936 - an arbitrary number )
na_values_april <- which(values(april) == 0.936)
# make those values NA values
values(april)[na_values_april] <- NA
#plot(april)
# May
may <- raster("data/ndvi/may.tiff")
na_values_may <- which(values(may) == 0.936)
values(may)[na_values_may] <- NA
#plot(may)
# June
june <- raster("data/ndvi/june.tiff")
na_values_june <- which(values(june) == 0.936)
values(june)[na_values_june] <- NA
#plot(june)
# July
july <- raster("data/ndvi/july.tiff")
na_values_july <- which(values(july) == 0.936)
values(july)[na_values_july] <- NA
#plot(july)
# August
august <- raster("data/ndvi/august.tiff")
na_values_august <- which(values(august) == 0.936)
values(august)[na_values_august] <- NA
#plot(august)
# Sept 
sept <- raster("data/ndvi/sept.tiff")
na_values_sept <- which(values(sept) == 0.936)
values(sept)[na_values_sept] <- NA
#plot(sept)

# combine all 6 months to make on ndvi layer
ndvi <- (april+may+june+july+august+sept)/6
plot(ndvi)
# reproject ndvi into correct projection
ndvi_reproj <- projectRaster(from = ndvi, crs = land_crs, method = "bilinear")
# resample ndvi to match the land cover dataset
ndvi_final <- raster::resample(ndvi_reproj, binary_land,  method = "bilinear")
# rescale ndvi between 0 and 1
ndvi_1 <- rescale0to1(ndvi_final)
#writeRaster(ndvi_1, 
 ##           filename = "output/03_ndvi.grd", 
   #         options = "INTERLEAVE=BAND", 
    #        overwrite = TRUE)
# ndvi_1 <- raster("output/03_ndvi.grd")
# plot final ndvi dataset
plot(ndvi_1)
# combine land and ndvi together to form a habitat layer
habitat <- (binary_land*ndvi_1)
# plot habitat layer
plot(habitat)

# combine with climatic predictions for Europe from worldwide model
climate <- raster::resample(europe_clim, habitat,  method = "bilinear")
plot(climate)
climate_1 <- rescale0to1(climate)
plot(climate_1)

hab_clim <- (habitat*climate_1)
plot(hab_clim)
writeRaster(hab_clim, 
            filename = "output/03_hab_clim.grd", 
            options = "INTERLEAVE=BAND", 
            overwrite = TRUE)

# overlay occurrences
occur <- read.csv("data/00_clean_occur.csv")
coordinates(occur) <- ~lon+lat
projection(occur) <- CRS('+proj=longlat +datum=WGS84')
occur2 <- spTransform(occur, crs(hab_clim))
occur_uk <- raster::crop(occur2, extent(hab_clim))

# plot
plot(hab_clim)
points(occur_uk, pch = 1, cex = 0.2)

xx <- extent(2854400, 3865000, 3006000, 4216700)
climate_plot <- raster::crop(climate_1, xx)
habitat_plot <- raster::crop(habitat, xx)
hab_clim_plot <- raster::crop(hab_clim, xx)
extent(climate_plot) <- alignExtent(xx, climate_plot)
extent(hab_clim_plot) <- alignExtent(xx, hab_clim_plot)


## 2.0 UK HABITAT SUITABILITY ####

# Get UK extent
uk_extent <- extent(c(-10.67, 1.87, 49.92, 59.48))

## 2.1 CORINE LAND COVER ####
# import land cover data
corine <- raster("data/02_corine_land.tif")
land_crs <- crs(corine)
# import climatic suitability projection of invaded model onto UK
uk_proj<- stack("BMSB/proj_invaded_onto_uk/proj_invaded_onto_uk_BMSB_ensemble.grd")
uk_proj <- uk_proj[[2]]
# crop climatic projection to same extent
uk_clim <- raster::crop(uk_proj, uk_extent)
# reproject climatic suitability into crs that matches land cover
uk_climatic <- projectRaster(from = uk_clim, crs = land_crs, method = "bilinear")

## 2.2 Binary Land Data ####
# Import raster file of binary land data (this file is the UK only)
binary_land <- raster("output/03_binary_land.grd")

## 2.3 NDVI Data ####
# import ndvi data
ndvi_1 <- raster("output/03_ndvi.grd")
# Crop to the extent of the UK
uk_ndvi <- raster::crop(ndvi_1, extent(uk_climatic))

# plot final ndvi dataset
plot(uk_ndvi)
plot(binary_land)

## 2.4 Binary Land * NDVI = HABITAT ####
# combine land and ndvi together to form a habitat layer
habitat <- (binary_land*uk_ndvi)
habitat <- rescale0to1(habitat)
# write habitat layer to raster
writeRaster(habitat, 
            filename = "output/03_habitat_only_uk.grd", 
            options = "INTERLEAVE=BAND", 
            overwrite = TRUE)
# plot habitat layer
plot(habitat)

# REIMPORT habitat layer
habitat <- raster("output/03_habitat_only_uk.grd")

## 2.5 HABITAT * CLIMATIC = HAB-CLIM ####
# Combine habitat with climatic predictions for Europe from worldwide model
climate <- raster::resample(uk_climatic, habitat,  method = "bilinear")
climate_1 <- rescale0to1(climate)
plot(climate_1)

writeRaster(climate_1,
            filename = "03_uk_clim.grd",
            options = "INTERLEAVE=BAND",
            overwrite = TRUE)

hab_clim <- (habitat*climate_1)
plot(hab_clim)
writeRaster(hab_clim, 
            filename = "03_hab_clim_uk.grd", 
            options = "INTERLEAVE=BAND", 
            overwrite = TRUE)


