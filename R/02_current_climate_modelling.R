## HEADER ####
## Florence Galliers
## Biomod2 SDM 
## Last Edited: 2021-07-22

## Install packages
library(biomod2)
library(blockCV)
library(dismo)
library(openxlsx)
library(maptools)
library(usdm)
library(tictoc)
library(tidyverse)
library(rasterVis)

## 0.0 Table of Contents ####
## 0.1 Read in Occurrences
## 0.2 Create extents for regions
## 0.3 Prepare predictor variables
## 0.4 Collinearity checks
## 1.0 EUROPE MODEL
##    1.1 Occurrences
##    1.2 Data Formatting
##    1.3 Pseudo-absence replication
##    1.4 Run the Model
##    1.5 Ensemble Modelling
##    1.6 Project model over UK Extent
##    1.7 Ensemble Forecasting
## 2.0 USA MODEL
##    2.1 Collinearity Checks
##    2.2 Data Formatting
##    2.3 Create pseudo-absence data
##    2.4 Run the Model
##    2.5 Ensemble Modelling
##    2.6 Project model over UK
##    2.7 Ensemble Forecasting
## 3.0 NATIVE MODEL
##    3.1 Collinearity Checks
##    3.2 Data Formatting
##    3.3 Pseudo-absence points
##    3.4 Run the model
##    3.5 Ensemble Modelling
##    3.6 Project model over UK
##    3.7 Ensemble Forecasting
## 4.0 WORLDWIDE MODEL
##    4.1 Occurrences
##    4.2 Data Formatting
##    4.3 Pseudo-absences
##    4.4 Run the model
##    4.5 Ensemble Modelling
##    4.6 Project models over uk
##    4.7 Ensemble Forecasting
## 5.0 INVADED RANGE MODEL
##    5.1 Occurrences
##    5.2 Data Formatting
##    5.3 Pseudo-absences
##    5.4 Run the model
##    5.5 Ensemble Modelling
##    5.6 Project models over UK
##    5.7 Ensemble Forecasting

# 0.1 Read in Occurrences ####
occur <- read.csv("data/00_clean_occur.csv")
# turn occurrences into spatial points data frame and add projection
coordinates(occur) <- ~lon+lat
projection(occur) <- CRS('+proj=longlat +datum=WGS84')

# 0.2 Create extents for regions ####
usa_extent <- extent(c(-128.2961, -55.3049, 20.79696, 55.83368))
europe_extent <- extent(c(-12, 45, 34, 72))
native_extent <- extent(c(92.54915, 149.632, 1.643559, 47.89202))
world_extent <- extent(c(-171.2854 , 171.4799, -65, 90.68569)) # world extent
## UK Extent
uk_extent <- extent(c(-10.67, 1.87, 49.92, 59.48))

# 0.3 Prepare predictor variables ####

## DATA CAN BE IMPORTED BELOW THIS
## THIS IS JUST FORMATTING THE FILES
## Import predictor file stack
path <- file.path("data/CHELSA")
files <- list.files(path, full.names=TRUE )
predictors <- stack(files)
predictors

# Upscale to 2.5 resolution from the baseline 0.5
preds_2.5_world <- aggregate(predictors, # stack to aggregate
                             fact = 5, # aggregate by a factor of 5
                             fun = mean) # use mean as the function
preds_2.5_world <- stack(preds_2.5_world)
world_predictors <- raster::crop(world_predictors, world_extent)
world_predictors <- stack(world_predictors)
# write world predictors at 2.5 resolution into a file 
writeRaster(world_predictors, 
            filename = "preds_2.5_world.grd", 
            options = "INTERLEAVE=BAND", 
            overwrite = TRUE)

# REIMPORT this file as the predictors to use for modelling
world_predictors <- stack("preds_2.5_world.grd")

# 0.4 Collinearity checks ####
# Reduce variables to those that are not correlated
# Use threshold of 5
low_vif_world <- vifstep(world_predictors, th = 5)
# exclude variables with multicollinearity from predictor set
world_preds_reduced <- exclude(world_predictors, low_vif_world)
# Write raster file with the reduced variables
writeRaster(world_preds_reduced, 
            filename = "world_preds_red_2.5.grd", 
            options = "INTERLEAVE=BAND", 
            overwrite = TRUE)

# REIMPORT raster file of reduced predictor variables
world_preds_reduced <- stack("data/world_preds_red_2.5.grd")

## 1.0 EUROPE MODEL ####
# 1.1 Occurrences ####
# Crop occurrences to Europe
europe_occur <- raster::crop(occur, europe_extent)
# Re import the predictors file, this is 2.5 resolution
# Crop to Europe extent
europe_preds_reduced <- raster::crop(world_preds_reduced, europe_extent)
europe_preds_reduced <- stack(europe_preds_reduced)
# Make occurrences have the same coordinate reference system as the predictors
europe_occur <- spTransform(europe_occur, crs(world_preds_reduced))

# 1.2 Data Formatting ####
# Format the data before it can be modeled
BMSB_data_europe <- BIOMOD_FormatingData(
  resp.var = europe_occur, # occurrence data
  expl.var = europe_preds_reduced, # predictors data
  resp.name = "BMSB",
  PA.nb.rep = 2, # no of sets of pseudo-absence points created
  PA.nb.absences = 10000, # Set the number of PA points to 10000
  PA.strategy = "random")

# 1.3 Pseudo-absence replication ####
# Get the pseudo-absence data from the formatted dataset
BMSB_all <- SpatialPointsDataFrame(coords = BMSB_data_europe@coord,
                                   data = data.frame(BMSB = BMSB_data_europe@data.species),
                                   proj = CRS(proj4string(europe_occur)))
# Make the pseudo-absences = 0, so presences will = 1
BMSB_all$BMSB[is.na(BMSB_all$BMSB)] <- 0

# use default settings for all of the models 
myBiomodOption <- BIOMOD_ModelingOptions()

# make data table of train/test data split
sb <- spatialBlock(speciesData = BMSB_all,
                   rasterLayer = europe_preds_reduced,
                   theRange = 150000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0,
                   seed = 11)

DataSplitTable <- sb$biomodTable

# 1.4 Run the model ####
# Model using 10 types of method
europe_model <- BIOMOD_Modeling(BMSB_data_europe,
                                models = c("GLM", "GBM", "GAM", "CTA", "ANN",
                                           "SRE", "FDA", "MARS", "RF",
                                           "MAXENT.Phillips"),
                                model.options = myBiomodOption,
                                DataSplitTable = DataSplitTable,
                                VarImport = 0, # one fold cross validation
                                models.eval.meth = c('TSS','ROC'),
                                modeling.id = "europe_model")

# Graph results
europe_model_results <- models_scores_graph(europe_model, by = "models", metrics = c("ROC", "TSS"))

# 1.5 Ensemble Modelling ####
biomodEM_europe <- BIOMOD_EnsembleModeling(modeling.output = europe_model,
                                           chosen.models = 'all', # all models
                                           em.by = 'all', # a total consensus model
                                           eval.metric = c('TSS'), # use TSS for evaluation
                                           eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                           models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                           prob.mean.weight = TRUE,
                                           prob.mean.weight.decay = 'proportional') 
# Get evaluations
get_evaluations(biomodEM_europe)

# 1.6 Project model over UK Extent ####
# UK predictors at 0.5 resolution
# write world predictors at 2.5 resolution into a file 
#writeRaster(uk, 
 #           filename = "uk_preds.grd", 
  #          options = "INTERLEAVE=BAND", 
   #         overwrite = TRUE)
# reimport this file as the predictors
uk <- stack("data/uk_preds.grd")
uk_preds <- exclude(uk, low_vif_world)

europe_onto_uk <- BIOMOD_Projection(modeling.output = europe_model,
                                    new.env = uk_preds,
                                    proj.name = 'europe_onto_uk',
                                    selected.models = "all",
                                    build.clamping.mask = FALSE)
# save predictions into an object
europe_onto_uk_proj <- get_predictions(europe_onto_uk)
europe_onto_uk_proj
# plot all predictions 
plot(europe_onto_uk_proj)

# 1.7 Ensemble Forecasting ####
EF_europe_onto_UK <- BIOMOD_EnsembleForecasting(projection.output = europe_onto_uk,
                                                EM.output = biomodEM_europe,
                                                total.consensus = TRUE)

ensemble_preds_europe <- get_predictions(EF_europe_onto_UK)
ensemble_preds_europe_weighted <- ensemble_preds_europe[[2]] 
plot(ensemble_preds_europe_weighted)

## 2.0 USA MODEL ####
usa_predictors <- stack("data/preds_2.5_usa.grd")
# crop occurences to just USA extent
usa_occur <- raster::crop(occur, usa_extent)
usa_occur <- spTransform(usa_occur, crs(usa_predictors))
# 2.1 Collinearity Checks ####
# reduce variables to those that are not correlated
usa_preds_reduced <- exclude(usa_predictors, low_vif_world)

# 2.2 Data Formatting ####
# Format the data before it can be modeled
BMSB_data_usa <- BIOMOD_FormatingData(
  resp.var = usa_occur, # occurrence data
  expl.var = usa_preds_reduced, # predictors data
  resp.name = "BMSB",
  PA.nb.rep = 2, # 2 different sets of pseudo-absence points created
  PA.nb.absences = 10000, # Set the number of PA points =  no. presences
  PA.strategy = "random"
)

# 2.3 Create pseudo-absence data ####
# Get the pseudo-absence data from the formatted dataset
BMSB_all_usa <- SpatialPointsDataFrame(coords = BMSB_data_usa@coord,
                                       data = data.frame(BMSB = BMSB_data_usa@data.species),
                                       proj = CRS(proj4string(usa_occur)))
# Make the pseudo-absences = 0, so presences will = 1
BMSB_all_usa$BMSB[is.na(BMSB_all_usa$BMSB)] <- 0

# make data table of train/test data split
sb_usa <- spatialBlock(speciesData = BMSB_all_usa,
                       rasterLayer = usa_preds_reduced,
                       theRange = 150000, # size of the blocks
                       k = 5,
                       selection = "random",
                       iteration = 100, # find evenly dispersed folds
                       biomod2Format = TRUE,
                       xOffset = 0, # shift the blocks horizontally
                       yOffset = 0,
                       seed = 11)

DataSplitTable_usa <- sb_usa$biomodTable

# 2.4 Run the Model ####
# Model using 10 types of method
model_usa <- BIOMOD_Modeling(BMSB_data_usa,
                             models = c("GLM", "GBM", "GAM", "CTA", "ANN",
                                        "SRE", "FDA", "MARS", "RF", "MAXENT.Phillips"),
                             model.options = myBiomodOption,
                             DataSplitTable = DataSplitTable_usa,
                             models.eval.meth = c('TSS','ROC'),
                             modeling.id = "usa"
)

# Evaluate all models that were run
model_usa
get_evaluations(model_usa)

# Graph results
model_usa_results <- models_scores_graph(model_usa, 
                                         by = "models", 
                                         metrics = c("ROC", "TSS"))

# 2.5 Ensemble Modelling ####
EM_usa <- BIOMOD_EnsembleModeling(modeling.output = model_usa,
                                  chosen.models = 'all', # all models
                                  em.by = 'all', # a total consensus model
                                  eval.metric = c('TSS'), # use TSS for evaluation
                                  eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                  models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                  prob.mean.weight = TRUE,
                                  prob.mean.weight.decay = 'proportional') 

EM_usa
get_evaluations(EM_usa)

# 2.6 Project model over UK ####
usa_onto_uk <- BIOMOD_Projection(modeling.output = model_usa,
                                          new.env = uk_preds,
                                          proj.name = 'usa_onto_uk',
                                          selected.models = "all",
                                          build.clamping.mask = FALSE)


usa_onto_uk_proj <- get_predictions(usa_onto_uk)
usa_onto_uk_proj
plot(usa_onto_uk_proj)

# 2.7 Ensemble Forcasting ####
EF_usa_onto_uk <- BIOMOD_EnsembleForecasting(projection.output = usa_onto_uk,
                                            EM.output = EM_usa,
                                            total.consensus = TRUE)

ensemble_preds_usa <- get_predictions(EF_usa_onto_uk)
ensemble_preds_usa_weighted <- ensemble_preds_usa[[2]] 
ensemble_preds_usa_mean <- ensemble_preds_usa[[1]]
plot(ensemble_preds_usa_weighted)
plot(ensemble_preds_usa_mean)

## 3.0 NATIVE MODEL ####
# crop occurences to just USA extent
native_occur <- raster::crop(occur, native_extent)
# Import file for the native predictors, 2.5 resolution
preds_2.5_nat <- raster::crop(world_preds_reduced, native_extent)
native_occur <- spTransform(native_occur, crs(preds_2.5_nat))
# 3.1 Collinearity Checks ####
# reduce variables to those that are not correlated
nat_preds_reduced <- stack(preds_2.5_nat)
# REIMPORT Native Variables 
nat_preds_reduced <- stack("data/preds_2.5_nat.grd")
# 3.2 Data Formatting ####
# Format the data before it can be modeled
nat_BMSB_data <- BIOMOD_FormatingData(
  resp.var = native_occur, # occurrence data
  expl.var = nat_preds_reduced, # predictors data
  resp.name = "BMSB",
  PA.nb.rep = 2, # 3 different sets of pseudo-absence points created
  PA.nb.absences = 10000, # Set the number of PA points =  no. presences
  PA.strategy = "random"
)

# use default settings for now 
myBiomodOption <- BIOMOD_ModelingOptions()

# 3.3 Pseudo-absence points ####
# Get the pseudo-absence data from the formatted dataset
BMSB_all_nat <- SpatialPointsDataFrame(coords = nat_BMSB_data@coord,
                                       data = data.frame(BMSB = nat_BMSB_data@data.species),
                                       proj = CRS(proj4string(native_occur)))
# Make the pseudo-absences = 0, so presences will = 1
BMSB_all_nat$BMSB[is.na(BMSB_all_nat$BMSB)] <- 0

# make data table of train/test data split
sb_nat <- spatialBlock(speciesData = BMSB_all_nat,
                       rasterLayer = nat_preds_reduced,
                       theRange = 150000, # size of the blocks
                       k = 5,
                       selection = "random",
                       iteration = 100, # find evenly dispersed folds
                       biomod2Format = TRUE,
                       xOffset = 0, # shift the blocks horizontally
                       yOffset = 0,
                       seed = 11)

DataSplitTable_nat <- sb_nat$biomodTable

# 3.4 Run the Model ####
# Model using 5 types of method
nat_model <- BIOMOD_Modeling(nat_BMSB_data,
                             models = c("GLM", "GBM", "GAM", "CTA", "ANN",
                                        "SRE", "FDA", "MARS", "RF", "MAXENT.Phillips"),
                             model.options = myBiomodOption,
                             DataSplitTable = DataSplitTable_nat,
                             models.eval.meth = c('TSS','ROC'),
                             modeling.id = "native_model"
)

# Evaluate all models that were run
nat_model
nat_model_eval <- get_evaluations(nat_model)
dimnames(nat_model_eval)

# Graph results
nat_model_results <- models_scores_graph(nat_model, 
                                         by = "models", 
                                         metrics = c("ROC", "TSS"))

# 3.5 Ensemble Modelling ####
EM_nat <- BIOMOD_EnsembleModeling(modeling.output = nat_model,
                                  chosen.models = 'all', # all models
                                  em.by = 'all', # a total consensus model
                                  eval.metric = c('TSS'), # use TSS for evaluation
                                  eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                  models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                  prob.mean.weight = TRUE,
                                  prob.mean.weight.decay = 'proportional') 


EM_nat
get_evaluations(EM_nat)

# 3.6 Project model over UK ####
nat_onto_uk <- BIOMOD_Projection(modeling.output = nat_model,
                                          new.env = uk_preds,
                                          proj.name = 'nat_onto_uk',
                                          selected.models = "all",
                                          build.clamping.mask = FALSE)

nat_onto_uk_proj <- get_predictions(nat_onto_uk)
nat_onto_uk_proj
plot(nat_onto_uk_proj)

# 3.7 Ensemble Forcasting ####
EF_nat_onto_uk <- BIOMOD_EnsembleForecasting(projection.output = nat_onto_uk,
                                            EM.output = EM_nat,
                                            total.consensus = TRUE)

nat_ensemble_preds <- get_predictions(EF_nat_onto_uk)
nat_ensemble_preds_weighted <- nat_ensemble_preds[[2]] 
nat_ensemble_preds_mean <- nat_ensemble_preds[[1]]
plot(nat_ensemble_preds_weighted)
plot(nat_ensemble_preds_mean)

## 4.0 WORLDWIDE MODEL ####

# 4.1 Occurrences ####
# Make occurrences have the same coordinate reference system as the predictors
# combine usa, europe and native occurrences to create world occurrences
world_occur <- rbind(usa_occur, europe_occur)
world_occur <- rbind(world_occur, native_occur)
# write into a csv file
write.csv(world_occur, "01_world_occur.csv")
# REIMPORT world occurrences from the .csv file
world_occur <- read.csv("data/01_world_occur")
# Import reduced set of world predictors
world_preds_reduced <- stack("data/world_preds_red_2.5.grd")

# 4.2 Data Formatting ####
# Format the data before it can be modeled
BMSB_data_world <- BIOMOD_FormatingData(
  resp.var = world_occur, # occurrence data
  expl.var = world_preds_reduced, # predictors data
  resp.name = "BMSB",
  PA.nb.rep = 2, # no of sets of pseudo-absence points created
  PA.nb.absences = 10000, # Set the number of PA points to 10000
  PA.strategy = "random")

# 4.3 Pseudo-absences ####
# Get the pseudo-absence data from the formatted dataset
world_BMSB_all <- SpatialPointsDataFrame(coords = BMSB_data_world@coord,
                                         data = data.frame(BMSB = BMSB_data_world@data.species),
                                         proj = CRS(proj4string(world_occur)))
# Make the pseudo-absences = 0, so presences will = 1
world_BMSB_all$BMSB[is.na(world_BMSB_all$BMSB)] <- 0

# make data table of train/test data split
world_sb <- spatialBlock(speciesData = world_BMSB_all,
                         rasterLayer = world_preds_reduced,
                         theRange = 150000, # size of the blocks
                         k = 5,
                         selection = "random",
                         iteration = 100, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, # shift the blocks horizontally
                         yOffset = 0,
                         seed = 11)

world_DataSplitTable <- world_sb$biomodTable

# use default settings for all of the models 
myBiomodOption <- BIOMOD_ModelingOptions()

# 4.4 Run the model ####
# Model using 10 types of method
world_model <- BIOMOD_Modeling(BMSB_data_world,
                               models = c("GLM", "GBM", "GAM", "CTA", "ANN",
                                          "SRE", "FDA", "MARS", "RF",
                                          "MAXENT.Phillips"),
                               model.options = myBiomodOption,
                               DataSplitTable = world_DataSplitTable,
                               VarImport = 0,
                               models.eval.meth = c('TSS','ROC'),
                               modeling.id = "world_model")

# Evaluate all models that were run
world_model
world_model_eval <- get_evaluations(world_model)
dimnames(world_model_eval)

# Graph results
world_model_graph <- models_scores_graph(world_model, by = "models", metrics = c("TSS", "ROC"))

# 4.5 Ensemble Modelling ####
EM_world <- BIOMOD_EnsembleModeling(modeling.output = world_model,
                                          chosen.models = 'all', # all models
                                          em.by = 'all', # a total consensus model
                                          eval.metric = c('TSS'), # use TSS for evaluation
                                          eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                          models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                          prob.mean.weight = TRUE,
                                          prob.mean.weight.decay = 'proportional') 

EM_world
get_evaluations(EM_world)

# 4.6 Project models over uk ####
uk <- stack("uk_preds.grd")
uk_preds_world <- exclude(uk, low_vif_world)

world_onto_uk <- BIOMOD_Projection(modeling.output = world_model,
                                            new.env = uk_preds_world,
                                            proj.name = 'world_onto_uk',
                                            selected.models = "all",
                                            build.clamping.mask = FALSE)
# save predictions into an object
world_onto_uk_proj <- get_predictions(world_onto_uk)
world_onto_uk_proj
# plot all predictions 
plot(world_onto_uk_proj)

# 4.7 Ensemble Forecasting ####
EF_world_onto_uk <- BIOMOD_EnsembleForecasting(projection.output = world_onto_uk,
                                              EM.output = EM_world,
                                              total.consensus = TRUE)

ensemble_preds_world <- get_predictions(EF_world_onto_uk)
ensemble_preds_world_weighted <- ensemble_preds_world[[2]] 
plot(ensemble_preds_world_weighted)
ensemble_preds_world_mean <- ensemble_preds_world[[1]]
plot(ensemble_preds_world_mean)

rasterVis::levelplot(ensemble_preds_world_weighted)
rasterVis::levelplot(ensemble_preds_europe_weighted)
rasterVis::levelplot(ensemble_preds_usa_weighted)


## 5.0 INVADED RANGE MODEL ####

# 5.1 Occurrences ####
# combine usa and europe occurrences to create invaded occurrences
occur <- read.csv("world-occur.csv")
occur <- occur[, -c(1,2,5)]
coordinates(occur) <- ~lon+lat
projection(occur) <- CRS('+proj=longlat +datum=WGS84')
usa_occur <- raster::crop(occur, usa_extent)
europe_occur <- raster::crop(occur, europe_extent)
nat_occur <- raster::crop(occur, native_extent)
invaded_occur <- rbind(usa_occur, europe_occur)
# Write csv file with invaded occurrences in
write.csv(invaded_occur, "data/02_invaded_occur.csv")

# REIMPORT invaded occurrences from csv file
invaded_occur <- read.csv("data/02_invaded_occur.csv")

# Import world predictors
world_preds_reduced <- stack("world_preds_red_2.5.grd")
invaded_occur <- spTransform(invaded_occur, crs(world_preds_reduced))

# 5.2 Data Formatting ####
# Format the data before it can be modeled
BMSB_data_invade <- BIOMOD_FormatingData(
  resp.var = invaded_occur, # occurrence data
  expl.var = world_preds_reduced, # predictors data
  resp.name = "BMSB",
  PA.nb.rep = 2, # no of sets of pseudo-absence points created
  PA.nb.absences = 10000, # Set the number of PA points to 10000
  PA.strategy = "random")

# 5.3 Pseudo-absences ####
# Get the pseudo-absence data from the formatted dataset
invade_BMSB_all <- SpatialPointsDataFrame(coords = BMSB_data_invade@coord,
                                         data = data.frame(BMSB = BMSB_data_invade@data.species),
                                         proj = CRS(proj4string(invaded_occur)))
# Make the pseudo-absences = 0, so presences will = 1
invade_BMSB_all$BMSB[is.na(invade_BMSB_all$BMSB)] <- 0

# make data table of train/test data split
invaded_sb <- spatialBlock(speciesData = invade_BMSB_all,
                         rasterLayer = world_preds_reduced,
                         theRange = 150000, # size of the blocks
                         k = 5,
                         selection = "random",
                         iteration = 100, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, # shift the blocks horizontally
                         yOffset = 0,
                         seed = 11)

invaded_DataSplitTable <- invaded_sb$biomodTable

# use default settings for all of the models 
myBiomodOption <- BIOMOD_ModelingOptions()

# 5.4 Run the model ####
# Model using 10 types of method
invaded_model <- BIOMOD_Modeling(BMSB_data_invade,
                               models = c("GLM", "GBM", "GAM", "CTA", "ANN",
                                          "SRE", "FDA", "MARS", "RF",
                                          "MAXENT.Phillips"),
                               model.options = myBiomodOption,
                               DataSplitTable = invaded_DataSplitTable,
                               VarImport = 0,
                               models.eval.meth = c('TSS','ROC'),
                               modeling.id = "invaded_model")

# Evaluate all models that were run
invaded_model
invaded_model_eval <- get_evaluations(invaded_model)
dimnames(invaded_model_eval)

# Graph results
invaded_model_graph <- models_scores_graph(invaded_model, by = "models", metrics = c("TSS", "ROC"))

# 5.5 Ensemble Modelling ####
EM_invaded <- BIOMOD_EnsembleModeling(modeling.output = invaded_model,
                                    chosen.models = 'all', # all models
                                    em.by = 'all', # a total consensus model
                                    eval.metric = c('TSS'), # use TSS for evaluation
                                    eval.metric.quality.threshold = c(0.7), # remove models that score below 0.7
                                    models.eval.meth = c('TSS','ROC'), # methods used for ensemble evaluation
                                    prob.mean.weight = TRUE,
                                    prob.mean.weight.decay = 'proportional') 

EM_invaded
get_evaluations(EM_invaded)

# 5.6 Project models over uk ####
uk <- stack("uk_preds.grd")
uk_preds_world <- exclude(uk, low_vif_world)

invaded_onto_uk <- BIOMOD_Projection(modeling.output = invaded_model,
                                   new.env = uk_preds_world,
                                   proj.name = 'invaded_onto_uk',
                                   selected.models = "all",
                                   build.clamping.mask = FALSE)
# save predictions into an object
invaded_onto_uk_proj <- get_predictions(invaded_onto_uk)
invaded_onto_uk_proj
# plot all predictions 
plot(invaded_onto_uk_proj)

# 5.7 Ensemble Forecasting ####
EF_invaded_onto_uk <- BIOMOD_EnsembleForecasting(projection.output = invaded_onto_uk,
                                               EM.output = EM_invaded,
                                               total.consensus = TRUE)

ensemble_preds_invaded <- get_predictions(EF_invaded_onto_uk)
ensemble_preds_invaded_weighted <- ensemble_preds_invaded[[2]] 
plot(ensemble_preds_invaded_weighted)
ensemble_preds_invaded_mean <- ensemble_preds_invaded[[1]]
plot(ensemble_preds_invaded_mean)




