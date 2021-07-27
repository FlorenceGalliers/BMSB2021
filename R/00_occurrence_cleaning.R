## HEADER ####
## Florence Galliers
## 00 Occurence Import and Cleaning 
## Last Edited: 2021-07-15

## Install Packages ####
library(biomod2)
library(dismo)
library(openxlsx)
library(CoordinateCleaner)
library(mapview)


## Contents ####
## 1.0 Download GBIF Occurrences
## 2.0 Import published occurrences
## 3.0 Data Cleaning
##    3.1 Coordinate Uncertainty
##    3.2 NA Values
##    3.3 Country Centroids and Duplicate Records
## 4.0 Save final file to .csv

## 1.0 Download GBIF Occurrences ####
# Download occurrences from GBIF for the brown marmorated stink bug 
# (Halyomorpha halys) using the {dismo} function gbif()
bmsb <- gbif("halyomorpha", "halys*",
             geo = T) # geo = T only downloads occurrences that have lat/lon
# Put these occurrences into xlsx file and save
write.xlsx(bmsb, "data/00_occurrences.xlsx")

# Read this file back in
bmsb <- read.xlsx(paste0(wd$data,"00_occurrences.xlsx"))
# Subset bmsb dataframe so it only contains required columns
bmsb <- bmsb[, c("country", 
                 "coordinateUncertaintyInMeters", 
                 "lat", 
                 "lon")]
# Re-name coordinateUncertaintyInMeters column to 'uncertainty'
names(bmsb) <- c("country", "uncertainty", "lat", "lon")
# View data frame
head(bmsb)

## 2.0 Import published occurrences ####
# Import file containing occurrences from published sources
extra <- read.xlsx("data/00_extra_occurrences.xlsx")
# Subset columns to only those that are required (to match GBIF columns)
extra <- extra[, c("country", "uncertainty", "lat", "lon")]
# Merge GBIF occurrences with extra occurrences
occur <- merge(bmsb, extra, 
               all.x = TRUE, # all GBIF occurrences
               all.y = TRUE) # all extra occurrences

## 3.0 Data Cleaning ####
# Check for NA values
summary(occur)
# There are 2903 NAs for lat and lon - this will need to be removed at some stage

## 3.1 Coordinate Uncertainty
# Remove occurrences with uncertainty greater than 5km (5000m)
uncertain <- which(occur$uncertainty>5000)
occur <- occur[-uncertain, ]

# Remove the uncertainty column from the occurrences data as it is no longer needed
occur <- occur[, c("country", "lat", "lon")]

## 3.2 NA Values
# Now remove the NAs
occur <- na.omit(occur)
# Left with 13608 observations

## 3.3 Country Centroids and Duplicate Records
# Overview of countries present in occurrences
table(occur$country)
# Use CoordinateCleaner package to remove points near country centroids
occur <- cc_cen(occur, lon = "lon", lat = "lat", verbose = TRUE) 
# removes 60 records
# Use same package to remove duplicate records that have identical lat/lon
occur <- cc_dupl(occur, lon = "lon", lat = "lat", species = NULL, additions = NULL) 
# removes 1577 records

# Now remove the country column
occur <- occur[, c("lon", "lat")]


## 4.0 Save final file to .csv ####
# Left with 11971 cleaned BMSB occurrence records ready to go.
# Save this to csv file
write.csv(occur, "data/00_clean_occur.csv")
