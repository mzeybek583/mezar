
# Cross sectioner

# Load libraries ----------------------------------------------------------
time <- proc.time()
library(lidR)
library(dbscan)
library(ggplot2)
library(lmodel2)

library(concaveman)
library(rgdal)
library(raster)

# Read data ---------------------------------------------------------------
path <- "E:/DR_Sonrasi_Projeler/Mezar/R_test/test4_full.las"

#las <- readLAS(files = "test_2.las")
las <- readLAS(files = path)
#plot(las)


# Ground Segmentation -----------------------------------------------------

mycsf <- csf(FALSE, 0.1,0.5, time_step = 0.65)
las <- classify_ground(las, mycsf)

#plot(las, color = "Classification")

#dtm = grid_terrain(las, res = 0.5, algorithm = tin())

#plot(dtm)

las <- normalize_height(las, tin())

#plot(las)

# Cross section -----------------------------------------------------------

las.crs <- filter_poi(las, Z>=0.3 & Z<=0.5)

#plot(las.crs)
writeLAS(las.crs, "las.crs.las")
