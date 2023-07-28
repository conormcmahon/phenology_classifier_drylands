
library(terra)
library(sf)
library(tidyverse)
library(here)

# NOTE - paths to file names will need to be changed for your use case. 
#   The source files are too large to be distributed with this github package, and in some cases are proprietary.



# Get extent coordinates to evenly sample a raster into chunks of equal size
#   Raster is a spatRaster (from the Terra package)
#   chunk_size is the width/height of each chunk, in pixels
# This code was sourced from Robert Hijmans on Stack Overflow
#  Source: https://stackoverflow.com/questions/71709795/how-to-iterate-crop-raster-in-r-using-a-moving-window
getRasterChunks <- function(rast_in, chunk_size)
{
  cells <- cellFromRowColCombine(rast_in, 
                                 seq(1, nrow(rast_in), by=chunk_size), 
                                 seq(1, ncol(rast_in), by=chunk_size))
  xy <- xyFromCell(rast_in, cells)
  rs <- res(rast_in)
  xy[,1] <- xy[,1] - rs[1]/2
  xy[,2] <- xy[,2] + rs[2]/2
  xy <- cbind(xy[,1], xy[,1] + chunk_size*rs[1], xy[,2] - chunk_size*rs[2], xy[,2])
  return(xy)
}


# ********************************************************************************
# First, for the San Pedro River: 
# ********************************************************************************

# Load in phenology data 
san_pedro_pheno <- rast(here::here("Phenology","san_pedro_2021_phenoseries.tif"))
# Load in DEM and reproject to 10 m from 1 m resolution
dem_san_pedro_2015 <- terra::rast("D:/SERDP/Huachuca/LiDAR/dem/dem_updated_datum.tif")
dem_san_pedro_2017 <- terra::rast("D:/SERDP/San_Pedro/LiDAR/dem/san_pedro_dem_mosaic.tif")
dem_san_pedro_2015_expanded <- terra::extend(dem_san_pedro_2015, terra::union(terra::ext(dem_san_pedro_2015), terra::ext(dem_san_pedro_2017)))
dem_san_pedro_2015_expanded[dem_san_pedro_2015_expanded == 0] <- NA
dem_san_pedro_2017_expanded <- terra::extend(dem_san_pedro_2017, terra::union(terra::ext(dem_san_pedro_2015), terra::ext(dem_san_pedro_2017)))
dem_san_pedro_2017_expanded[dem_san_pedro_2017_expanded == 0] <- NA
dem_san_pedro_1m <- sum(mean(dem_san_pedro_2015_expanded, dem_san_pedro_2017_expanded),
  (dem_san_pedro_2015_expanded * is.na(dem_san_pedro_2017_expanded)),
  (dem_san_pedro_2017_expanded * is.na(dem_san_pedro_2015_expanded)), na.rm=TRUE)
writeRaster(dem_san_pedro_1m, "D:/SERDP/San_Pedro/LiDAR_mosaic/san_pedro_dem_1m.tif", overwrite=TRUE)
# Reproject to 10 m resolution, by chunks (for memory reasons)
san_pedro_dem_chunks <- getRasterChunks(dem_san_pedro_1m, ncol(dem_san_pedro_1m) / 2)
reprojectDEMChunk <- function(chunk_index)
{
  print(paste("Beginning to work on chunk ", chunk_index, sep=""))
  dem_chunk <- crop(dem_san_pedro_1m, san_pedro_dem_chunks[chunk_index,])
  terra::project(dem_chunk, san_pedro_pheno, method="average", mask=TRUE, align=TRUE, filename = paste("D:/SERDP/San_Pedro/LiDAR_mosaic/dem_10m_chunks/chunk_",chunk_index,".tif",sep=""))
}
# Patch some holes in raster at tile edges
#   Fill holes with the mean value from a 5-pixel window
dem_san_pedro_10m <- terra::rast("D:/SERDP/San_Pedro/LiDAR_mosaic/san_pedro_dem_10m.tif")
#holes_area <- st_read("D:/SERDP/San_Pedro/LiDAR_mosaic/holes_to_clean.shp")
slope_san_pedro <- terra::terrain(dem_san_pedro_10m, v="slope", neighbors=8)
slope_san_pedro_holes <- slope_san_pedro #mask(crop(slope_san_pedro, holes_area), holes_area)
dem_san_pedro_holes <- dem_san_pedro_10m #mask(crop(dem_san_pedro_10m, holes_area), holes_area)
bad_areas <- slope_san_pedro_holes > 50
bad_areas <- (bad_areas + (dem_san_pedro_holes == 0)) > 0
dem_without_holes <- dem_san_pedro_holes
dem_without_holes[which(as.logical(as.vector(bad_areas)))] <- NA
neighborhood_mean = focal(dem_without_holes, w=11, fun="mean", na.rm=TRUE)
filled_raster <- sum((is.na(dem_without_holes) * neighborhood_mean),
                      dem_without_holes,
                      na.rm=TRUE)
writeRaster(filled_raster, "D:/SERDP/San_Pedro/LiDAR_mosaic/dem_10m_filled_final.tif", overwrite=TRUE)
writeRaster(terra::terrain(filled_raster, v="slope", neighbors=8), here::here("LiDAR","terrain_slope","san_pedro_slope.tif"), overwrite=TRUE)

srtm_dem <- terra::rast(here::here("SRTM","terrain_elevation","N31W111.hgt"))
srtm_dem_reproj <- terra::project(srtm_dem, san_pedro_pheno)
terra::writeRaster(srtm_dem_reproj, here::here("SRTM","terrain_elevation","san_pedro_srtm.tif"))
terra::writeRaster(terra::terrain(srtm_dem_reproj, v="slope", neighbors=8), here::here("SRTM","terrain_slope","san_pedro_srtm_slope.tif"))

# Load in CHM and reproject to 10 m from 1 m resolution
chm_mosaic_1m <- rast("D:/SERDP/San_Pedro/LiDAR_mosaic/san_pedro_lidar_mosaic.tif")
canopy_10m_max <- terra::project(chm_mosaic_1m, san_pedro_pheno, method="max", mask=TRUE, align=TRUE)
writeRaster(canopy_10m_max, "D:/SERDP/San_Pedro/LiDAR_mosaic/san_pedro_lidar_10m_max.tif", overwrite=TRUE)

# Load flowlines vector file with geometry
flowlines_san_pedro <- st_read("D:/SERDP/Huachuca/NHDPlus/NHDPLUS_H_1505_HU4_GDB.gdb", layer="NHDFlowline")
# Load feature containing stream order
vaa_san_pedro <- st_read("D:/SERDP/Huachuca/NHDPlus/NHDPLUS_H_1505_HU4_GDB.gdb", layer="NHDPlusFlowlineVAA")
# Combine flowline geometry and data
flowdata_san_pedro <- merge(flowlines_san_pedro, vaa_san_pedro, by="NHDPlusID")
flowdata_san_pedro <- st_transform(flowdata_san_pedro, st_crs(dem_san_pedro))
flowdata_san_pedro <- st_crop(flowdata_san_pedro, dem_san_pedro)
st_write(flowdata_san_pedro, here::here("hydrology","spr_flowlines_with_vaa.shp"))


# ********************************************************************************
# Next, for the Santa Ynez River: 
# ********************************************************************************

vandenberg_pheno <- terra::rast(here::here("phenology","vandenberg_2013_phenoseries.tif"))
vandenberg_dem_srtm <- terra::rast(here::here("SRTM","terrain_elevation","N34W121.hgt"))
vandenberg_dem_srtm_reproj <- terra::project(vandenberg_dem_srtm, vandenberg_pheno)
terra::writeRaster(vandenberg_dem_srtm_reproj, here::here("SRTM","terrain_elevation","vandenberg_srtm.tif"))
terra::writeRaster(terra::terrain(vandenberg_dem_srtm_reproj, v="slope", neighbors=8), here::here("SRTM","terrain_slope","vandenberg_srtm_slope.tif"))

# Load flowlines vector file with geometry
flowlines_vandenberg <- st_read("D:/SERDP/Vandenberg/hydrology/NHDPLUS_H_1806_HU4_GDB/NHDPLUS_H_1806_HU4_GDB.gdb", layer="NHDFlowline")
# Load feature containing stream order
vaa_vandenberg <- st_read("D:/SERDP/Vandenberg/hydrology/NHDPLUS_H_1806_HU4_GDB/NHDPLUS_H_1806_HU4_GDB.gdb", layer="NHDPlusFlowlineVAA")
# Combine flowline geometry and data
flowdata_vandenberg <- merge(flowlines_vandenberg, vaa_vandenberg, by="NHDPlusID")
flowdata_vandenberg <- st_transform(flowdata_vandenberg, st_crs(vandenberg_dem_srtm_reproj))
flowdata_vandenberg <- st_crop(flowdata_vandenberg, vandenberg_dem_srtm_reproj)
st_write(flowdata_vandenberg, here::here("hydrology","syr_flowlines_with_vaa.shp"))


# ********************************************************************************
# Next, for the Santa Margarita River: 
# ********************************************************************************

pendleton_pheno <- terra::rast(here::here("phenology","pendleton_2019_phenoseries.tif"))
pendleton_dem_srtm <- terra::rast(here::here("SRTM","terrain_elevation","N33W118.hgt"))
pendleton_dem_srtm_reproj <- terra::project(pendleton_dem_srtm, pendleton_pheno)
terra::writeRaster(pendleton_dem_srtm_reproj, here::here("SRTM","terrain_elevation","pendleton_srtm.tif"))
terra::writeRaster(terra::terrain(pendleton_dem_srtm_reproj, v="slope", neighbors=8), here::here("SRTM","terrain_slope","pendleton_srtm_slope.tif"))
pendleton_dem <- terra::rast("D:/SERDP/Pendleton/LiDAR/DEM/dem_mosaic_10m.tif")
pendleton_dem_reproj <- terra::project(pendleton_dem, pendleton_pheno)
terra::writeRaster(pendleton_dem_reproj, here::here("LiDAR","terrain_elevation","pendleton_dem.tif"))
terra::writeRaster(terra::terrain(pendleton_dem_reproj, v="slope", neighbors=8), here::here("LiDAR","terrain_slope","pendleton_slope.tif"))

pendleton_chm <- terra::rast("D:/SERDP/Pendleton/LiDAR/CHM/chm_10m.tif")
pendleton_chm_reproj <- terra::project(pendleton_chm, pendleton_pheno)
terra::writeRaster(pendleton_chm_reproj, here::here("LiDAR","canopy_height","pendleton_chm.tif"))

# Load flowlines vector file with geometry
flowlines_pendleton <- st_read("D:/SERDP/Pendleton/hydrology/NHDPLUS_H_1807_HU4_GDB/NHDPLUS_H_1807_HU4_GDB.gdb", layer="NHDFlowline")
# Load feature containing stream order
vaa_pendleton <- st_read("D:/SERDP/Pendleton/hydrology/NHDPLUS_H_1807_HU4_GDB/NHDPLUS_H_1807_HU4_GDB.gdb", layer="NHDPlusFlowlineVAA")
# Combine flowline geometry and data
flowdata_pendleton <- merge(flowlines_pendleton, vaa_pendleton, by="NHDPlusID")
flowdata_pendleton <- st_transform(flowdata_pendleton, st_crs(pendleton_dem_reproj))
flowdata_pendleton <- st_crop(flowdata_pendleton, pendleton_dem_reproj)
st_write(flowdata_pendleton, here::here("hydrology","smr_flowlines_with_vaa.shp"))
