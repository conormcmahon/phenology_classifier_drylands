
library(terra)
library(sf)
library(tidyverse)
library(RANN)

# This script generates a relative elevation and creek distance raster from an input NHDPlus flowlines model and DEM terrain model



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

# Generate a relative elevation raster
#   Output is a raster where each cell is height above nearest creek point
#   Inputs are a digital terrain model raster (elev_rast) and an SF feature set of linear features (flowlines) for creek positions
getRelativeElevation <- function(elev_rast, flowlines, raster_filename_prefix, order_thresh, stream_order_variable_name)
{
  print("")
  print("Beginning to process a new set of flowlines and terrain raster to generate relative elevations.")

  # Clip flowlines to extent of raster
  print(paste("  Working set of flowlines with ", nrow(flowlines), " features.", sep=""))
  print(paste("Starting out by cropping to the extent of the target raster: ", terra::ext(elev_rast), sep=""))
  flowlines <- st_crop(st_zm(flowlines), elev_rast)
  # Also remove any flowlines with zero length
  flowlines <- flowlines[as.numeric(st_length(flowlines)) > 0,]
  print(paste("  After cropping, flowlines has ", nrow(flowlines), " features", sep=""))
  
  # Remove flowlines below the order threshold
  flowlines <- flowlines[(as.numeric((as.data.frame(flowlines %>% 
                                                    dplyr::select(match(stream_order_variable_name,names(.))))[,1])) > 1), ]
  print(paste("  After filtering for stream order, flowlines has ", nrow(flowlines), " features", sep=""))
  
  # Rasterize flowlines to match grid of DEM
  print(paste("  Generating rasterized version of flowlines at input DEM resolution: ", terra::res(elev_rast)[1], " by ", terra::res(elev_rast)[1], sep=""))
  flowlines_raster <- terra::rasterize(flowlines, elev_rast, field=stream_order_variable_name)
  flowlines_elev <- terra::mask(elev_rast, flowlines_raster, maskvalues=NA, updatevalue=NA)
  flowlines_raster <- c(flowlines_raster, flowlines_elev)
  # Also, remove any flowline points which are outside the mapped extent of the LiDAR file. These will have NA or 0 for elevation.
  print("  Filtering flowlines raster to remove points with stream order 0 or outside of DEM extent.")
  flowlines_raster[sum(is.na(flowlines_raster)) > 0] <- NA   # Remove cases where order or elevation is NA
  flowlines_raster[(flowlines_raster[[2]] == 0)] <- NA       # Remove cases where elevation is NA
  print("Resulting raster: ")
  print(flowlines_raster)
  terra::writeRaster(flowlines_raster, "D:/SERDP/San_Pedro/LiDAR_mosaic/flowlines_raster_1m.tif", overwrite=TRUE)
  
  # Turn flowlines raster into a dataframe with XY coord and elevation
  #    Because the raster might be very large, it may be prohibitive to pull the whole thing into memory as a matrix at once
  #    So instead, we'll iterate over chunks of it separately
  chunk_extents <- getRasterChunks(flowlines_raster, 1e8)
  flowlines_df <- lapply(1:nrow(chunk_extents), 
                                function(i){ 
                                  print(paste("   Working on chunk ", i, "/", nrow(chunk_extents), " with extent ", chunk_extents[i,1], " ", chunk_extents[i,2], "..."))
                                  return( as.data.frame(crop(flowlines_raster, chunk_extents[i,]), xy=TRUE) )
                                })
  flowlines_df <- bind_rows(flowlines_df)
  
  # Turn dataframe into a SF object and build a KD Tree for searches on it
  #   https://search.r-project.org/CRAN/refmans/RANN/html/nn2.html
  #   Again, we'll need to break this into chunks for analysis on very large rasters
  #   We'll use the same chunk mapping as before: 
  print("Searching for neighbor creek points and getting relative elevation:")
  mapRelativeElevationChunk <- function(i)
  {
    print(paste("  Getting neighbors for chunk", i, "/", nrow(chunk_extents)))
    # Subset elevation raster and flowlines data to local chunk
    print("    Cropping raster...")
    elev_chunk <- crop(elev_rast, chunk_extents[i,])
    # Get indices of all empty pixels, so they can be skipped over later on
    empty_pixels <- which(as.logical((as.vector(is.na(elev_chunk)))))
    filled_pixels <- which(!as.logical((as.vector(is.na(elev_chunk)))))
    # Pull the raster into memory and get XY and elevation values in a dataframe
    print("    Building KD tree and searching for neighbors...")
    elev_df <- as.data.frame(elev_chunk, xy=TRUE)
    creek_neighbors <- nn2(data=flowlines_df[,c(1,2)], 
                           query=elev_df[,c(1,2)],
                           k=1,
                           searchtype="standard")
    print("    Getting elevation relative to creek neighbor.")
    rel_elev <- elev_df[,3] - flowlines_df[creek_neighbors[[1]],4]
    rel_elev_with_NA <- rep(NA,ncol(elev_chunk)*nrow(elev_chunk))
    rel_elev_with_NA[filled_pixels] <- rel_elev
    rel_elev_rast <- init(elev_chunk, fun=rel_elev_with_NA)
    terra::writeRaster(rel_elev_rast, paste(raster_filename_prefix, "/rel_elev_", i, ".tif", sep=""), overwrite=TRUE)
    print('    Getting horizontal distance to creek neighbor.')
    horz_dist_with_NA <- rep(NA,ncol(elev_chunk)*nrow(elev_chunk))
    horz_dist_with_NA[filled_pixels] <- creek_neighbors[[2]]
    creek_dist_rast <- init(elev_chunk, fun=horz_dist_with_NA)
    terra::writeRaster(creek_dist_rast, paste(raster_filename_prefix, "/creek_dist_", i, ".tif", sep=""), overwrite=TRUE)
  }
  lapply(1:nrow(chunk_extents), mapRelativeElevationChunk)
  
}


# ********************************************************************************
# San Pedro River: 
# ********************************************************************************

# First, using LiDAR...
# Load in elevation raster
dem_san_pedro <- terra::rast(here::here("LiDAR","terrain_elevation","san_pedro_dem.tif"))
# Load feature containing stream order
flowdata_san_pedro <- st_read(here::here("hydrology","spr_flowlines_with_vaa.shp"))
# Get relative elevation for San Pedro data, keeping only flowlines of at least 2nd order
san_pedro_relative_elevation <- getRelativeElevation(dem_san_pedro, flowdata_san_pedro, here::here("LiDAR","relative_elevation"), 2, "StrmOrd")

# Second, using SRTM...
# Load in elevation raster
dem_san_pedro <- terra::rast(here::here("SRTM","terrain_elevation","san_pedro_srtm.tif"))
# Load feature containing stream order
flowdata_san_pedro <- st_read(here::here("hydrology","spr_flowlines_with_vaa.shp"))
# Get relative elevation for San Pedro data, keeping only flowlines of at least 2nd order
san_pedro_relative_elevation <- getRelativeElevation(dem_san_pedro, flowdata_san_pedro, here::here("SRTM","relative_elevation"), 2, "StrmOrd")


# ********************************************************************************
# Santa Margarita River:  
# ********************************************************************************

# First, using LiDAR...
# Load in elevation raster
dem_pendleton <- terra::rast(here::here("LiDAR","terrain_elevation","pendleton_dem.tif"))
# Load feature containing stream order
flowdata_pendleton <- st_read(here::here("hydrology","smr_flowlines_with_vaa.shp"))
# Get relative elevation for San Pedro data, keeping only flowlines of at least 2nd order
pendleton_relative_elevation <- getRelativeElevation(dem_pendleton, flowdata_pendleton, here::here("LiDAR","relative_elevation"), 2, "StrmOrd")

# Second, using SRTM...
# Load in elevation raster
dem_pendleton <- terra::rast(here::here("SRTM","terrain_elevation","pendleton_srtm.tif"))
# Load feature containing stream order
flowdata_pendleton <- st_read(here::here("hydrology","smr_flowlines_with_vaa.shp"))
# Get relative elevation for San Pedro data, keeping only flowlines of at least 2nd order
pendleton_relative_elevation <- getRelativeElevation(dem_pendleton, flowdata_pendleton, here::here("SRTM","relative_elevation"), 2, "StrmOrd")


# ********************************************************************************
# Santa Ynez River: 
# ********************************************************************************

# Second, using SRTM...
# Load in elevation raster
dem_vandenberg <- terra::rast(here::here("SRTM","terrain_elevation","vandenberg_srtm.tif"))
# Load feature containing stream order
flowdata_vandenberg <- st_read(here::here("hydrology","syr_flowlines_with_vaa.shp"))
# Get relative elevation for San Pedro data, keeping only flowlines of at least 2nd order
vandenberg_relative_elevation <- getRelativeElevation(dem_vandenberg, flowdata_vandenberg, here::here("SRTM","relative_elevation"), 2, "StrmOrd")
