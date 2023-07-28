
library(sf)
library(terra)
library(caret)
library(tidyverse)

# Generic function to load a list of images and stack them into one image
loadImagery <- function(filename_list,
                        band_name_list)
{
  loadFile <- function(i)
  {
    print("")
    print(filename_list[[i]])
    print(band_name_list[[i]])
    new_img <- terra::rast(filename_list[[i]])
    names(new_img) <- band_name_list[[i]]
    return(new_img)
  }
  return(lapply(1:length(filename_list), loadFile))
}


# Function to extract raster data from under a vegetation polygon
getDataAtPolygon <- function(index, polygons, img, class_var)
{
  polygon <- polygons[index,]
  print(paste("Working on the following polygon with area ", st_area(polygon), sep=""))
  print(polygon)
  new_data <- as.data.frame(terra::extract(img, polygon))
  names(new_data) <- c("ID", names(img))
  new_data$class <- as.data.frame(polygon[,class_var])[,1]
  new_data$polygon_index <- index
  
  print("")
  return(new_data)
}

# Given a dataframe of class-labelled data, pull a random sample for training or validation
samplePoints <- function(class_label, class_var, veg_df, num_samples)
{
  set.seed(1)
  class_data <- veg_df[which(veg_df[,class_var] == class_label),]
  # Sample training data
  training_indices <- sample(1:nrow(class_data), num_samples)
  training <- class_data[training_indices,]
  training$type <- "training"
  # Sample validation data
  remaining_data <- class_data[-training_indices,]
  validation_indices <- sample(1:nrow(remaining_data), num_samples)
  validation <- class_data[validation_indices,]
  validation$type <- "validation"
  
  return(rbind(training, validation))
}
  
# Generate training data point samples from polygons of manually delineated vegetation
# This function generates random training / validation data, sampling without replacement from all polygons
# That means that a single polygon can contribute both training and validation data, but each pixel will be only one or the other
# this is deprecated and switched to below approach, which uses assigns each entire polygon to either training or validation
generateTrainingData <- function(veg_polygons, imagery, class_variable, class_list)
{
  print(paste("Input has ", nrow(veg_polygons), " polygons with a total area of ", round(sum(as.numeric(st_area(veg_polygons))),0), " m^2", sep=""))
  veg_polygons <- st_crop(veg_polygons, imagery)
  veg_polygons <- veg_polygons[as.numeric(st_area(veg_polygons)) > 0,] #as.numeric(st_area(veg_polygons)) > 0,]
  print(paste("After crop, ", nrow(veg_polygons), " polygons are left with a total area of ", round(sum(as.numeric(st_area(veg_polygons))),0), " m^2", sep=""))
  
  polygon_indices_with_target_classes <- which(as.data.frame(veg_polygons[,"class"])[,1] %in% class_list)
  veg_polygons <- veg_polygons[polygon_indices_with_target_classes,]
  veg_polygons[,class_variable] <- factor(as.data.frame(veg_polygons[,class_variable])[,1], 
                                          levels=class_list)
  
  print("")
  all_data <- bind_rows(lapply(1:nrow(veg_polygons), getDataAtPolygon, polygons=veg_polygons, img=imagery, class_var=class_variable))
  
  print("Labelled data include this many samples, by class:")
  print(all_data %>% group_by(class) %>% tally())
  min_labelled_points <- min((all_data %>% group_by(class) %>% tally())$n)
  points_per_class <- floor(min_labelled_points/2)
  print(paste("Minimum number of training points for a given class is ", min_labelled_points, ", so we'll use ", points_per_class, " points for each class in the training and test sets.", sep=""))
  
  set.seed(1)
  train_val <- bind_rows(lapply(class_list, FUN=samplePoints, class_var=class_variable, veg_df=all_data, num_samples=points_per_class))
  return(train_val)
}


# Generate Polygon Combinations
# max_points is the maximum number of pixels to include in training or validation datasets
# polygon_sizes is a dataframe with two columns:
#   polygon_index - index of the polygon
#   n - number of pixels in that polygon
generateCombinations <- function(max_points, polygon_sizes)
{
  i <- 0
  training_polygons <- c()
  training_points <- 0
  validation_polygons <- c()
  validation_points <- 0
  set.seed <- 1
  unassigned_polygons <- polygon_sizes$polygon_index
  print(unassigned_polygons)
  while((training_points < max_points) * (length(unassigned_polygons) > 0))
  {
    index_in_remainder <- sample(1:length(unassigned_polygons), 1)
    polygon_index <- unassigned_polygons[index_in_remainder]
    unassigned_polygons <- unassigned_polygons[-index_in_remainder]
    training_polygons <- c(training_polygons, polygon_index)
    training_points <- training_points + polygon_sizes[polygon_sizes$polygon_index == polygon_index,"n"]
    print(paste("Picked a new polygon for training data with index ",polygon_index," and ",polygon_sizes[polygon_sizes$polygon_index == polygon_index,"n"], " points.", sep=""))
    print(paste("New total number of training points for class is ", training_points,sep=""))
  }
  while((validation_points < max_points) * (length(unassigned_polygons) > 0))
  {
    index_in_remainder <- sample(1:length(unassigned_polygons), 1)
    polygon_index <- unassigned_polygons[index_in_remainder]
    unassigned_polygons <- unassigned_polygons[-index_in_remainder]
    validation_polygons <- c(validation_polygons, polygon_index)
    validation_points <- validation_points + polygon_sizes[polygon_sizes$polygon_index == polygon_index,"n"]
    print(paste("Picked a new polygon for validation data with index ",polygon_index," and ",polygon_sizes[polygon_sizes$polygon_index == polygon_index,"n"], " points.", sep=""))
    print(paste("New total number of validation points for class is ", validation_points,sep=""))
  }
  if(length(unassigned_polygons) > 0)
    num_remaining_pixels <- sum(polygon_sizes[polygon_sizes$polygon_index %in% unassigned_polygons,"n"])
  else 
    num_remaining_pixels <- 0
  print(paste("After selecting all polygons up to the number specified, ", length(unassigned_polygons), " polygons with ", num_remaining_pixels, " total pixels are still left over and will be added to the potential validation data.",sep=""))
  validation_polygons <- c(validation_polygons, unassigned_polygons)
  return(list(training_polygons, validation_polygons))
}

# This again samples training and validation data from polygons
# However, each polygon can only contribute training OR validation data
samplePointsByPolygon <- function(class_label, class_var, veg_df, num_samples)
{
  set.seed(1)
  class_data <- veg_df[which(veg_df[,class_var] == class_label),]
  # Choose polygons for training and test data
  polygon_list <- class_data %>% 
    group_by(polygon_index) %>% 
    tally()
  all_polygons_sorted <- generateCombinations(sum(polygon_list$n)/2, polygon_list)
  
  # Sample training data
  training_polygon_data <- class_data %>% 
    filter(polygon_index %in% all_polygons_sorted[[1]])
  training_indices <- sample(1:nrow(training_polygon_data), num_samples)
  training <- training_polygon_data[training_indices,]
  training$type <- "training"
  # Sample validation data
  validation_polygon_data <- class_data %>% 
    filter(polygon_index %in% all_polygons_sorted[[2]])
  if(nrow(validation_polygon_data) < num_samples)
    num_samples <- nrow(validation_polygon_data)
  validation_indices <- sample(1:nrow(validation_polygon_data), num_samples)
  validation <- validation_polygon_data[validation_indices,]
  validation$type <- "validation"
  
  return(rbind(training, validation))
}

# Generate training data point samples from polygons of manually delineated vegetation
# Unlike the above function, this one will take only training OR validation points from each polygon
generateTrainingDataSeparatePolygons <- function(veg_polygons, imagery, class_variable, class_list)
{
  print(paste("Input has ", nrow(veg_polygons), " polygons with a total area of ", round(sum(as.numeric(st_area(veg_polygons))),0), " m^2", sep=""))
  veg_polygons <- st_crop(veg_polygons, imagery)
  veg_polygons <- veg_polygons[as.numeric(st_area(veg_polygons)) > 0,] #as.numeric(st_area(veg_polygons)) > 0,]
  print(paste("After crop, ", nrow(veg_polygons), " polygons are left with a total area of ", round(sum(as.numeric(st_area(veg_polygons))),0), " m^2", sep=""))
  
  polygon_indices_with_target_classes <- which(as.data.frame(veg_polygons[,"class"])[,1] %in% class_list)
  veg_polygons <- veg_polygons[polygon_indices_with_target_classes,]
  veg_polygons[,class_variable] <- factor(as.data.frame(veg_polygons[,class_variable])[,1], 
                                          levels=class_list)
  
  print("")
  all_data <- bind_rows(lapply(1:nrow(veg_polygons), getDataAtPolygon, polygons=veg_polygons, img=imagery, class_var=class_variable))
  
  print("Labelled data include this many samples, by class:")
  print(all_data %>% group_by(class) %>% tally())
  min_labelled_points <- min((all_data %>% group_by(class) %>% tally())$n)
  points_per_class <- floor(min_labelled_points/2)
  print(paste("Minimum number of training points for a given class is ", min_labelled_points, ", so we'll use ", points_per_class, " points for each class in the training and test sets.", sep=""))
  
  
  train_val <- bind_rows(lapply(class_list, FUN=samplePointsByPolygon, class_var=class_variable, veg_df=all_data, num_samples=points_per_class))
  return(train_val)
}



predictRaster <- function(img, classifier, filename=NA)
{
  class_map <- terra::predict(object=img, model=classifier, na.rm=TRUE)
  if(!is.na(filename))
  {
    terra::writeRaster(class_map, filename, overwrite=TRUE)
  }
  return(class_map)
}

# Run a random forest classifier, print out a classified image, confusion matrix, and variable importance
#   data_stack_train and data_stack_valid should be dataframes with 'class' as a categorical variable and 'type' a character with 'training' or 'validation'
#   if data_stack_valid isn't specified, will default to using the same datastack for training and validation (but separated to training/validation by 'type')
#   predict_raster is a boolean - controls whether to predict an entire output map
#   if so, requires an input data raster (raster_to_predict) with bands that have the same names as variables in band_list and the data_stacks
#   raster_filename is the full path where the predicted raster will be written
classifyImage <- function(band_list, data_stack_train, data_stack_valid=data_stack_train, predict_raster=FALSE, raster_to_predict=NA, raster_filename=NA)
{
  rf_formula <- as.formula(paste("as.factor(class) ~ ", paste(band_list, collapse="+"), sep=""))
  
  modFit_rf <- train(rf_formula, 
                     method="rf", data = (data_stack_train %>% filter(type=="training") %>% drop_na()))
  
  # --- Validation  Model (Balanced) ---
  validation_results <- predict(modFit_rf, newdata=(data_stack_valid %>% filter(type=="validation") %>% drop_na()))
  # Confusion Matrix with Validation Data
  conf_mat <- confusionMatrix(validation_results, as.factor((data_stack_valid %>% filter(type=="validation") %>% drop_na())$class))
  # Variable Importance 
  var_imp <- varImp(modFit_rf)
  # Build class image, optionally write to disk
  if(predict_raster)
  {
    class_map <- predictRaster(raster_to_predict, modFit_rf, raster_filename)
  }
  else
    class_map <- NA
  
  # Output results
  print(conf_mat)
  print(var_imp)
  return( list(modFit_rf,
               conf_mat,
               var_imp,
               class_map) )
}
  

# Collapse all variable importance results into one coherent table
buildVarImpTable <- function(input_data)
{
  variable_names <- rownames(input_data[[1]]$importance)
  importances <- ((input_data[[1]])$importance)[,1]
  
  return(data.frame(variable = variable_names,
                    importance = importances,
                    test = input_data[[2]]))
}


# ********************************************************************************
# Vegetation polygon check-in
# ********************************************************************************
# run this code while building polygons to verify that enough training points are collected in each class
veg_polygons <- st_read(here::here("vegetation_map_inputs","huachuca_vegetation.shp"))
raster_resolution = 10
veg_polygons$area <- as.numeric(st_area(st_make_valid(veg_polygons)))
as.data.frame(veg_polygons) %>% 
  group_by(class) %>% 
  summarize(total_area = sum(area),
            mean_area = mean(area),
            num_polygons = n(),
            num_pixels = total_area/(raster_resolution^2))


# Lists of bands to use in different classifiers
band_list_all <- c("blue", "green", "red", "NIR", "SWIR1", "SWIR2", "month_1", "month_2", "month_3", "month_4", "month_5", "month_6", "month_7", "month_8", "month_9", "month_10", "month_11", "month_12", "chm_lidar", "dem_lidar", "slope_lidar", "rel_elev_lidar", "creek_dist_lidar")
band_list_no_lidar <- c("blue", "green", "red", "NIR", "SWIR1", "SWIR2", "month_1", "month_2", "month_3", "month_4", "month_5", "month_6", "month_7", "month_8", "month_9", "month_10", "month_11", "month_12", "dem_srtm", "slope_srtm", "rel_elev_srtm", "creek_dist_srtm")
band_list_spectra <- c("blue", "green", "red", "NIR", "SWIR1", "SWIR2")
band_list_pheno <- c("month_1", "month_2", "month_3", "month_4", "month_5", "month_6", "month_7", "month_8", "month_9", "month_10", "month_11", "month_12")
band_list_lidar <- c("chm_lidar", "dem_lidar", "slope_lidar", "rel_elev_lidar", "creek_dist_lidar")
band_list_srtm <- c("dem_srtm", "slope_srtm", "rel_elev_srtm", "creek_dist_srtm")

# ********************************************************************************
# San Pedro River River
# ********************************************************************************

# san_pedro_dem <- terra::project(san_pedro_imagery[[3]], san_pedro_imagery[[1]])
# san_pedro_slope <- terra::project(san_pedro_imagery[[4]], san_pedro_imagery[[1]])
# san_pedro_chm <- terra::project(san_pedro_imagery[[5]], san_pedro_imagery[[1]])
# spr_rel_elev <- terra::project(san_pedro_imagery[[6]], san_pedro_imagery[[1]])
# spr_creek_dist <- terra::project(san_pedro_imagery[[7]], san_pedro_imagery[[1]])
# writeRaster(san_pedro_dem, here::here("LiDAR","terrain_elevation","san_pedro_dem.tif"), overwrite=TRUE)
# writeRaster(san_pedro_slope, here::here("LiDAR","terrain_slope","san_pedro_slope.tif"), overwrite=TRUE)
# writeRaster(san_pedro_chm, here::here("LiDAR","canopy_height","san_pedro_chm.tif"), overwrite=TRUE)
# writeRaster(spr_rel_elev, here::here("LiDAR","relative_elevation","spr_rel_elev.tif"), overwrite=TRUE)
# writeRaster(spr_creek_dist, here::here("LiDAR","relative_elevation","spr_creek_dist.tif"), overwrite=TRUE)

san_pedro_imagery <- loadImagery(list(here::here("phenology","san_pedro_2021_phenoseries.tif"),
                                      here::here("Spectra","san_pedro_summer_spectra_2021.tif"),
                                      here::here("LiDAR","terrain_elevation","san_pedro_dem.tif"),
                                      here::here("LiDAR","terrain_slope","san_pedro_slope.tif"),
                                      here::here("LiDAR","canopy_height","san_pedro_chm.tif"),
                                      here::here("LiDAR","relative_elevation","spr_rel_elev.tif"),
                                      here::here("LiDAR","relative_elevation","spr_creek_dist.tif"),
                                      here::here("SRTM","terrain_elevation","san_pedro_srtm.tif"),
                                      here::here("SRTM","terrain_slope","san_pedro_srtm_slope.tif"),
                                      here::here("SRTM","relative_elevation","spr_rel_elev.tif"),
                                      here::here("SRTM","relative_elevation","spr_creek_dist.tif")),
                                 list(paste("month_",1:12,sep=""),
                                      c("blue","green","red","NIR","SWIR1","SWIR2","DOY"),
                                      "dem_lidar",
                                      "slope_lidar",
                                      "chm_lidar",
                                      "rel_elev_lidar",
                                      "creek_dist_lidar",
                                      "dem_srtm",
                                      "slope_srtm",
                                      "rel_elev_srtm",
                                      "creek_dist_srtm"))
san_pedro_imagery <- terra::rast(san_pedro_imagery)
san_pedro_polygons <- st_read(here::here("vegetation_map_inputs","huachuca_vegetation.shp"))
san_pedro_polygons[san_pedro_polygons$class %in% c("sycamore", "cottonwood"),]$class <- "riparian woodland"
san_pedro_polygons <- st_make_valid(san_pedro_polygons)
san_pedro_training <- generateTrainingDataSeparatePolygons(san_pedro_polygons, san_pedro_imagery, "class", c("riparian woodland","pine_oak","mesquite","shrubland","grassland","sparse","agriculture","urban"))
# Generate all classification results 
san_pedro_results_all <- classifyImage(band_list=band_list_all, data_stack_train=san_pedro_training, 
                                       predict_raster=TRUE, raster_to_predict=san_pedro_imagery, 
                                       raster_filename=here::here("class_outputs","san_pedro_all_data.tif"))
san_pedro_results_no_lidar <- classifyImage(band_list=band_list_no_lidar, data_stack_train=san_pedro_training, 
                                            predict_raster=TRUE, raster_to_predict=san_pedro_imagery, 
                                            raster_filename=here::here("class_outputs","san_pedro_no_lidar.tif"))
san_pedro_results_spectra <- classifyImage(band_list=band_list_spectra, data_stack_train=san_pedro_training, 
                                           predict_raster=TRUE, raster_to_predict=san_pedro_imagery, 
                                           raster_filename=here::here("class_outputs","san_pedro_spectra_only.tif"))
san_pedro_results_pheno <- classifyImage(band_list=band_list_pheno, data_stack_train=san_pedro_training, 
                                         predict_raster=TRUE, raster_to_predict=san_pedro_imagery, 
                                         raster_filename=here::here("class_outputs","san_pedro_phenology_only.tif"))
san_pedro_results_lidar <- classifyImage(band_list=band_list_lidar, data_stack_train=san_pedro_training, 
                                         predict_raster=TRUE, raster_to_predict=san_pedro_imagery, 
                                         raster_filename=here::here("class_outputs","san_pedro_lidar_only.tif"))
san_pedro_results_srtm <- classifyImage(band_list=band_list_srtm, data_stack_train=san_pedro_training, 
                                        predict_raster=TRUE, raster_to_predict=san_pedro_imagery, 
                                        raster_filename=here::here("class_outputs","san_pedro_srtm_only.tif"))
san_pedro_results_list <- list(san_pedro_results_all[1:3],
                               san_pedro_results_no_lidar[1:3],
                               san_pedro_results_spectra[1:3],
                               san_pedro_results_pheno[1:3],
                               san_pedro_results_lidar[1:3],
                               san_pedro_results_srtm[1:3])
save(san_pedro_results_list, file=here::here("random_forest_results","san_pedro.rdata"))

# Variable Importance Outputs
san_pedro_varimp <- bind_rows(lapply(list(list(san_pedro_results_all[[3]], "all"),
                                          list(san_pedro_results_no_lidar[[3]], "no_lidar"),
                                          list(san_pedro_results_spectra[[3]], "spectra"),
                                          list(san_pedro_results_pheno[[3]], "pheno"),
                                          list(san_pedro_results_lidar[[3]], "lidar"),
                                          list(san_pedro_results_srtm[[3]], "srtm")),
                                     buildVarImpTable)) %>% 
  pivot_wider(names_from=test, 
              values_from=importance)
View(san_pedro_varimp)




# ********************************************************************************
# Santa Margarita River
# ********************************************************************************

pendleton_imagery <- loadImagery(list(here::here("phenology","pendleton_2019_phenoseries.tif"),
                                      here::here("Spectra","pendleton_summer_spectra_2019.tif"),
                                      here::here("LiDAR","terrain_elevation","pendleton_dem.tif"),
                                      here::here("LiDAR","terrain_slope","pendleton_slope.tif"),
                                      here::here("LiDAR","canopy_height","pendleton_chm.tif"),
                                      here::here("LiDAR","relative_elevation","smr_rel_elev.tif"),
                                      here::here("LiDAR","relative_elevation","smr_creek_dist.tif"),
                                      here::here("SRTM","terrain_elevation","pendleton_srtm.tif"),
                                      here::here("SRTM","terrain_slope","pendleton_srtm_slope.tif"),
                                      here::here("SRTM","relative_elevation","smr_rel_elev.tif"),
                                      here::here("SRTM","relative_elevation","smr_creek_dist.tif")),
                                 list(paste("month_",1:12,sep=""),
                                      c("blue","green","red","NIR","SWIR1","SWIR2","DOY"),
                                      "dem_lidar",
                                      "slope_lidar",
                                      "chm_lidar",
                                      "rel_elev_lidar",
                                      "creek_dist_lidar",
                                      "dem_srtm",
                                      "slope_srtm",
                                      "rel_elev_srtm",
                                      "creek_dist_srtm"))
pendleton_imagery <- terra::rast(pendleton_imagery)
pendleton_polygons <- st_read(here::here("vegetation_map_inputs","pendleton_vegetation.shp"))
pendleton_training <- generateTrainingDataSeparatePolygons(pendleton_polygons, pendleton_imagery, "class", c("riparian_woodland","riparian_scrub","wetland","oak","chaparral","sage_scrub","grassland","urban","water","non_vegetated"))
# NOTE - removed eucalyptus from the class list, because it's rare on this landscape and difficult to generate enough training data for
# Generate all classification results 
pendleton_results_all <- classifyImage(band_list=band_list_all, data_stack_train=pendleton_training, 
                                       predict_raster=TRUE, raster_to_predict=pendleton_imagery, 
                                       raster_filename=here::here("class_outputs","pendleton_all_data.tif"))
pendleton_results_no_lidar <- classifyImage(band_list=band_list_no_lidar, data_stack_train=pendleton_training, 
                                            predict_raster=TRUE, raster_to_predict=pendleton_imagery, 
                                            raster_filename=here::here("class_outputs","pendleton_no_lidar.tif"))
pendleton_results_spectra <- classifyImage(band_list=band_list_spectra, data_stack_train=pendleton_training, 
                                           predict_raster=TRUE, raster_to_predict=pendleton_imagery, 
                                           raster_filename=here::here("class_outputs","pendleton_spectra_only.tif"))
pendleton_results_pheno <- classifyImage(band_list=band_list_pheno, data_stack_train=pendleton_training, 
                                         predict_raster=TRUE, raster_to_predict=pendleton_imagery, 
                                         raster_filename=here::here("class_outputs","pendleton_phenology_only.tif"))
pendleton_results_lidar <- classifyImage(band_list=band_list_lidar, data_stack_train=pendleton_training, 
                                         predict_raster=TRUE, raster_to_predict=pendleton_imagery, 
                                         raster_filename=here::here("class_outputs","pendleton_lidar_only.tif"))
pendleton_results_srtm <- classifyImage(band_list=band_list_srtm, data_stack_train=pendleton_training, 
                                        predict_raster=TRUE, raster_to_predict=pendleton_imagery, 
                                        raster_filename=here::here("class_outputs","pendleton_srtm_only.tif"))
pendleton_results_list <- list(pendleton_results_all[1:3],
                              pendleton_results_no_lidar[1:3],
                              pendleton_results_spectra[1:3],
                              pendleton_results_pheno[1:3],
                              pendleton_results_lidar[1:3],
                              pendleton_results_srtm[1:3])
save(pendleton_results_list, file=here::here("random_forest_results","Pendleton.rdata"))

pendleton_all_data_raster_masked <- terra::mask(x=pendleton_results_all[[4]], mask=(pendleton_imagery[[22]]==0), maskvalue=1)
terra::writeRaster(pendleton_all_data_raster_masked, here::here("class_outputs","pendleton_all_data_masked.tif"))

# Variable Importance Outputs
pendleton_varimp <- bind_rows(lapply(list(list(pendleton_results_all[[3]], "all"),
                                     list(pendleton_results_no_lidar[[3]], "no_lidar"),
                                     list(pendleton_results_spectra[[3]], "spectra"),
                                     list(pendleton_results_pheno[[3]], "pheno"),
                                     list(pendleton_results_lidar[[3]], "lidar"),
                                     list(pendleton_results_srtm[[3]], "srtm")),
                                buildVarImpTable)) %>% 
   pivot_wider(names_from=test, 
               values_from=importance)
View(pendleton_varimp)


# ********************************************************************************
# Santa Ynez River
# ********************************************************************************

vandenberg_imagery_2013 <- loadImagery(list(here::here("phenology","vandenberg_2013_phenoseries.tif"),
                                            here::here("Spectra","vandenberg_summer_spectra_2013.tif"),
                                            here::here("SRTM","terrain_elevation","vandenberg_srtm.tif"),
                                            here::here("SRTM","terrain_slope","vandenberg_srtm_slope.tif"),
                                            here::here("SRTM","relative_elevation","syr_rel_elev.tif"),
                                            here::here("SRTM","relative_elevation","syr_creek_dist.tif")),
                                      list(paste("month_",1:12,sep=""),
                                           c("blue","green","red","NIR","SWIR1","SWIR2","DOY"),
                                           "dem_srtm",
                                           "slope_srtm",
                                           "rel_elev_srtm",
                                           "creek_dist_srtm"))
vandenberg_imagery_2013 <- terra::rast(vandenberg_imagery_2013)
vandenberg_imagery_2015 <- loadImagery(list(here::here("phenology","vandenberg_2015_phenoseries.tif"),
                                            here::here("Spectra","vandenberg_summer_spectra_2015.tif"),
                                            here::here("SRTM","terrain_elevation","vandenberg_srtm.tif"),
                                            here::here("SRTM","terrain_slope","vandenberg_srtm_slope.tif"),
                                            here::here("SRTM","relative_elevation","syr_rel_elev.tif"),
                                            here::here("SRTM","relative_elevation","syr_creek_dist.tif")),
                                       list(paste("month_",1:12,sep=""),
                                            c("blue","green","red","NIR","SWIR1","SWIR2","DOY"),
                                            "dem_srtm",
                                            "slope_srtm",
                                            "rel_elev_srtm",
                                            "creek_dist_srtm"))
vandenberg_imagery_2015 <- terra::rast(vandenberg_imagery_2015)
vandenberg_imagery_2017 <- loadImagery(list(here::here("phenology","vandenberg_2017_phenoseries.tif"),
                                            here::here("Spectra","vandenberg_summer_spectra_2017.tif"),
                                            here::here("SRTM","terrain_elevation","vandenberg_srtm.tif"),
                                            here::here("SRTM","terrain_slope","vandenberg_srtm_slope.tif"),
                                            here::here("SRTM","relative_elevation","syr_rel_elev.tif"),
                                            here::here("SRTM","relative_elevation","syr_creek_dist.tif")),
                                       list(paste("month_",1:12,sep=""),
                                            c("blue","green","red","NIR","SWIR1","SWIR2","DOY"),
                                            "dem_srtm",
                                            "slope_srtm",
                                            "rel_elev_srtm",
                                            "creek_dist_srtm"))
vandenberg_imagery_2017 <- terra::rast(vandenberg_imagery_2017)
vandenberg_imagery_2019 <- loadImagery(list(here::here("phenology","vandenberg_2019_phenoseries.tif"),
                                            here::here("Spectra","vandenberg_summer_spectra_2019.tif"),
                                            here::here("SRTM","terrain_elevation","vandenberg_srtm.tif"),
                                            here::here("SRTM","terrain_slope","vandenberg_srtm_slope.tif"),
                                            here::here("SRTM","relative_elevation","syr_rel_elev.tif"),
                                            here::here("SRTM","relative_elevation","syr_creek_dist.tif")),
                                       list(paste("month_",1:12,sep=""),
                                            c("blue","green","red","NIR","SWIR1","SWIR2","DOY"),
                                            "dem_srtm",
                                            "slope_srtm",
                                            "rel_elev_srtm",
                                            "creek_dist_srtm"))
vandenberg_imagery_2019 <- terra::rast(vandenberg_imagery_2019)
# Generate training / validation data (common across all images by fixing seed)
vandenberg_polygons <- st_read(here::here("vegetation_map_inputs","vandenberg_vegetation.shp"))
set.seed(1)
vandenberg_training_2013 <- generateTrainingDataSeparatePolygons(vandenberg_polygons, vandenberg_imagery_2013, "class", c("riparian","oak","chaparral","sage_scrub","eucalyptus","annual","agriculture","urban","bare","water"))
set.seed(1)
vandenberg_training_2015 <- generateTrainingDataSeparatePolygons(vandenberg_polygons, vandenberg_imagery_2015, "class", c("riparian","oak","chaparral","sage_scrub","eucalyptus","annual","agriculture","urban","bare","water"))
set.seed(1)
vandenberg_training_2017 <- generateTrainingDataSeparatePolygons(vandenberg_polygons, vandenberg_imagery_2017, "class", c("riparian","oak","chaparral","sage_scrub","eucalyptus","annual","agriculture","urban","bare","water"))
set.seed(1)
vandenberg_training_2019 <- generateTrainingDataSeparatePolygons(vandenberg_polygons, vandenberg_imagery_2019, "class", c("riparian","oak","chaparral","sage_scrub","eucalyptus","annual","agriculture","urban","bare","water"))
# Generate all classification results - first, building a classifier in each year and testing it in that year
vandenberg_results_2013 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2013, data_stack_valid=vandenberg_training_2013,
                                         predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2013, 
                                         raster_filename=here::here("class_outputs","vandenberg_2013_test_2013_train.tif"))
vandenberg_results_2015 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2015, data_stack_valid=vandenberg_training_2015,
                                         predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2015, 
                                         raster_filename=here::here("class_outputs","vandenberg_2015_test_2015_train.tif"))
vandenberg_results_2017 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2017, data_stack_valid=vandenberg_training_2017,
                                         predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2017, 
                                         raster_filename=here::here("class_outputs","vandenberg_2017_test_2017_train.tif"))
vandenberg_results_2019 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2019, data_stack_valid=vandenberg_training_2019,
                                         predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2019, 
                                         raster_filename=here::here("class_outputs","vandenberg_2019_test_2019_train.tif"))
# Now, do the same thing, but train a classifier in a dry year and test it in a wet year. 
# Train in 2015, test in 2013
vandenberg_results_2013_2015 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2015, data_stack_valid=vandenberg_training_2013,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2013, 
                                              raster_filename=here::here("class_outputs","vandenberg_2013_test_2015_train.tif"))
# Train in 2017, test in 2013
vandenberg_results_2013_2017 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2017, data_stack_valid=vandenberg_training_2013,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2013, 
                                              raster_filename=here::here("class_outputs","vandenberg_2013_test_2017_train.tif"))
# Train in 2019, test in 2013
vandenberg_results_2013_2019 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2019, data_stack_valid=vandenberg_training_2013,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2013, 
                                              raster_filename=here::here("class_outputs","vandenberg_2013_test_2019_train.tif"))
 
# Train in 2013, test in 2015
vandenberg_results_2015_2013 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2013, data_stack_valid=vandenberg_training_2015,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2015, 
                                              raster_filename=here::here("class_outputs","vandenberg_2015_test_2013_train.tif"))
# Train in 2017, test in 2015
vandenberg_results_2015_2017 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2017, data_stack_valid=vandenberg_training_2015,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2015, 
                                              raster_filename=here::here("class_outputs","vandenberg_2015_test_2017_train.tif"))
# Train in 2019, test in 2015
vandenberg_results_2015_2019 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2019, data_stack_valid=vandenberg_training_2015,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2015, 
                                              raster_filename=here::here("class_outputs","vandenberg_2015_test_2019_train.tif"))

# Train in 2013, test in 2017
vandenberg_results_2017_2013 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2013, data_stack_valid=vandenberg_training_2017,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2017, 
                                              raster_filename=here::here("class_outputs","vandenberg_2017_test_2013_train.tif"))
# Train in 2015, test in 2017
vandenberg_results_2017_2015 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2015, data_stack_valid=vandenberg_training_2017,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2017, 
                                              raster_filename=here::here("class_outputs","vandenberg_2017_test_2015_train.tif"))
# Train in 2019, test in 2017
vandenberg_results_2017_2019 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2019, data_stack_valid=vandenberg_training_2017,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2017, 
                                              raster_filename=here::here("class_outputs","vandenberg_2017_test_2019_train.tif"))

# Train in 2013, test in 2019
vandenberg_results_2019_2013 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2013, data_stack_valid=vandenberg_training_2019,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2019, 
                                              raster_filename=here::here("class_outputs","vandenberg_2019_test_2013_train.tif"))
# Train in 2015, test in 2019
vandenberg_results_2019_2015 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2015, data_stack_valid=vandenberg_training_2019,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2019, 
                                              raster_filename=here::here("class_outputs","vandenberg_2019_test_2015_train.tif"))
# Train in 2017, test in 2019
vandenberg_results_2019_2017 <- classifyImage(band_list_no_lidar, data_stack_train=vandenberg_training_2017, data_stack_valid=vandenberg_training_2019,
                                              predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2019, 
                                              raster_filename=here::here("class_outputs","vandenberg_2019_test_2017_train.tif"))

# Finally, we'll train across two years to see how much that improves performance
# Train in 2013 and 2019, test in 2013 and 2019
vandenberg_results_2013_2019_2013_2019 <- classifyImage(band_list_no_lidar, data_stack_train=rbind(vandenberg_training_2013, vandenberg_training_2019), 
                                                        data_stack_valid=rbind(vandenberg_training_2013, vandenberg_training_2019),
                                                        predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2013,
                                                        raster_filename=here::here("class_outputs","vandenberg_2013_test_2013_2019_train.tif"))
vandenberg_2019_from_2013_2019_rast <- predictRaster(vandenberg_imagery_2019, vandenberg_results_2013_2019_2013_2019[[1]], 
                                                     here::here("class_outputs","vandenberg_2019_test_2013_2019_train.tif"))
# Train in 2013 and 2019, test in 2015 and 2017
vandenberg_results_2015_2017_2013_2019 <- classifyImage(band_list_no_lidar, data_stack_train=rbind(vandenberg_training_2013, vandenberg_training_2019), 
                                                         data_stack_valid=rbind(vandenberg_training_2015, vandenberg_training_2017),
                                                         predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2015, 
                                                         raster_filename=here::here("class_outputs","vandenberg_2015_test_2013_2019_train.tif"))
vandenberg_2017_from_2013_2019_rast <- predictRaster(vandenberg_imagery_2017, vandenberg_results_2015_2017_2013_2019[[1]], 
              here::here("class_outputs","vandenberg_2017_test_2013_2019_train.tif"))
# Train in 2015 and 2017, test in 2015 and 2017
vandenberg_results_2015_2017_2015_2017 <- classifyImage(band_list_no_lidar, data_stack_train=rbind(vandenberg_training_2015, vandenberg_training_2017), 
                                                        data_stack_valid=rbind(vandenberg_training_2015, vandenberg_training_2017),
                                                        predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2015,
                                                        raster_filename=here::here("class_outputs","vandenberg_2015_test_2015_2017_train.tif"))
vandenberg_2017_from_2015_2017_rast <- predictRaster(vandenberg_imagery_2017, vandenberg_results_2015_2017_2015_2017[[1]], 
                                                     here::here("class_outputs","vandenberg_2017_test_2015_2017_train.tif"))
# Train in 2015 and 2017, test in 2013 and 2019
vandenberg_results_2013_2019_2015_2017 <- classifyImage(band_list_no_lidar, data_stack_train=rbind(vandenberg_training_2015, vandenberg_training_2017), 
                                                         data_stack_valid=rbind(vandenberg_training_2013, vandenberg_training_2019),
                                                         predict_raster=TRUE, raster_to_predict=vandenberg_imagery_2013, 
                                                         raster_filename=here::here("class_outputs","vandenberg_2013_test_2015_2017_train.tif"))
vandenberg_2019_from_2015_2017_rast <- predictRaster(vandenberg_imagery_2019, vandenberg_results_2013_2019_2015_2017[[1]], 
              here::here("class_outputs","vandenberg_2019_test_2015_2017_train.tif"))


# Variable Importance Outputs
vandenberg_varimp <- bind_rows(lapply(list(list(vandenberg_results_2013[[3]], "2013"),
                                           list(vandenberg_results_2015[[3]], "2015"),
                                           list(vandenberg_results_2017[[3]], "2017"),
                                           list(vandenberg_results_2019[[3]], "2019")),
                                     buildVarImpTable)) %>% 
  pivot_wider(names_from=test, 
              values_from=importance)
View(vandenberg_varimp)

vandenberg_results_list <- list(vandenberg_results_2013,
                                vandenberg_results_2015,
                                vandenberg_results_2017,
                                vandenberg_results_2019,
                                vandenberg_results_2013_2015,
                                vandenberg_results_2013_2017,
                                vandenberg_results_2013_2019,
                                vandenberg_results_2015_2013,
                                vandenberg_results_2015_2017,
                                vandenberg_results_2015_2019,
                                vandenberg_results_2017_2013,
                                vandenberg_results_2017_2015,
                                vandenberg_results_2017_2019,
                                vandenberg_results_2019_2013,
                                vandenberg_results_2019_2015,
                                vandenberg_results_2019_2017,
                                vandenberg_results_2013_2019_2013_2019,
                                vandenberg_results_2015_2017_2013_2019,
                                vandenberg_results_2015_2017_2015_2017,
                                vandenberg_results_2013_2019_2015_2017)

save(vandenberg_results_list, here::here("random_forest_results","vandenberg.rdata"))






# Some visualizations of class differences in training data
pendleton_training_long <- pendleton_training %>% 
  drop_na() %>%
  pivot_longer(2:29, names_to="variable", values_to="value")
pendleton_phenology <- ggplot(pendleton_training_long %>% filter(variable %in% paste("month_", 1:12, sep=""))) + 
  facet_wrap(~class, nrow=1) + 
  geom_density_2d_filled(aes(x=as.numeric(substr(variable,7,10)), y=value), 
                         contour_var = "ndensity") + 
  geom_line(data=pendleton_training_long %>% 
              group_by(class,variable) %>% 
              summarize(median_NDVI = median(value)), 
            aes(x=as.numeric(substr(variable,7,10)),y=median_NDVI), 
            col="red") + 
  scale_y_continuous(limits=c(-0.2,1)) + 
  scale_x_continuous(breaks=c(4,8,12)) + 
  xlab("Month") + 
  ylab("NDVI") 
pendleton_phenology
ggsave(here::here("output_plots","pendleton_phenology.png"), pendleton_phenology, 
       width=14, height=3)

spectrum_bands <- c(mean(0.452,0.512),
                    mean(0.533,0.590),
                    mean(0.636,0.673),
                    mean(0.851,0.879),
                    mean(1.566,1.651),
                    mean(2.107,2.294))
names(spectrum_bands) <- c("blue","green","red","NIR","SWIR1","SWIR2")
pendleton_spectra <- pendleton_training_long %>% 
  filter(variable %in% names(spectrum_bands))
which(pendleton_spectra$variable == spectrum_bands)
pendleton_spectra$wavelength <- spectrum_bands[pendleton_spectra$variable]
pendleton_spectrum <- ggplot() + 
  facet_wrap(~class, nrow=1) + 
  geom_boxplot(data=pendleton_spectra,
               aes(x=wavelength, y=value/10000, group=wavelength)) + 
  geom_line(data=pendleton_spectra %>% 
              group_by(class,wavelength) %>% 
              summarize(median_reflectance = median(value)), 
            aes(x=wavelength,y=median_reflectance/10000), 
            col="red") + 
  xlab("Wavelength") + 
  ylab("Reflectance") +
  theme_bw()
pendleton_spectrum
ggsave(here::here("output_plots","pendleton_spectra.png"), pendleton_spectrum, 
       width=14, height=3)





# Some visualizations of class differences in training data
san_pedro_training_long <- san_pedro_training %>% 
  drop_na() %>%
  pivot_longer(2:29, names_to="variable", values_to="value")
san_pedro_phenology <- ggplot(san_pedro_training_long %>% filter(variable %in% paste("month_", 1:12, sep=""))) + 
  facet_wrap(~class, nrow=1) + 
  geom_density_2d_filled(aes(x=as.numeric(substr(variable,7,10)), y=value), 
                         contour_var = "ndensity") + 
  geom_line(data=san_pedro_training_long %>% 
              group_by(class,variable) %>% 
              summarize(median_NDVI = median(value)), 
            aes(x=as.numeric(substr(variable,7,10)),y=median_NDVI), 
            col="red") + 
  scale_y_continuous(limits=c(-0.2,1)) + 
  scale_x_continuous(breaks=c(4,8,12)) + 
  xlab("Month") + 
  ylab("NDVI")
san_pedro_phenology
ggsave(here::here("output_plots","san_pedro_phenology.png"), san_pedro_phenology, 
       width=12, height=3)

spectrum_bands <- c(mean(0.452,0.512),
                    mean(0.533,0.590),
                    mean(0.636,0.673),
                    mean(0.851,0.879),
                    mean(1.566,1.651),
                    mean(2.107,2.294))
names(spectrum_bands) <- c("blue","green","red","NIR","SWIR1","SWIR2")
san_pedro_spectra <- san_pedro_training_long %>% 
  filter(variable %in% names(spectrum_bands))
which(san_pedro_spectra$variable == spectrum_bands)
san_pedro_spectra$wavelength <- spectrum_bands[san_pedro_spectra$variable]
san_pedro_spectrum <- ggplot() + 
  facet_wrap(~class, nrow=1) + 
  geom_boxplot(data=san_pedro_spectra,
               aes(x=wavelength, y=value/10000, group=wavelength)) + 
  geom_line(data=san_pedro_spectra %>% 
              group_by(class,wavelength) %>% 
              summarize(median_reflectance = median(value)), 
            aes(x=wavelength,y=median_reflectance/10000), 
            col="red") + 
  xlab("Wavelength") + 
  ylab("Reflectance") +
  theme_bw()
san_pedro_spectrum
ggsave(here::here("output_plots","san_pedro_spectra.png"), san_pedro_spectrum, 
       width=12, height=3)




# Some visualizations of class differences in training data
# First, load the extra Vandeberg imagery (only generated later after first pass of classifier)
vandenberg_

vand_nv_classes <- c("riparian","oak","chaparral","sage_scrub","eucalyptus","annual")

vandenberg_pheno_long_2013 <- vandenberg_training_2013 %>% 
  pivot_longer(2:13, names_to="month_str", values_to="ndvi") %>%
  mutate(month = as.numeric(substr(month_str, 7, 10)),
         year = 2013,
         year_col = "red")
vandenberg_pheno_long_2015 <- vandenberg_training_2015 %>% 
  pivot_longer(2:13, names_to="month_str", values_to="ndvi") %>%
  mutate(month = as.numeric(substr(month_str, 7, 10)),
         year = 2015, 
         year_col = "orange")
vandenberg_pheno_long_2017 <- vandenberg_training_2017 %>% 
  pivot_longer(2:13, names_to="month_str", values_to="ndvi") %>%
  mutate(month = as.numeric(substr(month_str, 7, 10)),
         year = 2017,
         year_col = "green")
vandenberg_pheno_long_2019 <- vandenberg_training_2019 %>% 
  pivot_longer(2:13, names_to="month_str", values_to="ndvi") %>%
  mutate(month = as.numeric(substr(month_str, 7, 10)),
         year = 2019,
         year_col = "green3")

vandenberg_pheno_all <- rbind(vandenberg_pheno_long_2013,
                              vandenberg_pheno_long_2015,
                              vandenberg_pheno_long_2017,
                              vandenberg_pheno_long_2019)
vandenberg_pheno_summary <- vandenberg_pheno_all %>%
  group_by(class, month, year) %>%
  summarize(median_ndvi = median(ndvi, na.rm=TRUE),
            mean_ndvi = mean(ndvi, na.rm=TRUE),
            sd = sd(ndvi, na.rm=TRUE),
            se = sd/sqrt(n()))
vandenberg_phenology <- ggplot() + 
  facet_wrap(~class, nrow=1) + 
  geom_errorbar(data=vandenberg_pheno_summary %>% 
                  filter(class %in% vand_nv_classes), 
                aes(x=month, y=mean_ndvi, 
                    ymin=mean_ndvi-se, ymax=mean_ndvi+se,
                    group=year, col=as.character(year))) + 
  geom_line(data=vandenberg_pheno_summary %>% 
              filter(class %in% vand_nv_classes), 
            aes(x=month, y=mean_ndvi, 
                group=year, col=as.character(year)), linetype="dashed") + 
  scale_y_continuous(limits=c(0.2,0.85)) + 
  scale_x_continuous(breaks=c(4,8,12)) + 
  scale_color_manual(values=c("2013"="orange","2015"="red","2017"="green3","2019"="forestgreen")) + 
  xlab("Month") + 
  ylab("NDVI") +
  theme_bw()
vandenberg_phenology
ggsave(here::here("output_plots","vandenberg_phenology.png"), vandenberg_phenology, 
       width=10, height=3)

# Pixel-by-pixel comparison of NDVI in each month, by class
vandenberg_pheno_difference <- vandenberg_pheno_long_2013[,c(1:16,18)]
vandenberg_pheno_difference$ndvi_2013 <- vandenberg_pheno_long_2013$ndvi
vandenberg_pheno_difference$ndvi_2015 <- vandenberg_pheno_long_2015$ndvi
vandenberg_pheno_difference$ndvi_2017 <- vandenberg_pheno_long_2017$ndvi
vandenberg_pheno_difference$ndvi_2019 <- vandenberg_pheno_long_2019$ndvi
vandenberg_pheno_difference$avg_wet <- (vandenberg_pheno_difference$ndvi_2017+vandenberg_pheno_difference$ndvi_2019)/2
vandenberg_pheno_difference$avg_dry <- (vandenberg_pheno_difference$ndvi_2013+vandenberg_pheno_difference$ndvi_2015)/2

vand_pheno_ttest <- vandenberg_pheno_difference %>%
  group_by(class, month) %>%
  summarize(diff_dry_dry = mean(ndvi_2013-ndvi_2015),
            diff_wet_dry = mean(avg_wet-avg_dry),
            diff_wet_wet = mean(ndvi_2017-ndvi_2019),
            diff_diff_water = mean(c(abs(ndvi_2013-ndvi_2017),abs(ndvi_2013-ndvi_2019),abs(ndvi_2015-ndvi_2017),abs(ndvi_2015-ndvi_2019))),
            diff_same_water = mean(c(abs(ndvi_2013-ndvi_2015),abs(ndvi_2017-ndvi_2019))),
            diff_wet_dry_sd = sd(avg_wet-avg_dry),
            diff_wet_dry_se = diff_wet_dry_sd/n(),
            t_dry_dry = t.test(ndvi_2013, ndvi_2015)$p.value,
            t_wet_dry = t.test(avg_wet, avg_dry)$p.value,
            t_wet_wet = t.test(ndvi_2017, ndvi_2019)$p.value)

# Ignoring this for now... ABS of pairwise comparisons between similar or different water years
# For most classes, summer differences are greater for dissimilar water years than for similar water years

# ggplot(vand_pheno_ttest %>% filter(class %in% vand_nv_classes)) +
#   facet_wrap(~class, nrow=1) + 
#   geom_line(aes(x=month, y=diff_diff_water),
#             col="red") +  
#   geom_line(aes(x=month, y=diff_same_water),
#             col="blue") +
#   theme_bw()

# Comparison of greenness in wet vs. dry years at each pixel, by class
#  note that 2013 and 2015 were dry, 2017 and 2019 were wet
vand_wet_dry_diff <- ggplot(vand_pheno_ttest %>% filter(class %in% vand_nv_classes)) +
  facet_wrap(~class, nrow=1) + 
  geom_line(aes(x=month, y=diff_wet_dry),
            col="red") +
  geom_hline(yintercept=0, col="black", linetype="dashed") + 
  scale_x_continuous(breaks=c(4,8,12)) + 
  theme_bw() + 
  xlab("Month") + 
  ylab("δNDVI (Wet vs. Dry Years)")
ggsave(here::here("output_plots","vand_pheno_wet_dry_diff.png"), vand_wet_dry_diff,
       width=10, height=3)


# Bunch of complicated figures trying to further unpack annual differences
#  This became to complicated and I decided to leave this for later work 
# vand_pheno_r2 <- (vandenberg_pheno_difference %>% 
#                     group_by(class, month) %>%
#                     summarize(r_sqd_2013_2015 = summary(lm(ndvi_2013 ~ ndvi_2015))$r.squared,
#                               r_sqd_2013_2017 = summary(lm(ndvi_2013 ~ ndvi_2017))$r.squared,
#                               r_sqd_2013_2019 = summary(lm(ndvi_2013 ~ ndvi_2019))$r.squared,
#                               r_sqd_2015_2017 = summary(lm(ndvi_2015 ~ ndvi_2017))$r.squared,
#                               r_sqd_2015_2019 = summary(lm(ndvi_2015 ~ ndvi_2019))$r.squared,
#                               r_sqd_2017_2019 = summary(lm(ndvi_2017 ~ ndvi_2019))$r.squared))
# rmse <- function(vec1, vec2)
# {
#   sqrt(mean((vec1-vec2)^2))
# }
# vand_pheno_rmse <- (vandenberg_pheno_difference %>% 
#                       group_by(class, month) %>%
#                       summarize(rmse_2013_2015 = rmse(ndvi_2013, ndvi_2015),
#                                 rmse_2013_2017 = rmse(ndvi_2013, ndvi_2017),
#                                 rmse_2013_2019 = rmse(ndvi_2013, ndvi_2019),
#                                 rmse_2015_2017 = rmse(ndvi_2015, ndvi_2017),
#                                 rmse_2015_2019 = rmse(ndvi_2015, ndvi_2019),
#                                 rmse_2017_2019 = rmse(ndvi_2017, ndvi_2019)))
# # Plot RMSE across points for each class
# ggplot(vand_pheno_rmse %>% 
#          filter(class %in% vand_nv_classes)) + 
#   facet_wrap(~class) + 
#   geom_line(aes(x=month,
#                 y=rmse_2013_2015),
#             col="green1") + 
#   geom_line(aes(x=month,
#                 y=rmse_2013_2017),
#             col="maroon") + 
#   geom_line(aes(x=month,
#                 y=rmse_2013_2019),
#             col="red") + 
#   geom_line(aes(x=month,
#                 y=rmse_2015_2017),
#             col="orange") + 
#   geom_line(aes(x=month,
#                 y=rmse_2015_2019),
#             col="yellow") + 
#   geom_line(aes(x=month,
#                 y=rmse_2017_2019),
#             col="forestgreen") + 
#   theme_bw()
# # Plot R^2 across years for each class
# ggplot(vand_pheno_r2 %>% 
#          filter(class %in% vand_nv_classes)) + 
#   facet_wrap(~class) + 
#   geom_line(aes(x=month,
#                 y=r_sqd_2013_2015),
#             col="green1") + 
#   geom_line(aes(x=month,
#                 y=r_sqd_2013_2017),
#             col="maroon") + 
#   geom_line(aes(x=month,
#                 y=r_sqd_2013_2019),
#             col="red") + 
#   geom_line(aes(x=month,
#                 y=r_sqd_2015_2017),
#             col="orange") + 
#   geom_line(aes(x=month,
#                 y=r_sqd_2015_2019),
#             col="yellow") + 
#   geom_line(aes(x=month,
#                 y=r_sqd_2017_2019),
#             col="forestgreen") + 
#   theme_bw()
# 
# 
# ggplot(vand_pheno_rmse %>% 
#          filter(class %in% vand_nv_classes)) + 
#   facet_wrap(~class, nrow=1) + 
#   geom_line(aes(x=month,
#                 y=rmse_2013_2015),
#             col="blue") + 
#   geom_line(aes(x=month,
#                 y=rmse_2013_2017),
#             col="chocolate1") + 
#   geom_line(aes(x=month,
#                 y=rmse_2013_2019),
#             col="darkgoldenrod4") + 
#   geom_line(aes(x=month,
#                 y=rmse_2015_2017),
#             col="darkgoldenrod1") + 
#   geom_line(aes(x=month,
#                 y=rmse_2015_2019),
#             col="yellow3") + 
#   geom_line(aes(x=month,
#                 y=rmse_2017_2019),
#             col="cyan3") + 
#   theme_bw() + 
#   scale_x_continuous(breaks=c(4,8,12))


# Comparison between SRTM and LiDAR for terrain metrics
#   Linear regression
summary(lm(san_pedro_training$dem_lidar~san_pedro_training$dem_srtm))
summary(lm(san_pedro_training$slope_lidar~san_pedro_training$slope_srtm))
summary(lm(san_pedro_training$rel_elev_lidar~san_pedro_training$rel_elev_srtm))
summary(lm(pendleton_training$dem_lidar~pendleton_training$dem_srtm))
summary(lm(pendleton_training$slope_lidar~pendleton_training$slope_srtm))
summary(lm(pendleton_training$rel_elev_lidar~pendleton_training$rel_elev_srtm))
#   Scatterplot
plot(san_pedro_training$dem_lidar, san_pedro_training$dem_srtm)
plot(san_pedro_training$slope_lidar, san_pedro_training$slope_srtm)
plot(pendleton_training$dem_lidar, pendleton_training$dem_srtm)
plot(pendleton_training$slope_lidar, pendleton_training$slope_srtm)
plot(pendleton_training$rel_elev_lidar~pendleton_training$rel_elev_srtm)
plot(pendleton_training$rel_elev_lidar~pendleton_training$rel_elev_srtm)

# Comparisons of 3D metrics by class - Pendleton
#   Canopy Height
pendleton_chm_plot <- ggplot(pendleton_training) + 
  geom_boxplot(aes(x=class, y=chm_lidar/3.28084, group=class)) + 
  scale_y_continuous(limits=c(0,30)) + 
  ylab("Canopy Height (m)") + 
  xlab("Class") +
  theme_bw()
pendleton_chm_plot
ggsave(here::here("output_plots","pendleton_chm.png"), pendleton_chm_plot,
       width=14, height=3)
#   Relative Elevation
pendleton_rel_elev_plot <- ggplot(pendleton_training) + 
  geom_boxplot(aes(x=class, y=rel_elev_lidar/3.28084, group=class)) + 
  scale_y_continuous(limits=c(0,600)) + 
  ylab("Relative Elevation Over Creek (m)") + 
  xlab("Class") +
  theme_bw()
pendleton_rel_elev_plot
ggsave(here::here("output_plots","pendleton_rel_elev.png"), pendleton_rel_elev_plot,
       width=14, height=3)
pendleton_slope_plot <- ggplot(pendleton_training) + 
  geom_boxplot(aes(x=class, y=slope_lidar, group=class)) + 
  scale_y_continuous(limits=c(0,60)) + 
  ylab("Terrain Slope (m)") + 
  xlab("Class") +
  theme_bw()
pendleton_slope_plot
ggsave(here::here("output_plots","pendleton_slope.png"), pendleton_slope_plot,
       width=14, height=3)



# Comparisons of 3D metrics by class - San Pedro
#   Canopy Height
san_pedro_chm_plot <- ggplot(san_pedro_training) + 
  geom_boxplot(aes(x=class, y=chm_lidar, group=class)) + 
  scale_y_continuous(limits=c(0,30)) + 
  ylab("Canopy Height (m)") + 
  xlab("Class") +
  theme_bw()
san_pedro_chm_plot
ggsave(here::here("output_plots","san_pedro_chm.png"), san_pedro_chm_plot,
       width=12, height=3)
#   Relative Elevation
san_pedro_rel_elev_plot <- ggplot(san_pedro_training) + 
  geom_boxplot(aes(x=class, y=rel_elev_lidar, group=class)) + 
  scale_y_continuous(limits=c(0,600)) + 
  ylab("Relative Elevation Over Creek (m)") + 
  xlab("Class") +
  theme_bw()
san_pedro_rel_elev_plot
ggsave(here::here("output_plots","san_pedro_rel_elev.png"), san_pedro_rel_elev_plot,
       width=12, height=3)
san_pedro_slope_plot <- ggplot(san_pedro_training) + 
  geom_boxplot(aes(x=class, y=slope_lidar, group=class)) + 
  scale_y_continuous(limits=c(0,60)) + 
  ylab("Terrain Slope (m)") + 
  xlab("Class") +
  theme_bw()
san_pedro_slope_plot
ggsave(here::here("output_plots","san_pedro_slope.png"), san_pedro_slope_plot,
       width=12, height=3)




# Comparisons of year differences for Vandenberg
# First, compare overall accuracy between self-trained vs. off-trained years
t.test(c(.919, .874, .903, .896),c(.864, .853, .861, .843, .821, .829, .839, .829, .834, .833, .785, .873))
mean(c(.919, .874, .903, .896))
mean(c(.864, .853, .861, .843, .821, .829, .839, .829, .834, .833, .785, .873))
# Next, do the same for balanced riparian accuracy
t.test(c(.991, .989, .984, .974),c(.990, .991, .970, .981, .915, .931, .942, .964, .955, .970, .987, .964))
mean(c(.991, .989, .984, .974))
mean(c(.990, .991, .970, .981, .915, .931, .942, .964, .955, .970, .987, .964))
# Next, compare off-trained years which had similar water availability to off-trained years with different water availability 



acc_diff <- c(0.0549,
              0.0579,
              0.0859,
              0.0210,
              0.0310,
              0.0890,
              0.0740,
              0.0820,
              0.0300,
              0.0570,
              0.0670,
              0.0620,
              0.0000,
              0.0000,
              0.0000,
              0.0000)
water_diff <- c(19,
                379,
                335,
                19,
                360,
                316,
                379,
                360,
                44,
                335,
                316,
                44,
                0,
                0,
                0,
                0)

