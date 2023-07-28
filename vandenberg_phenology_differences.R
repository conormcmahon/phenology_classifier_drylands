
library(tidyverse)
library(sf)
library(terra)

# Load vegetation polygons
vandenberg_polygons <- st_read(here::here("vegetation_map_inputs","vandenberg_vegetation.shp"))

# Method to load and extract raster data at each polygon in each year
loadVandData <- function(year)
{
  # Load a list of images
  vandenberg_imagery <- loadImagery(list(here::here("phenology",paste("vandenberg_",year,"_phenoseries.tif",sep="")),
                                              here::here("Spectra",paste("vandenberg_summer_spectra_",year,".tif",sep="")),
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
  # Stack into one image
  vandenberg_imagery <- terra::rast(vandenberg_imagery)
  # Extract raster data at polygons
  vandenberg_training <- generateTrainingDataSeparatePolygons(vandenberg_polygons, vandenberg_imagery, "class", c("riparian","oak","chaparral","sage_scrub","eucalyptus","annual","agriculture","urban","bare","water"))
  # Convert to phenology format
  vandenberg_pheno_long <- vandenberg_training %>% 
    pivot_longer(2:13, names_to="month_str", values_to="ndvi") %>%
    mutate(month = as.numeric(substr(month_str, 7, 10)),
           year = year,
           year_col = "red")
  return(list(vandenberg_training, vandenberg_imagery, vandenberg_pheno_long))
}

# Load data for each year
vandenberg_data_2013 <- loadVandData(2013)
vandenberg_data_2014 <- loadVandData(2014)
vandenberg_data_2015 <- loadVandData(2015)
vandenberg_data_2016 <- loadVandData(2016)
vandenberg_data_2017 <- loadVandData(2017)
vandenberg_data_2018 <- loadVandData(2018)
vandenberg_data_2019 <- loadVandData(2019)
vandenberg_data_2020 <- loadVandData(2020)
vandenberg_data_2021 <- loadVandData(2021)
vandenberg_data_2022 <- loadVandData(2022)

# Stack into one dataframe
vandenberg_pheno <- rbind(vandenberg_data_2013[[3]],
                          vandenberg_data_2014[[3]],
                          vandenberg_data_2015[[3]],
                          vandenberg_data_2016[[3]],
                          vandenberg_data_2017[[3]],
                          vandenberg_data_2018[[3]],
                          vandenberg_data_2019[[3]],
                          vandenberg_data_2020[[3]],
                          vandenberg_data_2021[[3]],
                          vandenberg_data_2022[[3]])
# Set up average annual precipitation in previous water year for each phenology curve
#   Water data from Lompoc rain gauge: https://files.countyofsb.org/pwd/hydrology/historic%20data/rainfall/439mdd.pdf
precip_annual <- c(7.25,
                   7.20,
                   8.03,
                   11.69,
                   22.19,
                   8.59,
                   20.44,
                   12.97,
                   10.78,
                   10.19) * 25.4  # source in inches, output in mm
names(precip_annual) <- paste("year_20", 13:22, sep="")
vandenberg_pheno$precip <- precip_annual[paste("year_", vandenberg_pheno$year, sep="")]
# Summarize by vegetation class
vandenberg_pheno_summary <- vandenberg_pheno %>%
  group_by(class, month, year) %>%
  summarize(median_ndvi = median(ndvi, na.rm=TRUE),
            mean_ndvi = mean(ndvi, na.rm=TRUE),
            sd = sd(ndvi, na.rm=TRUE),
            se = sd/sqrt(n()),
            precip = mean(precip))
vandenberg_phenology <- ggplot() + 
  facet_wrap(~class, nrow=1) + 
  geom_errorbar(data=vandenberg_pheno_summary %>% 
                  filter(class %in% vand_nv_classes), 
                aes(x=month, y=mean_ndvi, 
                    ymin=mean_ndvi-se, ymax=mean_ndvi+se,
                    group=year, col=precip)) + 
  geom_line(data=vandenberg_pheno_summary %>% 
              filter(class %in% vand_nv_classes), 
            aes(x=month, y=mean_ndvi, 
                group=year, col=precip), linetype="dashed") + 
  scale_y_continuous(limits=c(0.2,0.85)) + 
  scale_x_continuous(breaks=c(4,8,12)) + 
  xlab("Month") + 
  ylab("NDVI") +
  theme_bw()
vandenberg_phenology
ggsave(here::here("output_plots","vandenberg_phenology_all_years.png"), vandenberg_phenology, 
       width=10, height=3)


