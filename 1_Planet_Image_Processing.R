##########################################################################
##########################################################################
#### Purpose: Planet imagery processing and vegetation index calculations
#### By: Jan Eitel and Grace Peven
#### Created: 4/27/2023
##########################################################################
##########################################################################


###########################################################################

files_per_image = 5   ### How many files are associate with a single image (e.g., can be four or five)
image_file_loc = 3    ### Which of the five files is the actual 4 band image
udm_file_loc = 4      ### which of files is the udm file that contains information about image quality
json_file_loc = 5     ### which of the files is the json file that contains information about image quality and acquisition time

### Band Order
#Blue = 1
#Green = 2
#Red = 3
#NIR = 4

require(jsonlite)
require(lubridate)
require(dplyr)
require(terra)
require(sp)
require(sf)
require(ggplot2)
require(writexl)
require(rstatix)
require(readxl)
library(tcltk)

################################################################
### 1.  Read in polygon file of spring and non-spring footprints 
################################################################

setwd()
spring_footprint = st_read(dsn = 'C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap2_Snow/Spring_Delineation/Shapefiles', layer = 'Spring_Control_Footprints_10th_pct', stringsAsFactors = T) ## sf package
View(spring_footprint) #### View attribute table

################################################################
### 2. Plot and check projection to make sure shapefile imported correctly 
################################################################

ggplot(data = spring_footprint) + 
geom_sf() +
theme_minimal() +
ggtitle("Spring Footprints")


### Check projection: Should be NAD 1983 UTM Zone 11N ###
st_crs(spring_footprint) 


################################################################
### 3. Import imagery 
################################################################

#### Iterate through each year of data and then bind at end of code

file_2017 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2018 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2019 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2020 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2021 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2022 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2023 <- tk_choose.files('C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap1_Satellite_Monitoring/PlanetScope_Images/', multi = TRUE)
file_2024 <- tk_choose.files('D:/PhD/PlanetScope/', multi = TRUE)

file_vector = merge(file_2017, file_2019, file_2019, file_2020, file_2021, file_2022, file_2023, file_2024)

### Check dimensions of file_vector (i.e., how many files were uploaded) ###

file_vector_dim = data.frame(file_vector)
file_vector_dim = dim(file_vector_dim)[1]
View(file_vector_dim)

length_extract_image = seq(3, file_vector_dim, 5)
length_extract_udm = seq(4, file_vector_dim, 5)         ### See info regarding udm here: https://developers.planet.com/docs/data/udm-2/
length_extract_json = seq(5, file_vector_dim, 5)

file_vector_image = file_vector[length_extract_image]
file_vector_udm = file_vector[length_extract_udm]
file_vector_json = file_vector[length_extract_json]


################################################################
### 3. Make Raster Brick
################################################################

ufc = rast(file_vector_image[1]) 
str(ufc)

###Plot raster and vector data ###

plotRGB(ufc, 3,2,1, stretch = 'hist')
plot(spring_footprint, add =TRUE, border = "red")   ### add polygon to raster image to check spatial alignment

### Make list of files a dataframe ###

file_vector_dim = data.frame(file_vector_image)
file_vector_dim = dim(file_vector_dim)[1]


################################################################################
### Loop to derive vegetation indices and metadata from all images in folder ###
################################################################################

output = c()

for(i in 1:file_vector_dim){                 

ufc = rast(file_vector_image[i])


ufc = ufc*0.01      #### Scaling factor - see pagge 


ndvi = (ufc[[4]]-ufc[[3]])/(ufc[[4]]+ ufc[[3]]) # (NIR-Red)/(NIR+Red) 

##########################################################
################### Extract statistics ###################
##########################################################

### NDVI stats
ndvi_mean = extract(ndvi, spring_footprint, fun = mean, na.rm =TRUE)
ndvi_mean$ID = NULL #### the output adds an ID column which I don't need. This removes the ID column and leaves me with one column of the extracted values
colnames(ndvi_mean) = c("ndvi_mean")

ndvi_max  = extract(ndvi, spring_footprint, fun = max, na.rm =TRUE)
ndvi_max$ID = NULL #### the output adds an ID column which I don't need. This removes the ID column and leaves me with one column of the extracted values
colnames(ndvi_max) = c("ndvi_max")

ndvi_min= extract(ndvi, spring_footprint, fun = min, na.rm =TRUE)
ndvi_min$ID = NULL #### the output adds an ID column which I don't need. This removes the ID column and leaves me with one column of the extracted values
colnames(ndvi_min) = c("ndvi_min")

ndvi_sd= extract(ndvi, spring_footprint, fun = sd, na.rm =TRUE)
ndvi_sd$ID = NULL #### the output adds an ID column which I don't need. This removes the ID column and leaves me with one column of the extracted values
colnames(ndvi_sd) = c("ndvi_sd")

### Pixel based quality information from udm file (see https://developers.planet.com/docs/data/udm-2/

udm = rast(file_vector_udm[i])
clear = udm[[1]]
clear= clear== 1
clear = extract(clear, spring_footprint,fun = mean,na.rm=TRUE) 
clear$ID = NULL
colnames(clear) = c("clear")

snow = udm[[2]]
snow= snow== 0
snow = extract(snow, spring_footprint,fun = mean,na.rm=TRUE)
snow$ID = NULL
colnames(snow) = c("snow")

shadow = udm[[3]]
shadow= shadow== 0
shadow = extract(shadow, spring_footprint,fun = mean,na.rm=TRUE)
shadow$ID = NULL
colnames(shadow) = c("shadow")

light_haze = udm[[4]]
light_haze= light_haze== 0
light_haze  = extract(light_haze, spring_footprint,fun = mean,na.rm=TRUE)
light_haze$ID = NULL
colnames(light_haze) = c("light_haze")

heavy_haze = udm[[5]]
heavy_haze= heavy_haze== 0
heavy_haze  = extract(heavy_haze, spring_footprint,fun = mean,na.rm=TRUE)
heavy_haze$ID = NULL
colnames(heavy_haze) = c("heavy_haze")

cloud = udm[[6]]
cloud= cloud== 0
cloud  = extract(cloud, spring_footprint,fun = mean,na.rm=TRUE)
cloud$ID = NULL
colnames(cloud) = c("cloud")

confidence = udm[[7]]
confidence= confidence> 50
confidence  = extract(confidence, spring_footprint,fun = mean,na.rm=TRUE) 
confidence$ID = NULL
colnames(confidence) = c("confidence")

### Metadata from json file. Note: will write metadata for the entire image downloaded, not extracted polygon.

metadata_info = fromJSON(file_vector_json[i], flatten=TRUE)

acquisition_time = metadata_info$properties$acquired
acquisition_time = ymd_hms(acquisition_time)

view_angle = metadata_info$properties$view_angle
cloud_percentage = metadata_info$properties$cloud_percent
heavy_haze_percentage = metadata_info$properties$heavy_haze_percent
pixel_resolution = metadata_info$properties$pixel_resolution
quality_category = metadata_info$properties$quality_category
snow_ice_percentage = metadata_info$properties$snow_ice_percent



out = cbind.data.frame(spring_footprint, acquisition_time, ndvi_mean, ndvi_min, ndvi_max, ndvi_sd,clear, snow, shadow, light_haze, heavy_haze, cloud, confidence, view_angle, cloud_percentage, heavy_haze_percentage, pixel_resolution, quality_category, snow_ice_percentage)
output = rbind(output, out)

} 



## Check number of rows and make sure all images were processed ##
View(output)


### Add Year and DOY to output table ###

as_datetime(output$acquisition_time)
ymd_hms(output$acquisition_time, tz = "UTC")
output$DOY = yday(output$acquisition_time)
output$Year = year(output$acquisition_time)


##################################
## Export Output file to csv ##
##################################


write.csv(output, "C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap2_Snow/Chapter 2 R/NDVI_timeseries/Spring_VIs_ALL_2017_2024.csv") 

###################################################################################
  
