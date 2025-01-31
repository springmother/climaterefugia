######################################################################
#### Purpose: Prepare temperature and precipitation data for analysis
#### By: Grace Peven
#### Date Created: 10/14/2024
######################################################################

# load libraries

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lubridate)
require(DHARMa)
require(lme4)
require(mgcv)
require(tidygam)

###############################################################################
# PRECIPITATION DATA PREP 
# sourced from PRISM
###############################################################################
setwd("")  # set working directory
precip = read.csv("PRISM_ppt_800m_2016_2024_GP.csv")
str(precip)

# Make updates to imported dataset

precip<- precip %>%
  filter(SpringName != "" & !is.na(SpringName)) # blank rows were importing

# Update spring names that didn't import correctly

precip$SpringName = ifelse(precip$SpringName == "Drinking Bul", "Drinking Bull Snake", precip$SpringName)
precip$SpringName = ifelse(precip$SpringName == "Nervous Grou", "Nervous Grouse", precip$SpringName)
precip$SpringName = ifelse(precip$SpringName == "North Pionee", "North Pioneer", precip$SpringName)
precip$SpringName = ifelse(precip$SpringName == "Mountain Ald", "Mountain Alder", precip$SpringName)

### Convert to cm and mm ###
precip$precip_cm = precip$ppt_in*2.54
precip$precip_mm = precip$precip_cm*10


# Date updates
precip$Date = ym(precip$Date)
precip$Year = year(precip$Date)



# Function to calculate Water Year
calculate_water_year <- function(date) {
  ifelse(month(date) >= 10, year(date) + 1, year(date))
}

# Apply function to create the "Water Year" column
precip$Water_Year <- sapply(precip$Date, calculate_water_year)


### Add month labels

precip$Month = month(precip$Date,label = TRUE, abbr = TRUE)

##############################################################################
# Calculate cumulative monthly and seasonal precip
##############################################################################
# remove water year 2025 data
precip = precip %>% filter(Water_Year!= 2025)
precip = precip %>% filter(Water_Year!= 2016)
#cumulative precip - can't do this with monthly data

precip = precip%>%
  arrange(Water_Year, Month) %>%
  group_by(Water_Year, SpringName) %>%
  mutate(precip_cumulative_cm = cumsum(precip_cm))

# cumulative precip by June
precip <- precip %>%
  group_by(Water_Year, SpringName) %>%
  mutate(June_cumulative_precip_cm = precip_cumulative_cm[Month == "Jun"]) %>%
  ungroup()

# total June precip
precip <- precip %>%
  group_by(Water_Year, SpringName) %>%
  mutate(June_tot_precip_cm = sum(precip_cm[Month == "Jun"])) %>%
  ungroup()



#seasonal
precip = precip %>%
  mutate(
    Season = case_when(
      Month %in% c("Jan", "Feb", "Mar") ~ "Winter",
      Month %in% c("Apr", "May", "Jun") ~ "Spring",
      Month %in% c("Jul", "Aug", "Sep") ~ "Summer",
      Month %in% c("Oct", "Nov", "Dec") ~ "Fall"
    )
  )


precip = precip %>% group_by(SpringName, Water_Year, Season) %>%
  mutate (seasonal_tot_precip_cm = sum(precip_cm))

# Calculate total water year precip
precip = precip %>% group_by(SpringName, Water_Year) %>%
  mutate (wy_tot_precip_cm = sum(precip_cm))


##########################################

# quick visual check 
  ggplot(precip, aes(Water_Year, wy_tot_precip_cm))+
  geom_col()+
  facet_wrap(~SpringName)








##############################################################################
##############################################################################
# Merge PRISM precip and temp data
##############################################################################
##############################################################################

# load temp data if not already loaded
springs_temp_data = read.csv(file = "springs_temp_data.csv")


precip <- precip %>%
  filter(SpringName != "Fireweed")

precip = precip %>% rename(DOY = jday)

# merge
springs_climate_all_prism = merge(precip, springs_temp_data, by = c("SpringName", "Year", "Month"), all.x = TRUE)

# remove duplicate fields and clean up table
springs_climate_all_prism = springs_climate_all_prism %>% select(-Latitude.y, -Longitude.y, -TWRS_elev_diff, -Date.x, -Aspect,   -temp_lapse, -Season.y, -MaxTempC, -MinTempC, -Day, -Elevation_ft)
springs_climate_all_prism = springs_climate_all_prism %>% rename(PRISM_precip_cm = precip_cm, PRISM_wy_tot_precip_cm = wy_tot_precip_cm, PRISM_seasonal_tot_precip_cm = seasonal_tot_precip_cm)
springs_climate_all_prism = springs_climate_all_prism %>% rename(Latitude = Latitude.x, Longitude = Longitude.x, Season = Season.x, Date = Date.y)


# save
write.csv(springs_climate_all_prism, file = "springs_climate_all_prism.csv")



####################################################################
#### Calculate Climate Water Balance
####################################################################

library(SPEI)

# merge precip and average monthly temp 

springs_spei = springs_temp_spei %>%
  merge(precip, by = c("SpringName", "Year", "Month"))

springs_spei = unique(springs_spei)
springs_spei = springs_spei %>% select(SpringName, Latitude, Longitude, Year, Water_Year, Month, Season, precip_mm, wy_tot_precip_cm, TempC_avg_sum_month)
springs_spei$Month <- match(springs_spei$Month, month.abb) # change back to number

springs_spei$PET <- thornthwaite(springs_spei$TempC_avg_sum_month, lat = 45) # estimate PET using Thornthwaite method
springs_spei$CWB_mm <- springs_spei$precip_mm - springs_spei$PET

springs_spei$CWB_mm = format(springs_spei$CWB_mm, scientific = FALSE)
springs_spei$CWB_mm = as.numeric(springs_spei$CWB_mm)
springs_spei$CWB_mm = round(springs_spei$CWB_mm, 2)

springs_CWB = springs_spei





# ## TWRS average annual climate
# 
# TWRS_cwb = read.csv("C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Climate/gridMET/gridmet_TWRS_1981_2023_temp_precip.csv")
# 
# TWRS_cwb$daily_avg_T = (TWRS_cwb$tmin_degC +TWRS_cwb$tmax_degC)/2
# 
# TWRS_cwb = TWRS_cwb %>% group_by(Year, Month) %>% summarize(monthly_avg_T = mean(daily_avg_T), monthly_precip_mm = sum(precip_mm))
# 
# TWRS_cwb$PET = thornthwaite(TWRS_cwb$monthly_avg_T, lat = 45)
# TWRS_cwb$CWB_mm = TWRS_cwb$monthly_precip_mm - TWRS_cwb$PET
# 
# TWRS_cwb_annual = TWRS_cwb %>% group_by(Year) %>% summarize(CWB_tot = sum(CWB_mm), precip_tot = sum(monthly_precip_mm))
# 
# TWRS_cwb$Month = month(TWRS_cwb$Month, label = TRUE, abbr = TRUE)
# 
# TWRS_cwb = TWRS_cwb %>%
#   mutate(
#     Season = case_when(
#       Month %in% c("Jan", "Feb", "Mar") ~ "Winter",
#       Month %in% c("Apr", "May", "Jun") ~ "Spring",
#       Month %in% c("Jul", "Aug", "Sep") ~ "Summer",
#       Month %in% c("Oct", "Nov", "Dec") ~ "Fall"
#     )
#   )
# 
# TWRS_seasonal = TWRS_cwb %>% group_by(WY, Season) %>% summarize(precip_tot = sum(monthly_precip_mm))
###############################################################################






################################################################################
### MODIS: Calculate  snow disappearance date (SDD) 
################################################################################
# MODIS attributes documentation: https://nsidc.org/sites/default/files/mod10a1-and-myd10a1-v5-global-and-local-attributes.pdf

# 1. filter out data with relevant QA flags

MODIS = read.csv("NDSI-Springs-V2-MOD10A1-061-results.csv")
View(MODIS)
MODIS$Date = as.Date(MODIS$Date,format = "%Y-%m-%d")

# remove any duplicates and category field since pixel is large enough to capture both spring and control
MODIS = MODIS %>%
  distinct(ID, Date, .keep_all = TRUE)

MODIS = MODIS %>% select(-Category, -Latitude, -Longitude)

# QA Flags

MODIS$GP_QA = ""
MODIS <- MODIS %>%
  mutate(GP_QA = ifelse(MOD10A1_061_NDSI_Snow_Cover  %in% c(250, 201, 200, 255) ,"Unusable", GP_QA))

# Create new snow cover field and remove unusable data

MODIS = MODIS %>% filter(GP_QA != "Unusable") %>%
  mutate(NDSI_Snow_Cover_PCT = MOD10A1_061_NDSI_Snow_Cover)

# Add water year and doy 

calculate_dowy <- function(date) {
  # Define the water year start based on the current date's year
  wy_start <- as.Date(paste0(ifelse(month(date) >= 10, year(date), year(date) - 1), "-10-01"))
  # Calculate Day of Water Year
  as.numeric(difftime(date, wy_start, units = "days")) + 1
}

# Apply function to create the "Day of Water Year" column
MODIS$doy_wy <- sapply(MODIS$Date, calculate_dowy)


# Function to calculate Water Year
calculate_water_year <- function(date) {
  ifelse(month(date) >= 10, year(date) + 1, year(date))
}

# Apply function to create the "Water Year" column
MODIS$Water_Year <- sapply(MODIS$Date, calculate_water_year)

# Limit to water years 2016 - 2024

MODIS = MODIS %>% filter(Water_Year %in% c(2017,2018,2019,2020,2021,2022,2023,2024) )

# visualize one spring 

MODIS %>% filter(ID == "Artemis") %>%
  ggplot(aes(doy_wy, NDSI_Snow_Cover_PCT))+ geom_point()+geom_line()+facet_wrap(~Water_Year)


####################################################################################
# Many studies use a threshold of 0.4 to describe snow vs. snow free (snow >0.4, no snow <0.4) 
# snowmelt timing could be determined by fitting a curve to the MODIS data and extracting the doy when NDSI = 40
####################################################################################


# 1. fit gam curve to data and extract snowmelt timing

spring_name_vector = unique(MODIS$ID)
spring_year_vector = unique(MODIS$Water_Year)


MODIS_snowmelt_summary = c() 
NDSI_time_series = c()

for(i in 1: length(spring_name_vector)) {
  
  extract = which(MODIS$ID == spring_name_vector[i])
  sub1 = MODIS[extract,]
  
  for(ii in 1:length(spring_year_vector)) {
    
    extract = which(sub1$Water_Year == spring_year_vector[ii])
    sub2 = sub1[extract,]
    
    

model_ndsi = gam(NDSI_Snow_Cover_PCT ~ s(doy_wy, bs = "cs", k =10, id = 1), data = sub2) # issue with k value but should be fixed with updated dataset


# predict values based on final model
p1_ndsi = data.frame(predict_gam(model_ndsi, length_out = 150))   

colnames(p1_ndsi)[2] = c("pred_ndsi")                           # change column name so we don't confuse it with observed value

# calculate doy of max snow cover
p1_ndsi = p1_ndsi %>% 
mutate(wy_doy_max_NDSI = doy_wy[which.max(pred_ndsi)])

# snowmelt timing: the first day that fractional snowcover is less than 40%

snowmelt_timing = p1_ndsi$doy_wy[which.min(ifelse(p1_ndsi$doy_wy > p1_ndsi$wy_doy_max_NDSI, p1_ndsi$pred_ndsi<40, Inf))]

SpringName = spring_name_vector[i]
Water_Year = spring_year_vector[ii]

out = cbind.data.frame(SpringName, Water_Year, snowmelt_timing)

MODIS_snowmelt_summary = rbind(MODIS_snowmelt_summary, out)


# time series

p1_ndsi$Water_Year = Water_Year
p1_ndsi$SpringName = SpringName

NDSI_time_series = rbind(NDSI_time_series, p1_ndsi)

  }
}





MODIS_snowmelt_summary$Water_Year = factor(MODIS_snowmelt_summary$Water_Year, levels = c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))

save(MODIS_snowmelt_summary, file ="MODIS_snowmelt_summary.RData" )
save(NDSI_time_series, file ="NDSI_time_series.RData" )


write.csv(MODIS_snowmelt_summary, "MODIS_snowmelt_summary.csv" )
write.csv(NDSI_time_series, "NDSI_time_series.csv" )



###############################################################################
# visually check results
###############################################################################


NDSI_time_series$Water_Year = as.factor(NDSI_time_series$Water_Year)

NDSI_time_series %>% filter(SpringName == "Wolf Skull")%>%
ggplot(aes(doy_wy, pred_ndsi, color = Water_Year))+ geom_point()

MODIS_snowmelt_summary$Water_Year = as.factor(MODIS_snowmelt_summary$Water_Year)

MODIS_snowmelt_summary %>% filter(SpringName == "Wolf Skull")%>%
ggplot( aes(Water_Year,snowmelt_timing))+geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


ggplot(MODIS_snowmelt_summary, aes(Water_Year,snowmelt_timing))+geom_boxplot()



###########################################################################
### Merge final climate file ###
###########################################################################

# Note: This is an earlier version of the final climate file prior to landing on 
#       using CWB

springs_climate_all_prism$Water_Year = as.factor(springs_climate_all_prism$Water_Year)

MODIS_snowmelt_summary= MODIS_snowmelt_summary %>% 
   mutate(SpringName = ifelse(SpringName == "Dunce9", "Dunce_9", SpringName), 
  SpringName = ifelse(SpringName == "Goat8", "Goat_8", SpringName), 
  SpringName = ifelse(SpringName == "Cougar37", "Cougar_37", SpringName), 
  SpringName = ifelse(SpringName == "Cougar39", "Cougar_39", SpringName), 
  SpringName = ifelse(SpringName == "Dunce4", "Dunce_4", SpringName))

springs_climate_final = springs_climate_all_prism%>% 
  right_join(MODIS_snowmelt_summary, by = c("Water_Year", "SpringName"))

# clean up table
springs_climate_final = springs_climate_final %>% 
  rename(NDSI_snowmelt_doy = snowmelt_timing)


write.csv(springs_climate_final, file = "springs_climate_final.csv")




################################################################################
#### Merge CWB with MODIS
################################################################################ 

springs_CWB$Water_Year = as.factor(springs_CWB$Water_Year)
springs_CWB = springs_CWB%>% 
  full_join(MODIS_snowmelt_summary, by = c("Water_Year", "SpringName"))

# clean up table
springs_CWB = springs_CWB %>% 
  rename(NDSI_snowmelt_doy = snowmelt_timing)

write.csv(springs_CWB, file = "springs_CWB.csv")













###############################################################################
## ready for 4_climate_phen_data_exploration.R ##
###############################################################################












###########################################################################
## temperature graph for Appendix
##########################################################################

temp_plot = twrs_temp %>% filter(Year != 2016) %>%
  ggplot() +
  geom_line(aes(x = jday, y = MaxTempC, color = "Max Temperature"), lwd = 0.5) +
  geom_line(aes(x = jday, y = MinTempC, color = "Min Temperature"), lwd = 0.5) +
  facet_wrap(~Year, nrow = 4) +
  ylab("Temperature (Celsius)") +
  xlab("Julian Day") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 11, family = "serif", color = "black"),
    axis.text.x = element_text(size = 11, family = "serif", color = "black"),
    axis.title.y = element_text(size = 11, family = "serif", color = "black"),
    axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
    legend.position = "bottom", 
    legend.text = element_text(size = 11, family = "serif", color = "black"), 
    legend.title = element_text(size = 11, family = "serif", color = "black")
  ) +
  scale_color_manual(
    name = "", 
    values = c("Max Temperature" = "tomato2", "Min Temperature" = "steelblue"),
    labels = c("Max Temperature" = "Max Temp", "Min Temperature" = "Min Temp")
  )

  ggsave(filename = "C:/Users/gpeven/OneDrive - University of Idaho/Springs Research/Data/Chap2_Snow/FIGURES/temp_plots.tif", plot = temp_plot,  width = 6, height = 5, units = "in",dpi = 500)
  
