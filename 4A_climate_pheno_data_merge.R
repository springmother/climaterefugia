######################################################################
#### Purpose: Climate and Phenology Data Merge 
#### By: Grace Peven
#### Date Created: 10/16/2024
######################################################################

# load libraries

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lubridate)
require(lme4)
require(glmmTMB)
library(RColorBrewer)
library(corrplot)
library(performance)
library(moments)

######################################################################
# 1.  Load and prep phenology data
######################################################################

setwd()  # set working directory

output_summary = read.csv("phenometrics_summary.csv") # reads in as output_summary object
output_summary_2024 = read.csv("phenometrics_summary_2024.csv") # processed 2024 separately when data came available

output_summary_2024$year = as.factor(output_summary_2024$year) 

# join all years
phen_springs = output_summary %>% full_join(output_summary_2024, by = c("spring_name", "year", "site_type", "max_ndvi", "max_ndvi_doy", "amplitude_ndvi", "dur_80th_ndvi", "sos_doy_ndvi", "eos_doy_ndvi", "gsl_ndvi"))


# round timing metrics to closest jday so we can extract climate info for those days

phen_springs$max_ndvi_doy = round(phen_springs$max_ndvi_doy)
phen_springs$sos_doy_ndvi = round(phen_springs$sos_doy_ndvi)
phen_springs$eos_doy_ndvi = round(phen_springs$eos_doy_ndvi)
phen_springs$gsl_ndvi = round(phen_springs$gsl_ndvi)

# only keep relevant fields 
phen_springs = phen_springs %>% 
  select(spring_name, site_type, year, max_ndvi, max_ndvi_doy, amplitude_ndvi, dur_80th_ndvi, sos_doy_ndvi, eos_doy_ndvi, gsl_ndvi)

# rename to match climate field
phen_springs = phen_springs %>%
  rename(SpringName = spring_name, Water_Year = year)





###################################################################
### 2. Summazrize monthly, seasonal, and annual CWB
###################################################################

springs_CWB = read.csv("springs_CWB.csv")

# filter to relevant years
springs_CWB = springs_CWB %>% filter(Water_Year %in% c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))

# calculate total CWB per water year
springs_CWB = springs_CWB %>%
  group_by(SpringName, Water_Year) %>%
  mutate(CWB_tot_mm = sum(CWB_mm, na.rm = TRUE))

# Make sure there are no duplicates
springs_CWB_annual = springs_CWB %>% 
  distinct(SpringName, Water_Year, CWB_tot_mm, snowmelt_timing)

# Calculate seasonal total CWB per water year
springs_CWB_annual= springs_CWB %>%
  group_by(SpringName, Water_Year, Season) %>%
  mutate(seasonal_tot_CWB = sum(CWB_mm, na.rm = TRUE)) 

# Calculate monthly total CWB per water year
springs_CWB_annual_month= springs_CWB_annual %>%
  group_by(SpringName, Water_Year, Month) %>%
  mutate(monthly_tot_CWB = sum(CWB_mm, na.rm = TRUE))
  
# Make sure no duplicates exist
springs_CWB_annual = springs_CWB_annual %>%
  distinct(SpringName, Water_Year, Season, seasonal_tot_CWB, CWB_tot_mm, snowmelt_timing) %>%
  pivot_wider(names_from = Season, values_from = seasonal_tot_CWB)

# clean up field names
springs_CWB_annual = springs_CWB_annual %>% 
  rename(Winter_CWB_mm = Winter, Fall_CWB_mm = Fall, Spring_CWB_mm = Spring, Summer_CWB_mm = Summer)

# change month labels to abbreviated text
springs_CWB_annual_month$Month = month(springs_CWB_annual_month$Month, label = TRUE,abbr = TRUE)
  
# Change structure of monthly CWB columns
springs_CWB_annual_month = springs_CWB_annual_month %>%
  distinct(SpringName, Water_Year, monthly_tot_CWB, CWB_tot_mm, snowmelt_timing) %>%
  pivot_wider(names_from = Month, values_from = monthly_tot_CWB)
  
springs_CWB_annual_month = springs_CWB_annual_month %>% 
  rename(Jan_CWB_mm = Jan, Feb_CWB_mm = Feb, Mar_CWB_mm = Mar, Apr_CWB_mm = Apr, May_CWB_mm = May, Jun_CWB_mm = Jun, Jul_CWB_mm = Jul, Aug_CWB_mm = Aug, Sep_CWB_mm = Sep, Oct_CWB_mm = Oct, Nov_CWB_mm = Nov, Dec_CWB_mm = Dec)


# Merge annual and monthly CWB metrics
springs_CWB_annual = merge(springs_CWB_annual_month, springs_CWB_annual, by = c("SpringName", "Water_Year", "CWB_tot_mm", "snowmelt_timing"))



##################################################################
# 3. Join Phenology and Climate tables 
###################################################################

all_climate_phen = phen_springs %>%
  full_join(springs_CWB_annual, by = c("SpringName", "Water_Year"), relationship = "many-to-one")


all_climate_phen$site_type <- as.character(all_climate_phen$site_type)
all_climate_phen$site_type <- ifelse(all_climate_phen$site_type == "Matrix", "Non-Spring", all_climate_phen$site_type)
all_climate_phen$site_type = factor(all_climate_phen$site_type, levels = c("Spring", "Non-Spring"))


################################################################### 
## 4. Add microclimate data to dataframe
###################################################################

springs_topo = read.csv("all_sites_topo_vars.csv")

springs_topo$site_type <- as.character(springs_topo$site_type)
springs_topo$site_type <- ifelse(springs_topo$site_type == "Control", "Non-Spring", springs_topo$site_type)

#### Note: The Wind Exposure Index (WEI) was calculated in ArcGIS, but estimates for wind are unreliable in the study area
####       given the complex and steep canyon terrain, so was not included in final analysis

all_climate_phen <- merge(all_climate_phen, springs_topo, by = c("SpringName", "site_type"), all.x = TRUE)

all_climate_phen = all_climate_phen %>% select(-Slope_rad, -Aspect_rad, -Lat_rad, -Aspect_fold) # remove unnecessary columns


# readjust levels
all_climate_phen$Water_Year = factor(all_climate_phen$Water_Year, levels = c(2021,2017,2018,2019,2020,2022,2023,2024))

# make sure SpringName is a factor
all_climate_phen$SpringName = as.factor(all_climate_phen$SpringName)

cwb_phen = all_climate_phen
cwb_phen  = cwb_phen  %>% filter(Water_Year %in% c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))

# Save final file
write.csv(cwb_phen, "cwb_phen.csv")


######################################################################
### Ready for 4B_data_exploration.R
######################################################################






##########################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################################


##########################################################################
# After running 4B_data_exploration.R, discovered 2 outliers. 
# Came back to this script and removed outliers and re-ran lines 
# 146 - 150 with final cwb phen file
##########################################################################

# Even though the outliers are only for the control sites I removed both sites
# to maintain even sample sizes for each year/spring

cwb_phen <- cwb_phen%>%
  filter(!(SpringName == "Gun Shot" & Water_Year == "2017"))
cwb_phen <- cwb_phen %>%
  filter(!(SpringName == "Snail" & Water_Year == "2019"))





