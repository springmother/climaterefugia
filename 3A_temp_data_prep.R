######################################################################
#### Purpose: Prepare temperature data for analysis
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


###############################################################################
# TEMPERATURE DATA PREP 
###############################################################################
# 2016 - 2024 twrs observed temperature data

setwd() # set working directory

twrs_temp = read.csv("TWRS_weather_combined_GP.csv")
str(twrs_temp)
twrs_temp$date <- paste(twrs_temp$Month, twrs_temp$Day, twrs_temp$Year, sep = "/")
twrs_temp$date = mdy(twrs_temp$date)
twrs_temp$jday = yday(twrs_temp$date)

twrs_temp <- twrs_temp %>%
  select(-X) # I don't know how this column got in 

##############################################################
## Yellow Pine NOAA weather station
##############################################################
#bring in closest weather station data to fill in gaps

YP_noaa = read.csv("NOAA_YellowPine.csv")
YP_noaa$DATE = as.Date(YP_noaa$DATE,format = "%m/%d/%Y")
YP_noaa$Year = year(YP_noaa$DATE)
YP_noaa$jday = yday(YP_noaa$DATE)


# limit to study timeframe
YP_noaa = YP_noaa %>% filter(Year  %in% c(2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024))

# convert to Celsius
YP_noaa$TMAX_C = (YP_noaa$TMAX - 32)* (5/9)
YP_noaa$TMIN_C = (YP_noaa$TMIN - 32)* (5/9)


####################################################################
## Explore statistical relationship between YP and TWRS temp data
####################################################################

# join YP and TWRS temp data keeping all rows from YP 
twrs_yp_temp<- twrs_temp %>%
  right_join(YP_noaa, by = c("jday", "Year"))


ggplot(twrs_yp_temp, aes(TMAX_C, MaxTempC))+
  geom_point()+
  ylab("TWRS observed Temp")+
  xlab("YP observed Temp")+
  geom_smooth(method = "lm")

##############################################################################
# linear regression models to predict TWRS temp
##############################################################################

temp_relate_tmax = lm(MaxTempC ~TMAX_C, data = twrs_yp_temp)
summary(temp_relate_tmax)
plot(temp_relate_tmax)


# R-squared = 0.81, RSE = 5.179
# formula:  TWRS temp = 0.916649*YPtemp + 0.813148


temp_relate_tmin = lm(MinTempC ~TMIN_C, data = twrs_yp_temp)
summary(temp_relate_tmin)
plot(temp_relate_tmin)

# R-squared = 0.86, RSE = 3.13
# formula: TWRS temp = 0.973809*YPtemp +2.758162

################################################################################
# Supplement the linear relationship between YP weather station temperature data 
# and TWRS temp when TWRS daily values were missing
################################################################################

twrs_yp_temp <-twrs_yp_temp %>%
  mutate(
    MaxTempC = ifelse(is.na(MaxTempC), (0.916649*TMAX_C+0.813148), MaxTempC),
    MinTempC = ifelse(is.na(MinTempC), (0.973809*TMIN_C+2.758162), MinTempC)
  )

# populate date field with all dates
twrs_yp_temp <-twrs_yp_temp %>%
  mutate(
    date = DATE)

# remove fields from YP weather station
twrs_yp_temp = twrs_yp_temp%>%
  select(-names(YP_noaa))

# add year and jday field back in 
twrs_yp_temp$Year = year(twrs_yp_temp$date)
twrs_yp_temp$jday = yday(twrs_yp_temp$date)


#############################################################################
twrs_temp = twrs_yp_temp # final processed temperature data
#############################################################################


###############################################################################
# Calculate temperature for each spring using the mean environmental lapse rate 
# +/- 6.5 C/1000 meters
###############################################################################

springs_elev = read.csv("all_sites_topo_vars.csv")
springs_elev = springs_elev %>% filter(SpringName != "Fireweed") 

# add field with elevation difference relative to twrs
springs_elev$TWRS_elev_diff = (springs_elev$Elevation_m - 1176)

# calculate temperature difference relative to Taylor for each spring
springs_elev$temp_lapse = (-6.5/1000)*springs_elev$TWRS_elev_diff

# replicate all temperature data for each unique spring name
springs_elev_unique = springs_elev %>% filter(site_type == "Spring")
springs_temp_data <- expand_grid(Site = springs_elev_unique$SpringName, twrs_temp)
springs_temp_data = springs_temp_data %>% rename(SpringName = Site)

# combine springs elevation with temp table
springs_temp_data <- springs_temp_data %>%
  left_join(springs_elev_unique, by = "SpringName")

# calculate new temperature values based on MELR for each spring
springs_temp_data$MaxTempC_lapse = springs_temp_data$MaxTempC + springs_temp_data$temp_lapse
springs_temp_data$MinTempC_lapse = springs_temp_data$MinTempC + springs_temp_data$temp_lapse

# plot example spring
springs_temp_data %>% filter(SpringName == "Hole") %>%
  ggplot(aes(date, MaxTempC_lapse, color = SpringName)) +
  geom_point()+
  geom_line()


###############################################################################
# Calculate Growing Degree Days
###############################################################################

springs_temp_data$GDD = ((springs_temp_data$MaxTempC_lapse + springs_temp_data$MinTempC_lapse)/2) - 5

springs_temp_data = springs_temp_data %>%
  mutate(GDD = ifelse(GDD <0, 0, GDD))

springs_temp_data = springs_temp_data %>%
  arrange(Year, jday) %>%
  group_by(Year, SpringName) %>%
  mutate(GDD_cumulative = cumsum(GDD))


springs_temp_data$Month = month(springs_temp_data$date, abbr= TRUE, label = TRUE)

#seasonal
springs_temp_data <- springs_temp_data %>%
  mutate(
    Season = case_when(
      Month %in% c("Jan", "Feb", "Mar") ~ "Winter",
      Month %in% c("Apr", "May", "Jun") ~ "Spring",
      Month %in% c("Jul", "Aug", "Sep") ~ "Summer",
      Month %in% c("Oct", "Nov", "Dec") ~ "Fall"
    )
  )


############################################################
# calculate cumulative GDD by season, month, and year
############################################################

springs_temp_data = springs_temp_data %>% group_by(SpringName, Year, Month) %>%
  mutate (monthly_GDD = sum(GDD))

springs_temp_data = springs_temp_data %>% group_by(SpringName, Year, Season) %>%
  mutate (seasonal_GDD = sum(GDD))

springs_temp_data = springs_temp_data %>% group_by(SpringName, Year) %>%
  mutate (annual_GDD = sum(GDD))

############################################################
# calculate GDD on a monthly basis
############################################################

springs_monthly_GDD = springs_temp_data %>% group_by(SpringName, Year, Month) %>%
  select(monthly_GDD)

springs_monthly_GDD = unique(springs_monthly_GDD)

springs_monthly_GDD = springs_monthly_GDD%>%
  pivot_wider(names_from = Month, values_from = monthly_GDD)

springs_monthly_GDD = springs_monthly_GDD %>% 
  select(-Jan, -Feb, -Mar, -Oct, -Nov, -Dec)

springs_monthly_GDD = springs_monthly_GDD %>%
  rename(Apr_GDD = Apr, May_GDD = May, Jun_GDD = Jun, Jul_GDD = Jul, Aug_GDD = Aug, Sep_GDD = Sep)

springs_temp_data = merge(springs_temp_data, springs_monthly_GDD, by = c("Year", "SpringName"))

############################################################
#rename to match precip data for next script
############################################################

springs_temp_data= springs_temp_data %>% rename(DOY = jday, Date = date)

############################################################
# save temperature data
############################################################
write.csv(springs_temp_data, "springs_temp_data.csv")


#################################################
## Calculate Monthly average temperature for SPEI 
#################################################

springs_temp_data$Temp_daily_avg = (springs_temp_data$MaxTempC_lapse +springs_temp_data$MinTempC_lapse)/2

springs_temp_data$Day = day(springs_temp_data$Date)

springs_temp_data = springs_temp_data %>%
  group_by(SpringName, Year, Month) %>%
  mutate(days_month = n_distinct(Day))

springs_temp_spei = springs_temp_data %>% 
  group_by(SpringName, Month, Year) %>%
  mutate(TempC_avg_sum_month = sum(Temp_daily_avg)/days_month)%>%
  select(SpringName, Year, Month, TempC_avg_sum_month)


write.csv(springs_temp_spei, "springs_temp_spei.csv")


###############################################################
## READY FOR 3B_climate_data_prep.R
###############################################################




