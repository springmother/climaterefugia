######################################################################
#### Purpose: Climate data visualization for Big Creek study area
####          to identify interannual climate differences
#### By: Grace Peven
#### Date Created: 08/09/2024
######################################################################

# load libraries

library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lubridate)



## the following dataset was downloaded from climate toolbox on 8/9/2024 for a 4 km grid cell 
## centered at TWRS. 

###############################################################################
# 1. Load in data and prep
##############################################################################
setwd()
climate1 = read.csv("gridMET_temp_precip_vpd_TWRS_1991_2024.csv")
climate1$Date = ymd(climate1$Date)
climate1$Year = factor(climate1$Year) # make year a factor



## Date recognition
climate1$Month= month(climate1$Month, label = TRUE, abbr = TRUE)
View(climate1)


# # Function to calculate Day of Water Year
calculate_dowy <- function(date) {
  # Define the water year start based on the current date's year
  wy_start <- as.Date(paste0(ifelse(month(date) >= 10, year(date), year(date) - 1), "-10-01"))
  # Calculate Day of Water Year
  as.numeric(difftime(date, wy_start, units = "days")) + 1
}

# Apply function to create the "Day of Water Year" column
climate1$doy_wy <- sapply(climate1$Date, calculate_dowy)


# Function to calculate Water Year
calculate_water_year <- function(date) {
  ifelse(month(date) >= 10, year(date) + 1, year(date))
}

# Apply function to create the "Water Year" column
climate1$Water_Year <- sapply(climate1$Date, calculate_water_year)


# filter out irrelevant dates
climate1=
climate1 %>%
  filter(!Water_Year %in% 1991:2016)



##########################################################
# 2. Calculate cumulative precipitation 
###########################################################

climate1 = climate1 %>% 
  arrange(Water_Year, doy_wy) %>%
  group_by(Water_Year) %>%
mutate(precip_cumulative_mm = cumsum(precip_mm)) 


### Import 30 year average
normals_30 = read.csv("gridMET_TWRS_30yr.csv")

normals_30$Water_Year= "1991-2020"
climate1$Water_Year = as.character(climate1$Water_Year)

climate1 = bind_rows(normals_30, climate1)



View(climate1)

climate1 <- climate1 %>%
  mutate(
    Season = case_when(
      Month %in% c("Jan", "Feb", "Mar") ~ "Winter",
      Month %in% c("Apr", "May", "Jun") ~ "Spring",
      Month %in% c("Jul", "Aug", "Sep") ~ "Summer",
      Month %in% c("Oct", "Nov", "Dec") ~ "Fall"
    )
  )

climate1$Year = factor(climate1$Year)

##########################################################
# 3. Calculate cumulative seasonal precipitation and plot
###########################################################


seasonal_precip = climate1 %>% group_by(Water_Year, Season) %>%
  summarize (total_precip = sum(precip_mm))

seasonal_precip2 = seasonal_precip %>% filter(Season == "Fall") %>%
  mutate(diff_avg = total_precip - 97)

seasonal_precip3 = seasonal_precip %>% filter(Season == "Winter") %>%
  mutate(diff_avg = total_precip - 83)

seasonal_precip4 = seasonal_precip %>% filter(Season == "Spring") %>%
  mutate(diff_avg = total_precip - 136)

seasonal_precip5 = seasonal_precip %>% filter(Season == "Summer") %>%
  mutate(diff_avg = total_precip - 65)

seasonal_precip = rbind(seasonal_precip2, seasonal_precip3, seasonal_precip4, seasonal_precip5)

 sp_diff = seasonal_precip %>% filter(Water_Year != "1991-2020")%>%
  ggplot(aes(x = Water_Year, y = diff_avg, fill = Season))+
  geom_col(position = position_dodge(), width = 0.7, color = "black")+
  scale_fill_manual(values = c( "tan3","springgreen3", "skyblue2",  "snow3"))+
  theme_minimal()+
  ylab("Seasonal precipitation difference from 30-year (1991-2020) average (mm)")+
  xlab("Water Year")+
  scale_y_continuous(breaks = seq(-80, 80, by = 10))+
  geom_hline(yintercept = 0, linetype = "longdash", color = "black", linewidth = 1) +
  
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 11, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 11, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 11, family = "serif", color = "black"))

 
 ggsave(filename = "seasonal_precip.tiff", plot = sp_diff,  width = 8, height = 6, units = "in",dpi = 500)
 
 

View(seasonal_precip)

seasonal_precip_wy = seasonal_precip %>% filter(Water_Year != "1991-2020")%>%
ggplot(aes(x = Water_Year, y = total_precip, fill = Season))+
  geom_col(position = position_dodge(), width = 0.7, color = "black")+
  scale_fill_manual(values = c( "tan3","springgreen3", "skyblue2",  "snow3"))+
  theme_minimal()+
  ylab("Total Seasonal Precipitation (mm)")+
  xlab("Water Year")+
  scale_y_continuous(breaks = seq(0, 175, by = 25))+
  
  geom_hline(yintercept = 83, linetype = "longdash", color = "snow4", linewidth = 1.5) +

  geom_hline(yintercept = 136, linetype = "longdash", color = "springgreen4", linewidth = 1.5) +

  geom_hline(yintercept = 65, linetype = "longdash", color = "skyblue3", linewidth = 1.5) +

  geom_hline(yintercept = 97, linetype = "longdash", color = "tan4", linewidth = 1.5) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))

ggsave(filename = "seasonal_precip.tiff", plot = seasonal_precip_wy,  width = 8, height = 5, units = "in",dpi = 500)


##########################################################
# 4. Calculate annual precipitation 
###########################################################

tot_precip = climate1 %>% group_by(Year) %>%
  summarize (total_precip = sum(precip_mm))

View(tot_precip)

ggplot(tot_precip, aes(x = Year, y = total_precip)) +
  geom_col(fill = "lightblue") +
  theme_minimal()+
  ylab("Total Annual Precipitation (mm)")

# save climate objects

save(climate1, file = "all_climate.RData")
save(seasonal_precip, file = "seasonal_precip.RData")
save(tot_precip, file = "tot_precip.RData")


############################################################################


##########################################################
# 5. Summarize drought data
###########################################################


drought = read.csv("SPEI_TWRS_2024.csv")
drought$Date = mdy(drought$Date)
drought$Year = year(drought$Date)
drought$Month = month(drought$Date)


drought$Month= month(drought$Month, label = TRUE, abbr = TRUE)

drought <- drought %>%
  mutate(
    Season = case_when(
      Month %in% c("Jan", "Feb", "Mar") ~ "Winter",
      Month %in% c("Apr", "May", "Jun") ~ "Spring",
      Month %in% c("Jul", "Aug", "Sep") ~ "Summer",
      Month %in% c("Oct", "Nov", "Dec") ~ "Fall"
    )
  )

drought_seasonal = drought %>%
  group_by(Year, Season)%>%
  summarize(spei90 = mean(SPEI))

drought_seasonal$Season = factor(drought_seasonal$Season, levels = c("Winter", "Spring", "Summer", "Fall"))


drought_seasonal=
  drought_seasonal %>%
  filter(!Year %in% 1980:2016)

### Monthly SPEI

drought_month = drought %>%
  group_by(Year, Month)%>%
  summarize(PDSI = mean(pdsi), spei90 = mean(spei90d), spi90 = mean(spi90d))

drought_month=
  drought_month %>%
  filter(!Year %in% 1991:2016)

## boxplots
drought=
  drought%>%
  filter(!Year %in% 1991:2016)

ggplot(drought, aes(x = Year, group = Year, y = SPEI))+
  geom_boxplot(position = position_dodge(), width = 0.7, color = "black")+
  scale_x_continuous(breaks = seq(min(drought$Year), max(drought$Year), by = 1))+
  theme_minimal()+
  ylab("90-day SPEI")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) 

#### Seasonal SPEI 

seasonal_spei = ggplot(drought_seasonal, aes(x = Year, y = spei90, fill = Season))+
  geom_col(position = position_dodge(), width = 0.7, color = "black")+
  scale_fill_manual(values = c("snow3", "springgreen3", "skyblue2", "tan3"))+
  scale_x_continuous(breaks = seq(min(drought_seasonal$Year), max(drought_seasonal$Year), by = 1))+
  theme_minimal()+
  ylab("90-day SPEI")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))


ggsave(filename = "seasonal_spei.tiff", plot = seasonal_spei,  width = 7, height = 5, units = "in",dpi = 500)


### annual drought 

drought_annual = drought %>%
  group_by(Year)%>%
  summarize(PDSI = mean(pdsi), spei90 = mean(spei90d), spei30 = mean(spei30d), spei1y = mean(spei1y))

drought_annual=
  drought_annual %>%
  filter(!Year %in% 1991:2016)

ggplot(drought_annual, aes(x = Year, y = spei90))+
  geom_col(position = position_dodge(), width = 0.7)+
  scale_x_continuous(breaks = seq(min(drought_annual$Year), max(drought_annual$Year), by = 1))+
  theme_minimal()\

#############################################################################
#### Visualize cumulative precipitation trends ####
#############################################################################

## gridmET data

climate1$Month = factor(climate1$Month, levels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")) 
                                      
## 2017
c1 = climate1 %>% filter(Water_Year %in% c("1991-2020", "2017"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2017" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2017")+
  guides(color = guide_legend(title = NULL))
  

## 2018
c2 = climate1 %>% filter(Water_Year %in% c("1991-2020", "2018"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2018" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2018")+
  guides(color = guide_legend(title = NULL))


## 2019
c3 = climate1 %>% filter(Water_Year %in% c("1991-2020", "2019"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2019" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2019")+
  guides(color = guide_legend(title = NULL))

## 2020
c4 = climate1 %>% filter(Water_Year %in% c("1991-2020", "2020"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2020" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2020")+
  guides(color = guide_legend(title = NULL))

## 2021
c5 = climate1 %>% filter(Water_Year %in% c("1991-2020", "2021"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2021" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2021")+
  guides(color = guide_legend(title = NULL))

## 2022
c6 = climate1 %>% filter(Water_Year %in% c("1991-2020", "2022"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2022" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2022")+
  guides(color = guide_legend(title = NULL))


## 2023
c7= climate1 %>% filter(Water_Year %in% c("1991-2020", "2023"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2023" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2023")+
  guides(color = guide_legend(title = NULL))


## 2024
c8= climate1 %>% filter(Water_Year %in% c("1991-2020", "2024"))%>%
  ggplot(aes(x = Month, y = precip_cumulative_mm, color = Water_Year))+
  ylab("")+
  xlab("")+
  geom_smooth(aes(group = Water_Year), se= FALSE)+
  theme_bw()+
  ylim(0, 500)+
  scale_color_manual(values = c("1991-2020" = "blue", "2024" = "tan3"),
                     labels = c("1991-2020" = "30-year average (1991-2020)"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, family = "serif", color = "black"), 
        axis.text.y = element_text(size = 9, family = "serif", color = "black"), 
        axis.title.y = element_text(size = 11, family = "serif", color = "black"),
        axis.title.x = element_text(size = 11, family = "serif", color = "black"), 
        plot.title = element_text(size = 12, family = "serif", color = "black"), 
        legend.text = element_text(size = 10, family = "serif", color = "black"))+
  ggtitle("2024")+
  guides(color = guide_legend(title = NULL))

## make gridded precip plots

precip_plots = grid.arrange(c1, c2, c3, c4, c5, c6, c7, c8,  ncol = 2, left = text_grob("Cumulative precipitation (mm)", rot = 90, vjust = 1, family =  "serif",size = 12))

ggsave(filename = "precip_plots.tiff", plot = precip_plots,  width = 7, height = 8, units = "in",dpi = 500)







