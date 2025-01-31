######################################################################
#### Purpose: Climate and Phenology Data Exploration
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
library(RColorBrewer)
library(corrplot)


#################################################
setwd()
all_climate_phen = read.csv("cwb_phen.csv")


##############################################################################
## 5. Data Exploration ##
##############################################################################

# split into springs and non springs for data exploration

springs_climate_phen = all_climate_phen %>% filter(site_type == "Spring")
control_climate_phen = all_climate_phen %>% filter(site_type == "Non-Spring")

# 1. Look for outliers

#############################################################################

#### SPRINGS
ggplot(springs_climate_phen, aes(x = eos_doy_ndvi, y = reorder(SpringName, eos_doy_ndvi), color = Water_Year)) +
  geom_point(size = 3) +
  ggtitle("SPRINGS: End of Season DOY NDVI by Spring") +
  xlab("EOS DOY NDVI") +
  ylab("Spring Name") +
  theme_bw()

ggplot(springs_climate_phen, aes(x = CWB_tot_mm, y = reorder(SpringName, CWB_tot_mm), color = Water_Year)) +
  geom_point(size = 3) +
  ggtitle("SPRINGS: CWB") +
  xlab("CWB") +
  ylab("Spring Name") +
  theme_bw()


ggplot(springs_climate_phen, aes(x = snowmelt_timing, y = reorder(SpringName, snowmelt_timing), color = Water_Year)) +
  geom_point(size = 3) +
  ggtitle("SPRINGS: SDD") +
  xlab("NDSI SDD") +
  ylab("Spring Name") +
  theme_bw()

# no obvious outliers


#### CONTROLS
ggplot(control_climate_phen, aes(x = eos_doy_ndvi, y = reorder(SpringName, eos_doy_ndvi), color = Water_Year)) +
  geom_point(size = 3) +
  ggtitle("NON-SPRINGS: End of Season DOY NDVI by Spring") +
  xlab("EOS DOY NDVI") +
  ylab("Spring Name") +
  theme_bw()


# Gun Shot 2017, Snail 2019 are outliers

ggplot(control_climate_phen, aes(x = CWB_tot_mm, y = reorder(SpringName, CWB_tot_mm), color = Water_Year)) +
  geom_point(size = 3) +
  ggtitle("NON-SPRINGS: CWB") +
  xlab("CWB") +
  ylab("Spring Name") +
  theme_bw()


ggplot(control_climate_phen, aes(x = snowmelt_timing, y = reorder(SpringName, snowmelt_timing), color = Water_Year)) +
  geom_point(size = 3) +
  ggtitle("NON-SPRINGS: SDD") +
  xlab("NDSI SDD") +
  ylab("Spring Name") +
  theme_bw()

#############################################################################
# 3. Check for normal distribution

# springs
ggplot(springs_climate_phen, aes(x = eos_doy_ndvi)) +
  geom_histogram(stat="count") +
  theme_bw() +
  ylab("Frequency") +
  xlab("EOS timing") +
  ggtitle("SPRINGS: EOS timing")

# note: appears to be slightly left-skewed

# tests for normality

kurtosis(springs_climate_phen$eos_doy_ndvi) # kurtosis = 3.34, meaning the data contain a sharp peak

skewness(springs_climate_phen$eos_doy_ndvi) # skewness = -0.654, well within the normal range

shapiro.test(springs_climate_phen$eos_doy_ndvi) # suggests non-normality with p-value but W is very close to 1

# note: getting mixed results for normality, seems like it's slighlty non-normal but very very close

# control sites 
ggplot(control_climate_phen, aes(x = eos_doy_ndvi)) +
  geom_histogram(stat="count") +
  theme_bw() +
  ylab("Frequency") +
  xlab("EOS timing") +
  ggtitle("CONTROL: EOS timing")

# note: appears normal with some high values on the right tail

# check normality

kurtosis(control_climate_phen$eos_doy_ndvi) # kurtosis = 2.58, meaning the data indicate a slight sharp peak, but okay

skewness(control_climate_phen$eos_doy_ndvi) # skewness = 0.657, well within the normal range

shapiro.test(control_climate_phen$eos_doy_ndvi) # suggests non-normality with p-value but W is very close to 1

#############################################################################
# 4. Collinearity between predictor variables 
#############################################################################

###########################################
# SPRINGS
###########################################

## Correlation between all variables
cwb_phen_springs = cwb_phen %>% filter(site_type == "Spring")

cor_climate_springs = cor(cwb_phen[, c("eos_doy_ndvi", "TWI","WEI", "HLI","snowmelt_timing",  "CWB_tot_mm", 
                                       "Fall_CWB_mm", "Spring_CWB_mm", "Summer_CWB_mm", "Winter_CWB_mm" 
)], use = "na.or.complete", method = "spearman")

# visualize correlation matrix

col1 <- colorRampPalette(brewer.pal(9,"BrBG"))
corrplot(cor_climate_springs,method = "square", order = "FPC", tl.col = "black", tl.cex = 0.75, sig.level = 0.05, insig = "pch", pch.cex = 1, col = col1(100))

# variables with spearman's rank => 0.5 
View(cor_climate_springs)

# 1. eos_precip_cumulative_mm & sos_precip_cumulative_mm
# 2. eos_precip_cumulative_mm & seasonal_tot_precip_Winter
# 3. eos_precip_cumulative_mm & seasonal_tot_precip_Spring
# 4. eos_precip_cumulative_mm & seasonal_tot_precip_Fall
# 5. sos_precip_cumulative_mm & seasonal_tot_precip_Winter
# 6. sos_precip_cumulative_mm & seasonal_tot_precip_Spring
# 7. sos_precip_cumulative_mm & seasonal_tot_precip_Fall
# 8. sos_precip_cumulative_mm & seasonal_GDD_Spring
# 9. eos_GDD_cumulative & seasonal_GDD_Spring
# 10. eos_GDD_cumulative & seasonal_GDD_Summer
# 11. sos_GDD_cumulative & seasonal_GDD_Spring
# 12. seasonal_tot_precip_Fall & seasonal_tot_precip_Winter
# 13. seasonal_GDD_Summer & seasonal_GDD_Winter
# 14. seasonal_GDD_Summer & seasonal_GDD_Spring
# 15. seasonal_GDD_Fall & seasonal_GDD_Spring


##################################################
# NON SPRING
##################################################
## Correlation between all variables
cwb_phen_control = cwb_phen %>% filter(site_type == "Non-Spring")

cor_climate_control = cor(cwb_phen_control[, c("eos_doy_ndvi", "TWI","WEI", "HLI","snowmelt_timing",  "CWB_tot_mm", 
                                               "Fall_CWB_mm", "Spring_CWB_mm", "Summer_CWB_mm", "Winter_CWB_mm" 
)], use = "na.or.complete", method = "spearman")

#visualize correlation matrix

col1 <- colorRampPalette(brewer.pal(9,"BrBG"))
corrplot(cor_climate_control,method = "square", order = "FPC", tl.col = "black", tl.cex = 0.75, sig.level = 0.05, insig = "pch", pch.cex = 1, col = col1(100))


View(cor_climate_control)

#############################################################################
# 5. Dependent and independent variable relationships 
#############################################################################


# Explore relationships between dependent and independent variables

# precip vars

ggplot(all_climate_phen, aes(CWB_tot_mm, eos_doy_ndvi)) +
  geom_point()+
  ggtitle("Total annual precip") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~site_type)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))



ggplot(all_climate_phen, aes(PRISM_seasonal_tot_precip_cm_Winter, eos_doy_ndvi)) +
  geom_point(color = "blue")+
  ggtitle("Total Winter precip") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~site_type)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))

ggplot(all_climate_phen, aes(PRISM_seasonal_tot_precip_cm_Spring, eos_doy_ndvi)) +
  geom_point(color = "blue")+
  ggtitle("Total Spring precip") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~site_type)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))


ggplot(all_climate_phen, aes(PRISM_seasonal_tot_precip_cm_Summer, eos_doy_ndvi)) +
  geom_point(color = "blue")+
  ggtitle("Total Summer precip") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~site_type)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))

# temp vars

ggplot(all_climate_phen, aes(annual_GDD, eos_doy_ndvi)) +
  geom_point(color = "red4")+
  ggtitle("Total Annual GDDs") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~site_type)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))



ggplot(all_climate_phen, aes(seasonal_GDD_Winter, eos_doy_ndvi)) +
  geom_point(color = "red4")+
  ggtitle("Total Winter GDDs") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  facet_wrap(~site_type)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))

ggplot(all_climate_phen, aes(seasonal_GDD_Spring, eos_doy_ndvi)) +
  geom_point(color = "red4")+
  ggtitle("Total Spring GDDs") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  facet_wrap(~site_type)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))

ggplot(all_climate_phen, aes(seasonal_GDD_Summer, eos_doy_ndvi)) +
  geom_point(color = "red4")+
  ggtitle("Total Summer GDDs") + 
  xlab("")+ 
  ylab("")+
  theme_bw()+
  facet_wrap(~site_type)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(axis.title = element_text(family = "serif",size = 11), 
        title = element_text(family = "serif", size = 11))


# appears to be strong relationships between annual GDD, spring GDD, and summer GDD, 
# total annual precip, precip by EOS day, winter precip,and spring precip


######################################################################
## eos timing and SDD
#####################################################################

cwb_phen %>% filter(site_type == "Spring" & SpringName %in% c("Salix", "Artemis", "Cave", "West Goat","Goat Ridge")) %>%
  ggplot(aes(SpringName, max_ndvi_doy))+
  geom_boxplot()


NDSI_time_series %>% filter(SpringName == "Cave" & Water_Year == "2017") %>%
  ggplot(aes(doy_wy, pred_ndsi))+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab("MODIS Fractional Snow Cover")+
  xlab("Day of Water Year")+
  annotate("text", x = 175, y = 40, label =  paste("Snow Disappearance Date"),
           size = 10, color = "black", hjust = 0, fontface = "bold", family = "serif")+
  theme(axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13))











