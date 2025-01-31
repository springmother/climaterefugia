######################################################################
#### Purpose: Climate to Phenology Regression
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
require(DHARMa)
require(lme4)
library(performance)
library(moments)
library(car)
library(effectsize)



#############################################################################
## Model Fitting
#############################################################################
# 1. start with annual metrics of temp, precip, and snow
# 2. explore interactions vs. no interactions with site type 
# 3. assess fit using anova, AIC comparison
# 4. add microclimate variables to model
# 5. seasonal precip variables

#############################################################################


setwd()
all_climate_phen = read.csv("cwb_phen.csv")




#####################################
### 1. Climate Water Balance Models
#####################################

# remove water years with no phen data
all_climate_phen = all_climate_phen %>% filter(Water_Year %in% c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))


########################################################################################
# fit linear mixed effects models with different predictors and RE structure and compare
#######################################################################################


# Annual CWB with spring as RE
M1 = lmer(eos_doy_ndvi ~ CWB_tot_mm*site_type + HLI*site_type + TWI*site_type + Water_Year+(1|SpringName), data = all_climate_phen)
standardize_parameters(M1) 

# Spring CWB with spring as RE
M1_seasonal = lmer(eos_doy_ndvi ~ Spring_CWB_mm*site_type +snowmelt_timing*site_type+ HLI*site_type + TWI*site_type +(1|SpringName), data = all_climate_phen)
standardize_parameters(M1_seasonal) 

# compare
anova(M1, M1_seasonal) # Annual has better fit


############################################################################################################
## Cycle through seasons by changing the "Season_CWB_mm" predictor variable and compare all to annual metric
#############################################################################################################


# annual CWB with water year and spring as RE
M1_re = lmer(eos_doy_ndvi ~ CWB_tot_mm*site_type + HLI*site_type + TWI*site_type + (1|SpringName)+(1|Water_Year), data = all_climate_phen)
standardize_parameters(M1_re)
standardise(M1_re)

anova(M1_re,M1) # M1_re has best fit



##########################################
# SDD MODEL
##########################################


M2 = lmer(eos_doy_ndvi ~ snowmelt_timing*site_type  +(1|Water_Year), data = all_climate_phen)
standardize_parameters(M2) 

M2_re = lmer(eos_doy_ndvi ~ snowmelt_timing*site_type +(1|SpringName) +(1|Water_Year), data = all_climate_phen)
standardize_parameters(M2_re) 

anova(M2_re, M2) # M2_re has best fit similar to cwb model




################################################################
FINAL_MODEL = M1_re
FINAL_MODEL_snow = M2_re
standardize_parameters(FINAL_MODEL)
standardise(FINAL_MODEL)
################################################################



#############################################################################
# MODEL VALIDATION
#############################################################################
# DHARMA

#1. create DHARMa object of residuals  
simulationOutput <- simulateResiduals(fittedModel = FINAL_MODEL_snow)

#2. Plot scaled residuals
plot(simulationOutput)

# no overdispersion, distribution is correct, res vs. predicted quantiles deviate significantly

#3. Plot residuals against all predictors
### Residual issues 
plotResiduals(simulationOutput, all_climate_phen$annual_GDD)
plotResiduals(simulationOutput, all_climate_phen$PRISM_wy_tot_precip_cm)
plotResiduals(simulationOutput, all_climate_phen$Water_Year)
plotResiduals(simulationOutput, all_climate_phen$HLI)
plotResiduals(simulationOutput, all_climate_phen$TWI)
plotResiduals(simulationOutput, all_climate_phen$CWB_tot_mm)
plotResiduals(simulationOutput, all_climate_phen$snowmelt_timing)


par(mfrow = c(2,2))
plot(FINAL_MODEL_snow)
qqnorm(resid(FINAL_MODEL_snow))
qqline(resid(FINAL_MODEL_snow))


#4. Goodness-of-fit tests
testUniformity(simulationOutput) # no problems
testDispersion(simulationOutput) # no problems
testOutliers(simulationOutput)   # no problems
testZeroInflation(simulationOutput) # looks good, don't need to use ZI


# calculating x, y positions per Spring_Name
groupLocations = aggregate(all_climate_phen, list(all_climate_phen$SpringName), mean)

# calculating residuals per group
res2 = recalculateResiduals(simulationOutput, group = all_climate_phen$SpringName)
# running the spatial test on grouped residuals
testSpatialAutocorrelation(res2, x=groupLocations$Longitude, y=groupLocations$Latitude)
# need to add in lat/long data

### Distribution and normality of model residuals 
hist(resid(M1_re)) 
kurtosis(resid(M7_all_year))
skewness(resid(M7_all_year))








####################################################################
### Code to pull model results for table 
####################################################################

effectsize(FINAL_MODEL)
standardize_parameters(FINAL_MODEL)

# creates nice output table in Word
sjPlot::tab_model(standardize(FINAL_MODEL), show.re.var = TRUE,  file = "FINAL_MODEL_results.doc")
sjPlot::tab_model(standardize(FINAL_MODEL_snow), show.re.var = TRUE,  file = "FINAL_MODEL_snow_results.doc")


#############################################################################
## Visualize results 
############################################################################

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggeffects)
library(RColorBrewer)
library(scales)

####################################
### Marginal Effect Plots ###
####################################


### Climate effects 


CWB_effect = plot_model(FINAL_MODEL, type = "pred", terms = c("CWB_tot_mm", "site_type"), colors = c("skyblue", "tan3"), legend.title = "Site",
           line.size = 1.5, show.data = TRUE, jitter = 0.75, dot.size = 1.5)+
  xlab("Climate Water Balance (mm)")+
  ylab("End of growing season day of year")+
  ggtitle("a. Predicted CWB effect")+
  theme_bw()+
  theme(
    legend.text = element_text(color = "black", size = 11, family = "serif"), 
    axis.text = element_text(color = "black", size = 12, family = "serif"), 
    axis.title = element_text(color = "black", size = 12, family = "serif"), 
    legend.title = element_text(color = "black", size = 12, family = "serif", face = "bold"), 
    axis.text.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.text.y = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.y = element_text(size = 11, color = "black", family = "serif"),  
    plot.title = element_text(size = 13, color = "black", face = "bold", family = "serif"))+
  annotate("text", x = -420, y = 265, label =  paste("Spring coeff: 0.49"),
          size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")+
  annotate("text", x = -420, y = 170, label =  paste("Non-Spring coeff: 0.83"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")



snowmelt_effect = plot_model(FINAL_MODEL_snow, type = "pred", terms = c("snowmelt_timing", "site_type"), colors = c("skyblue", "tan3"), legend.title = "Site",
                        line.size = 1.5, show.data = TRUE, jitter = 0.75, dot.size = 1.5)+
  xlab("Snow Disapperance Date (day of water year)")+
  ylab("End of growing season day of year")+
  ggtitle("b. Predicted SDD effect")+
  theme_bw()+
  ylim(200,280)+
  theme(
    legend.text = element_text(color = "black", size = 11, family = "serif"), 
    axis.text = element_text(color = "black", size = 12, family = "serif"), 
    axis.title = element_text(color = "black", size = 12, family = "serif"), 
    legend.title = element_text(color = "black", size = 12, family = "serif", face = "bold"), 
    axis.text.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.text.y = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.y = element_text(size = 11, color = "black", family = "serif"),  
    plot.title = element_text(size = 13, color = "black", face = "bold", family = "serif"))+
  annotate("text", x = 87, y = 271, label =  paste("Spring coeff: -0.11"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")+
  annotate("text", x = 87, y = 240, label =  paste("Non-Spring coeff: 0.17"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")

climate_effects = grid.arrange(CWB_effect, snowmelt_effect)

ggsave(filename = "marginal_climate_effects.tif", plot = climate_effects,  width = 7, height = 6, units = "in",dpi = 500)


#######################################################################
### Topographic effects
#######################################################################


TWI_effect = plot_model(FINAL_MODEL, type = "pred", terms = c("TWI", "site_type"), colors = c("skyblue", "tan3"), legend.title = "Site",
                        line.size = 1.5, show.data = TRUE, jitter = 0.75, dot.size = 1.5)+
  xlab("Topographic Wetness Index (TWI)")+
  ylab("End of growing season day of year")+
  ggtitle("a. Predicted TWI effect")+
  theme_bw()+
  ylim(210,300)+
  theme(
    legend.text = element_text(color = "black", size = 11, family = "serif"), 
    axis.text = element_text(color = "black", size = 12, family = "serif"), 
    axis.title = element_text(color = "black", size = 12, family = "serif"), 
    legend.title = element_text(color = "black", size = 12, family = "serif", face = "bold"), 
    axis.text.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.text.y = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.y = element_text(size = 11, color = "black", family = "serif"),  
    plot.title = element_text(size = 13, color = "black", face = "bold", family = "serif"))+
  annotate("text", x = 0.07, y = 270, label =  paste("Spring coeff: -0.07"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")+
  annotate("text", x = 0.07, y = 237, label =  paste("Non-Spring coeff: 0.11"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")


HLI_effect = plot_model(FINAL_MODEL, type = "pred", terms = c("HLI", "site_type"), legend.title = "Site",colors = c("skyblue", "tan3"),
                        line.size = 1.5, show.data = TRUE, jitter = 0.75, dot.size = 1.5)+
  xlab("Heat Load Index (HLI)")+
  ylab("End of growing season day of year")+
  ggtitle("b. Predicted HLI effect")+
  xlim(0,1)+
  ylim(210,300)+
  theme_bw()+
  theme(
    legend.text = element_text(color = "black", size = 11, family = "serif"), 
    axis.text = element_text(color = "black", size = 12, family = "serif"), 
    axis.title = element_text(color = "black", size = 12, family = "serif"), 
    legend.title = element_text(color = "black", size = 12, family = "serif", face = "bold"), 
    axis.text.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.text.y = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.x = element_text(size = 11, color = "black", family = "serif"),  
    axis.title.y = element_text(size = 11, color = "black", family = "serif"),  
    plot.title = element_text(size = 13, color = "black", face = "bold", family = "serif"))+
  annotate("text", x = 0.01, y = 268, label =  paste("Spring coeff: 0.04"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")+
  annotate("text", x = 0.01, y = 230, label =  paste("Non-Spring coeff: -0.14"),
           size = 3.5, color = "black", hjust = 0, fontface = "bold", family = "serif")

topo_effects = grid.arrange(TWI_effect, HLI_effect)

ggsave(filename = "marginal_topo_effects.tif", plot = topo_effects,  width = 7, height = 5, units = "in",dpi = 500)




