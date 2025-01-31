#####################################################################
##### Purpose: Comparing phenophases for springs and non-springs
##### Date Created: 09/16/2024
##### By: Jan Eitel & Grace Peven
####################################################################


require(lubridate)
require(tidygam)
require(mgcv)
require(dplyr)
require(ggplot2)
require(lme4)
library(emmeans)
require(grid)
require(gridExtra)
require(DHARMa)   


##############################################################
## 1. Load data and prep
##############################################################
setwd()
output_summary = read.csv(file = "phenometrics_summary.csv") 
output_summary_2024 = read.csv(file = "phenometrics_summary_2024.csv") # 2024 was processed separately when data came available

output_summary_2024$year= as.factor(output_summary_2024$year)
output_summary$year= as.factor(output_summary$year)

# join tables
output_summary = output_summary %>% full_join(output_summary_2024, by = c("spring_name", "year", "site_type", "max_ndvi", "max_ndvi_doy", "amplitude_ndvi", "dur_80th_ndvi", "sos_doy_ndvi", "eos_doy_ndvi", "gsl_ndvi"))


# round timing metrics to closest jday so we can extract climate info for those days
output_summary$max_ndvi_doy = round(output_summary$max_ndvi_doy)
output_summary$sos_doy_ndvi = round(output_summary$sos_doy_ndvi)
output_summary$eos_doy_ndvi = round(output_summary$eos_doy_ndvi)
output_summary$gsl_ndvi = round(output_summary$gsl_ndvi)

# only keep relevant fields 
output_summary = output_summary %>% 
  select(spring_name, site_type, year, max_ndvi, max_ndvi_doy, amplitude_ndvi, dur_80th_ndvi, sos_doy_ndvi, eos_doy_ndvi, gsl_ndvi)


# calculate IQR

phenometrics_iqr = output_summary %>%
  group_by(site_type, year)%>%
  summarize(maxndvi_iqr = IQR(max_ndvi), maxndvi_doy_iqr = IQR(max_ndvi_doy), 
            amplitude_ndvi_iqr = IQR(amplitude_ndvi), dur80th_ndvi_iqr = IQR(dur_80th_ndvi), 
            sos_doy_ndvi_iqr = IQR(sos_doy_ndvi), eos_doy_ndvi_iqr = IQR(eos_doy_ndvi),gsl_ndvi_iqr = IQR(gsl_ndvi),
            maxndvi_sd = sd(max_ndvi), maxndvi_doy_sd = sd(max_ndvi_doy), 
            amplitude_ndvi_sd = sd(amplitude_ndvi), dur80th_ndvi_sd = sd(dur_80th_ndvi), 
            sos_doy_ndvi_sd = sd(sos_doy_ndvi), eos_doy_ndvi_sd = sd(eos_doy_ndvi),gsl_ndvi_sd = sd(gsl_ndvi))


### make sure these variables are factors

output_summary$site_type = factor(output_summary$site_type, levels = c("Spring", "Matrix")) 

output_summary$spring_name = as.factor(output_summary$spring_name)

# remove outliers (identified in 4_climate_pheno_daa_exploration.R)
output_summary<- output_summary %>%
  filter(!(spring_name == "Gun Shot" & year == "2017"))

output_summary<- output_summary%>%
  filter(!(spring_name == "Snail" & year == "2019"))



###############################################################################
## 2. ANOVA and post-hoc tests
###############################################################################

##############################################################################
# max ndvi
##############################################################################


phen_model1 = lmer(max_ndvi ~ site_type + (1|spring_name), data = output_summary)
summary(phen_model1)
emmeans(phen_model1, pairwise ~ site_type, adjust = "tukey") 


#conclusion: springs have a statistically significant higher mean than non-springs sites (>0.138 greater average max ndvi)

###########################################
### Check model assumptions of lmer model: 
###########################################

### Simulate residuals

sim_res = simulateResiduals(phen_model1)
plot(sim_res)

### Plot residuals vs. fitted values

plot(fitted(phen_model1), residuals(phen_model1), 
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted values")
abline(h = 0, col = "red")


### Normality of residuals (Q-Q plot)

qqnorm(residuals(phen_model1))
qqline(residuals(phen_model1), col = "red")


##############################################################################
# max ndvi doy
##############################################################################

phen_model2 = lmer(max_ndvi_doy ~ site_type + (1|spring_name), data = output_summary)
summary(phen_model2)
emmeans(phen_model2, pairwise ~ site_type, adjust = "tukey") 

#conclusion: springs have a statistically significant later max ndvi timing than non-springs sites (9.75 days later than non-springs sites)

###########################################
### Check model assumptions of lmer model: 
###########################################

sim_res = simulateResiduals(phen_model2)
plot(sim_res)

### Plot residuals vs. fitted values

plot(fitted(phen_model2), residuals(phen_model2), 
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted values")
abline(h = 0, col = "red")


### Normality of residuals (Q-Q plot)

qqnorm(residuals(phen_model2))
qqline(residuals(phen_model2), col = "red")


##############################################################################
# start of season julian day
##############################################################################

phen_model3 = lmer(sos_doy_ndvi ~ site_type + (1|spring_name), data = output_summary)
summary(phen_model3)
emmeans(phen_model3, pairwise ~ site_type, adjust = "tukey") 

#conclusion: springs have a statistically significant later sos timing than non-springs sites (6.23 days later than non-springs sites)

###########################################
### Check model assumptions of lmer model: 
###########################################

sim_res = simulateResiduals(phen_model3)
plot(sim_res)

### Plot residuals vs. fitted values

plot(fitted(phen_model3), residuals(phen_model3), 
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted values")
abline(h = 0, col = "red")


### Normality of residuals (Q-Q plot)

qqnorm(residuals(phen_model3))
qqline(residuals(phen_model3), col = "red")


##############################################################################
# end of season julian day
##############################################################################

phen_model4 = lmer(eos_doy_ndvi ~ site_type + (1|spring_name), data = output_summary)
summary(phen_model4)
emmeans(phen_model4, pairwise ~ site_type, adjust = "tukey") 

#conclusion: springs have a statistically significant later eos timing than non-springs sites (21.6 days (!!!) later than non-springs sites)

###########################################
### Check model assumptions of lmer model: 
###########################################

sim_res = simulateResiduals(phen_model4)
plot(sim_res)

### Plot residuals vs. fitted values

plot(fitted(phen_model4), residuals(phen_model4), 
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted values")
abline(h = 0, col = "red")


### Normality of residuals (Q-Q plot)

qqnorm(residuals(phen_model4))
qqline(residuals(phen_model4), col = "red")


##############################################################################
# growing season length 
##############################################################################


phen_model5 = lmer(gsl_ndvi ~ site_type + (1|spring_name), data = output_summary)
summary(phen_model5)
emmeans(phen_model5, pairwise ~ site_type, adjust = "tukey") 

#conclusion: springs have a statistically significant longer growing seasons than non-springs sites (15.4 days longer than non-springs sites)



##############################################################################
# amplitude 
##############################################################################


phen_model6 = lmer(amplitude_ndvi ~ site_type + (1|spring_name), data = output_summary)
summary(phen_model6)
emmeans(phen_model6, pairwise ~ site_type, adjust = "tukey") 

#conclusion: springs have a statistically significant larger seasonal amplitude than non-springs sites (>0.06 larger than non-springs sites)
#            but I interpret this as not that large or important ecologically


###########################################
### Check model assumptions of lmer model: 
###########################################

sim_res = simulateResiduals(phen_model5)
plot(sim_res)

### Plot residuals vs. fitted values

plot(fitted(phen_model5), residuals(phen_model5), 
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted values")
abline(h = 0, col = "red")


### Normality of residuals (Q-Q plot)

qqnorm(residuals(phen_model5))
qqline(residuals(phen_model5), col = "red")





#############################################################################
## 3. T-test IQR:: Comparing springs to non-springs sites
#############################################################################

phenometrics_iqr$year = factor(phenometrics_iqr$year, levels = c("2021", "2017", "2018", "2019", "2020", "2022", "2023", "2024")) 
phenometrics_iqr$site_type = as.factor(phenometrics_iqr$site_type)

# originally I was using aov(), but since I'm only comparing two groups a paired t-test makes more sense

# phenometrics_iqr_springs =phenometrics_iqr %>% filter(site_type == "Spring")

# max ndvi iqr

t.test(maxndvi_iqr~site_type, paired = TRUE, data = phenometrics_iqr)

#conclusion: the non-springs variability (0.05) is about 2X higher compared to the spring (0.1)

##############################################################################
# max ndvi doy iqr 

t.test(maxndvi_doy_iqr~site_type, paired = TRUE, data = phenometrics_iqr)

#conclusion: springs have a slightly higher IQR than non-springs sites (>2.54) but it is not statistically significant

##############################################################################
# start of season julian day iqr

t.test(sos_doy_ndvi_iqr~site_type, paired = TRUE, data = phenometrics_iqr)

#conclusion: springs have slightly lower variability (1 day) but not significant

##############################################################################
# end of season julian day iqr


t.test(eos_doy_ndvi_iqr~site_type, paired = TRUE, data = phenometrics_iqr)

#conclusion: springs have significantly lower variability (25.2 days (!!)) than non-springs sites

##############################################################################
# growing season length iqr


t.test(gsl_ndvi_iqr~site_type, paired = TRUE, data = phenometrics_iqr)

#conclusion: springs have a statistically significant lower variability (19 days)

##############################################################################
# amplitude 
phen_model12 = aov (amplitude_ndvi ~ site_type , data = phenometrics_iqr)
summary(phen_model12)
emmeans(phen_model12, pairwise ~ site_type, adjust = "tukey") 

#conclusion: springs have a statistically significant larger seasonal amplitude than non-springs sites (>0.06 larger than non-springs sites)
#            but I interpret this as not that large or important ecologically













##############################################################################
## 4. Comparing annual differences at springs
##############################################################################




output_summary$year = factor(output_summary$year, levels = c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")) 
output_summary_springs =output_summary %>% filter(site_type == "Spring")

# max ndvi 
phen_model13 = lmer(max_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model13)
emmeans(phen_model13, pairwise ~ year, adjust = "tukey") 
tidy(phen_model13)

# conclusion: 2017,2018,2019,2020 are not significantly different from each other but (2017-2020) are significantly 
#             different than (2021-2023) with an increasing max ndvi trend in later years. no sig difference between 
#             2022 and 2023 but sig differences between 2021 and all years.
          

# max ndvi doy
phen_model14 = lmer(max_ndvi_doy ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model14)
emmeans(phen_model14, pairwise ~ year, adjust = "tukey") 
tidy(phen_model14)

# conclusion: fairly large variability between springs (sd = 14.5). could be worth doing a follow-up 
#             test to understand what's driving the large variability. 
#             2021 has the earliest max ndvi timing while 2022 the latest. 2021
#             is significantly earlier than 2017, 2019,  2020, 2022, 2023 but not 2018


# sos doy
phen_model15 = lmer(sos_doy_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model15)
emmeans(phen_model15, pairwise ~ year, adjust = "tukey") 

# conclusion: less variability than max ndvi timing, but still decent sd (8 days) between springs
#             variability across years is fairly small, with 2017 being the earliest. 


# eos doy
phen_model16 = lmer(eos_doy_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model16)
emmeans(phen_model16, pairwise ~ year, adjust = "tukey") 

# conclusion: large variability between springs (20 days). 2017 has the earliest end of season timing
#             while 2022 has the latest. 2021 is middle of the road


# gsl 
phen_model17 = lmer(gsl_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model17)
emmeans(phen_model17, pairwise ~ year, adjust = "tukey") 

# conclusion: 2022 has significantly longer growing season than other years. 

##############################################################################
# thoughts on interannual phenology differences at springs:
# interesting that 2021 has large max ndvi value and earlier timing, could be that warmer temps
# encourage earlier/steeper plant growth. There may be a trend towards a longer growing season
# over time. 



##############################################################################

#############################################################################
## Comparing annual differences at non-springs
############################################################################
output_summary$year = factor(output_summary$year, levels = c("2017", "2018", "2019", "2020", "2021", "2022", "2023")) 
output_summary_non-springs =output_summary %>% filter(site_type == "Matrix")

# max ndvi 
phen_model18 = lmer(max_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model18)
emmeans(phen_model18, pairwise ~ year, adjust = "tukey") 

# conclusion: 2017,2018,2019,2020 are not significantly different from each other but (2017-2020) are significantly 
#             different than (2021-2023) with an increasing max ndvi trend in later years. no sig difference between 
#             2022 and 2023 but sig differences between 2021 and all years.


# max ndvi doy
phen_model19 = lmer(max_ndvi_doy ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model19)
emmeans(phen_model19, pairwise ~ year, adjust = "tukey") 

# conclusion: fairly large variability between springs (sd = 14.5). could be worth doing a follow-up 
#             test to understand what's driving the large variability. 
#             2021 has the earliest max ndvi timing while 2022 the latest. 2021
#             is significantly earlier than 2017, 2019,  2020, 2022, 2023 but not 2018


# sos doy
phen_model20 = lmer(sos_doy_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model20)
emmeans(phen_model20, pairwise ~ year, adjust = "tukey") 

# conclusion: less variability than max ndvi timing, but still decent sd (8 days) between springs
#             variability across years is fairly small, with 2017 being the earliest. 


# eos doy
phen_model21 = lmer(eos_doy_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model21)
emmeans(phen_model21, pairwise ~ year, adjust = "tukey") 

# conclusion: large variability between springs (20 days). 2017 has the earliest end of season timing
#             while 2022 has the latest. 2021 is middle of the road


# gsl 
phen_model22= lmer(gsl_ndvi ~ year + (1|spring_name), data = output_summary_springs)
summary(phen_model22)
emmeans(phen_model22, pairwise ~ year, adjust = "tukey") 

# conclusion: 2022 has significantly longer growing season than other years. 






##############################################################################
# Visualize differences for paper
##############################################################################                                   
phen_springs$site_type <- factor(phen_springs$site_type, 
                                   levels = c("Spring", "Matrix"), 
                                   labels = c("Spring", "Non-Spring"))


p1 = ggplot(phen_springs, aes(x = site_type, y = max_ndvi, fill = site_type)) +
  geom_boxplot() +  
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("Maximum NDVI") +  
  xlab("Site") +  
  ylab("NDVI") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(size = 11, color = "black"),  
    plot.title = element_text(size = 11, color = "black", hjust = 0.5, face = "bold"),    
    legend.position = "none"               
  )

p2 = ggplot(phen_springs, aes(x = site_type, y = max_ndvi_doy, fill = site_type)) +
  geom_boxplot() +  
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("Maximum NDVI timing") +  
  xlab("Site") +  
  ylab("Day of year") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(size = 11, color = "black"),  
    plot.title = element_text(size = 11, color = "black", hjust = 0.5, face = "bold"),    
    legend.position = "none"               
  )


p3 = ggplot(phen_springs, aes(x = site_type, y = sos_doy_ndvi, fill = site_type)) +
  geom_boxplot() +  
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("Start of season (SOS) timing") +  
  xlab("Site") +  
  ylab("Day of year") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(size = 11, color = "black"),  
    plot.title = element_text(size = 11, color = "black", hjust = 0.5, face = "bold"),    
    legend.position = "none"               
  )

p4 = ggplot(phen_springs, aes(x = site_type, y = eos_doy_ndvi, fill = site_type)) +
  geom_boxplot() +  
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("End of season (EOS) timing") +  
  xlab("Site") +  
  ylab("Day of year") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(size = 11, color = "black"),  
    plot.title = element_text(size = 11, color = "black", hjust = 0.5, face = "bold"),    
    legend.position = "none"           
  )


p5 = ggplot(phen_springs, aes(x = site_type, y = gsl_ndvi, fill = site_type)) +
  geom_boxplot() +  
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("Annual growing season length") +  
  xlab("Site") +  
  ylab("Number of days") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(size = 11, color = "black"),  
    plot.title = element_text(size = 11, color = "black", hjust = 0.5, face = "bold"),    
    legend.position = "none"           
  )
  


p6 = ggplot(phen_springs, aes(x = site_type, y = amplitude_ndvi, fill = site_type)) +
  geom_boxplot() +  
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("Annual NDVI amplitude") +  
  xlab("Site") +  
  ylab("NDVI") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(size = 11, color = "black"),  
    plot.title = element_text(size = 11, color = "black", hjust = 0.5, face = "bold"),    
    legend.position = "none"         
  )
                                   

phen_all = grid.arrange(p2, p3, p4, p5, ncol = 2)
ggsave(filename = "phen_all.tif", plot = phen_all,  width = 8, height = 7, units = "in",dpi = 500)
           

#############################################################################
## Interannual differences 
#############################################################################

eos_interannual = ggplot(phen_springs, aes(year, eos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Non-Spring" = "tan3"), labels = c("Spring" = "Spring", "Non-Spring" = "Non-Spring")) +  
  theme_bw() +
  ggtitle("End of growing season timing") +  
  xlab("Year") +  
  ylab("Day of year") +
  labs(fill = "Site")+
  theme(
    axis.text.x = element_text(size = 11, color = "black"),  
    axis.text.y = element_text(size = 11, color = "black"),  
    axis.title.x = element_text(size = 12, color = "black"),  
    axis.title.y = element_text(size = 12, color = "black"),  
    plot.title = element_text(size = 12, color = "black", hjust = 0.5, face = "bold"),    
    legend.title = element_text(size = 12, color = "black"),  
    legend.text = element_text(size = 12, color = "black")               
  )


ggsave(filename = "eos_interannual.tif", plot = eos_interannual,  width = 8, height = 5, units = "in",dpi = 500)





############################################################################
## further exploration of spring eos variability
###########################################################################
output_summary %>% filter(site_type == "Spring") %>%
  ggplot(aes(eos_doy_ndvi, y = reorder(spring_name, eos_doy_ndvi)))+
  geom_boxplot()






