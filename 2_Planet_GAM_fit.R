#####################################################################
##### Purpose: Fit GAM to NDVI time series data
##### Date Created: 09/16/2024
##### By: Jan Eitel & Grace Peven
####################################################################


require(lubridate)
require(tidygam)
require(mgcv)
require(dplyr)
require(ggplot2)


setwd()
springs_all = read.csv("Spring_VIs_ALL_2017_2024.csv") 


############################################################
## 1. Prepare dataset for gam fitting
############################################################


extract = which(springs_all$DOY <=300)          ### Select Julian Day Length
springs_all = springs_all[extract,]                                   

springs_all = filter(springs_all, SpringName !="Fireweed")         ### Remove Fireweed spring. After running through loop it's clear that the spring and matrix veg is not disinguishable.

spring_name_vector = unique(springs_all$SpringName)
spring_type_vector = unique(springs_all$Type)
spring_year_vector = unique(springs_all$Year)


output_summary = c()        # create blank object for summary results
output_time_series_ndvi = c()     # create blank object for all predicted results

output_aic = c()           # create blank object for final model stats

for(i in 1: length(spring_name_vector)) {
  
  extract = which(springs_all$SpringName == spring_name_vector[i])
  sub1 = springs_all[extract,]
  
  for(ii in 1:length(spring_type_vector)) {
    
    extract = which(sub1$Type == spring_type_vector[ii])
    sub2 = sub1[extract,]

    for(iii in 1:length(spring_year_vector)) {

      extract = which(sub2$Year == spring_year_vector[iii])
      sub3 = sub2[extract,]

      extract = which(!is.na(sub3$ndvi_mean))  # removes NA values
      sub3 = sub3[extract,]

     
      
  #############################################################################
  ## Model fitting and selection 
  #############################################################################
      
  # fit NDVI gam models with varying k-values 
    
        model5_ndvi = gam(ndvi_mean ~ s(DOY, bs = "cs", k = 5, id = 1), data = sub3)
        model6_ndvi = gam(ndvi_mean ~ s(DOY, bs = "cs", k = 6, id = 1), data = sub3)
        model7_ndvi = gam(ndvi_mean ~ s(DOY, bs = "cs", k = 7, id = 1), data = sub3)
        model8_ndvi = gam(ndvi_mean ~ s(DOY, bs = "cs", k = 8, id = 1), data = sub3)
        model9_ndvi = gam(ndvi_mean ~ s(DOY, bs = "cs", k = 9, id = 1), data = sub3)
        model10_ndvi = gam(ndvi_mean ~ s(DOY, bs = "cs", k =10, id = 1), data = sub3)
        
        model_list = list(model5_ndvi, model6_ndvi, model7_ndvi, model8_ndvi, model9_ndvi, model10_ndvi)
        
        aic = data.frame(AIC(model5_ndvi, model6_ndvi, model7_ndvi, model8_ndvi, model9_ndvi, model10_ndvi))
    
        best_aic = which(aic$AIC == min(aic$AIC))
        
        final_model_ndvi = model_list[[best_aic]]
        
        
 
      # predict values based on final model
      p1_ndvi = data.frame(predict_gam(final_model_ndvi, length_out = 175))   
      
      colnames(p1_ndvi)[2] = c("pred_ndvi")                           # change column name so we don't confuse it with observed value
      
 
  #############################################################################
  ## Extract metrics from predicted values    
  #############################################################################    
      # ndvi 
      max_ndvi = max(p1_ndvi$pred_ndvi, na.rm = TRUE) # seasonal max
      extract = which(p1_ndvi$pred_ndvi == max_ndvi) 
      
      max_ndvi_doy = p1_ndvi$DOY[extract]             # doy of seasonal max
      
      amplitude_ndvi = (quantile(p1_ndvi$pred_ndvi, 0.9) - quantile(p1_ndvi$pred_ndvi, 0.1))
      

      # => 80th percentile duration  
      dur_80th_ndvi = which(p1_ndvi$pred_ndvi >= quantile(p1_ndvi$pred_ndvi, 0.8))
      dur_80th_ndvi_start = p1_ndvi$DOY[dur_80th_ndvi[1]]
      dur_80th_ndvi_end = p1_ndvi$DOY[dur_80th_ndvi[length(dur_80th_ndvi)]]   
      dur_80th_ndvi = dur_80th_ndvi_end - dur_80th_ndvi_start
      
      
      
      # sos, eos, and gsl
    
      threshold_ndvi <- 0.5 * amplitude_ndvi       
      
      sos_doy_ndvi <- p1_ndvi$DOY[min(which(p1_ndvi$pred_ndvi >= max_ndvi - threshold_ndvi))]

      eos_doy_ndvi <- p1_ndvi$DOY[max(which(p1_ndvi$pred_ndvi >= max_ndvi -threshold_ndvi))]
      
      gsl_ndvi = eos_doy_ndvi - sos_doy_ndvi
  
  
  ############################################################################ 
      
    # create spring name, type, and year columns
      
      spring_name = spring_name_vector[i]
      site_type = spring_type_vector[ii]
      year = spring_year_vector[iii]
      
  ############################################################################
    # bind together summary metric columns 
      
      out = cbind.data.frame(spring_name,  site_type , year, max_ndvi, max_ndvi_doy, amplitude_ndvi, dur_80th_ndvi, sos_doy_ndvi, eos_doy_ndvi, gsl_ndvi)
                             
  
      output_summary = rbind(output_summary, out) # final summary table
  #############################################################################
      # predicted time series table
      
       # add year, spring name, and site type to predicted values table
         p1_ndvi$year = year
         p1_ndvi$spring_name = spring_name
         p1_ndvi$type = site_type

         
         
    # bind together all time series columns     
        
        output_time_series_ndvi = rbind(output_time_series_ndvi, p1_ndvi)


   #############################################################################       
     # extract final model stats  
         
        aic = AIC(final_model_ndvi)
        model_names = c("ndvi")

        out_aic = cbind.data.frame(spring_name, site_type, year, model_names, aic)
       output_aic = rbind(output_aic, out_aic)
         

    }
    
  }
}





# save R objects
save(output_summary_2024, file = "R_objects/phenometrics_summary_2024.RData")
save(output_summary, file = "R_objects/phenometrics_summary.RData")
save(output_aic, file = "R_objects/gam_models_aic.RData")
save(output_time_series_ndvi, file = "R_objects/ndvi_timeseries_gam.RData")

# save csv files
write.csv(output_summary_2024, file = "CSVs/phenometrics_summary_2024.csv")
write.csv(output_summary, file = "CSVs/phenometrics_summary.csv")
write.csv(output_time_series_ndvi, file = "CSVs/ndvi_timeseries_gam.csv")






  
################################################
## Visual quality check 
###############################################

## fitted
output_time_series_ndvi$year = as.factor(output_time_series_ndvi$year)


output_time_series_ndvi%>% filter(year == 2017  )%>% #cycle through various springs and indices to check
  ggplot(aes(DOY, pred_ndvi, color = type))+
  geom_line(lwd = 1.5)+
  facet_wrap(~year)+
  facet_wrap(~spring_name)


output_time_series_evi %>% filter(year == 2017 )%>% 
  ggplot(aes(DOY, pred_evi, color = type))+
  geom_line(lwd = 1.5)+
  ylim(0,2)+
  facet_wrap(~spring_name)


## summary metrics 

output_summary$year = as.factor(output_summary$year)
output_summary$site_type = as.factor(output_summary$site_type)

################################################################################
# 80th percentile duration metric 
################################################################################
library(gridExtra)

plot1 = ggplot(output_summary, aes(site_type,  dur_80th_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Duration =>80th percentile NDVI")+
  xlab("Site type")+
  ylab("Julian day duration =>80th percentile")

plot2 = ggplot(output_summary, aes(site_type,  dur_80th_evi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Duration =>80th percentile evi")+
  xlab("Site type")+
  ylab("Julian day duration =>80th percentile")

grid.arrange(plot1, plot2, ncol = 2)

# by year
ggplot(output_summary, aes(site_type,  dur_80th_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Duration =>80th percentile NDVI")+
  xlab("Site type")+
  ylab("Julian day duration =>80th percentile")+
  facet_wrap(~year)

ggplot(output_summary, aes(site_type,  dur_80th_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Duration =>80th percentile NDVI")+
  xlab("Site type")+
  ylab("Julian day duration =>80th percentile")+
  facet_wrap(~year)


# by spring
ggplot(output_summary, aes(site_type,  dur_80th_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Duration =>80th percentile NDVI")+
  xlab("Site type")+
  ylab("Julian day duration =>80th percentile")+
  facet_wrap(~spring_name)

################################################################################
# start of season 
################################################################################

# all data
plot1 = ggplot(output_summary, aes(site_type,  sos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Start of Season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")



# by year
ggplot(output_summary, aes(site_type,  sos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Start of Season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~year)

# by spring
ggplot(output_summary, aes(site_type,  sos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Start of season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~spring_name)

################################################################################
# end of season
################################################################################

# all data
plot1 = ggplot(output_summary, aes(site_type,  eos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("End of Season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")



# by year
ggplot(output_summary, aes(site_type,  eos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("End of Season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~year)

# by spring
ggplot(output_summary, aes(site_type,  eos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("End of season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~spring_name)

################################################################################
# Growing Season Length
################################################################################

# all data
plot1 = ggplot(output_summary, aes(site_type,  gsl_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Growing season length (NDVI)")+
  xlab("Site type")+
  ylab("Number of days")



# by year
ggplot(output_summary, aes(site_type,gsl_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Growing season length (NDVI)")+
  xlab("Site type")+
  ylab("Number of days")+
  facet_wrap(~year)

# by spring
ggplot(output_summary, aes(site_type,gsl_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Growing season length (NDVI)")+
  xlab("Site type")+
  ylab("Number of days")+
  facet_wrap(~spring_name)
################################################################################
# max NDVI DOY
################################################################################

# all data
plot1 = ggplot(output_summary, aes(site_type,  max_ndvi_doy, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Max DOY (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")



grid.arrange(plot1, plot2, ncol = 2)

# by year
ggplot(output_summary, aes(site_type,  max_ndvi_doy, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Max DOY (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~year)

# by spring
ggplot(output_summary, aes(site_type,  eos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("End of season (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~spring_name)

################################################################################
# max NDVI
################################################################################

# all data
plot1 = ggplot(output_summary, aes(site_type,  max_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Seasonal max NDVI")+
  xlab("Site type")+
  ylab("NDVI")



# by year
ggplot(output_summary, aes(site_type, max_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Seasonal max NDVI")+
  xlab("Site type")+
  ylab("NDVI")+
  facet_wrap(~year)

# by spring
ggplot(output_summary, aes(site_type, max_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Seasonal max NDVI")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~spring_name)



################################################################################
# Seasonal amplitude
################################################################################

# all data
plot1 = ggplot(output_summary, aes(site_type,  amplitude_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Seasonal amplitude (NDVI)")+
  xlab("Site type")+
  ylab("NDVI")



# by year
ggplot(output_summary, aes(site_type,   amplitude_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Seasonal amplitude (NDVI)")+
  xlab("Site type")+
  ylab("NDVI")+
  facet_wrap(~year)

# by spring
ggplot(output_summary, aes(site_type,   amplitude_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("Seasonal amplitude (NDVI)")+
  xlab("Site type")+
  ylab("Julian day")+
  facet_wrap(~spring_name)


ggplot(output_summary, aes(year,  sos_doy_ndvi, fill = site_type))+
  geom_boxplot()+
  scale_fill_manual(values = c("Spring" = "skyblue", "Matrix" = "tan3"))+
  theme_bw()+
  ggtitle("End of Season timing")+
  xlab("Year")+
  ylab("EOS julian day")







  
  
  
  
  
  
  
  

                            
                            