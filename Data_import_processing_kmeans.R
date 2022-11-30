##### setup #####
library(tidyverse)
library(lubridate)
library(zoo)
library(imputeTS)
library(ggpubr)
library(signal)
library(strucchange)
library(xts)
library(dygraphs)
library(RcppRoll)
library(factoextra)
library(rioja)
library(dataRetrieval)
library(hydrostats)
library(devtools)
#devtools::install_github("markwh/streamstats")
#library(streamstats)

########### Setting up details for this script #############
# set smoothing details 
swindow <- 192

# option to write daily data as .csv
# if this is set as 1, then it will write csvs for all years and all sites, any other number and it won't 
writecsv <- 1 

#Cut off to decide how many groups there should be for each year. 
#It uses the first number of groups that has a difference between the SSE and broken stick model is less than this value
b_cutoff <- 0.5

#Option to decide if you want to use baseflow or total Q in cluster analysis
# set to 1 if using baseflow, any other number will use regular Q
use_base <- 1

# generating a custom theme to get rid of the ugly ggplot defaults 
theme_cust <- function(base_size = 11, base_family = "") {
  theme_classic() %+replace%
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text = element_text(color = "black")
    )
}

#loading a csv that has basin characteristics pulled from the NWIS gage stream stats page (Unuk and Kennicott are missing the full dataset)
gage_stat <- read.csv("gage_stats.csv") %>%
  mutate(area_m3 = Area_mi2*2.59e6)

# NWIS data download organization function. 
data_org <- function(gage_data){
  gage_data %>%
  select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
         Q_m3s = Q*0.028316847,
         Q_mm_d = (Q_m3s/gage_stat$area_m3[ gage_stat$Site_id == site])*1000, #normalizing Q by wtrshd area and converting to mm 
         period = NA,
         Q_filled = na_interpolation(Q_mm_d, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_mm_d,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')
}

# Define function to convert time series data to daily 
# second part of this function calculates base flow and adds it to the daily dataset 
daily<- function(merged){
  merged %>%
    mutate(datetime_daily = cut(datetime, 'day')) %>%
    group_by(datetime_daily) %>% 
    summarise(SC_mean = mean(SC_filled),
              Q_mean = mean(Q_filled)) %>%
    na.omit() %>%
    mutate(Q_mean = Q_mean*60*60*24)
}

base <- function(mean_dat){ 
    mean_dat %>%  
    rename(Date = datetime_daily,Q= Q_mean) %>%
    select(Date,Q) %>%
    baseflows( ts = "daily")
}



# Define function to normalize ts data and convert to an xts object 
daily_norm <- function(daily_D){
  daily_D %>%
  mutate(datetime_daily = ymd(datetime_daily))
        SC_norm <- ((daily_D$SC_mean)-min(na.omit(daily_D$SC_mean)))/(max(na.omit(daily_D$SC_mean))-min(na.omit(daily_D$SC_mean)))
       if (use_base == 1){
       Q_norm <- ((daily_D$base-min(na.omit(daily_D$base)))/(max(na.omit(daily_D$base))-min(na.omit(daily_D$base)))) 
        }else{
       Q_norm <- ((daily_D$Q_mean-min(na.omit(daily_D$Q_mean)))/(max(na.omit(daily_D$Q_mean))-min(na.omit(daily_D$Q_mean))))
        }
  daily_D <- cbind(daily_D,SC_norm,Q_norm)
  daily_D %>%
  select(datetime_daily, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_daily))
}

# define function to find membership breaks
membership_breaks <- function(memb, daily){
  for (k  in 1:length(unique(memb))) {
    temp <- daily$datetime_daily[which(memb==k)]
    if (k== 1){
      p_time <- as.POSIXct(tail(temp,n=1)) 
    } else{
      p_time [k] <- tail(temp, n=1)}
  }
  return(p_time)
}

##### nwis data #####
### loading gauge data
### add column of q in cms
### add index column
### interpolate to fill NAs in '_filled' ts columns
### add smoothed data to '_smoothed' ts columns
### NOTE: data does not need to be trimmed, just set bounds as start/end in NWIS call

#### Wolverine ########## 
## set bounds for all years
bounds16 <- as.POSIXct(c('04/28/2016 00:00:00','10/06/2016 23:45:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds17 <- as.POSIXct(c('04/20/2017 14:00:00','10/17/2017 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds18 <- as.POSIXct(c('05/22/2018 14:00:00','10/17/2018 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds19 <- as.POSIXct(c('04/30/2019 14:00:00','10/12/2019 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds20 <- as.POSIXct(c('04/18/2020 14:00:00','10/18/2020 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds21 <- as.POSIXct(c('06/04/2021 10:00:00','10/18/2021 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds22 <- as.POSIXct(c('05/13/2022 10:00:00','10/11/2022 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

## set nwis details for Wolverine 
site <- '15236900'
params <- c('00060', '00095') #discharge and EC respectively 

## 2016
gauge_data16 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = as.Date(bounds16[1]), endDate = as.Date(bounds16[2]))
gauge_data16 <- data_org(gauge_data16)

## 2017 
gauge_data17 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = as.Date(bounds17[1]), endDate = as.Date(bounds17[2])) 
gauge_data17 <- data_org(gauge_data17)

## 2018
gauge_data18 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = as.Date(bounds18[1]), endDate = as.Date(bounds18[2])) 
gauge_data18 <- data_org(gauge_data18)

## 2019
gauge_data19 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = as.Date(bounds19[1]), endDate = as.Date(bounds19[2])) 
gauge_data19 <- data_org(gauge_data19)

## 2020
gauge_data20 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using max and min dates from other years for now
                             startDate = as.Date(bounds20[1]), endDate =as.Date(bounds20[2])) 
gauge_data20 <- data_org(gauge_data20)

## 2021
gauge_data21 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using max and min dates from other years for now
                             startDate = as.Date(bounds21[1]), endDate =as.Date(bounds21[2])) 
gauge_data21 <- data_org(gauge_data21)

## 2021
gauge_data22 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using max and min dates from other years for now
                             startDate = as.Date(bounds22[1]), endDate =as.Date(bounds22[2])) 
gauge_data22 <- data_org(gauge_data22)


##### weather data #####
## load in weather data
# loading the 15 min met data from the 990m weather station. Data from Baker et al., 2019 (USGS data publication) 
# data found at: 'https://alaska.usgs.gov/data/glacier/benchmarkGlacierMassBalance/benchmarkGlacier_weather/benchmarkGlacier_weather.zip'
full_wx <- read.csv('data/wolverine990_15min_LVL2.csv') %>% 
  mutate(local_time = ymd_hm(local_time, tz = 'America/Anchorage'))

## subset and rename each year
#2016
precip16 <- full_wx %>%
  dplyr::filter(local_time >= bounds16[1],
                local_time <= bounds16[2]) %>%
  select(precip_int = TPGIncremental,
         datetime = local_time, 
         temp = site_temp)

# 2017
precip17 <- full_wx %>%
  dplyr::filter(local_time >= bounds17[1],
                local_time <= bounds17[2]) %>%
  select(precip_int = TPGIncremental,
         datetime = local_time, 
         temp = site_temp)

# 2018
precip18 <- full_wx %>%
  dplyr::filter(local_time >= min(gauge_data18$datetime),
                local_time <= max(gauge_data18$datetime)) %>%
  select(precip_int = TPGIncremental,
         datetime = local_time, 
         temp = site_temp)

# 2019 
precip19 <- full_wx %>%
  dplyr::filter(local_time >= min(gauge_data19$datetime),
                local_time <= max(gauge_data19$datetime)) %>%
  select(precip_int = TPGIncremental,
         datetime = local_time, 
         temp = site_temp)

# 2020 
precip20 <- full_wx %>%
  dplyr::filter(local_time >= min(gauge_data20$datetime),
                local_time <= max(gauge_data20$datetime)) %>%
  select(precip_int = TPGIncremental,
         datetime = local_time, 
         temp = site_temp)

# 2021 
precip21 <- full_wx %>%
  dplyr::filter(local_time >= min(gauge_data21$datetime),
                local_time <= max(gauge_data21$datetime)) %>%
  select(precip_int = TPGIncremental,
         datetime = local_time, 
         temp = site_temp)


##### Working with the time series data #####
 
##### Depth-constrained clustering (to constrain clusters in time) #####
# daily intervals - averaging the Q and EC

### 2016
## prep data
# merge nwis and met data
merged16 <- gauge_data16 %>%
  full_join(.,precip16, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged16$precip_int[is.nan(merged16$precip_int)] <- NA
merged16$precip_int <- na_interpolation(merged16$precip_int, maxgap = Inf, option = 'linear')

# create daily ts 
daily16 <- daily(merged16)

#calculate baseflow and add it to the daily dataset
base16 <-base(daily16)
daily16<-cbind(daily16,base = base16$bf)

# normalize data, convert to xts
norm16 <- daily_norm(daily16)
  
## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(norm16[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))


# find membership breaks
memb_break16 <- membership_breaks(memb, daily16)
memb_break16

## visualize results
# c-q plot
ggplot(daily16, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()
 
# q ts
ggplot(daily16, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()
 
### 2017
## prep data
# merge nwis and met data
merged17 <- gauge_data17 %>%
  full_join(.,precip17, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
  merged17$precip_int[is.nan(merged17$precip_int)] <- NA
  merged17$precip_int <- na_interpolation(merged17$precip_int, maxgap = Inf, option = 'linear')
 
# create daily ts
  daily17 <- daily(merged17)
  
#calculate baseflow and add it to the daily dataset
  base17 <-base(daily17)
  daily17<-cbind(daily17,base = base17$bf)
  
# normalize data, convert to xts
  norm17 <- daily_norm(daily17)
 
## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(norm17[,-1], method = "euclidean")
 
# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")
 
# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
memb_break17 <- membership_breaks(memb, daily17)
memb_break17
 
## visualize results
# c-q plot
ggplot(daily17, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q ( mm d"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()
 
# q ts
ggplot(daily17, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

ggplot(daily17, aes(x= ymd(datetime_daily), y= SC_mean/Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("EC/Q")))+
  xlab(NULL) + 
  theme_cust()

### 2018
## prep data
# merge nwis and met data
merged18 <- gauge_data18 %>%
  full_join(.,precip18, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged18$precip_int[is.nan(merged18$precip_int)] <- NA
merged18$precip_int <- na_interpolation(merged18$precip_int, maxgap = Inf, option = 'linear')

# create daily ts
daily18 <- daily(merged18)

#calculate baseflow and add it to the daily dataset
base18 <-base(daily18)
daily18<-cbind(daily18,base = base18$bf)

# normalize data, convert to xts
norm18 <- daily_norm(daily18)

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(norm18[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
memb_break18 <- membership_breaks(memb, daily18)
memb_break18

## visualize results
# c-q plot
ggplot(daily18, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(daily18, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

ggplot(daily18, aes(x= ymd(datetime_daily), y= SC_mean/Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("EC/Q")))+
  xlab(NULL) + 
  theme_cust()


### 2019
## prep data
# merge nwis and met data
merged19 <- gauge_data19 %>%
  full_join(.,precip19, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged19$precip_int[is.nan(merged19$precip_int)] <- NA
merged19$precip_int <- na_interpolation(merged19$precip_int, maxgap = Inf, option = 'linear')

# create daily ts
daily19 <- daily(merged19)

#calculate baseflow and add it to the daily dataset
base19 <-base(daily19)
daily19<-cbind(daily19,base = base19$bf)

# normalize data, convert to xts
norm19 <- daily_norm(daily19)

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(norm19[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <- bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,k= min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
memb_break19 <- membership_breaks(memb, daily19)
memb_break19

## visualize results
# c-q plot
ggplot(daily19, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("ln Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("ln EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(daily19, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

ggplot(daily19, aes(x= ymd(datetime_daily), y= SC_mean/Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("EC/Q")))+
  xlab(NULL) + 
  theme_cust()

### 2020
## prep data
# merge nwis and met data
merged20 <- gauge_data20 %>%
  full_join(.,precip20, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged20$precip_int[is.nan(merged20$precip_int)] <- NA
merged20$precip_int <- na_interpolation(merged20$precip_int, maxgap = Inf, option = 'linear')

# create daily ts
daily20 <- daily(merged20)

#calculate baseflow and add it to the daily dataset
base20 <-base(daily20)
daily20<-cbind(daily20,base = base20$bf)

# normalize data, convert to xts
norm20 <- daily_norm(daily20)

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(norm20[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <- bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
memb_break20 <- membership_breaks(memb, daily20)
memb_break20

## visualize results
# c-q plot
ggplot(daily20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(daily20, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

ggplot(daily20, aes(x= ymd(datetime_daily), y= SC_mean/Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("EC/Q")))+
  xlab(NULL) + 
  theme_cust()

### 2021
## prep data
# merge nwis and met data
merged21 <- gauge_data21 %>%
  full_join(.,precip21, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged21$precip_int[is.nan(merged21$precip_int)] <- NA
#merged21$precip_int <- na_interpolation(merged21$precip_int, maxgap = Inf, option = 'linear')

# create daily ts
daily21 <- daily(merged21)

#calculate baseflow and add it to the daily dataset
base21 <-base(daily21)
daily21<-cbind(daily21,base = base21$bf)

# normalize data, convert to xts
norm21 <- daily_norm(daily21)

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <- bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
#memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 
memb <- cutree(dep_clust,5) 

# find membership breaks
memb_break21 <- membership_breaks(memb, daily21)
memb_break21

## visualize results
# c-q plot
ggplot(daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(daily21, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

### 2022
## prep data
# merge nwis and met data
merged22 <- gauge_data22 %>%
  full_join(.,precip21, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged21$precip_int[is.nan(merged21$precip_int)] <- NA
#merged21$precip_int <- na_interpolation(merged21$precip_int, maxgap = Inf, option = 'linear')

# create daily ts
daily22 <- daily(merged22)

#calculate baseflow and add it to the daily dataset
base22 <-base(daily22)
daily22<-cbind(daily22,base = base22$bf)

# normalize data, convert to xts
norm22 <- daily_norm(daily22)

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(norm22[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <- bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
#memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 
memb <- cutree(dep_clust,5) 

# find membership breaks
memb_break22 <- membership_breaks(memb, daily22)
memb_break22

## visualize results
# c-q plot
ggplot(daily22, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(daily22, aes(x= ymd(datetime_daily), y= base)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


##### visualize membership breaks #####
max_length <- max(c(length(memb_break16), length(memb_break17), length(memb_break18), length(memb_break19), length(memb_break20), length(memb_break21)))
memb_break_df <- data.frame(col1 = c(as.integer(strftime(memb_break16, '%j')),                 # Create data frame with unequal vectors
                            rep(NA, max_length - length(memb_break16))),
                   col2 = c(as.integer(strftime(memb_break17, '%j')),
                            rep(NA, max_length - length(memb_break17))),
                   col3 = c(as.integer(strftime(memb_break18, '%j')),
                            rep(NA, max_length - length(memb_break18))),
                   col4 = c(as.integer(strftime(memb_break19, '%j')),
                            rep(NA, max_length - length(memb_break19))),
                   col5 = c(as.integer(strftime(memb_break20, '%j')),
                            rep(NA, max_length - length(memb_break20))),
                   col6 = c(as.integer(strftime(memb_break21, '%j')),
                            rep(NA, max_length - length(memb_break21))))
colnames(memb_break_df)<-c(2016,2017,2018,2019,2020,2021)
# memb_break_df$row_num <- seq.int(nrow(memb_break_df))
memb_break_df <- tibble(memb_break_df)%>%
  rowid_to_column('break_num') %>%
  pivot_longer(cols = -break_num, names_to = 'year', values_to = 'doy') %>%
  mutate(year = year)

ggplot(memb_break_df, aes(x = doy, y = year, color = as.factor(break_num))) +
 geom_point()

if (writecsv == 1){ 
  write.csv(daily16,"daily_16.csv", row.names = FALSE)
  write.csv(daily17,"daily_17.csv", row.names = FALSE)
  write.csv(daily18,"daily_18.csv", row.names = FALSE)
  write.csv(daily19,"daily_19.csv", row.names = FALSE)
  write.csv(daily20,"daily_20.csv", row.names = FALSE)
  write.csv(daily21,"daily_21.csv", row.names = FALSE)
  write.csv(daily22,"daily_22.csv", row.names = FALSE)
}

########### Pulling in NWIS data from other sites ##############
##### nwis data #####
## set nwis details

#### Taku #####
site <- '15041200' #Taku R. Nr Juneau
params <- c('00060', '00095') # discharge and SC

### 2019 ####
bounds <- as.POSIXct(c('05/06/2019 00:00:00','10/01/2019 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Taku_data19 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, 
                             startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Taku_data19 <- data_org(Taku_data19)

### 2020 ####
bounds <- as.POSIXct(c('04/15/2020 00:00:00','10/01/2020 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Taku_data20 <- readNWISdata(sites = site, service = 'iv', parameterCd = params,
                            startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Taku_data20 <- data_org(Taku_data20)

### 2021 ####
bounds <- as.POSIXct(c('04/28/2021 00:00:00','10/01/2021 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Taku_data21 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, 
                            startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Taku_data21 <- data_org(Taku_data21)

ggplot(Taku_data19, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled*10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Taku R. 2019')+
  theme_cust()

ggplot(Taku_data20, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled*10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Taku R. 2020')+
  theme_cust()

ggplot(Taku_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled*10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Taku R. 2021')+
  theme_cust()

### Clustering
## prep data

# create daily ts
T_daily19 <- daily(Taku_data19)

# normalize data, convert to xts
T_norm19 <- daily_norm(T_daily19) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(T_norm19[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
T_memb_break19 <- membership_breaks(memb, T_daily19)
T_memb_break19

## visualize results
# c-q plot
ggplot(T_daily19, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(T_daily19, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

########## 2020 ############
## prep data
# create daily ts
T_daily20 <- daily(Taku_data20) 

# normalize data, convert to xts
T_norm20 <- daily_norm(T_daily20) 

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(T_norm20[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))


# find membership breaks
T_memb_break20 <- membership_breaks(memb, T_daily20)
T_memb_break20

## visualize results
# c-q plot
ggplot(T_daily20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(T_daily20, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()



########## 2021 ############
## prep data
# create daily ts
T_daily21 <- daily(Taku_data21)

# normalize data, convert to xts
T_norm21 <- daily_norm(T_daily21) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(T_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))


# find membership breaks
T_memb_break21 <- membership_breaks(memb, T_daily21)
T_memb_break21

## visualize results
# c-q plot
ggplot(T_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(T_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

# Writing Taku data to .csv files
if (writecsv == 1){
  write.csv(T_daily19,"Taku_hr19.csv", row.names = FALSE)
  write.csv(T_daily20,"Taku_hr20.csv", row.names = FALSE)
  write.csv(T_daily21,"Taku_hr21.csv", row.names = FALSE)
}

############ Stikine ####################
site <- '15024800' #Stikine R. Nr Wrangell
params <- c('00060', '00095') # discharge and SC

### 2020 ####
bounds <- as.POSIXct(c('04/15/2020 00:00:00','10/01/2020 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Stikine_data20 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, 
                            startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Stikine_data20 <- data_org(Stikine_data20)

ggplot(Stikine_data20, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled*10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Stikine R. 2020')+
  theme_cust()

########## 2020 ############
## prep data
# create daily ts
S_daily20 <- daily(Stikine_data20) 

# normalize data, convert to xts
S_norm20 <- daily_norm(S_daily20) 

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(S_norm20[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
S_memb_break20 <- membership_breaks(memb, S_daily20)
S_memb_break20

## visualize results
# c-q plot
ggplot(S_daily20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(S_daily20, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv == 1){
  write.csv(S_daily20,"Stikine_hr20.csv", row.names = FALSE)
}

############ Unuk ####################
site <- '15015595' #Unuk R. Bl Blue R. Nr Wrangell
params <- c('00060', '00095') # discharge and SC

### 2021 ####
bounds <- as.POSIXct(c('04/15/2021 00:00:00','10/01/2021 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Unuk_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, 
                           startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Unuk_data21<- data_org(Unuk_data21) 

ggplot(Unuk_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled*10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Unuk R. 2021')+
  theme_cust()

########## 2021 ############
## prep data
# create daily ts
U_daily21 <- daily(Unuk_data21) 

# normalize data, convert to xts
U_norm21 <- daily_norm(U_daily21) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(U_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
U_memb_break21 <- membership_breaks(memb, U_daily21)
U_memb_break21

## visualize results
# c-q plot
ggplot(U_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(U_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv == 1){
  write.csv(U_daily21,"Unuk_hr21.csv", row.names = FALSE)
}

############ Salmon #################### Too much missing data? Both years
site <- '15008000' #Salmon R. Nr Hyder
params <- c('00060', '00095') # discharge and SC
### 2021 ####
bounds <- as.POSIXct(c('04/15/2020 00:00:00','10/01/2020 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Salmon_data20<- readNWISdata(sites = site, service = 'iv', parameterCd = params,
                             startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Salmon_data20<- data_org(Salmon_data20)

ggplot(Salmon_data20, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled*10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Salmon R. 2020')+
  theme_cust()

## prep data
# create daily ts
Sa_daily20 <- daily(Salmon_data20) 

# normalize data, convert to xts
Sa_norm20 <- daily_norm(Sa_daily20) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Sa_norm20[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Sa_memb_break20 <- membership_breaks(memb, Sa_daily20)

## visualize results
# c-q plot
ggplot(Sa_daily20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sa_daily20, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()



### 2021 ####
bounds <- as.POSIXct(c('05/13/2021 00:00:00','09/18/2021 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Salmon_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                           startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) 
Salmon_data21<- data_org(Salmon_data21)

ggplot(Salmon_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Salmon R. 2021')+
  theme_cust()

## prep data
# create daily ts
Sa_daily21 <- daily(Salmon_data21) 

# normalize data, convert to xts
Sa_norm21 <- daily_norm(Sa_daily21)

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Sa_norm21[,-1], method = "euclidean")

# clustering
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,k=5) 

# find membership breaks
Sa_memb_break21 <- membership_breaks(memb, Sa_daily21)

## visualize results
# c-q plot
ggplot(Sa_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sa_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv==1) {
  write.csv(Sa_daily20,"Salmon_hr20.csv", row.names = FALSE)
  write.csv(Sa_daily21,"Salmon_hr21.csv", row.names = FALSE)
}

############ Kennicott ####################
site <- '15209700' #Kennicott R. at McCarthy
params <- c('00060', '00095') # discharge and SC
### 2020 ####
bounds <- as.POSIXct(c('05/24/2020 00:00:00','09/04/2020 00:00:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

Kennicott_data20<- readNWISdata(sites = site, service = 'iv', parameterCd = params,
                             startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2]))
Kennicott_data20<- data_org(Kennicott_data20)

ggplot(Kennicott_data20, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_filled), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Kennicott R. 2020')+
  theme_cust()

########## 2020 ############
## prep data
# create daily ts
K_daily20 <- daily(Kennicott_data20)

# normalize data, convert to xts
K_norm20 <- daily_norm(K_daily20)

## run depth constrained clustering
# distance matrix
dep_con_dist <- dist(K_norm20[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
K_memb_break20 <- membership_breaks(memb, K_daily20)

## visualize results
# c-q plot
ggplot(K_daily20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(K_daily20, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv == 1){
  write.csv(K_daily20,"Kennicott_hr20.csv", row.names = FALSE)
}

############### Hobo Data #########################

site <- '15272502' #Glacier C. at Alyeska 
params <- '00060' # discharge 

# load EC data from Hobo Loggers
Glacier21 <- read.csv('data/Glacier_C.csv')  
Glacier21$datetime <- as.POSIXct(Glacier21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Glacier21$datetime <- as.POSIXct(format(Glacier21$datetime), tz = "America/Anchorage")

Glacier22 <- read.csv('data/Glacier_C22.csv')  
Glacier22$datetime <- as.POSIXct(Glacier22$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Glacier22$datetime <- as.POSIXct(format(Glacier22$datetime), tz = "America/Anchorage")


bounds <- as.POSIXct(c(first(Glacier21$datetime),last(Glacier21$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Glacier_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, 
        startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
        select(datetime = dateTime, Q = X_00060_00000) %>%
        mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Glacier21, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
         mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

bounds <- as.POSIXct(c(first(Glacier22$datetime),last(Glacier22$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Glacier_data22<- readNWISdata(sites = site, service = 'iv', parameterCd = params, 
                              startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Glacier22, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

ggplot(Glacier_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Glacier C. 2021')+
  theme_cust()

ggplot(Glacier_data22, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Glacier C. 2022')+
  theme_cust()

## prep data
# create daily ts
Gl_daily21 <- daily(Glacier_data21) 

# normalize data, convert to xts
Gl_norm21 <- daily_norm(Gl_daily21) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Gl_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Gl_memb_break21 <- membership_breaks(memb, Gl_daily21)


## visualize results
# c-q plot
ggplot(Gl_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Gl_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()
### 2022 #####3
## prep data
# create daily ts
Gl_daily22 <- daily(Glacier_data22) 

# normalize data, convert to xts
Gl_norm22 <- daily_norm(Gl_daily22) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Gl_norm22[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Gl_memb_break22 <- membership_breaks(memb, Gl_daily22)


## visualize results
# c-q plot
ggplot(Gl_daily22, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Gl_daily22, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

####### Snow R. #########
site <- '15243900' #Snow R. nr Seward AK 
params <- '00060' # discharge 

Snow21 <- read.csv('data/Snow_R.csv')
Snow21$datetime <- as.POSIXct(Snow21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Snow21$datetime <- as.POSIXct(format(Snow21$datetime), tz = "America/Anchorage")

Snow22 <- read.csv('data/Snow_R22.csv')
Snow22$datetime <- as.POSIXct(Snow22$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Snow22$datetime <- as.POSIXct(format(Snow22$datetime), tz = "America/Anchorage")

bounds <- as.POSIXct(c(first(Snow21$datetime),last(Snow21$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Snow_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params,
                              startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Snow21, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

bounds <- as.POSIXct(c(first(Snow22$datetime),last(Snow22$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Snow_data22<- readNWISdata(sites = site, service = 'iv', parameterCd = params,
                           startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Snow22, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

ggplot(Snow_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Snow R. 2021')+
  theme_cust()

ggplot(Snow_data22, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Snow R. 2022')+
  theme_cust()


## prep data
# create daily ts
Sn_daily21 <- daily(Snow_data21)

# normalize data, convert to xts
Sn_norm21 <- daily_norm(Sn_daily21) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Sn_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Sn_memb_break21 <- membership_breaks(memb, Sn_daily21)

## visualize results
# c-q plot
ggplot(Sn_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sn_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

## prep data
# create daily ts
Sn_daily22 <- daily(Snow_data22)

# normalize data, convert to xts
Sn_norm22 <- daily_norm(Sn_daily22) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Sn_norm22[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Sn_memb_break22 <- membership_breaks(memb, Sn_daily22)

## visualize results
# c-q plot
ggplot(Sn_daily22, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sn_daily22, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

####### Knik R. #########
knik21 <- read.csv('data/Knick.csv')
knik21$datetime <- as.POSIXct(knik21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
knik21$datetime <- as.POSIXct(format(knik21$datetime), tz = "America/Anchorage")

site <- '15281000' #Knik R. nr palmer AK 
params <- '00060' # discharge and SC
bounds <- as.POSIXct(c(first(Snow21$datetime),last(Snow21$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Knik_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                           startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,knik21, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

ggplot(Knik_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Knik R. 2021')+
  theme_cust()

## prep data
# create daily ts
kn_daily21 <- daily(Knik_data21)

# normalize data, convert to xts
Kn_norm21 <- daily_norm(kn_daily21)

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(Kn_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Kn_memb_break21 <- membership_breaks(memb, kn_daily21)

## visualize results
# c-q plot
ggplot(kn_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(kn_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


####### Kenai R. #########
Kenai21 <- read.csv('data/Kenai_R.csv')
Kenai21$datetime <- as.POSIXct(Kenai21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Kenai21$datetime <- as.POSIXct(format(Kenai21$datetime), tz = "America/Anchorage")

site <- '15258000' #Kenai R. at Cooper Landing AK 
params <- '00060' # discharge 
bounds <- as.POSIXct(c(first(Kenai21$datetime),last(Kenai21$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Kenai_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                           startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Kenai21, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

ggplot(Kenai_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Kenai R. 2021')+
  theme_cust()

## prep data
# create daily ts
ken_daily21 <- daily(Kenai_data21) 

# normalize data, convert to xts
ken_norm21 <- daily_norm(ken_daily21)

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(ken_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset(bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff)))

# find membership breaks
Ken_memb_break21 <- membership_breaks(memb, ken_daily21)

## visualize results
# c-q plot
ggplot(ken_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(ken_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


####### Sixmile C. #########
site <- '15271000' #Sixmile C. nr Hope AK 
params <- '00060' # discharge 

Sixmi21 <- read.csv('data/Sixmile_C.csv')
Sixmi21$datetime <- as.POSIXct(Sixmi21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Sixmi21$datetime <- as.POSIXct(format(Sixmi21$datetime), tz = "America/Anchorage")

Sixmi22 <- read.csv('data/Sixmile_C22.csv')
Sixmi22$datetime <- as.POSIXct(Sixmi22$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Sixmi22$datetime <- as.POSIXct(format(Sixmi22$datetime), tz = "America/Anchorage")

bounds <- as.POSIXct(c(first(Sixmi21$datetime),last(Sixmi21$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Sixmi_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                            startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Sixmi21, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

bounds <- as.POSIXct(c(first(Sixmi22$datetime),last(Sixmi22$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

Sixmi_data22<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                            startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,Sixmi22, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

ggplot(Sixmi_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = Sp_cond), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Sixmi C. 2021')+
  theme_cust()

ggplot(Sixmi_data22, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = Sp_cond), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Sixmi C. 2022')+
  theme_cust()

## prep data
# create daily ts
sx_daily21 <- daily(Sixmi_data21)

# normalize data, convert to xts
sx_norm21 <- daily_norm(sx_daily21) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(sx_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
Sxm_memb_break21 <- membership_breaks(memb, sx_daily21)

## visualize results
# c-q plot
ggplot(sx_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(sx_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

## prep data
# create daily ts
sx_daily22 <- daily(Sixmi_data22)

# normalize data, convert to xts
sx_norm22 <- daily_norm(sx_daily22) 

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(sx_norm22[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
Sxm_memb_break22 <- membership_breaks(memb, sx_daily22)

## visualize results
# c-q plot
ggplot(sx_daily22, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(sx_daily22, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


####### Little Susitna C. #########

site <- '15290000' #L Susitna R. nr Palmer AK 
params <- '00060' # discharge 

LSus21 <- read.csv('data/L_Sus.csv')
LSus21$datetime <- as.POSIXct(LSus21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
LSus21$datetime <- as.POSIXct(format(LSus21$datetime), tz = "America/Anchorage")

LSus22 <- read.csv('data/L_Sus22.csv')
LSus22$datetime <- as.POSIXct(LSus22$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
LSus22$datetime <- as.POSIXct(format(LSus22$datetime), tz = "America/Anchorage")

bounds <- as.POSIXct(c(first(LSus21$datetime),last(LSus21$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

LSus_data21<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                            startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,LSus21, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

bounds <- as.POSIXct(c(first(LSus22$datetime),last(LSus22$datetime)),format="%m/%d/%y %H:%M", TZ = "America/Anchorage")

LSus_data22<- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                           startDate = as.Date(bounds[1]), endDate = as.Date(bounds[2])) %>%
  select(datetime = dateTime, Q = X_00060_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'))%>%
  full_join(.,LSus22, by = 'datetime') %>%
  select(datetime, Q, Sp_cond, temp)%>%
  mutate(Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(Sp_cond, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(Sp_cond,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')

ggplot(LSus_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Little Susitna. 2021')+
  theme_cust()

ggplot(LSus_data22, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Little Susitna. 2022')+
  theme_cust()

## prep data
# create daily ts
ls_daily21 <- daily(LSus_data21)

# normalize data, convert to xts
ls_norm21 <- daily_norm(ls_daily21)

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(ls_norm21[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
LSmemb_break21 <- membership_breaks(memb, ls_daily21)


## visualize results
# c-q plot
ggplot(ls_daily21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(ls_daily21, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

## prep data
# create daily ts
ls_daily22 <- daily(LSus_data22)

# normalize data, convert to xts
ls_norm22 <- daily_norm(ls_daily22)

## run depth constrained clusterin
# distance matrix
dep_con_dist <- dist(ls_norm22[,-1], method = "euclidean")

# clustering
dep_clust <- chclust(dep_con_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstk <-bstick(dep_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep_clust,min(subset( bstk$nGroups, bstk$dispersion-bstk$bstick < b_cutoff))) 

# find membership breaks
LSmemb_break21 <- membership_breaks(memb, ls_daily22)


## visualize results
# c-q plot
ggplot(ls_daily22, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(ls_daily22, aes(x= ymd(datetime_daily), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


write.csv(ls_daily21,"l_sus_hr21.csv", row.names = FALSE)

##### visualize membership breaks #####
max_length <- max(c(length(T_memb_break19), length(T_memb_break20), length(T_memb_break21), length(S_memb_break20), 
                    length(U_memb_break21), length(Sa_memb_break20), length(Sa_memb_break21), length(K_memb_break20), 
                    length(Gl_memb_break21), length(Sn_memb_break21), length(Kn_memb_break21), length(Ken_memb_break21), length(Sxm_memb_break21)))

memb_break_dfall <- data.frame(col1 = c(as.integer(strftime(T_memb_break19, '%j')),                 # Create data frame with unequal vectors
                                     rep(NA, max_length - length(T_memb_break19))),
                            col2 = c(as.integer(strftime(T_memb_break20, '%j')),
                                     rep(NA, max_length - length(T_memb_break20))),
                            col3 = c(as.integer(strftime(T_memb_break21, '%j')),
                                     rep(NA, max_length - length(T_memb_break21))),
                            col4 = c(as.integer(strftime(S_memb_break20, '%j')),
                                     rep(NA, max_length - length(S_memb_break20))),
                            col5 = c(as.integer(strftime(U_memb_break21, '%j')),
                                     rep(NA, max_length - length(U_memb_break21))),
                            col6 = c(as.integer(strftime(Sa_memb_break20, '%j')),
                                     rep(NA, max_length - length(Sa_memb_break20))),
                            col7 = c(as.integer(strftime(Sa_memb_break21, '%j')),
                                     rep(NA, max_length - length(Sa_memb_break21))),
                            col8 = c(as.integer(strftime(K_memb_break20, '%j')),
                                     rep(NA, max_length - length(K_memb_break20))),
                            col9 = c(as.integer(strftime(Gl_memb_break21, '%j')),
                                     rep(NA, max_length - length(Gl_memb_break21))),
                            col10 = c(as.integer(strftime(Sn_memb_break21, '%j')),
                                     rep(NA, max_length - length(Sn_memb_break21))),
                            col11 = c(as.integer(strftime(Kn_memb_break21, '%j')),
                                     rep(NA, max_length - length(Kn_memb_break21))),
                            col12 = c(as.integer(strftime(Ken_memb_break21, '%j')),
                                     rep(NA, max_length - length(Ken_memb_break21))),
                            col13 = c(as.integer(strftime(Sxm_memb_break21, '%j')),
                                     rep(NA, max_length - length(Sxm_memb_break21))))
colnames(memb_break_dfall)<-c("Taku19","Taku20","Taku21","Stikine20","Unuk21","Salmon20","Salmon21","Kennicott20","Glacier21","snow21","Knick21", "Kenai21","Sixmile21")

memb_break_dfall <- tibble(memb_break_dfall)%>%
  rowid_to_column('break_num') %>%
  pivot_longer(cols = -break_num, names_to = 'year', values_to = 'doy') %>%
  mutate(year = year)

ggplot(memb_break_dfall, aes(x = doy, y = year, color = as.factor(break_num))) +
  geom_point()

write.csv(memb_break_df,"combinedUSGS_clusters.csv", row.names = FALSE)

if (writecsv ==1){
  write.csv(Gl_daily21,"Glacier_hr21.csv", row.names = FALSE)
  write.csv(Sn_daily21,"Snow_hr21.csv", row.names = FALSE)
  write.csv(kn_daily21,"Knik_hr21.csv", row.names = FALSE)
  write.csv(ken_daily21,"Kenai_hr21.csv", row.names = FALSE)
  write.csv(sx_daily21,"Sixmi_hr21.csv", row.names = FALSE)
  write.csv(ls_daily21,"l_sus_hr21.csv", row.names = FALSE)
}

############## Eran Hood's Data ####################
CoweeQ18 <- read.csv('data/Cowee 15 minute 2018 cms.csv', header=TRUE, skip=2) %>%
  rename(Q_m3s = Value)

CoweeQ19 <- read.csv('data/Cowee Cr 15 min discharge cms 2019.csv', header=TRUE, skip=2) %>%
  rename(Q_m3s = Value)

CoweeEC18 <-read.csv('data/CoweeCr_EC_2018.csv') %>%
  rename(datetime = Date) 
CoweeEC18$datetime = mdy_hms(paste(CoweeEC18$datetime, CoweeEC18$Time))
CoweeEC18 <- CoweeEC18[-c(2)]

CoweeEC19 <-read.csv('data/CoweeCr_EC_2019.csv') %>%
  rename(datetime = Date) 
CoweeEC19$datetime = mdy_hms(paste(CoweeEC19$datetime, CoweeEC19$Time))
CoweeEC19 <- CoweeEC19[-c(2)]

HerbertQ18 <- read.csv('data/Herber 15 min Q cfs 2018.csv', header=TRUE, skip=2) %>%
rename(Q_m3s = Value)
HerbertQ18$Q_m3s <- HerbertQ18$Q_m3s/35.315

HerbertQ19 <- read.csv('data/Herbert 15 min discharge cfs 2019.csv', header=TRUE, skip=2) %>%
rename(Q_m3s = Value)
HerbertQ19$Q_m3s <- HerbertQ19$Q_m3s/35.315
  
HerbertEC18 <-read.csv('data/HerbertR_EC_2018.csv')%>%
  rename(datetime = Date) 
HerbertEC18$datetime = mdy_hms(paste(HerbertEC18$datetime, HerbertEC18$Time))
HerbertEC18 <- HerbertEC18[-c(2)]
HerbertEC18$datetime <- round_date(HerbertEC18$datetime, "hour")

HerbertEC19 <-read.csv('data/HerbertR_EC_2019.csv')%>%
  rename(datetime = Date) 
HerbertEC19$datetime = mdy_hms(paste(HerbertEC19$datetime, HerbertEC19$Time))
HerbertEC19 <- HerbertEC19[-c(2)]
HerbertEC19$datetime <- round_date(HerbertEC19$datetime, "hour")
  
################################ Extra Code ############################## 

 ##### Ca time series #########
 # ca_slope<- 0.1655
 # ca_int <- 1.6576
 # 
 # daily_data17$Ca_ts <- (ca_slope* daily_data17$SC_filled + ca_int) *1000* 86400* daily_data17$Q_m3sdaily /1e6
 # 
 # ggplot(new_df17, aes(x= date , y= daily_data17$Ca_ts)) +geom_point(aes(color= factor(memb)))+
 #   #scale_color_discrete(drop=FALSE)+
 #   scale_colour_brewer(palette = "Dark2")+
 #   ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
 #   xlab(NULL) + 
 #   theme_cust()