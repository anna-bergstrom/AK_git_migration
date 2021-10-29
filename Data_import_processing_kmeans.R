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

# load and configure all data sets 
#full_data16 <- read.csv('data/WG_chemsitry_2016.csv') #2016 chemistry - right now this is not used 
#full_data17 <- read.csv('data/WG_chem_2017.csv') #2017 chemistry - right now this is not used 

## set bounds for 2016 and 2017
bounds16 <- as.POSIXct(c('04/28/2016 00:00:00','10/06/2016 23:45:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")
bounds17 <- as.POSIXct(c('04/20/2017 14:00:00','10/17/2017 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage")

##### nwis data #####
## set nwis details
site <- '15236900'
params <- c('00060', '00095')

## set smoothing details 
swindow <- 192

### loading Wolverine gauge data from 2016-2020
### add column of q in cms
### add index column
### interpolate to fill NAs in '_filled' ts columns
### add smoothed data to '_smoothed' ts columns
### NOTE: data does not need to be trimmed, just set bounds as start/end in NWIS call
## 2016
gauge_data16 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = as.Date(bounds16[1]), endDate = as.Date(bounds16[2])) %>%
   select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
   mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
          Q_m3s = Q*0.028316847,
          period = NA,
          Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
          SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
          Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
          SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
   rowid_to_column('count')
# assign manual breaks
breaks16 <- c(mdy_hms('05-19-2016 23:45:00', tz = 'America/Anchorage'),
              mdy_hms('06-13-2016 23:45:00', tz = 'America/Anchorage'),
              mdy_hms('07-16-2016 23:45:00', tz = 'America/Anchorage'))
gauge_data16$period[gauge_data16$datetime < breaks16[1]] <- "p1"
gauge_data16$period[gauge_data16$datetime >= breaks16[1] & gauge_data16$datetime < breaks16[2]] <- "p2"
gauge_data16$period[gauge_data16$datetime >= breaks16[2] & gauge_data16$datetime < breaks16[3]] <- "p3"
gauge_data16$period[gauge_data16$datetime >= breaks16[3]] <- "p4"

## 2017 
gauge_data17 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = as.Date(bounds17[1]), endDate = as.Date(bounds17[2])) %>%
   select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
   mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
          Q_m3s = Q*0.028316847,
          period = NA,
          Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
          SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
          Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
          SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
   rowid_to_column('count')
# assign manual breaks
breaks17 <- c(mdy_hms('05/12/2017 23:45:00', tz = 'America/Anchorage'),
              mdy_hms('06/22/2017 23:45:00', tz = 'America/Anchorage'),
              mdy_hms('08/19/2017 23:45:00', tz = 'America/Anchorage'))
gauge_data17$period[gauge_data17$datetime < breaks17[1]] <- "p1"
gauge_data17$period[gauge_data17$datetime >= breaks17[1] & gauge_data17$datetime < breaks17[2] ] <- "p2"
gauge_data17$period[gauge_data17$datetime >= breaks17[2] & gauge_data17$datetime < breaks17[3] ] <- "p3"
gauge_data17$period[gauge_data17$datetime >= breaks17[3]] <- "p4"

## 2018
gauge_data18 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = '2018-05-22', endDate = '2018-10-17') %>%
   select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
   mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
          Q_m3s = Q*0.028316847,
          Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
          SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
          Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
          SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
   rowid_to_column('count')

## 2019
gauge_data19 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using start/end date from original for now
                             startDate = '2019-04-30', endDate = '2019-10-12') %>%
   select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
   mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
          Q_m3s = Q*0.028316847,
          Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
          SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
          Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
          SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
   rowid_to_column('count')

## 2020
gauge_data20 <- readNWISdata(sites = site, service = 'iv', parameterCd = params, # using max and min dates from other years for now
                             startDate = '2020-04-18', endDate = '2020-10-18') %>%
   select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
   mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
         Q_m3s = Q*0.028316847,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
   rowid_to_column('count')

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

# 2019 (data only partially available)
 precip19 <- full_wx %>%
   dplyr::filter(local_time >= min(gauge_data19$datetime),
                 local_time <= max(gauge_data19$datetime)) %>%
   select(precip_int = TPGIncremental,
          datetime = local_time, 
          temp = site_temp)

# 2020 (data currently not available)
 precip20 <- full_wx %>%
   dplyr::filter(local_time >= min(gauge_data20$datetime),
                 local_time <= max(gauge_data20$datetime)) %>%
   select(precip_int = TPGIncremental,
          datetime = local_time, 
          temp = site_temp)

##### Working with the time series data #####

### Recreating Figure 2 plots from the geochem paper (in review as of Aug. '21)
# generating a custom theme to get rid of the ugly ggplot defaults 
theme_cust <- function(base_size = 11, base_family = "") {
  theme_classic() %+replace%
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text = element_text(color = "black")
    )
}

# 2016 EC-Q plot
p16<-ggplot(gauge_data16, aes(x= Q_m3s, y= SC)) +geom_point(aes(col=period), size=2 )+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
 scale_colour_brewer(palette = "Dark2")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
 theme(legend.position="None")
p16

# 2017 EC-Q plot
p17<- ggplot(gauge_data17, aes(x= Q_m3s, y= SC)) +geom_point(aes(col=period), size=2 )+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_brewer(palette = "Dark2")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position= c(0.9,0.6))+
  theme(legend.title = element_blank())
p17

# 2018 EC-Q plot
p18<- ggplot(gauge_data18, aes(x= Q_m3s, y= SC)) +geom_point(aes(color= count))+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_gradient(low= "green", high="blue")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position="None") 
p18

# 2019 EC-Q plot
p19<- ggplot(gauge_data19, aes(x= Q_m3s, y= SC)) +geom_point(aes(color= count))+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_gradient( name = NULL, breaks= c(130, 11865),labels=c("May","Sep"),low= "green", high="blue" )+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position=c(0.89,0.6))
p19

# 2020 EC-Q plot
p20<- ggplot(gauge_data20, aes(x= Q_m3s, y= SC)) +geom_point(aes(color= count))+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_gradient(low= "green", high="blue" )+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position="None")
p20

# putting all plots in a multipanel figure
figure <- ggarrange(p16, p17, p18,p19,p20,
                    labels = c("2016", "2017", "2018", "2019", "2020"),
                    ncol = 3, nrow = 2, hjust= -1.6, vjust = 1.8)
figure

### Time series plots
# 2016
ggplot(gauge_data16, aes(x = datetime, y = Q_smoothed)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed/10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = '2016')+
  theme_cust()

# 2017
ggplot(gauge_data17, aes(x = datetime, y = Q_smoothed)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed/10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = '2017')+
  theme_cust()

# 2018
ggplot(gauge_data18, aes(x = datetime, y = Q_smoothed)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed/10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = '2018')+
  theme_cust()

# 2019
ggplot(gauge_data19, aes(x = datetime, y = Q_smoothed)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed/10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = '2019')+
  theme_cust()

# 2020
ggplot(gauge_data20, aes(x = datetime, y = Q_smoothed)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed/10), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = '2020')+
  theme_cust()
 
##### Depth-constrained clustering (to constrain clusters in time) #####
#Tried running it with 15 min data and it was taking too long (may try overnight, but not super practical)
#So moving to try it with 1-hour intervals - averaging the Q and EC, and summing precip

### 2016
## prep data
# merge nwis and met data
merged16 <- gauge_data16 %>%
  full_join(.,precip16, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged16$precip_int[is.nan(merged16$precip_int)] <- NA
merged16$precip_int <- na_interpolation(merged16$precip_int, maxgap = Inf, option = 'linear')

# create hourly ts
hourly16 <- merged16 %>%
  mutate(datetime_hourly = cut(datetime, 'hour')) %>% # create hour aggregated datetime
  group_by(datetime_hourly) %>% 
  summarise(SC_mean = mean(SC_filled),
            Q_mean = mean(Q_filled),
            precip_sum = sum(precip_int)) %>%
  na.omit()

# normalize data, convert to xts
norm16 <- hourly16 %>%
  mutate(datetime_hourly = ymd_hms(datetime_hourly), 
         precip_norm = (1-max(na.omit(precip_sum))-precip_sum)/(max(na.omit(precip_sum))-min(na.omit(precip_sum))),
         SC_norm = (1-max(na.omit(SC_mean))-SC_mean)/(max(na.omit(SC_mean))-min(na.omit(SC_mean))),
         Q_norm = (1-max(na.omit(Q_mean))-Q_mean)/(max(na.omit(Q_mean))-min(na.omit(Q_mean)))) %>%
  select(datetime_hourly, precip_norm, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_hourly))

## run depth constrained clusterin
# distance matrix
dep_con16_dist <- dist(norm16[,-1], method = "euclidean")

# clustering
dep16_clust <- chclust(dep_con16_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstick(dep16_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep16_clust,k=5) 

# define function to find membership breaks
membership_breaks <- function(memb, hourly){
  for (k  in 1:length(unique(memb))) {
    temp <- hourly$datetime_hourly[which(memb==k)]
    if (k== 1){
      p_time <- as.POSIXct(tail(temp,n=1)) 
    } else{
      p_time [k] <- tail(temp, n=1)}
  }
  return(p_time)
}

# find membership breaks
memb_break16 <- membership_breaks(memb, hourly16)
memb_break16

## visualize results
# c-q plot
ggplot(hourly16, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Dark2")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()
 
# q ts
ggplot(hourly16, aes(x= ymd_hms(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Dark2")+
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
 
# create hourly ts
hourly17 <- merged17 %>%
  mutate(datetime_hourly = cut(datetime, 'hour')) %>% # create hour aggregated datetime
  group_by(datetime_hourly) %>% 
  summarise(SC_mean = mean(SC_filled),
             Q_mean = mean(Q_filled),
             precip_sum = sum(precip_int)) %>%
   na.omit()
 
# normalize data, convert to xts
norm17 <- hourly17 %>%
  mutate(datetime_hourly = ymd_hms(datetime_hourly), 
          precip_norm = (1-max(na.omit(precip_sum))-precip_sum)/(max(na.omit(precip_sum))-min(na.omit(precip_sum))),
          SC_norm = (1-max(na.omit(SC_mean))-SC_mean)/(max(na.omit(SC_mean))-min(na.omit(SC_mean))),
          Q_norm = (1-max(na.omit(Q_mean))-Q_mean)/(max(na.omit(Q_mean))-min(na.omit(Q_mean)))) %>%
  select(datetime_hourly, precip_norm, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_hourly))
 
## run depth constrained clustering
# distance matrix
dep_con17_dist <- dist(norm17[,-1], method = "euclidean")
 
# clustering
dep17_clust <- chclust(dep_con17_dist, method = "coniss")
 
# broken stick plot showing reduction in sse over n of clusters
bstick(dep17_clust, ng=20, plot=TRUE)
 
# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep17_clust,k=5)

# find membership breaks
memb_break17 <- membership_breaks(memb, hourly17)
memb_break17
 
## visualize results
# c-q plot
ggplot(hourly17, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Dark2")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()
 
# q ts
ggplot(hourly17, aes(x= ymd_hms(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Dark2")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
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

# create hourly ts
hourly18 <- merged18 %>%
  mutate(datetime_hourly = cut(datetime, 'hour')) %>% # create hour aggregated datetime
  group_by(datetime_hourly) %>% 
  summarise(SC_mean = mean(SC_filled),
            Q_mean = mean(Q_filled),
            precip_sum = sum(precip_int)) %>%
  na.omit()

# normalize data, convert to xts
norm18 <- hourly18 %>%
  mutate(datetime_hourly = ymd_hms(datetime_hourly), 
         precip_norm = (1-max(na.omit(precip_sum))-precip_sum)/(max(na.omit(precip_sum))-min(na.omit(precip_sum))),
         SC_norm = (1-max(na.omit(SC_mean))-SC_mean)/(max(na.omit(SC_mean))-min(na.omit(SC_mean))),
         Q_norm = (1-max(na.omit(Q_mean))-Q_mean)/(max(na.omit(Q_mean))-min(na.omit(Q_mean)))) %>%
  select(datetime_hourly, precip_norm, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_hourly))

## run depth constrained clustering
# distance matrix
dep_con18_dist <- dist(norm18[,-1], method = "euclidean")

# clustering
dep18_clust <- chclust(dep_con18_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstick(dep18_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep18_clust,k=4)

# find membership breaks
memb_break18 <- membership_breaks(memb, hourly18)
memb_break18

## visualize results
# c-q plot
ggplot(hourly18, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Dark2")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(hourly18, aes(x= ymd_hms(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Dark2")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
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

# create hourly ts
hourly19 <- merged19 %>%
  mutate(datetime_hourly = cut(datetime, 'hour')) %>% # create hour aggregated datetime
  group_by(datetime_hourly) %>% 
  summarise(SC_mean = mean(SC_filled),
            Q_mean = mean(Q_filled),
            precip_sum = sum(precip_int)) %>%
  na.omit()

# normalize data, convert to xts
norm19 <- hourly19 %>%
  mutate(datetime_hourly = ymd_hms(datetime_hourly), 
         precip_norm = (1-max(na.omit(precip_sum))-precip_sum)/(max(na.omit(precip_sum))-min(na.omit(precip_sum))),
         SC_norm = (1-max(na.omit(SC_mean))-SC_mean)/(max(na.omit(SC_mean))-min(na.omit(SC_mean))),
         Q_norm = (1-max(na.omit(Q_mean))-Q_mean)/(max(na.omit(Q_mean))-min(na.omit(Q_mean)))) %>%
  select(datetime_hourly, precip_norm, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_hourly))

## run depth constrained clusterin
# distance matrix
dep_con19_dist <- dist(norm19[,-1], method = "euclidean")

# clustering
dep19_clust <- chclust(dep_con19_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstick(dep19_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep19_clust,k=5) 

# define function to find membership breaks
membership_breaks <- function(memb, hourly){
  for (k  in 1:length(unique(memb))) {
    temp <- hourly$datetime_hourly[which(memb==k)]
    if (k== 1){
      p_time <- as.POSIXct(tail(temp,n=1)) 
    } else{
      p_time [k] <- tail(temp, n=1)}
  }
  return(p_time)
}

# find membership breaks
memb_break19 <- membership_breaks(memb, hourly19)
memb_break19

## visualize results
# c-q plot
ggplot(hourly19, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Dark2")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(hourly19, aes(x= ymd_hms(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Dark2")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
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

# create hourly ts
hourly20 <- merged20 %>%
  mutate(datetime_hourly = cut(datetime, 'hour')) %>% # create hour aggregated datetime
  group_by(datetime_hourly) %>% 
  summarise(SC_mean = mean(SC_filled),
            Q_mean = mean(Q_filled),
            precip_sum = sum(precip_int)) %>%
  na.omit()

# normalize data, convert to xts
norm20 <- hourly20 %>%
  mutate(datetime_hourly = ymd_hms(datetime_hourly), 
         precip_norm = (1-max(na.omit(precip_sum))-precip_sum)/(max(na.omit(precip_sum))-min(na.omit(precip_sum))),
         SC_norm = (1-max(na.omit(SC_mean))-SC_mean)/(max(na.omit(SC_mean))-min(na.omit(SC_mean))),
         Q_norm = (1-max(na.omit(Q_mean))-Q_mean)/(max(na.omit(Q_mean))-min(na.omit(Q_mean)))) %>%
  select(datetime_hourly, precip_norm, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_hourly))

## run depth constrained clusterin
# distance matrix
dep_con20_dist <- dist(norm20[,-1], method = "euclidean")

# clustering
dep20_clust <- chclust(dep_con20_dist, method = "coniss")

# broken stick plot showing reduction in sse over n of clusters
bstick(dep20_clust, ng=20, plot=TRUE)

# cut and assign membership
# set the number of clusters off at 5 and defining the membership of those 5 clusters
memb <- cutree(dep20_clust,k=5) 

# define function to find membership breaks
membership_breaks <- function(memb, hourly){
  for (k  in 1:length(unique(memb))) {
    temp <- hourly$datetime_hourly[which(memb==k)]
    if (k== 1){
      p_time <- as.POSIXct(tail(temp,n=1)) 
    } else{
      p_time [k] <- tail(temp, n=1)}
  }
  return(p_time)
}

# find membership breaks
memb_break20 <- membership_breaks(memb, hourly20)
memb_break20

## visualize results
# c-q plot
ggplot(hourly20, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Dark2")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(hourly20, aes(x= ymd_hms(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Dark2")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

##### visualize membership breaks #####
memb_break_df <- tibble('2016' = as.integer(strftime(memb_break16, '%j')),
                        '2017' = as.integer(strftime(memb_break17, '%j')),
                        '2018' = as.integer(strftime(memb_break18, '%j')), 
                        '2019' = as.integer(strftime(memb_break19, '%j')),
                        '2020' = as.integer(strftime(memb_break20, '%j')))%>%
  rowid_to_column('break_num') %>%
  pivot_longer(cols = -break_num, names_to = 'year', values_to = 'doy') %>%
  mutate(year = year)

ggplot(memb_break_df, aes(x = doy, y = year, color = as.factor(break_num))) +
  geom_point()

################################ Extra Code ############################## 
 

 
 ##### Ca time series #########
 # ca_slope<- 0.1655
 # ca_int <- 1.6576
 # 
 # hourly_data17$Ca_ts <- (ca_slope* hourly_data17$SC_filled + ca_int) *1000* 86400* hourly_data17$Q_m3shourly /1e6
 # 
 # ggplot(new_df17, aes(x= date , y= hourly_data17$Ca_ts)) +geom_point(aes(color= factor(memb)))+
 #   #scale_color_discrete(drop=FALSE)+
 #   scale_colour_brewer(palette = "Dark2")+
 #   ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
 #   xlab(NULL) + 
 #   theme_cust()