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

########### Setting up details for this script #############
# set smoothing details 
swindow <- 192

# option to write hourly data as .csv
# if this is set as 1, then it will write csvs for all years and all sites, any other number and it won't 
writecsv <- 0 

b_cutoff <- 0.5

# generating a custom theme to get rid of the ugly ggplot defaults 
theme_cust <- function(base_size = 11, base_family = "") {
  theme_classic() %+replace%
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text = element_text(color = "black")
    )
}

# NWIS data download organization function. 
data_org <- function(gage_data){
  gage_data %>%
  select(datetime = dateTime, Q = X_00060_00000, SC = X_00095_00000) %>%
  mutate(datetime = with_tz(datetime, tz = 'America/Anchorage'),
         Q_m3s = Q*0.028316847,
         period = NA,
         Q_filled = na_interpolation(Q_m3s, option = 'linear', maxgap = Inf),
         SC_filled = na_interpolation(SC, option = 'linear', maxgap = Inf),
         Q_smoothed = rollapply(Q_m3s,swindow,mean, na.rm = TRUE, fill = NA),
         SC_smoothed = rollapply(SC,swindow,mean, na.rm = TRUE, fill = NA)) %>%
  rowid_to_column('count')
}

# Define function to convert time series data to hourly
hourly<- function(merged){
  merged %>%
    mutate(datetime_hourly = cut(datetime, 'day')) %>%
    group_by(datetime_hourly) %>% 
    summarise(SC_mean = mean(SC_filled),
              Q_mean = mean(Q_filled)) %>%
    na.omit()
}

# Define function to normalize ts data and convert to an xts object 
hourly_norm <- function(hourly_D){
  hourly_D %>%
  mutate(datetime_hourly = ymd(datetime_hourly),   
       SC_norm = ((SC_mean)-min(na.omit(SC_mean)))/(max(na.omit(SC_mean))-min(na.omit(SC_mean))),
       Q_norm = (Q_mean-min(na.omit(Q_mean)))/(max(na.omit(Q_mean))-min(na.omit(Q_mean)))) %>%
  select(datetime_hourly, SC_norm, Q_norm) %>%
  xts(.[,-1], order.by = as.POSIXct(.$datetime_hourly))
}

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

##### Working with the time series data #####
 
##### Depth-constrained clustering (to constrain clusters in time) #####
# 1-hour intervals - averaging the Q and EC

### 2016
## prep data
# merge nwis and met data
merged16 <- gauge_data16 %>%
  full_join(.,precip16, by = 'datetime') %>%
  select(datetime, Q_filled, SC_filled, precip_int, temp)
merged16$precip_int[is.nan(merged16$precip_int)] <- NA
merged16$precip_int <- na_interpolation(merged16$precip_int, maxgap = Inf, option = 'linear')

# create hourly ts
hourly16 <- hourly(merged16)

# normalize data, convert to xts
norm16 <- hourly_norm(hourly16)
  
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
ggplot(hourly16, aes(x= ymd(datetime_hourly), y= Q_mean)) +
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
  hourly17 <- hourly(merged17)
  
# normalize data, convert to xts
  norm17 <- hourly_norm(hourly17)
 
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
ggplot(hourly17, aes(x= ymd(datetime_hourly), y= Q_mean)) +
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
hourly18 <- hourly(merged18)

# normalize data, convert to xts
norm18 <- hourly_norm(hourly18)

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
memb_break18 <- membership_breaks(memb, hourly18)
memb_break18

## visualize results
# c-q plot
ggplot(hourly18, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(hourly18, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
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
hourly19 <- hourly(merged19)

# normalize data, convert to xts
norm19 <- hourly_norm(hourly19)

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
memb_break19 <- membership_breaks(memb, hourly19)
memb_break19

## visualize results
# c-q plot
ggplot(hourly19, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(hourly19, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
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
hourly20 <- hourly(merged20)

# normalize data, convert to xts
norm20 <- hourly_norm(hourly20)

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
memb_break20 <- membership_breaks(memb, hourly20)
memb_break20

## visualize results
# c-q plot
ggplot(hourly20, aes(x=Q_mean , y= SC_mean)) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(hourly20, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

##### visualize membership breaks #####
max_length <- max(c(length(memb_break16), length(memb_break17), length(memb_break18), length(memb_break19), length(memb_break20)))
memb_break_df <- data.frame(col1 = c(as.integer(strftime(memb_break16, '%j')),                 # Create data frame with unequal vectors
                            rep(NA, max_length - length(memb_break16))),
                   col2 = c(as.integer(strftime(memb_break17, '%j')),
                            rep(NA, max_length - length(memb_break17))),
                   col3 = c(as.integer(strftime(memb_break18, '%j')),
                            rep(NA, max_length - length(memb_break18))),
                   col4 = c(as.integer(strftime(memb_break19, '%j')),
                            rep(NA, max_length - length(memb_break19))),
                   col5 = c(as.integer(strftime(memb_break20, '%j')),
                            rep(NA, max_length - length(memb_break20))))
colnames(memb_break_df)<-c(2016,2017,2018,2019,2020)
# memb_break_df$row_num <- seq.int(nrow(memb_break_df))
memb_break_df <- tibble(memb_break_df)%>%
  rowid_to_column('break_num') %>%
  pivot_longer(cols = -break_num, names_to = 'year', values_to = 'doy') %>%
  mutate(year = year)

ggplot(memb_break_df, aes(x = doy, y = year, color = as.factor(break_num))) +
 geom_point()

if (writecsv == 1){ 
  write.csv(hourly16,"hourly_16.csv", row.names = FALSE)
  write.csv(hourly17,"hourly_17.csv", row.names = FALSE)
  write.csv(hourly18,"hourly_18.csv", row.names = FALSE)
  write.csv(hourly19,"hourly_19.csv", row.names = FALSE)
  write.csv(hourly20,"hourly_20.csv", row.names = FALSE)
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

# create hourly ts
T_hourly19 <- hourly(Taku_data19)

# normalize data, convert to xts
T_norm19 <- hourly_norm(T_hourly19) 

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
T_memb_break19 <- membership_breaks(memb, T_hourly19)
T_memb_break19

## visualize results
# c-q plot
ggplot(T_hourly19, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(T_hourly19, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

########## 2020 ############
## prep data
# create hourly ts
T_hourly20 <- hourly(Taku_data20) 

# normalize data, convert to xts
T_norm20 <- hourly_norm(T_hourly20) 

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
T_memb_break20 <- membership_breaks(memb, T_hourly20)
T_memb_break20

## visualize results
# c-q plot
ggplot(T_hourly20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(T_hourly20, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()



########## 2021 ############
## prep data
# create hourly ts
T_hourly21 <- hourly(Taku_data21)

# normalize data, convert to xts
T_norm21 <- hourly_norm(T_hourly21) 

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
T_memb_break21 <- membership_breaks(memb, T_hourly21)
T_memb_break21

## visualize results
# c-q plot
ggplot(T_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(T_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

# Writing Taku data to .csv files
if (writecsv == 1){
  write.csv(T_hourly19,"Taku_hr19.csv", row.names = FALSE)
  write.csv(T_hourly20,"Taku_hr20.csv", row.names = FALSE)
  write.csv(T_hourly21,"Taku_hr21.csv", row.names = FALSE)
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
# create hourly ts
S_hourly20 <- hourly(Stikine_data20) 

# normalize data, convert to xts
S_norm20 <- hourly_norm(S_hourly20) 

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
S_memb_break20 <- membership_breaks(memb, S_hourly20)
S_memb_break20

## visualize results
# c-q plot
ggplot(S_hourly20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(S_hourly20, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv == 1){
  write.csv(S_hourly20,"Stikine_hr20.csv", row.names = FALSE)
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
# create hourly ts
U_hourly21 <- hourly(Unuk_data21) 

# normalize data, convert to xts
U_norm21 <- hourly_norm(U_hourly21) 

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
U_memb_break21 <- membership_breaks(memb, U_hourly21)
U_memb_break21

## visualize results
# c-q plot
ggplot(U_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(U_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv == 1){
  write.csv(U_hourly21,"Unuk_hr21.csv", row.names = FALSE)
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
# create hourly ts
Sa_hourly20 <- hourly(Salmon_data20) 

# normalize data, convert to xts
Sa_norm20 <- hourly_norm(Sa_hourly20) 

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
Sa_memb_break20 <- membership_breaks(memb, Sa_hourly20)

## visualize results
# c-q plot
ggplot(Sa_hourly20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sa_hourly20, aes(x= ymd(datetime_hourly), y= Q_mean)) +
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
# create hourly ts
Sa_hourly21 <- hourly(Salmon_data21) 

# normalize data, convert to xts
Sa_norm21 <- hourly_norm(Sa_hourly21)

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
Sa_memb_break21 <- membership_breaks(memb, Sa_hourly21)

## visualize results
# c-q plot
ggplot(Sa_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sa_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv==1) {
  write.csv(Sa_hourly20,"Salmon_hr20.csv", row.names = FALSE)
  write.csv(Sa_hourly21,"Salmon_hr21.csv", row.names = FALSE)
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
# create hourly ts
K_hourly20 <- hourly(Kennicott_data20)

# normalize data, convert to xts
K_norm20 <- hourly_norm(K_hourly20)

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
K_memb_break20 <- membership_breaks(memb, K_hourly20)

## visualize results
# c-q plot
ggplot(K_hourly20, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(K_hourly20, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

if (writecsv == 1){
  write.csv(K_hourly20,"Kennicott_hr20.csv", row.names = FALSE)
}

############### Hobo Data #########################

# load EC data from Hobo Loggers
Glacier21 <- read.csv('data/Glacier_C.csv')  
Glacier21$datetime <- as.POSIXct(Glacier21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Glacier21$datetime <- as.POSIXct(format(Glacier21$datetime), tz = "America/Anchorage")

site <- '15272502' #Glacier C. at Alyeska 
params <- '00060' # discharge 
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

ggplot(Glacier_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Glacier C. 2021')+
  theme_cust()

## prep data
# create hourly ts
Gl_hourly21 <- hourly(Glacier_data21) 

# normalize data, convert to xts
Gl_norm21 <- hourly_norm(Gl_hourly21) 

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
Gl_memb_break21 <- membership_breaks(memb, Gl_hourly21)


## visualize results
# c-q plot
ggplot(Gl_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Gl_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

####### Snow R. #########
Snow21 <- read.csv('data/Snow_R.csv')
Snow21$datetime <- as.POSIXct(Snow21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Snow21$datetime <- as.POSIXct(format(Snow21$datetime), tz = "America/Anchorage")

site <- '15243900' #Snow R. nr Seward AK 
params <- '00060' # discharge 
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

ggplot(Snow_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Snow R. 2021')+
  theme_cust()

## prep data
# create hourly ts
Sn_hourly21 <- hourly(Snow_data21)

# normalize data, convert to xts
Sn_norm21 <- hourly_norm(Sn_hourly21) 

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
Sn_memb_break21 <- membership_breaks(memb, Sn_hourly21)

## visualize results
# c-q plot
ggplot(Sn_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(Sn_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
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
# create hourly ts
kn_hourly21 <- hourly(Knik_data21)

# normalize data, convert to xts
Kn_norm21 <- hourly_norm(kn_hourly21)

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
Kn_memb_break21 <- membership_breaks(memb, kn_hourly21)

## visualize results
# c-q plot
ggplot(kn_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(kn_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
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
# create hourly ts
ken_hourly21 <- hourly(Kenai_data21) 

# normalize data, convert to xts
ken_norm21 <- hourly_norm(ken_hourly21)

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
Ken_memb_break21 <- membership_breaks(memb, ken_hourly21)

## visualize results
# c-q plot
ggplot(ken_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(ken_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


####### Sixmile C. #########
Sixmi21 <- read.csv('data/Sixmile_C.csv')
Sixmi21$datetime <- as.POSIXct(Sixmi21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
Sixmi21$datetime <- as.POSIXct(format(Sixmi21$datetime), tz = "America/Anchorage")

site <- '15271000' #Sixmile C. nr Hope AK 
params <- '00060' # discharge 
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

ggplot(Sixmi_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = Sp_cond), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Sixmi C. 2021')+
  theme_cust()

## prep data
# create hourly ts
sx_hourly21 <- hourly(Sixmi_data21)

# normalize data, convert to xts
sx_norm21 <- hourly_norm(sx_hourly21) 

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
Sxm_memb_break21 <- membership_breaks(memb, sx_hourly21)

## visualize results
# c-q plot
ggplot(sx_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(sx_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()


####### Little Susitna C. #########
LSus21 <- read.csv('data/L_Sus.csv')
LSus21$datetime <- as.POSIXct(LSus21$datetime, format="%m/%d/%y %H:%M", tz = "America/Anchorage")
LSus21$datetime <- as.POSIXct(format(LSus21$datetime), tz = "America/Anchorage")

site <- '15271000' #Sixmile C. nr Hope AK 
params <- '00060' # discharge 
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

ggplot(LSus_data21, aes(x = datetime, y = Q_m3s)) + 
  geom_line(color = "blue2") +
  geom_line(aes(y = SC_smoothed), color = "darkorchid2") +
  xlab(NULL) + 
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  labs(title = 'Sixmi C. 2021')+
  theme_cust()

## prep data
# create hourly ts
ls_hourly21 <- hourly(LSus_data21)

# normalize data, convert to xts
ls_norm21 <- hourly_norm(ls_hourly21)

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
memb_break20 <- membership_breaks(memb, ls_hourly21)
memb_break20

## visualize results
# c-q plot
ggplot(ls_hourly21, aes(x=log(Q_mean) , y= log(SC_mean))) +
  geom_point(aes(color= factor(memb)))+
  scale_colour_brewer(palette = "Paired")+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
  theme_cust()

# q ts
ggplot(ls_hourly21, aes(x= ymd(datetime_hourly), y= Q_mean)) +
  geom_point(aes(color= factor(memb)))+
  #scale_color_discrete(drop=FALSE)+
  scale_colour_brewer(palette = "Paired")+
  ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  xlab(NULL) + 
  theme_cust()

write.csv(ls_hourly21,"l_sus_hr21.csv", row.names = FALSE)

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
  write.csv(Gl_hourly21,"Glacier_hr21.csv", row.names = FALSE)
  write.csv(Sn_hourly21,"Snow_hr21.csv", row.names = FALSE)
  write.csv(kn_hourly21,"Knik_hr21.csv", row.names = FALSE)
  write.csv(ken_hourly21,"Kenai_hr21.csv", row.names = FALSE)
  write.csv(sx_hourly21,"Sixmi_hr21.csv", row.names = FALSE)
  write.csv(ls_hourly21,"l_sus_hr21.csv", row.names = FALSE)
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
 # hourly_data17$Ca_ts <- (ca_slope* hourly_data17$SC_filled + ca_int) *1000* 86400* hourly_data17$Q_m3shourly /1e6
 # 
 # ggplot(new_df17, aes(x= date , y= hourly_data17$Ca_ts)) +geom_point(aes(color= factor(memb)))+
 #   #scale_color_discrete(drop=FALSE)+
 #   scale_colour_brewer(palette = "Dark2")+
 #   ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
 #   xlab(NULL) + 
 #   theme_cust()