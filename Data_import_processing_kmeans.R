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


rm(list= ls())

# load and configure all datasets 
full_data16 = read_csv('data/WG_chemsitry_2016.csv')
full_data17 = read_csv('data/WG_chem_2017.csv')

full_wx = read_csv('data/wolverine990_15min_LVL2_1619.csv')
full_wx$local_time<- as.POSIXct(full_wx$local_time,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M")

gauge_data16 = read_delim('data/wolvQ_16season_m.txt',"\t")
gauge_data17 = read_delim('data/wolvQ_17season_m.txt',"\t")
gauge_data16$datetime<- as.POSIXct(gauge_data16$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")
gauge_data17$datetime<- as.POSIXct(gauge_data17$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")

# loading 2018 and 2019 gauge data
gauge_data18 = read_csv('data/wolvQ_18season_m.csv')
gauge_data19 = read_csv('data/wolvQ_19season_m.csv')
gauge_data18$datetime<- as.POSIXct(gauge_data18$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")
gauge_data19$datetime<- as.POSIXct(gauge_data19$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")

###### Working with the timeseries data ###############
# Creating a variable for the manually identified breaks in the 16 and 17 datasets
dates16 <- c('05-19-2016 23:45:00','06-13-2016 23:45:00','07-16-2016 23:45:00')
breaks16 <- data.frame("datetime" = as.POSIXct(dates16,TZ = "America/Anchorage", format="%m-%d-%Y %H:%M:%S"))
breaks17 <- data.frame("datetime" =as.POSIXct(c('05/12/2017 23:45:00','06/22/2017 23:45:00','08/19/2017 23:45:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage"))



#Calculating q in m3/s
gauge_data16$Q_m3s <- gauge_data16$Q_cfs * 0.028316847 
gauge_data17$Q_m3s <- gauge_data17$Q_cfs * 0.028316847 
gauge_data18$Q_m3s <- gauge_data18$Q * 0.028316847 
gauge_data19$Q_m3s <- gauge_data19$Q * 0.028316847 

# Adding a column designatng the periods in 16 and 17
gauge_data16$period <-NA
gauge_data16$period[gauge_data16$datetime < breaks16[1,]] <- "p1"
gauge_data16$period[gauge_data16$datetime >= breaks16[1,] & gauge_data16$datetime < breaks16[2,] ] <- "p2"
gauge_data16$period[gauge_data16$datetime >= breaks16[2,] & gauge_data16$datetime < breaks16[3,] ] <- "p3"
gauge_data16$period[gauge_data16$datetime >= breaks16[3,]] <- "p4"

gauge_data17$period <-NA
gauge_data17$period[gauge_data17$datetime < breaks17[1,]] <- "p1"
gauge_data17$period[gauge_data17$datetime >= breaks17[1,] & gauge_data17$datetime < breaks17[2,] ] <- "p2"
gauge_data17$period[gauge_data17$datetime >= breaks17[2,] & gauge_data17$datetime < breaks17[3,] ] <- "p3"
gauge_data17$period[gauge_data17$datetime >= breaks17[3,]] <- "p4"

# Recreating Figure 2 plots
# generating a custom theme to get rid of the shitty ggplot defaults 
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

# 2017 EC-Q plot
p17<- ggplot(gauge_data17, aes(x= Q_m3s, y= SC)) +geom_point(aes(col=period), size=2 )+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_brewer(palette = "Dark2")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position= c(0.9,0.6))+
  theme(legend.title = element_blank())



gauge_data18$count <- seq(1,length(gauge_data18$SC))
gauge_data19$count <- seq(1,length(gauge_data19$SC))

# 2018 EC-Q plot
p18<- ggplot(gauge_data18, aes(x= Q_m3s, y= SC)) +geom_point(aes(color= count))+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_gradient(low= "green", high="blue")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position="None") 

# 2019 EC-Q plot
p19<- ggplot(gauge_data19, aes(x= Q_m3s, y= SC)) +geom_point(aes(color= count))+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_gradient( name = NULL, breaks= c(130, 11865),labels=c("May","Sep"),low= "green", high="blue" )+
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  theme_cust()+
  theme(legend.position=c(0.89,0.6))
 


#putting all plots in a multipanel figure
figure <- ggarrange(p16, p17, p18,p19,
                    labels = c("2016", "2017", "2018", "2019"),
                    ncol = 2, nrow = 2, hjust= -1.6, vjust = 1.8)
figure


# Trimming the 2016 and 2017 timeseries data to the window where we consistently have Q and SC. 
bounds16 <- data.frame("datetime" =as.POSIXct(c('04/28/2016 00:00:00','10/06/2016 23:45:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage"))
bounds17 <- data.frame("datetime" =as.POSIXct(c('04/20/2017 14:00:00','10/17/2017 11:30:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage"))
gauge_data16<-gauge_data16[ which(gauge_data16$datetime >= bounds16[1,] & gauge_data16$datetime <= bounds16[2,]), ]
gauge_data17<-gauge_data17[ which(gauge_data17$datetime >= bounds17[1,] & gauge_data17$datetime <= bounds17[2,]), ]




# Creating a moving average for timeseries plots and analysis
swindow <- 192
smoothed_16 <- data.frame("Q_m3s" = rollapply(gauge_data16$Q_m3s,swindow,mean, na.rm = TRUE, fill = NA))
smoothed_16$SC <- rollapply(gauge_data16$SC,swindow,mean, na.rm = TRUE, fill = NA)

smoothed_17 <- data.frame("Q_m3s" = rollapply(gauge_data17$Q_m3s,swindow,mean, na.rm = TRUE, fill = NA))
smoothed_17$SC <- rollapply(gauge_data17$SC,swindow,mean, na.rm = TRUE, fill = NA)

smoothed_18 <- data.frame("Q_m3s" = rollapply(gauge_data18$Q_m3s,swindow,mean, na.rm = TRUE, fill = NA))
smoothed_18$SC <- rollapply(gauge_data18$SC,swindow,mean, na.rm = TRUE, fill = NA)

smoothed_19 <- data.frame("Q_m3s" = rollapply(gauge_data19$Q_m3s,swindow,mean, na.rm = TRUE, fill = NA))
smoothed_19$SC <- rollapply(gauge_data19$SC,swindow,mean, na.rm = TRUE, fill = NA)

# Time series plots
#2016
ggplot(data = NULL, aes(x = gauge_data16$datetime, y=smoothed_16$Q_m3s)) +geom_line(color = "blue2") +
  geom_line(aes(x = gauge_data16$datetime, y=smoothed_16$SC/10), color = "darkorchid2") +
  xlab(NULL) + ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  theme_cust()

#2017
ggplot(data = NULL, aes(x = gauge_data17$datetime, y=smoothed_17$Q_m3s)) +geom_line(color = "blue2") +
  geom_line(aes(x = gauge_data17$datetime, y=smoothed_17$SC/10), color = "darkorchid2") +
  xlab(NULL) + ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  theme_cust()

#2018
ggplot(data = NULL, aes(x = gauge_data18$datetime, y=smoothed_18$Q_m3s)) +geom_line(color = "blue2") +
  geom_line(aes(x = gauge_data18$datetime, y=smoothed_18$SC/10), color = "darkorchid2") +
  xlab(NULL) + ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  theme_cust()

#2019
ggplot(data = NULL, aes(x = gauge_data19$datetime, y=smoothed_19$Q_m3s)) +geom_line(color = "blue2") +
  geom_line(aes(x = gauge_data19$datetime, y=smoothed_19$SC/10), color = "darkorchid2") +
  xlab(NULL) + ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
  theme_cust()


# Interpolating to fill NAs 
gauge_data16$SC_filled<- na_interpolation(gauge_data16$SC, option = "linear", maxgap = Inf)
gauge_data17$SC_filled<- na_interpolation(gauge_data17$SC, option = "linear", maxgap = Inf)
gauge_data17$Qm3s_filled<- na_interpolation(gauge_data17$Q_m3s, option = "linear", maxgap = Inf)

gauge_data18$SC_filled<- na_interpolation(gauge_data18$SC, option = "linear", maxgap = Inf)
gauge_data18$Qm3s_filled<- na_interpolation(gauge_data18$Q_m3s, option = "linear", maxgap = Inf)

gauge_data19$SC_filled<- na_interpolation(gauge_data19$SC, option = "linear", maxgap = Inf)

###### Sub-setting weather data ############

#2016
precip16 <- data.frame(full_wx$Precip_Weighing_Incremental[full_wx$local_time >= bounds16[1,] & full_wx$local_time <= bounds16[2,]])
precip16$datetime <- full_wx$local_time[full_wx$local_time >= bounds16[1,] & full_wx$local_time <= bounds16[2,]]
names(precip16)[1] <- "precip16int"

# 2017
precip17 <- data.frame(full_wx$Precip_Weighing_Incremental[full_wx$local_time >= bounds17[1,] & full_wx$local_time <= bounds17[2,]])
precip17$datetime <- full_wx$local_time[full_wx$local_time >= bounds17[1,] & full_wx$local_time <= bounds17[2,]]
names(precip17)[1] <- "precip17int"
merged17<- merge.xts(xts(gauge_data17$Q_m3s, as.POSIXct(gauge_data17$datetime)),xts(gauge_data17$SC_filled, as.POSIXct(gauge_data17$datetime)),xts(precip17$precip17int, as.POSIXct(precip17$datetime)),fill=NA)
names(merged17)<- c("Q_m3s","SC_filled", "precip17")  
merged17$SC_filled<- na_interpolation(merged17$SC_filled, option = "linear", maxgap = Inf)
merged17$Q_m3s<- na_interpolation(merged17$Q_m3s, option = "linear", maxgap = Inf)

# first, normalizing all the data
gauge_data16$SC_norm <- 1-(max(gauge_data16$SC_filled)-gauge_data16$SC_filled)/(max(gauge_data16$SC_filled)-min(gauge_data16$SC_filled))
gauge_data16$Q_norm <- 1-(max(gauge_data16$Q_cfs)-gauge_data16$Q_cfs)/(max(gauge_data16$Q_cfs)-min(gauge_data16$Q_cfs))
gauge_data16$integer <- seq(from = 1, to= length(gauge_data16$Q_cfs), by=1)
gauge_data16$int_as_prop <- 1-(max(gauge_data16$integer)-gauge_data16$integer)/(max(gauge_data16$integer)-min(gauge_data16$integer))

gauge_data17$SC_norm <- 1-(max(gauge_data17$SC_filled)-gauge_data17$SC_filled)/(max(gauge_data17$SC_filled)-min(gauge_data17$SC_filled))
gauge_data17$Q_norm <- 1-(max(gauge_data17$Q_cfs)-gauge_data17$Q_cfs)/(max(gauge_data17$Q_cfs)-min(gauge_data17$Q_cfs))
gauge_data17$integer <- seq(from = 1, to= length(gauge_data17$Q_cfs), by=1)
gauge_data17$int_as_prop <- 1-(max(gauge_data17$integer)-gauge_data17$integer)/(max(gauge_data17$integer)-min(gauge_data17$integer))

# Using smoothing filters 

# Loess
 
# pretty darned righteous. no unwanted bumps. no head/tail anomalies. hard to
# optimize (lots of parameters), but maybe search can be automated.

 
 loess_mod <- loess(SC_filled ~ integer, data = gauge_data16, span = 0.02)
 
 gauge_data16$SC_loess <- predict(loess_mod, newdata = gauge_data16)
 
 
 loess_modQ <- loess(Q_m3s ~ integer, data = gauge_data16, span = 0.02)
 
 gauge_data16$Q_loess <- predict(loess_modQ, newdata = gauge_data16)
 
 all_xts <- xts(gauge_data16 %>%
                  select(c('SC_loess', 'SC_filled','Q_m3s','Q_loess')),
                order.by=gauge_data16$datetime)
 dygraph(all_xts)
 
 loess_mod17 <- loess(SC_filled ~ integer, data = gauge_data17, span = 0.02)
 
 gauge_data17$SC_loess <- predict(loess_mod17, newdata = gauge_data17)
 
 
 loess_modQ17 <- loess(Q_m3s ~ integer, data = gauge_data17, span = 0.02)
 
 gauge_data17$Q_loess <- predict(loess_modQ17, newdata = gauge_data17)
 
 all_xts <- xts(gauge_data17 %>%
                  select(c('SC_loess', 'SC_filled','Q_m3s','Q_loess')),
                order.by=gauge_data17$datetime)
 dygraph(all_xts)
 
 #normalizing smoothed data
 gauge_data16$SC_loessnorm <- 1-(max(gauge_data16$SC_loess)-gauge_data16$SC_loess)/(max(gauge_data16$SC_loess)-min(gauge_data16$SC_loess))
 gauge_data16$Q_loessnorm <- 1-(max(gauge_data16$Q_loess)-gauge_data16$Q_loess)/(max(gauge_data16$Q_loess)-min(gauge_data16$Q_loess))

 
 gauge_data17$SC_loessnorm <- 1-(max(gauge_data17$SC_loess)-gauge_data17$SC_loess)/(max(gauge_data17$SC_loess)-min(gauge_data17$SC_loess))
 gauge_data17$Q_loessnorm <- 1-(max(gauge_data17$Q_loess)-gauge_data17$Q_loess)/(max(gauge_data17$Q_loess)-min(gauge_data17$Q_loess))
 
 #ggplot(data = NULL, aes(x = gauge_data16$datetime, y=gauge_data16$Q_m3s)) +geom_line(color = "blue2") +
 #  geom_line(aes(x = gauge_data16$datetime, y= filtered$SC_loess/10), color = "darkorchid2") +
 #  xlab(NULL) + ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
 #  theme_cust()
 
 #normalizing precip data and adding it to the gauge data frame
 is.nan.data.frame <- function(x)
   do.call(cbind, lapply(x, is.nan))
 precip16[is.nan(precip16)] <- 0
 
 gauge_data16$Precip_norm <- 1-(max(precip16$precip16int)-precip16$precip16int)/(max(precip16$precip16int)-min(precip16$precip16int))
 
 ###### testing k-means #######
 #setting up matrix for cluster analysis
 kmeans16 <- gauge_data16[ ,c("SC_loessnorm", "Q_loessnorm", "int_as_prop","Precip_norm")]
 kmeans17 <- gauge_data17[ ,c("SC_loessnorm", "Q_loessnorm", "int_as_prop")]
 
 #Determining optimum number of clusters, elbow method - total within sum of square
 fviz_nbclust(kmeans16, kmeans, method = "wss") +
   geom_vline(xintercept = 4, linetype = 2)
 
 fviz_nbclust(kmeans17, kmeans, method = "wss") +
   geom_vline(xintercept = 4, linetype = 2)
 
 # testing k-means on 2016 - just smoothed data with no
 
 kmeans16ans <- kmeans(kmeans16,4,nstart=25)
 
 ggplot(gauge_data16, aes(x= datetime, y= SC_filled)) +geom_point(aes(color= kmeans16ans$cluster))+
   ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
   xlab(NULL) + 
   theme_cust()


 ggplot(gauge_data16, aes(x= Q_m3s, y= SC_filled)) +geom_point(aes(color= kmeans16ans$cluster))+
   ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
   xlab(NULL) + 
   theme_cust()

 # testing k-means on 2017
 kmeans17 <- gauge_data17[ ,c("SC_loessnorm", "Q_loessnorm", "int_as_prop")]
 kmeans17ans <- kmeans(kmeans17,5,nstart=25)
 
 ggplot(gauge_data17, aes(x= datetime, y= Q_cfs)) +geom_point(aes(color= kmeans17ans$cluster))+
   ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
   xlab(NULL) + 
   theme_cust()
 
 
 ggplot(gauge_data17, aes(x= Q_loess, y= SC_loess)) +geom_point(aes(color= kmeans17ans$cluster))+
   ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
   xlab(NULL) + 
   theme_cust()
 
 ######### Messing with depth-constrained clustering (to constrain clusters in time)##############
 
#Tried running it with 15 min data and it was taking too long (may try overnight, but not super practical)
 #So moving to try it with 1-hour intervals - averaging the Q and EC, and summing precip
#2016 
SC.xts <- xts(gauge_data16$SC_filled, as.POSIXct(gauge_data16$datetime)) 
hourly_data <- data.frame(period.apply(SC.xts,endpoints(SC.xts,"hours"),mean))
names(hourly_data)[1] <- "SC_filledhourly"
Q.xts <- xts(gauge_data16$Q_m3s, as.POSIXct(gauge_data16$datetime)) 
hourly_data$Q_m3shourly <- period.apply(Q.xts,endpoints(Q.xts,"hours"),mean)
precip.xts <- xts(precip16, as.POSIXct(gauge_data16$datetime)) 
hourly_data$precip_hourlysum <- period.sum(precip.xts[,1],endpoints(precip.xts,"hours"))

hourly_norm<- data.frame(1-(max(hourly_data$precip_hourlysum)-hourly_data$precip_hourlysum)/(max(hourly_data$precip_hourlysum)-min(hourly_data$precip_hourlysum)))
names(hourly_norm)[1] <- "Precip_hr_norm"
hourly_norm$SC_hr_norm<- data.frame(1-(max(hourly_data$SC_filledhourly)-hourly_data$SC_filledhourly)/(max(hourly_data$SC_filledhourly)-min(hourly_data$SC_filledhourly)))
hourly_norm$Q_hr_norm<- data.frame(1-(max(hourly_data$Q_m3shourly)-hourly_data$Q_m3shourly)/(max(hourly_data$Q_m3shourly)-min(hourly_data$Q_m3shourly)))

 #dep_con16 <- gauge_data16[ ,c("SC_loessnorm", "Q_loessnorm","Precip_norm")]
 #dep_con16_dist <- dist(dep_con16, method = "euclidean")
 
dep_con16_dist <- dist(hourly_norm, method = "euclidean")
 dep16_clust <- chclust(dep_con16_dist, method = "coniss")
 bstick(dep16_clust, ng=20, plot=TRUE)
 memb<-cutree(dep16_clust,k=5)
 #cent<- NULL
 #for (k in 1:5){cent<- rbind(cent,colMeans(hourly_norm[memb==k, ,drop = FALSE]))}
 #dep16clust5 <- chclust(dist(cent, method = "euclidean"), method = "coniss",members = table(memb))
 
 
 
 ggplot(hourly_data, aes(x=Q_m3shourly , y= SC_filledhourly)) +geom_point(aes(color= factor(memb)))+
    scale_colour_brewer(palette = "Dark2")+
   xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
   ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
   theme_cust()
 
 new_df = data.frame(new_df = index(hourly_data$Q_m3shourly), coredata(hourly_data$Q_m3shourly))
 colnames(new_df) = c("date","Q_m3shourly")
 new_df$SC_filledhourly <- hourly_data$SC_filledhourly
 
 ggplot(new_df, aes(x= date , y= Q_m3shourly)) +geom_point(aes(color= factor(memb)))+
   #scale_color_discrete(drop=FALSE)+
   scale_colour_brewer(palette = "Dark2")+
   ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
   xlab(NULL) + 
   theme_cust()
 
 for (k  in 1:length(unique(memb))) {
   temp <- new_df$date[which(memb==k)]
   if (k== 1){
     p_time <- as.POSIXct(tail(temp,n=1)) 
   } else{
     p_time [k] <- tail(temp, n=1)}
 }
 print(p_time)
 

 ## Plotting actual tree (takes a while and not particularly useful)
 # plot(dep16_clust, labels = NULL, hang = 0.1, axes = TRUE, xvar=1:(length(dep16_clust$height)+1), xlim=NULL, ylim=NULL,  x.rev = FALSE, y.rev=FALSE, horiz=FALSE)

 ### 2017 ###
 hourly_data17 <- data.frame(period.apply(merged17$SC_filled,endpoints(merged17,"hours"),mean))
 hourly_data17$Q_m3shourly <- period.apply(merged17$Q_m3s,endpoints(merged17,"hours"),mean)
 hourly_data17$precip_hourlysum <- period.sum(merged17$precip17,endpoints(merged17,"hours"))
 
 hourly_norm<- data.frame(1-(max(hourly_data17$precip_hourlysum)-hourly_data17$precip_hourlysum)/(max(hourly_data17$precip_hourlysum)-min(hourly_data17$precip_hourlysum)))
 names(hourly_norm)[1] <- "Precip_hr_norm"
 hourly_norm$SC_hr_norm<- data.frame(1-(max(hourly_data17$SC_filled)-hourly_data17$SC_filled)/(max(hourly_data17$SC_filled)-min(hourly_data17$SC_filled)))
 hourly_norm$Q_hr_norm<- data.frame(1-(max(hourly_data17$Q_m3shourly)-hourly_data17$Q_m3shourly)/(max(hourly_data17$Q_m3shourly)-min(hourly_data17$Q_m3shourly)))
 
 dep_con17_dist <- dist(hourly_norm, method = "euclidean")
 dep17_clust <- chclust(dep_con17_dist, method = "coniss")
 bstick(dep17_clust, ng=20, plot=TRUE)
 memb<-cutree(dep17_clust,k=7)

 
 ggplot(hourly_data17, aes(x=Q_m3shourly , y= SC_filled)) +geom_point(aes(color= factor(memb)))+
   scale_colour_brewer(palette = "Dark2")+
   xlab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
   ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")"))) + 
   theme_cust()
 
 new_df17 = data.frame(new_df = index(hourly_data17$Q_m3shourly), coredata(hourly_data17$Q_m3shourly))
 colnames(new_df17) = c("date","Q_m3shourly")
 new_df$SC_filledhourly <- hourly_data17$SC_filledhourly
 
 ggplot(new_df17, aes(x= date , y= Q_m3shourly)) +geom_point(aes(color= factor(memb)))+
   #scale_color_discrete(drop=FALSE)+
   scale_colour_brewer(palette = "Dark2")+
   ylab(expression(paste("Q (m"^"3","s"^"-1", ")")))+
   xlab(NULL) + 
   theme_cust()

 for (k  in 1:length(unique(memb))) {
   temp <- new_df17$date[which(memb==k)]
   if (k== 1){
  p_time <- as.POSIXct(tail(temp,n=1)) 
  } else{
   p_time [k] <- tail(temp, n=1)}
 }
 print(p_time)
 
 
 