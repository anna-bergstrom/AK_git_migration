library(tidyverse)
library(lubridate)
library(zoo)
library(imputeTS)
library(ggpubr)
library(signal)
library(strucchange)


rm(list= ls())

# load and configure all datasets 
full_data16 = read_csv('data/WG_chemsitry_2016.csv')
full_data17 = read_csv('data/WG_chem_2017.csv')

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


# Using Signal Package to smooth data with "butter" filter


gauge_data16$SC_filter<- signal::filter(signal::butter(2,1/5,type='low'),gauge_data16$SC_filled) 
