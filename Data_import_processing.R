library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)

rm(list= ls())

# load all datasets
full_data16 = read_csv('data/WG_chemsitry_2016.csv')
full_data17 = read_csv('data/WG_chem_2017.csv')

gauge_data16 = read_delim('data/wolvQ_16season_m.txt',"\t")
gauge_data17 = read_delim('data/wolvQ_17season_m.txt',"\t")
gauge_data16$datetime<- as.POSIXct(gauge_data16$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")
gauge_data17$datetime<- as.POSIXct(gauge_data17$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")

# Creating a variable for the manually identified breaks in the 16 and 17 datasets
dates16 <- c('05-19-2016 23:45:00','06-13-2016 23:45:00','07-16-2016 23:45:00')
breaks16 <- data.frame("datetime" = as.POSIXct(dates16,TZ = "America/Anchorage", format="%m-%d-%Y %H:%M:%S"))
breaks17 <- data.frame("datetime" =as.POSIXct(c('05/12/2017 23:45:00','06/22/2017 23:45:00','08/19/2017 23:45:00'), format="%m/%d/%Y %H:%M:%S", TZ = "America/Anchorage"))


#Calculating q in m3/s
gauge_data16$Q_m3s <- gauge_data16$Q_cfs * 0.028316847 
gauge_data17$Q_m3s <- gauge_data17$Q_cfs * 0.028316847 

# Adding a column designatng the periods in 16 and 17
gauge_data16$period <-NA
gauge_data16$period[gauge_data16$datetime < breaks16[1,]] <- "p1"
gauge_data16$period[gauge_data16$datetime >= breaks16[1,] & gauge_data16$datetime < breaks16[2,] ] <- "p2"
gauge_data16$period[gauge_data16$datetime >= breaks16[2,] & gauge_data16$datetime < breaks16[3,] ] <- "p3"
gauge_data16$period[gauge_data16$datetime >= breaks16[3,]] <- "p4"

gauge_data17$period <-NA
gauge_data17$period[gauge_data17$datetime < breaks17[1,]] <- 1
gauge_data17$period[gauge_data17$datetime >= breaks17[1,] & gauge_data17$datetime < breaks16[2,] ] <- 2
gauge_data17$period[gauge_data17$datetime >= breaks17[2,] & gauge_data17$datetime < breaks16[3,] ] <- 3
gauge_data17$period[gauge_data17$datetime >= breaks17[3,]] <- 4

# Recreating Figure 2 plots
CQ16 <- ggplot(gauge_data16, aes(x= Q_m3s, y= SC)) +geom_point(aes(col=period, order= -order), size=3 )
CQ16 + xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))
CQ16 + theme(legend.position="None")
CQ16 + scale_colour_brewer(palette = "Dark2")
