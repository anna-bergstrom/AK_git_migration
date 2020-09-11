library(tidyverse)
library(lubridate)
library(zoo)


rm(list= ls())

# load and configure all datasets 
full_data16 = read_csv('data/WG_chemsitry_2016.csv')
full_data17 = read_csv('data/WG_chem_2017.csv')

gauge_data16 = read_delim('data/wolvQ_16season_m.txt',"\t")
gauge_data17 = read_delim('data/wolvQ_17season_m.txt',"\t")
gauge_data16$datetime<- as.POSIXct(gauge_data16$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")
gauge_data17$datetime<- as.POSIXct(gauge_data17$datetime,TZ = "America/Anchorage", format="%m/%d/%Y %H:%M:%S")

###### Working with the timeseries data ###############
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
ggplot(gauge_data16, aes(x= Q_m3s, y= SC)) +geom_point(aes(col=period), size=2 )+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
 scale_colour_brewer(palette = "Dark2")+
 theme_cust()+
 theme(legend.position="None")

# 2017 EC-Q plot
ggplot(gauge_data17, aes(x= Q_m3s, y= SC)) +geom_point(aes(col=period), size=2 )+
  xlab(expression(paste("Q (m"^"3","s"^"-1", ")"))) + ylab(expression(paste("EC (" ,  mu,  "S cm"^"-1", ")")))+
  scale_colour_brewer(palette = "Dark2")+
  theme_cust()+
  theme(legend.position="None") 

# Creating a moving average for timeseries plots and analysis
swindow <- 192
smoothed_16 <- data.frame("Q_m3s" = rollapply(gauge_data16$Q_m3s,swindow,mean, na.rm = TRUE, fill = NA))
smoothed_16$SC <- rollapply(gauge_data16$SC,swindow,mean, na.rm = TRUE, fill = NA)

smoothed_17 <- data.frame("Q_m3s" = rollapply(gauge_data17$Q_m3s,swindow,mean, na.rm = TRUE, fill = NA))
smoothed_17$SC <- rollapply(gauge_data17$SC,swindow,mean, na.rm = TRUE, fill = NA)

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
