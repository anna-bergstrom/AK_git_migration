library(tidyverse)
library(lubridate)


# load all datasets
full_data16 = read_csv('data/WG_chemsitry_2016.csv')
full_data17 = read_csv('data/WG_chem_2017.csv')

gauge_data16 = read_delim('data/wolvQ_16season_m.txt',"\t")
gauge_data17 = read_delim('data/wolvQ_17season_m.txt',"\t")

# Creating a variable for the manually identified breaks in the 16 and 17 datasets
breaks16 <- mdy_hms(c('05/19/2016 23:45:00','06/13/2016 23:45:00','07/16/2016 23:45:00'))
breaks17 <- mdy_hms(c('05/12/2017 23:45:00','06/22/2017 23:45:00','08/19/2017 23:45:00'))

