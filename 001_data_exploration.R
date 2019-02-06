# Williams et al 2019 
# Update Feb 2019

##### Data preperation and exploration

library(tidyverse)

# read in data
dat <- read_csv("Data/Blackwood-data.csv")

# dat1 = only samples where data for prey IS available
# dat2 = prey NOT included as a factor
# dat3 = 'plugging' NAs in predictor variable purely for data exploration




dat <- as.data.frame(dat)
class(dat)
glimpse(dat)

as.numeric(dat$Eggs)
hist(dat$Eggs)







# model concept
y1:3 ~ . + (1|Depth/Site) + (1|Month)