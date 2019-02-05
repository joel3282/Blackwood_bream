# Williams et al 2019 
# Update Feb 2019

##### Data preperation and exploration

library(tidyverse)

dat <- read_csv("Data/Blackwood-data.csv")
dat <- as.data.frame(dat)
class(dat)
glimpse(dat)

dat %>%
  as.numeric(Eggs) %>% 
  dat %>% 
  hist(Eggs)
