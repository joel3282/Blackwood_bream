# Williams et al 2019 
# Update Feb 2019

rm(list=ls())


##### Data preperation and exploration

library(tidyverse)



# read in data
dat <- read_csv("Data/Blackwood-data.csv")

# dat1 = only samples where data for prey IS available

# dat2 = 'plugging' NAs in predictor variable purely for data exploration
# dat3 = prey NOT included as a factor


# dat1 exploration
dat1 <- as.data.frame(dat)
glimpse(dat1)

plot(dat1$Temp, dat1$Eggs)
plot(dat1$Sal, dat1$Eggs)
plot(dat1$DO, dat1$Eggs)

plot(dat1$Temp, dat1$Feeding)
plot(dat1$Sal, dat1$Feeding)
plot(dat1$DO, dat1$Feeding)

plot(dat1$Plank_vol, dat1$Feeding)
plot(dat1$Naupli, dat1$Feeding)

hist(dat$Eggs)
hist(dat1$Yolk_sac)
hist(dat1$Feeding)

# dat2 exploration
dat2 <- dat1
str(dat2)
colSums(is.na(dat2))  

set.seed(12345678)

Extract <- dat2 %>% sele %>% 
  rfImputeUnsupervised() %>% 
  within({
    Turb <- round(Turb)
    Naupli <- round(Naupli)
  })
dat2[.,names(Extract)] <- Extract
