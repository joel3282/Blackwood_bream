#####
# Williams et al 2019

# Full subsets analysis as per Fischer et al 2016

# Script information----
# This script is designed to work with long format data - where response variables are stacked one upon each other (see http://tidyr.tidyverse.org/)
# There are random factors
# We have used a Tweedie error distribution to account for the high occurence of zero values in the dataset.
# We have implemented the ramdom effects and Tweedie error distribution using the mgcv() package

# librarys----

library(tidyr)
library(dplyr)
library(forcats)
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
library(doParallel)
library(gamm4)
library(RCurl)#needed to download data from GitHub
library(doSNOW)
library(ggplot2)
library(devtools)
# install package----
# devtools::install_github("beckyfisher/FSSgam_package") #run once
library(FSSgam)
?full.subsets.gam

rm(list=ls())

# GAMM models----

#name<-"snapper"


# Read in data----

dat<-read.csv("Data/Blackwood-data.csv")

# Need for two datasets as there are NAs for prey data due to sampling not being collected.

# dat1 = only samples where data for prey IS available
# dat2 = prey NOT included as a factor

%>%
  select(-no3_m,-po4_m, -sst1961_20)%>% # drop these as they have NA's
  filter(dst2ramp<150000)%>%
  mutate(ausbath=-ausbath)%>%
  filter(ausbath<(50))%>% #tim added this - might need to re-run
  filter(depth<50)%>%
  filter(dst2water<80000)%>%
  filter(relief<40)%>%
  filter(dst2mland<100000)%>%
  filter(no3_sd<2.5)%>%
  filter(t_sd<2.5)%>%
  filter(t_m<26)%>%
  filter(win_sst_m>12)%>% #Tim removed outliers
  # filter(win_sst_m<23.5)%>% #Tim removed outliers
  
  # dplyr::select(Taxa,response,latitude,longitude,year,time,location,status,ausbath,surface,slope,
  #        relief,t_m,t_sd,no3_m,no3_sd,po4_m,po4_sd,dst2water,dst2shelf,dst2coast,dst2ramp,dst2road,
  #        dist2town,dst2townc,day,dst2reef,depth,aspect)%>%
  na.omit()%>%
  mutate(year=as.factor(year))%>% ## Code year as a factor - to include as a random effect
  glimpse()


glimpse(dat)
summary(dat$win_sst_m)
summary(dat$response)



# Check out the depth distribution----
ggplot(data=dat%>%filter(response>0),aes(x=depth, colour=status))+
  geom_histogram()+
  facet_grid(state~.,scales="free")

# Check out the wsst distribution----
ggplot(data=dat,aes(x=win_sst_m,y=response, colour=status))+
  geom_point()+
  geom_smooth(method="gam", formula= y ~ s(x, bs = "cs",k=3))+
  facet_grid(taxa~.,scales="free")


# Check for NA's----
extra_NA<-dat %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise_all(funs(sum(is.na(.)))) 
extra_NA


#### Check distribution of the response
unique(dat$taxa)

## legal

par(mfrow=c(1,1))

legal<-dat%>%
  filter(taxa=="legal")%>%
  glimpse()
hist(legal$response)

### sublegal

sub.legal<-dat%>%
  filter(taxa=="sub.legal")%>%
  glimpse()
hist(sub.legal$response)

### largest

largest<-dat%>%
  filter(taxa=="largest")%>%
  glimpse()
hist(largest$response)


# Remeber - to include the ratio----



# Set predictor (i.e. continous) variables----
glimpse(dat)
pred.vars=c("depth","relief","t_sd","no3_sd","dst2shelf","dst2water","dst2ramp","dst2reef","dst2coast","win_sst_m")
# We have excluded those that were strongly correlated - could try to bring some back in?
# "t_m",


test.dat<-dat

# Check for correlation of predictor variables- remove anything highly correlated (>0.90)---
round(cor(test.dat[,pred.vars]),2)  ## There is some error here when applied to Bio-oracle covariates

# nothing is highly correlated


# Plot of likely transformations
plot.new()

par(mfrow=c(3,2))
for (i in pred.vars) {
  x<-test.dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

# Review of individual predictors - we have to make sure they have an even distribution---
#If the data are squewed to low numbers try sqrt>log or if squewed to high numbers try ^2 of ^3
# Is data is clumped into distances - of a feature of interest - try as a linear predictor

## Remove potential outlier
## This is done when I bring in the data

## Transform predictors if necessary - or put them as linear?
# dst2ramp - could be better as a tranform?


## However - I don't like to transfor predictor as the transformation is tricky to interpret

# Run the full subset model selection----
# Set directory for the model outputs----
setwd(model.out)
dir()


# Check to make sure Response vector has not more than 85% zeros---

unique.vars=unique(as.character(dat$taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<0.9){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
#there are not many legal????
unique.vars.use

# unique.vars.use<-c("legal","sub.legal") #have removed largest - to make this quicker
write.csv(unique.vars.use,file=paste(name,"unique.vars.use.csv",sep = "_"))


# Set variables needed for the FSS function-
glimpse(dat)
resp.vars=unique.vars.use
use.dat=dat
factor.vars=c("status") # Status as a Factor with two levels
# linear.vars=c("dst2coast","dst2reef") #Need to put back in
# cyclic.vars=c("aspect") #would be good to bring in ,"time"_ need to put aspect back in
null.vars=c("year","state") #,"day"
pred.vars=c("depth","relief","t_sd","no3_sd","dst2shelf","dst2water","dst2ramp","dst2reef","dst2coast","win_sst_m") #Need to many put back in
#pred.vars=c("depth", "relief","win_sst_m")


?full.subsets.gam


# Loop through the FSS function for each Taxa----
out.all=list()
var.imp=list()
for(i in 1:length(resp.vars)){use.dat=dat[which(dat$taxa==resp.vars[i]),]

Model1=gam(response~s(relief,k=3,bs="cr")+s(year,state,bs="re"),
           family=tw(), data=use.dat)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars,
                             pred.vars.fact=factor.vars,
                             factor.smooth.interactions = list(
                               fact.vars=c("status"),
                               cont.vars=c("depth","win_sst_m")), #"dst2ramp",
                             # linear.vars=c("dst2reef","dst2coast")),
                             # linear.vars = linear.vars,
                             # cyclic.vars=cyclic.vars,
                             # factor.factor.interactions = T,
                             k=3,
                             max.predictors=4,
                             
                             null.terms="s(year,state,bs='re')")

out.list=fit.model.set(model.set,
                       max.models=10000,
                       save.model.fits=T, #to save on memory
                       parallel=T)


names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
out.i=mod.table[which(mod.table$delta.AICc<=4),]
out.all=c(out.all,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) 

# plot the best models
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text=resp.vars[i],outer=F)}
  dev.off()
}
}

# Model fits and importance scores---
names(out.all)=resp.vars
names(var.imp)=resp.vars
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
# write.csv(all.mod.fits[,-2],file=paste(name,"all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"all.var.imp.csv",sep="_"))
write.csv(all.mod.fits,file=paste(name,"all.mod.fits.csv",sep="_"))

# write.csv(mod.table,"mod.table.csv")



