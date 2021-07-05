#wd
#setwd("C:/Users/nickg/OneDrive/Desktop/P")
install.packages("readr")
require(readr)
data<-read.csv("data/nest_defense_test_final.csv")

View(data)


library(lme4)
library(performance)
library(MASS)
library(arm)
library(broom)
library(tidyverse)
require(moments)
require(e1071)

#summary stats

var(data$MinDisRaw)
sd(data$MinDisRaw)
mean(data$MinDisRaw)
median(data$MinDisRaw)

summary(data$MinDisRaw) #600m was from a female (during laying) that was perched on the ground and calling
summary(data$Dives)
summary(data$Score)

#data exploration----

# Outliers Y
boxplot(data$MinDisRaw)
boxplot(data$Dives)
boxplot(data$Score)

require(lattice)
Z <- cbind(data$Dives, data$Min, data$Score, data$MinDisRaw*-1)

colnames(Z) <- c("Stoops", "Min", "Score", "MinDisRaw")

dotplot(as.matrix(Z), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        jitter.data=TRUE,
        xlab = "Value of the variable",
        ylab = "Order of the data from text file")

#homogentiy Y
  
  #MinDis
boxplot(MinDisRaw~NestStage, data=data)
boxplot(MinDisRaw~Year, data=data)
boxplot(MinDisRaw~Sex, data= data)

require(ggplot2)
Fig3=ggplot(data,aes(x=Sex, y=MinDisRaw, fill = factor(NestStage)))+
  geom_boxplot()+ 
  labs( x=c("", "", ""),y= "Minimum distance to observer (m)")+
  theme_classic()+
  
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())+ 
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig3+ facet_grid(~Year) + scale_y_continuous(limits = c(0,100))

  #Stoops
boxplot(Dives~NestStage, data=data)
boxplot(Dives~Year, data=data)
boxplot(Dives~Sex, data=data)



Fig4=ggplot(data,aes(x=Sex, y=Dives, fill = factor(NestStage)))+
  geom_boxplot()+ 
  labs( x=c("", "", ""),y= "Number of Stoops")+
  theme_classic()+
  
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())+ 
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig4+ facet_grid(~Year)


  #Score
boxplot(Score~NestStage, data=data) # the scoring limits the amount of variation possible b/n observations 
boxplot(Score~Year, data=data)
boxplot(Score~Sex, data=data)



Fig4=ggplot(data,aes(x=Sex, y=Score, fill = factor(NestStage)))+
  geom_boxplot()+ 
  labs( x=c("", "", ""),y= "Nest Defense Score")+
  theme_classic()+
  
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())+ 
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig4+ facet_grid(~Year)

#Normality 

  #Min Dis
hist(data$MinDisRaw)

library(lattice)
histogram( ~ MinDisRaw | NestStage , type = "count",
           xlab = "Min Distance (m)",
           ylab = "Frequency",
           nint=30,layout=c(1,3),
           strip.left = strip.custom(bg = 'white'),
           strip = F,
           col.line = "black", col = "white",
           scales = list(x = list(relation = "same"),
                         y = list(relation = "same"),
                         draw = TRUE), data=data)


  #Stoops

histogram( ~ Dives | NestStage , type = "count",
           xlab = "Number of Stoops",
           ylab = "Frequency",
           nint=30,layout=c(1,3),
           strip.left = strip.custom(bg = 'white'),
           strip = F,
           col.line = "black", col = "white",
           scales = list(x = list(relation = "same"),
                         y = list(relation = "same"),
                         draw = TRUE), data=data)


  #Score
hist(data$Score)

histogram( ~ Score | NestStage , type = "count",
           xlab = "Nest Defense Score",
           ylab = "Frequency",
           nint=30,layout=c(1,3),
           strip.left = strip.custom(bg = 'white'),
           strip = F,
           col.line = "black", col = "white",
           scales = list(x = list(relation = "same"),
                         y = list(relation = "same"),
                         draw = TRUE), data=data)

#Transform Min Dis----


#Transform
data$Min<- -1*(log10(max(data$MinDis + 1)-data$MinDis))
hist(data$Min)
# Check histograms for normality after transformation
histogram( ~ Min| NestStage , type = "count",
           xlab = "Nest Defense Score",
           ylab = "Frequency",
           nint=30,layout=c(1,3),
           strip.left = strip.custom(bg = 'white'),
           strip = F,
           col.line = "black", col = "white",
           scales = list(x = list(relation = "same"),
                         y = list(relation = "same"),
                         draw = TRUE), data=data)


#backtransform to original Min distance value
data$Min1<- -1*(10^(max(data$Min) -data$Min)-0.9)
hist(data$Min1)#backtransformed histogram
hist(data$MinDis)#original data no transformation

# Year as a factor
data100$Year2<-as.factor(data100$Year2)


#Full model----

data100<-subset(data, data$MinDisRaw <= 200)

hist(data100$MinDisRaw)
data100$MinSq<-sqrt(data100$MinR+0.01)
hist(data100$MinSq)
require(glmmTMB)

hist(data100$MinSq)

boxplot(data100$MinR~ data100$ID)

# require(highstat)
MyVar <- c("Mass.std", "Ivore.std", "Song.std", "DoYr", "Day")
FD3$E3 <- E3
MyMultipanel.ggp2(Z = FD3, 
                  varx = MyVar, 
                  vary = "E3", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)



require(lme4)
library(glmmTMB)
install.packages("glmmTMB")
install.packages('TMB', type = 'source')
m4a<-glmmTMB(MinR~ Sex + NestStage + Year2 +(1|ID) ,
         data=data100, family=Gamma("log"))

m4<-lmer(MinR~ -1+Sex + NestStage + Year2 +(1|ID) + (1|ID_Series),
            data=data100)
m4<-lmer(MinSq~ -1+Sex + NestStage + Year2 +(1|ID) + (1|ID_Series) + (1|SiteID),
         data=data100)

#summary
summary(m4)
summary(m4b)


table(data$MinDisRaw)
table(data$Score)
table(data$Score, data$MinDisRaw)
table(data$Score, data$MinDisRaw, data$NestStage)
      

#diagnostic plots----
plot(m4)


# Assumption 1-linearity-------

#Min distance values vs residuals
plot(m4b, MinSq ~ resid(.)) #shows linear relationship

#Min distance values vs fitted
plot(m4b, MinSq ~ fitted(.), abline = c(0,1)) #shows linear relationship


#Assumption 2-homogeneity of variance---------------

#residual vs fitted
plot(m4, type = c("p", "smooth"), col.line = 1) #homogeneity, no patterns&randomly spaced (line almost horizontal)






# Get residuals and fitted values of the lmm
E1 <- resid(m4b)
F1 <- fitted(m4b)


# Plot residuals versus fitted values.
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = 1)
abline(v=0)

#Assumption 3- Residuals normally distributed----

#hist of residuals
hist(resid(m4b)) #residuals normally distributed

#random effects histograms
require(MCMCglmm)
require(arm)
hist(ranef(m4)$ID)#normally distributed
hist(ranef(m4)$ID_Series)#normally distributed
hist(ranef(m4)$SiteID) #normally distributed
#qq plot 
lattice::qqmath(m4) #normality of errors

#simulated parameters + repeatability----
require(dplyr)
library(arm)
library(MCMCglmm)
smod<-sim(m4,1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))
##Between individual variance
bID<-smod@ranef$ID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##ID variance posterior distribution
require(MCMCglmm)
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)
##Between individual variance, ID_Series
bID2<-smod@ranef$ID_Series[,,1]
bvar2<-as.vector(apply(bID2, 1, var)) ##ID_Series variance posterior distribution
require(MCMCglmm)
bvar2<-as.mcmc(bvar2)
posterior.mode(bvar2 )## mode of the distribution
HPDinterval(bvar2)
##Between individual variance, NestID
bID3<-smod@ranef$SiteID[,,1]
bvar3<-as.vector(apply(bID3, 1, var)) ##ID_Series variance posterior distribution
require(MCMCglmm)
bvar3<-as.mcmc(bvar3)
posterior.mode(bvar3 )## mode of the distribution
HPDinterval(bvar3)
###residual variance
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)

### Long-term repeatability
r<-bvar/(rvar+bvar +bvar2 +bvar3 )
posterior.mode(r)
HPDinterval(r)  ##repeatability 

###Short-term
r1<-(bvar +bvar2 + bvar3 )/(rvar+bvar +bvar2 +bvar3)
posterior.mode(r1)
HPDinterval(r1)  ##repeatability

#backtransform----
(10^(max(data$Min)  +1.30395)-0.9) # FEMALE MIN DIS IN METERS

(10^(max(data$Min)  +1.23289)-0.9) # MALE MIN DIS IN METERS

(4.93^2)


#plot----
require(ggplot2)

Fig2=ggplot(data,aes(x=factor(Sex), y=MinDisRaw))+
  geom_boxplot()+
  labs( x="",y= "Min Distance (m)")+
   theme_classic()+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig2 + theme_classic()


Fig22=ggplot(data,aes(x=factor(Sex), y=Dives))+
  geom_boxplot()+
  labs( x="",y= "Number of Stoops")+
  theme_classic()+
  
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())+ 
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig22 + theme_classic()

Fig3=ggplot(data,aes(x=Sex, y=MinDisRaw, fill = factor(NestStage)))+
   geom_boxplot()+ 
  labs( x=c("", "", ""),y= "Minimum distance to observer (m)")+
  theme_classic()+
  
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())+ 
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig3+ facet_grid(~Year)



ggplot(aes(x, y), data = d) +
  geom_point()


Fig33=ggplot(data,aes(x=Sex, y=Dives, fill = factor(NestStage)))+
  geom_boxplot()+ 
  labs( x=c("", "", ""),y= "Number of Stoops")+
  theme_classic()+
  
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())+ 
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 14))

Fig33+ facet_grid(~Year)


