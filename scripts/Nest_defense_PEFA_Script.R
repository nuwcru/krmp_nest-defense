
# add some code here


#Aggression Scoring PCA----
setwd("C:/Users/Nick/Desktop/PEFA/P")
data<-read.csv("nest_defense_pca.csv")


#packages
install.packages("missMDA")
install.packages("ggplot2")
install.packages("ggfortify")

library(ggplot2)
library(ggfortify)
library(missMDA)
library(lme4)



#PCA with three variables FID, Min Dis, Dives
mypr<-prcomp(data[5:7], center= TRUE, scale = TRUE)
summary(mypr)
print(mypr)
plot(mypr, type= "l")
biplot(mypr, scale = 0)


# Extract PC Scores
str(mypr)
mypr$x
dataPC<-cbind(data,mypr$x[,1:2])


# Eigenvalues/Scree plot
ev<-mypr$sdev^2
ev #eigenvalues
ev1<-round(ev/sum(ev)*100,1)
barplot(ev1, main= "Scree Plot", xlab= "Principal Compent", ylab= "Percent Variation")


# Plot with ggfortify
autoplot(mypr, loadings = TRUE, loadings.label = TRUE,
         data = dataPC, repel = TRUE)




# State dependence and long/short term repeatabilty----
setwd("C:/Users/Nick/Desktop/PEFA/P")
data<-read.csv("nest_defense_test_final.csv")
library(lme4)
library(performance)
library(MASS)
library(arm)
library(broom)
library(tidyverse)
require(moments)


hist(data$MinDisRaw)
hist(data$MinDis)


#Transform
data$Min<- -1*(log10(max(data$MinDis + 1)-data$MinDis))
hist(data$Min)



# Year as a factor
data$Year2<-as.factor(data$Year2)


#Repeatability of aggression model

require(lme4)

m3<-lmer(Min ~ -1+Sex+ NestStage+  + Year2 +(1|ID) +(1|ID_Series),
         data=data)
summary(m3)
hist(resid(m3))
plot(m3)
plot(resid(m3))



#simulated parameters + repeatability

library(arm)
library(MCMCglmm)
smod<-sim(m3,1000)
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
###residual variance
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)

### Long-term repeatability
r<-bvar/(rvar+bvar +bvar2)
posterior.mode(r)
HPDinterval(r)  ##repeatability 

###Short-term
r1<-(bvar +bvar2)/(rvar+bvar +bvar2)
posterior.mode(r1)
HPDinterval(r1)  ##repeatability





#### Dis-assortative mating----

#load data
data<-read.csv("assortative_mating_final.csv")

#option 1 transformation
data$Male<- -1*(log10(max(data$Male + 1)-data$Male))
hist(data$Male)


data$Female<- -1*(log10(max(data$Female + 1)-data$Female))
hist(data$Female)


##Model 1 Female + Male for 2018
require(MCMCglmm)



#subset for 2018
data2018<-subset(data, data$Year=="2018")

#create matrix for prior
phen.var<-matrix(c(var(data2018 $Male,na.rm=TRUE),0,0,
                   var(data2018 $Female,na.rm=TRUE)),2,2)

# Prior specification for a bivariate
# model with a single random effect
prior2.1<-list(G=list(G1=list(V=phen.var/2,n=2)),
               R=list(V=phen.var/2,n=2))


#Dis-assortative model 2018
mod.12 <- MCMCglmm(cbind(Male, Female) ~ (trait-1),  
                   random = ~us(trait):NestID ,
                   rcov = ~us(trait):units, 
                   family = c("gaussian", "gaussian"),
                   data=data2018, 
                   prior = prior2.1, 
                   verbose = TRUE,
                   nitt=590000,thin=500,burnin=90000)

plot(mod.12)
summary(mod.12)


posteriors<-as.mcmc(mod.12$VCV)
posterior.mode(posteriors)

# correlations
pair.correlation2.1<-posteriors[,"traitFemale:traitMale.NestID"]/
  sqrt(posteriors[,"traitFemale:traitFemale.NestID"]*
         posteriors[,"traitMale:traitMale.NestID"])

# correlations (residuals)
residual.correlation2.1<-posteriors[,"traitFemale:traitMale.units"]/
  sqrt(posteriors[,"traitFemale:traitFemale.units"]*
         posteriors[,"traitMale:traitMale.units"])


# 2018 Correlations as mode with CI

posterior.mode(pair.correlation2.1)
HPDinterval(pair.correlation2.1)

posterior.mode(residual.correlation2.1)
HPDinterval(residual.correlation2.1)

#bayesian p-value for intercept since 95% CrI overlapped zero-Among nest
#the proportion of estimates that are >0 if the point estimate is negative, or where the point estimate is positive, the proportion of estimates that is < 0
countsI=ifelse(pair.correlation2.1<0, 1,0) #number of positive estimates for 2nd covariance estimate (i.e., the among-individual correlation). For the within-indivdiual correlation, you would use [,4]
sum(countsI)
countsI2=ifelse(pair.correlation2.1>0, 1,0 ) #number of negative estimates for 1st fixed effect (i.e., intercept)
sum(countsI2)
###bayesia p = sum(counts)/#simultations
sum(countsI)/(sum(countsI)+sum(countsI2))


#bayesian p-value for intercept since 95% CrI overlapped zero-Among nest
#the proportion of estimates that are >0 if the point estimate is negative, or where the point estimate is positive, the proportion of estimates that is < 0
countsI=ifelse(residual.correlation2.1<0, 1,0) #number of positive estimates for 2nd covariance estimate (i.e., the among-individual correlation). For the within-indivdiual correlation, you would use [,4]
sum(countsI)
countsI2=ifelse(residual.correlation2.1>0, 1,0 ) #number of negative estimates for 1st fixed effect (i.e., intercept)
sum(countsI2)
###bayesia p = sum(counts)/#simultations
sum(countsI)/(sum(countsI)+sum(countsI2))


#Dis-assortative 2019

require(MCMCglmm)

#subset for 2019
data2019<-subset(data, data$Year=="2019")


#create matrix for prior


phen.var1<-matrix(c(var(data2019 $Male,na.rm=TRUE),0,0,
                    var(data2019 $Female,na.rm=TRUE)),2,2)

# prior specification for a bivariate
# model with a single random effect
prior2.12<-list(G=list(G1=list(V=phen.var1/2,n=2)),
                R=list(V=phen.var/2,n=2))
#Dis-assortative 2019 model

mod.22 <- MCMCglmm(cbind(Male, Female) ~ (trait-1),  
                   random = ~us(trait):NestID ,
                   rcov = ~us(trait):units, 
                   family = c("gaussian", "gaussian"),
                   data=data2019, 
                   prior =prior2.12, 
                   verbose = TRUE,
                   nitt=590000,thin=500,burnin=90000)

plot(mod.22)
summary(mod.22)


#####correlations for figures-Run these in order to run graph below
posteriors212<-as.mcmc(mod.22$VCV)


# correlations
pair.correlation2.12<-posteriors212[,"traitFemale:traitMale.NestID"]/
  sqrt(posteriors212[,"traitFemale:traitFemale.NestID"]*
         posteriors212[,"traitMale:traitMale.NestID"])

# correlations (residuals)
residual.correlation2.12<-posteriors212[,"traitFemale:traitMale.units"]/
  sqrt(posteriors212[,"traitFemale:traitFemale.units"]*
         posteriors212[,"traitMale:traitMale.units"])





# 2019 Correlations as mode with CI

posterior.mode(pair.correlation2.12)
HPDinterval(pair.correlation2.12)

posterior.mode(residual.correlation2.12)
HPDinterval(residual.correlation2.12)





#Within-nest p-value
#bayesian p-value for intercept since 95% CrI overlapped zero-Among nest
#the proportion of estimates that are >0 if the point estimate is negative, or where the point estimate is positive, the proportion of estimates that is < 0
countsI=ifelse(residual.correlation2.12<0, 1,0) #number of positive estimates for 2nd covariance estimate (i.e., the among-individual correlation). For the within-indivdiual correlation, you would use [,4]
sum(countsI)
countsI2=ifelse(residual.correlation2.12>0, 1,0 ) #number of negative estimates for 1st fixed effect (i.e., intercept)
sum(countsI2)
###bayesia p = sum(counts)/#simultations
sum(countsI)/(sum(countsI)+sum(countsI2))


#between nest p-value
#bayesian p-value for intercept since 95% CrI overlapped zero-Among nest
#the proportion of estimates that are >0 if the point estimate is negative, or where the point estimate is positive, the proportion of estimates that is < 0
countsI=ifelse(pair.correlation2.12>0, 1,0) #number of positive estimates for 2nd covariance estimate (i.e., the among-individual correlation). For the within-indivdiual correlation, you would use [,4]
sum(countsI)
countsI2=ifelse(pair.correlation2.12<0, 1,0 ) #number of negative estimates for 1st fixed effect (i.e., intercept)
sum(countsI2)
###bayesia p = sum(counts)/#simultations
sum(countsI)/(sum(countsI)+sum(countsI2))



### create data frame and within-subject center for plot


#create data frame subtracting subjects mean from each observation (i.e., cw2= within-subject centered)
require(dplyr)

d2018 <- data2018 %>% group_by(NestID) %>% 
  mutate(
    male_mean = mean(Male, na.rm = T),
    female_mean= mean(Female, na.rm = T)
  ) %>% ungroup() %>% 
  mutate(
    male_c2 = scale(Male, center = T, scale = F),
    female_c2 = scale(Female, center = T, scale = F),
    
    male_cw2 = Male - male_mean,
    female_cw2 = Female - female_mean,
    
    male_cb2 = male_c2 - male_cw2,
    female_cb2 = female_c2 - female_cw2
  )

d2019 <- data2019 %>% group_by(NestID) %>% 
  mutate(
    male_mean = mean(Male, na.rm = T),
    female_mean= mean(Female, na.rm = T)
  ) %>% ungroup() %>% 
  mutate(
    male_c2 = scale(Male, center = T, scale = F),
    female_c2 = scale(Female, center = T, scale = F),
    
    male_cw2 = Male - male_mean,
    female_cw2 = Female - female_mean,
    
    male_cb2 = male_c2 - male_cw2,
    female_cb2 = female_c2 - female_cw2
  )


#load package
require(ggplot2)




# Calculate regression lines from the model fit (co)variances
between_slope2018 <- mean(mod.12$VCV[,"traitMale:traitFemale.NestID"]/
                            mod.12$VCV[,"traitMale:traitMale.NestID"])
between_slope2019 <- mean(mod.22$VCV[,"traitMale:traitFemale.NestID"]/
                            mod.22$VCV[,"traitMale:traitMale.NestID"])
within_slope2018 <- mean(mod.12$VCV[,"traitMale:traitFemale.units"]/
                           mod.12$VCV[,"traitMale:traitMale.units"])
within_slope2019 <- mean(mod.22$VCV[,"traitMale:traitFemale.units"]/
                           mod.22$VCV[,"traitMale:traitMale.units"])



#Within-nest plot 2018

Within_nest2018<-ggplot(d2018, aes(x = male_cw2, y = female_cw2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept =0 , slope = within_slope2018,size=1)+
  labs(x = "",
       y = "") +  scale_x_continuous(limits = c(-1.5,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic()



Within_nest2018





#within nest plot 2019


Within_nest2019<-ggplot(d2019, aes(x = male_cw2, y = female_cw2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept =0 , slope = within_slope2019,size=1)+
  labs(x = "",
       y = "") +  scale_x_continuous(limits = c(-1,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic()



Within_nest2019


#Between-nest plot 2018

Between_nest2018<-ggplot(d2018, aes(x = male_cb2, y = female_cb2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept = 0, slope =between_slope2018, size=1)+
  labs(x = " ",
       y = " ") +  scale_x_continuous(limits = c(-1.5,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic()


Between_nest2018








#Between-nest plot 2019

Between_nest2019<-ggplot(d2019, aes(x = male_cb2, y = female_cb2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept = 0, slope =between_slope2019, size=1)+
  labs(x = " ",
       y = " ") +  scale_x_continuous(limits = c(-1.5,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic()


Between_nest2019


#combine plots
require(gridExtra)
require(cowplot)


g4<-grid.arrange(Between_nest2018, top= "Between-nest, 2018")

g5<-grid.arrange(Between_nest2019, top="Between-nest, 2019")

g6<-grid.arrange(Within_nest2018, top="Within-nest, 2018")

g7<-grid.arrange(Within_nest2019, top="Within-nest, 2019")


#plot dis-assortative

final_plot<-grid.arrange(g4,g6,g5,g7, left= "Female nest defense (Log-tranformed)", 
                         bottom= "Male nest defense (Log-transformed)")






#Fitness models----

# Probablity of successfully fledging
#load data
setwd("C:/Users/Nick/Desktop/PEFA/P")
data<-read.csv("Fitness_FT_final.csv")

#Transform
data$Aggression<- -1*(log10(max(data$MinDis + 1)-data$MinDis))
hist(data$Aggression)



# Year as a factor
data$Year2<-as.factor(data$Year2)


require(lme4)

#Probality of sucessfully fledging fit as a binary response
mod.31<-glmer(Fledge_pro~ -1+Sex:Year2:Aggression  + (1|SiteID) + (1|ID)  , data= data, family= binomial)

hist(resid(mod.31))
summary(mod.31)



# Simulate to get effect sizes of parameters
library(arm)
require(MCMCglmm)
smod<-sim(mod.31, 1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))
##Between Site variance, SiteID
bID<-smod@ranef$SiteID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
require(MCMCglmm)
bvar<-as.mcmc(bvar)
posterior.mode(bvar)## mode of the distribution
HPDinterval(bvar)
##Between Site variance, ID
bID2<-smod@ranef$ID[,,1]
bvar2<-as.vector(apply(bID2, 1, var)) ##between individual variance posterior distribution
require(MCMCglmm)
bvar2<-as.mcmc(bvar2)
posterior.mode(bvar2)## mode of the distribution
HPDinterval(bvar2)
###residual variance is equal to 1 



# Number fledged model

#load data
data<-read.csv("Fitness_FT_final.csv")

#Transform
data$Aggression<- -1*(log10(max(data$MinDis + 1)-data$MinDis))
hist(data$Aggression)



# Year as a factor
data$Year2<-as.factor(data$Year2)


hist(data$NumberFledged)
data<-subset(data, data$NumberFledged>0)   


mod.10<-lmer(NumberFledged~ -1+Sex:Year2:Aggression + (1|SiteID), data=data)

hist(resid(mod.10))
summary(mod.10)




# Simulate to get effect sizes of parameters
library(arm)
require(MCMCglmm)
smod<-sim(mod.10, 1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))
##Between Site variance, SiteID
bID<-smod@ranef$SiteID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
require(MCMCglmm)
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)
###residual variance  
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)





#Chick Mass model


#load data
data<-read.csv("Fitness_FT_final.csv")

#Transform
data$Aggression<- -1*(log10(max(data$MinDis + 1)-data$MinDis))
hist(data$Aggression)



# Year as a factor
data$Year2<-as.factor(data$Year2)




hist(data$MassAvg)
data<-subset(data, data$MassAvg>0)   


mod.10<-lmer(MassAvg~ -1+Sex:Year2:Aggression + (1|SiteID)  , data=data)

hist(resid(mod.10))
summary(mod.10)

# Simulate to get effect sizes of parameters
library(arm)
require(MCMCglmm)
smod<-sim(mod.10, 1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))
##Between Site variance, SiteID
bID<-smod@ranef$SiteID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
require(MCMCglmm)
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)
###residual variance is equal to 1 
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)


