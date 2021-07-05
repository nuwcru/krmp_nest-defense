#format data
#wd
setwd("C:/Users/nickg/OneDrive/Desktop/P")
data<-read.csv("nest_defense_test_final.csv")

library(tidyr)
library(dplyr)


#load data
data<-read.csv("Fitness_FT_final.csv")

data1<-as_tibble(data)

d2<-filter(data1, Sex=="Male")
d2<-select(d2, NestID, MinDisRaw,MinDis, Year, Date, SiteID, NestStage)          
d2<-rename(d2, male_raw=MinDisRaw,Male= MinDis, NestID=NestID)

d3<-filter(data1, Sex=="Female")
d3<-select(d3, NestID, MinDisRaw,MinDis, Year, Date, SiteID, NestStage)          
d3<-rename(d3, female_raw=MinDisRaw,Female= MinDis, NestID=NestID)

d5<-full_join(d2,d3, key=NestID)
d7<-subset(d5, Male_Aggression>-1, Female_Aggression>-1)
data<-na.omit(d5)

library(writexl)

write_xlsx(x = data, path = "assortative_mating_final.xlsx", col_names = TRUE)




#load data
data<-read.csv("assortative_mating_final.csv")

#option 1 transformation
data$Male<- -1*(log10(max(data$Male + 1)-data$Male))
hist(data$Male)

data$Female<- -1*(log10(max(data$Female + 1)-data$Female))
hist(data$Female)







##Model 1... Female + Male for 2018----
require(MCMCglmm)




#subset for 2018
data2018<-subset(data, data$Year=="2018")

#create matrix for prior
phen.var<-matrix(c(var(data2018 $Male,na.rm=TRUE),0,0,
                   var(data2018 $Female,na.rm=TRUE)),2,2)

# an example prior specification for a bivariate
# model with a single random effect
prior2.1<-list(G=list(G1=list(V=phen.var/2,n=2)),
               R=list(V=phen.var/2,n=2))

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






#####correlations for figures-Run these in order to run graph below
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


# 2018 Correlations as mode with CI----

posterior.mode(pair.correlation2.1)
HPDinterval(pair.correlation2.1)

posterior.mode(residual.correlation2.1)
HPDinterval(residual.correlation2.1)



###script to quanitify means and quantiles from actual values rather than extrapolating (as in previous code)
###rounding to 2 decimals
re<-round(apply(mod.12$Sol,2,mean),2)%/%
round(apply(mod.12$Sol,2, quantile, c(0.025, 0.975)),2)

round(apply(mod.12$VCV,2,mean),2)
round(apply(mod.12$VCV,2, quantile, c(0.025, 0.975)),2)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is the within indivdiual covariance matrix
c1 <- posterior.cor(mod.12$VCV[,1:4])#Between nest
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)

c2 <- posterior.cor(mod.12$VCV[,5:8]) #Within nest
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)


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









###Model 2..Male + Female (2019)----


require(MCMCglmm)

#subset for 2019
data2019<-subset(data, data$Year=="2019")


#create matrix for prior


phen.var1<-matrix(c(var(data2019 $Male,na.rm=TRUE),0,0,
                   var(data2019 $Female,na.rm=TRUE)),2,2)

# an example prior specification for a bivariate
# model with a single random effect
prior2.12<-list(G=list(G1=list(V=phen.var1/2,n=2)),
               R=list(V=phen.var/2,n=2))


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





# 2019 Correlations as mode with CI----

posterior.mode(pair.correlation2.12)
HPDinterval(pair.correlation2.12)

posterior.mode(residual.correlation2.12)
HPDinterval(residual.correlation2.12)



###script to quanitify means and quantiles from actual values rather than extrapolating (as in previous code)
###rounding to 2 decimals
re<-round(apply(mod.22$Sol,2,mean),2)%/%
  round(apply(mod.22$Sol,2, quantile, c(0.025, 0.975)),2)

round(apply(mod.22$VCV,2,mean),2)
round(apply(mod.22$VCV,2, quantile, c(0.025, 0.975)),2)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is the within indivdiual covariance matrix
c1 <- posterior.cor(mod.22$VCV[,1:4])#Between nest
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)

c2 <- posterior.cor(mod.22$VCV[,5:8]) #Within nest
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)



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






 #Save plot
 ggsave("Correlation_All.jpg", plot = last_plot(), device= "jpeg", units="mm", dpi= "retina", limitsize = TRUE)
 

 
### create data frame and within-subject center for plot----


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



data <- data %>% group_by(NestID) %>% 
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


#Within-nest plot  2019
coef(lm(female_cw2 ~ male_cw2, data = d2019))

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
coef(lm(female_cw2 ~ male_cw2, data = d2018))
Within_nest2018<-ggplot(d2018, aes(x = male_cw2, y = female_cw2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept =0 , slope = within_slope2018,size=1)+
  labs(x = "",
       y = "") +  scale_x_continuous(limits = c(-1.5,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic() +  ggtitle("2018") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))

  

Within_nest2018





#within nest plot 2019

coef(lm(female_cw2 ~ male_cw2, data = d2019))
Within_nest2019<-ggplot(d2019, aes(x = male_cw2, y = female_cw2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept =0 , slope = within_slope2019,size=1)+
  labs(x = "",
       y = "") +  scale_x_continuous(limits = c(-1,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic() +  ggtitle("2019") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))



Within_nest2019


wnp<-grid.arrange(Within_nest2018, Within_nest2019, ncol=2)

#Between-nest plot 2018




coef(lm(male_cb2 ~ female_cb2, data = d2018))



Between_nest2018<-ggplot(d2018, aes(x = male_cb2, y = female_cb2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept = 0, slope =between_slope2018, size=1)+
  labs(x = " ",
       y = " ") +  scale_x_continuous(limits = c(-1.5,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic() + ggtitle("2018") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))



Between_nest2018








#Between-nest plot 2019
coef(lm(female_mean ~ male_mean, data = d2019))

Between_nest2019<-ggplot(d2019, aes(x = male_cb2, y = female_cb2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept = 0, slope =between_slope2019, size=1)+
  labs(x = " ",
       y = " ") +  scale_x_continuous(limits = c(-1.5,1.25)) + 
  scale_y_continuous( limits = c(-1,1)) +
  theme_classic() +  ggtitle("2019") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))



Between_nest2019



bnp<-grid.arrange(Between_nest2018, Between_nest2019, ncol=2)

#combine plots
require(gridExtra)
require(cowplot)
g1<-plot_grid(Between_nest2019, Between_nest2018,Within_nest2019, Within_nest2018,
               ncol=2, align='vh', axis='l')


g4<-grid.arrange(Between_nest2018, top= "Between-nest, 2018")

g4
g5<-grid.arrange(Between_nest2019, top="Between-nest, 2019")

g6<-grid.arrange(Within_nest2018, top="Within-nest, 2018")

g7<-grid.arrange(Within_nest2019, top="Within-nest, 2019")


#plot----

final_plot<-grid.arrange(g4,g6,g5,g7, left= "Female nest defense (Log-tranformed)", 
             bottom= "Male nest defense (Log-transformed)")



#save plot
ggsave("Within_Between_centered.jpg", plot = final_plot, device= "jpeg", units="mm", dpi= "retina")
 








require(tidyr)

df_mcmc_cors <- data_frame(Traits = c("Between-nest, 2018", 
                                      "Between-nest, 2019"),
                           Estimate = c(mean(pair.correlation2.1),
                                        mean(pair.correlation2.12)),
                           Lower = c(HPDinterval(pair.correlation2.1)[,"lower"],
                                     HPDinterval(pair.correlation2.12)[,"lower"]),
                           Upper = c(HPDinterval(pair.correlation2.1)[,"upper"],
                                     HPDinterval(pair.correlation2.12)[,"upper"]))

require(ggplot2)
g1<-ggplot(df_mcmc_cors, aes(x = Traits, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) +
  geom_hline(yintercept = 0,
             alpha = 0.3) +
  scale_x_discrete(limits = c("Between-nest, 2018",
                              "Between-nest, 2019")) +
  labs(x = "",
       y = "Correlation (Estimate +/- 95% CIs)", title= "") +
  ylim(-1,1) +
  coord_flip() +
  theme_classic()

g1


df_mcmc_cors1 <- data_frame(Traits = c("Within-nest, 2018", 
                                      "Within-nest, 2019"),
                           Estimate = c(mean(residual.correlation2.1),
                                        mean(residual.correlation2.12)),
                           Lower = c(HPDinterval(residual.correlation2.1)[,"lower"],
                                     HPDinterval(residual.correlation2.12)[,"lower"]),
                           Upper = c(HPDinterval(residual.correlation2.1)[,"upper"],
                                     HPDinterval(residual.correlation2.12)[,"upper"]))

require(ggplot2)
g2<-ggplot(df_mcmc_cors1, aes(x = Traits, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) +
  geom_hline(yintercept = 0,
             alpha = 0.3) +
  scale_x_discrete(limits = c("Within-nest, 2018",
                              "Within-nest, 2019")) +
  labs(x = "",
       y = "Correlation (Estimate +/- 95% CIs)", title = "") +
  ylim(-1,1) +
  coord_flip() + 
  theme_classic() 


require(gridExtra)
grid.arrange(g1,g2, nrow = 1)


