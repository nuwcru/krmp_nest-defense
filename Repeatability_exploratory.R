#data
data<-read.csv("data/nest_defense_test_final.csv")

#Transform
data$Min<- -1*(log10(max(data$MinRN + 1)-data$MinRN))
hist(data$Min)

#exploratory view of stoops & Min Dis
#prior (inverse gamma prior) for 1 random effect
prior.bivar<-list(R=list(V=diag(2),nu=1.002),G=list(G1=list(V=diag(2), nu=1.002)))


hist(data$Dives)
hist(data$MinR)
require(MCMCglmm)
#model 1
yz.bivarC.mcmc<-
  MCMCglmm(cbind(Dives,MinR)~(trait-1) ,
           random=~us(trait):ID, rcov=~us(trait):units ,
           family=c("poisson","poisson"), prior=prior.bivar,
           nitt=1300000,thin=1000,burnin=300000, data=data,verbose=TRUE)


summary(yz.bivarC.mcmc)
plot(yz.bivarC.mcmc)
autocorr.diag(yz.bivarC.mcmc$Sol)
autocorr.diag(yz.bivarC.mcmc$VCV)
#Repeatability
rep.Stoop<-
  yz.bivarC.mcmc$VCV[,1]/ (yz.bivarC.mcmc$VCV[,1]+ yz.bivarC.mcmc$VCV[,5])
rep.Min<-
  yz.bivarC.mcmc$VCV[,4]/ (yz.bivarC.mcmc$VCV[,4]+ yz.bivarC.mcmc$VCV[,8])

#posterior mode for repeatability of each response:
posterior.mode(rep.Stoop)
posterior.mode(rep.Min)



#95% credibility intervals for repeatabilities: 
HPDinterval(rep.Stoop)
HPDinterval(rep.Min)

#correlation
corr_SM <- yz.bivarC.mcmc$VCV[,"traitDives:traitMinR.ID"]/
  (sqrt(yz.bivarC.mcmc$VCV[,"traitDives:traitDives.ID"])*
     sqrt(yz.bivarC.mcmc$VCV[,"traitMinR:traitMinR.ID"]))
plot(corr_SM)
mean(corr_SM)
HPDinterval(corr_SM)


#Model 3 selection gradients--------------------------------
d18<-subset(data, data$Year2=="0")
#relative fitness
d18$rfit <- d18$NumberFledged/mean(d18$NumberFledged)

d19<-subset(data, data$Year2=="1")
#relative fitness
d19$rfit <- d19$NumberFledged/mean(d19$NumberFledged)


#prior0 = list(R = list(V = diag(c(1,0.0001),2,2), nu = 0.002, fix = 2),
                        # G = list(G1 = list(V = diag(2), nu = 2,
                                     #       alpha.mu = rep(0,2),
             

#prior




da88<-as.data.frame(yz.bivarC.mcmc$VCV)
da89<-as.data.frame(yz.bivarC.mcmc$Sol)
hist(data$MinR)



prior0 <- list(R = list(R1 = list(V = diag(2), nu = 2, covu = TRUE), R2 = list(V = diag(1),
                                                                               nu = 0.002)))

#2018 selection gradient (fledge-yes/no)----
modbi <- MCMCglmm(Fledged ~ MinR - 1  , 
                  random=~us(at.level(Sex,"Male")):ID + us(at.level(Sex,"Female")):ID,
                rcov=~us(at.level(Sex, "Male"):MinR):ID+us(at.level(Sex, "Female"):MinR):ID,
                  ,
                  data = d18, prior = prior0, nitt = 109000,
                  burnin = 9000, thin = 100, verbose = TRUE)


summary(modbi)
plot(modbi)
#selection gradient
est18<-posterior.mode(modbi$VCV[,2]/modbi$VCV[,1]*sqrt(modbi$VCV[,1])/modbi$Sol[,2])

est18
#2019 selection gradient (Fledge-yes/no)----


d19$Sex<-as.factor(d19$Sex)
modbi1 <- MCMCglmm(cbind(Min, Fledged) ~ (trait- 1), 
                  random = ~ us(at.level(Sex, "Male"):trait):ID + 
                    us(at.level(Sex, "Female"):trait):ID, 
                  rcov= ~ us(at.level(Sex, "Male"):trait):units + 
                    us(at.level(Sex, "Female"):trait):units,
                  family= c("gaussian", "categorial"),
                  data = d19, prior = priorb, nitt = 109000,
                  burnin = 9000, thin = 100, verbose = TRUE)
summary(modbi1)
plot(modbi1)








prior0 <- list(R=list(V=diag(2), nu=0.002, fix=2), G=list(G1=list(V=diag(2), nu=0.002),G1=list(V=diag(2), nu=0.002)))
prior0$R$V[2,2]<-0.0001
modbi <- MCMCglmm(cbind(MinR, Fledged) ~ (trait-1), 
                  random = ~us(at.level(Sex, "Male"):trait):ID + us(at.level(Sex, "Female"):trait):ID,
                  rcov = ~idh(trait-1):units, 
                  family = c("gaussian", "categorical"), 
                  data=d19, prior = prior0,nitt=100500, thin=100, burnin=500)

summary(modbi)

est<-posterior.mode(modbi$VCV[,2]/modbi$VCV[,1]*sqrt(modbi$VCV[,1])/modbi$Sol[,2])

est


#selection gradient
#among-individual correaltion b/n Fledge:Min / among-individual correaltion in Min:Min * sqrt among-individual correaltion in Min:Min/ population level fledge value
#Fledged:Min.ID/ Min:Min.ID * sqrt Min:Min.ID/ traitFledge
est19<-posterior.mode(modbi1$VCV[,2]/modbi1$VCV[,1]*sqrt(modbi1$VCV[,1])/modbi1$Sol[,2])
est19

d99<-as.data.frame(modbi1$VCV)

#2018 selection gradient (Number Fledged)----
#subset for nests that didn't fledge at least 1 nestling
d18<-subset(d18, d18$Fledged=="1")



#2019 selection gradient (Number Fledged)----
#subset for nests that didn't fledge at least 1 nestling
d19<-subset(d19, d19$Fledged=="1")
