#data
data<-read.csv("data/nest_defense_test_final.csv")

#exploratory view of stoops & Min Dis
#prior (inverse gamma prior) for 1 random effect
prior.bivar<-list(R=list(V=diag(2),nu=1.002),G=list(G1=list(V=diag(2), nu=1.002)))


hist(data$Dives)
hist(data$MinR)
require(MCMCglmm)
#model 1- bivarite model 
yz.bivarC.mcmc<-
  MCMCglmm(cbind(Dives,MinR)~(trait-1) ,
           random=~us(trait):ID, rcov=~us(trait):units ,
           family=c("poisson","poisson"), prior=prior.bivar,
           nitt=1300000,thin=1000,burnin=300000, data=data,verbose=TRUE)


summary(yz.bivarC.mcmc)
plot(yz.bivarC.mcmc)

# no autocorrelation since <0.1
autocorr.diag(yz.bivarC.mcmc$Sol)
autocorr.diag(yz.bivarC.mcmc$VCV)

#Repeatability for each mesure of nest defense
# among-individual variance/ (within-individual variance (i.e.,residual) + among-individual variance)
rep.Stoop<-
  yz.bivarC.mcmc$VCV[,1]/ (yz.bivarC.mcmc$VCV[,1]+ yz.bivarC.mcmc$VCV[,5])

rep.Min<-
  yz.bivarC.mcmc$VCV[,4]/ (yz.bivarC.mcmc$VCV[,4]+ yz.bivarC.mcmc$VCV[,8])

# Variance explained by ID (among-individual)
Stoop_am<-yz.bivarC.mcmc$VCV[,1] #stoops
posterior.mode(Stoop_am)

Min_am<-yz.bivarC.mcmc$VCV[,4] #min distance
posterior.mode(Min_am)

# Variance explained by residual (within-individual)
Stoop_res<-yz.bivarC.mcmc$VCV[,5] #stoops
posterior.mode(Stoop_res)

Min_res<-yz.bivarC.mcmc$VCV[,8] #min distance
posterior.mode(Min_res)



#posterior mode for repeatability of each response:
posterior.mode(rep.Stoop) #stoop
posterior.mode(rep.Min) #min distance

#95% credibility intervals for repeatabilities: 
HPDinterval(rep.Stoop) #stoop
HPDinterval(rep.Min) #min distance

#correlation b/n both variables
#highly negatively correlated suggesting peregrines that stoop more often approach closer to the observer
corr_SM <- yz.bivarC.mcmc$VCV[,"traitDives:traitMinR.ID"]/
  (sqrt(yz.bivarC.mcmc$VCV[,"traitDives:traitDives.ID"])*
     sqrt(yz.bivarC.mcmc$VCV[,"traitMinR:traitMinR.ID"]))


plot(corr_SM)
mean(corr_SM) #correlation
HPDinterval(corr_SM)


