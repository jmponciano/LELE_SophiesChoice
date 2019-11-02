
# **********************************************************************************************************************

# This R code generates a data set from the linear regression through origin model. These observed data can be used to obtain different predictive densities, estimated, Classical, Bayesian and Evidential for comparison. 

# We can change the priors to see how the densities look like for different priors.
# We can change the new covariate value (Xnew) to see how it affects the predictive densities. Farther from the observed data, wider should be the prediction limits. 

# Date: November 1, 2019

# ***************************************************************************************************************************


# True parameters
par(mfrow=c(1,1))
beta.true = 2
N = 8               # Sample size
X1 = rnorm(N,0,1)    # These are the covariates for the observed data.
Xnew = 2             # This is the covariate value for the new observation (to be predicted). 

# True density for the new observation is given by dnorm(Ynew,Xnew*beta.true,1) plotted as a function of Ynew. The appropriate limits for evaluation are (Xnew*beta.true - 4,Xnew*beta.true + 4 )

LL = Xnew*beta.true - 4
UL = Xnew*beta.true + 4
Ynew = seq(LL,UL,0.1)
True.den = cbind(Ynew,dnorm(Ynew,Xnew*beta.true),1)

Y1 = rnorm(N,X1*beta.true,1)    # Observed responses
beta.hat = sum(X1*Y1)/sum(X1^2)   # This is the usual MLE of the slope parameter.


# Estimated predictive density and true density

Est.den = cbind(Ynew,dnorm(Ynew,Xnew*beta.hat),1)

# Predictive density integrated over the unconditional sampling distribution

sigma.unc = sqrt(1 + (Xnew^2)/(N-2))
UncondPred.den = cbind(Ynew,dnorm(Ynew,Xnew*beta.hat,sigma.unc))

# Predictive density integrated over the conditional sampling distribution

sigma.c = sqrt(1 + (Xnew^2)/sum(X1^2))
CondPred.den = cbind(Ynew,dnorm(Ynew,Xnew*beta.hat,sigma.c))

# Predictive likelihood based predictive density

Pred.den = matrix(0,length(Ynew),2)
for (i in 1:length(Ynew)){
Y.new = Ynew[i]
X1.new = c(X1,Xnew)
Y1.new = c(Y1,Y.new)
beta.new = sum(X1.new*Y1.new)/sum(X1.new^2)

norm.predlike = exp(sum(dnorm(Y1.new,X1.new*beta.new,log=T)) - sum(dnorm(Y1,X1*beta.hat,log=T)))
Pred.den[i,]=cbind(Y.new,norm.predlike)

}

# Figure 2: Different classical predictive densities and evidential predictive density along with the true density of the new observation. This is based on the same observed data. 

plot(Est.den,lty=2,col=2, type="l",xlab="",ylab="",xlim=c(LL,UL),ylim=c(0,max(1.5*True.den[,2])),main="Comparing predictive densities")
par(new=T)
plot(UncondPred.den,lty=3,col=3, type="l",xlab="",ylab="",xlim=c(LL,UL),ylim=c(0,max(1.5*True.den[,2])),main="Comparing predictive densities")
par(new=T)
plot(CondPred.den,lty=4,col=4, type="l",xlab="",ylab="",xlim=c(LL,UL),ylim=c(0,max(1.5*True.den[,2])),main="Comparing predictive densities")
par(new=T)
plot(Pred.den,col=5, type="l",lty = 5,xlab="",ylab="",xlim=c(LL,UL),ylim=c(0,max(1.5*True.den[,2])),main="Comparing predictive densities")
par(new=T)
plot(True.den,type="l",xlab="Ynew",ylab="Predictive density",xlim=c(LL,UL),ylim=c(0,max(1.5*True.den[,2])),main="Comparing predictive densities")
legend(0.5,0.6,legend=c("True","Estimated","Unconditional","Conditional","Evidential"),lty=1:5,col=1:5,cex=0.8)

# Observations:

# It might be worthwhile normalizing the predictive likelihood to get prediction intervals. 
# The information adjustment might reduce the bias as for profile likelihood but at the cost of invariance. Which is more important?
# The overlap between the true and evidential predictive density is the largest



# ****************************************************************************************************

# Bayesian predictive density using different priors

# Prior predictive distribution and sensitivity to the prior

# *****************************************************************************************************


B=100000
X = 2

beta.true = 2
Y.true = rnorm(B, X*beta.true,1)

beta0 = rnorm(B,0,1)
Y1.prior = rnorm(B,X*beta0,1)   # This is the mixed density for Y, mixed over the prior.

beta1 = runif(B,-5,5)
Y2.prior = rnorm(B,X*beta1,1)

LL.prior = min(c(Y1.prior,Y2.prior))-1
UU.prior = max(c(Y1.prior,Y2.prior))-1

# Figure 5: Prior predictive distributions for different priors

plot(density(Y.true),xlim=c(LL.prior,UU.prior),ylim=c(0,0.5),main="Prior predictive density for different priors",xlab="")
par(new=T)
plot(density(Y1.prior),xlim=c(LL.prior,UU.prior),ylim=c(0,0.5),main="",xlab="",col=2,lty=2)
par(new=T)
plot(density(Y2.prior),xlim=c(LL.prior,UU.prior),ylim=c(0,0.5),main="",xlab="Ynew",col=3,lty=3)
legend(-10,0.4,legend=c("True","Prior 1: N(0,1)","Prior 2: U(-5,5)"),lty=1:3,col=1:3,cex=1)


# It is obvious that the prior predictive density is highly sensitive the choice of the prior. No surprise there. Their coverage properties can be quite bad.

# ************************************************************************************************
# Now let us write a JAGS program to get the posterior predictive density.

library(dclone)

# We will use different priors to see how it affects the post-predictive density.

reg.fn = function(){
	for (i in 1:N){
		Y1[i] ~ dnorm(mu[i],1)
		mu[i] <- X1[i]*beta0
	}
	beta0 ~ dnorm(0,1)
}

dat= list(X1=X1,Y1=Y1,N=N)

reg.fit = jags.fit(dat,"beta0",reg.fn)
summary(reg.fit)

# Now compute the post-predictive density and plot on the true density
Bayeshat.1 = as.matrix(reg.fit)

Y.B1 = rnorm(length(Bayeshat.1),X*Bayeshat.1,1)


# Change the prior and get post-predictive density

reg.fn = function(){
	for (i in 1:N){
		Y1[i] ~ dnorm(mu[i],1)
		mu[i] <- X1[i]*beta0
	}
	beta0 ~ dunif(-5,5)
}

dat= list(X1=X1,Y1=Y1,N=N)

reg.fit2 = jags.fit(dat,"beta0",reg.fn)
summary(reg.fit2)

# Now compute the post-predictive density and plot on the true density
Bayeshat.2 = as.matrix(reg.fit2)
Y.B2 = rnorm(length(Bayeshat.2),X*Bayeshat.2,1)

# Figure 4: Plot the posterior distributions corresponding to different priors.  

plot(density(Bayeshat.1),xlab="",xlim=range(c(Bayeshat.1,Bayeshat.2)),ylim=c(0,1.5),main="",col=1,lty=1)
par(new=T)
plot(density(Bayeshat.2),xlim=range(c(Bayeshat.1,Bayeshat.2)),ylim=c(0,1.5),xlab="Slope parameter",col=2,lty=2,main="Posterior distributions of slope under two different priors")
legend(0.5,1.4,legend=c("Prior 1: N(0,1)","Prior 2: U(-5,5)"),lty=1:2,col=1:2,cex=0.7)



# Figure 6: Effect of prior on post-predictive densities under different priors but same observed data. 

plot(density(Y.B1),xlim=c(-1,10),ylim=c(0,0.5),main="",xlab="",col=2,lty=2)
par(new=T)
plot(density(Y.B2),xlim=c(-1,10),ylim=c(0,0.5),main="",xlab="",col=3,lty=3)
par(new=T)
plot(density(Y.true),xlim=c(-1,10),ylim=c(0,0.5),main="Post-predictive density",xlab="Ynew",lty=1,col=1)
legend(0.1,0.5,legend=c("True","Prior 1: N(0,1)","Prior 2: U(-5,5)"),lty=1:3,col=1:3,cex=0.7)


# Plot of the prior, likelihood and posterior distribution for the above data. You can change the prior and the JAGS program to see the effect of the choice of the prior on the posterior. 

B = 100000
beta.prior = runif(B,-2,3)
beta.prior = sort(beta.prior)
like.fn = matrix(0,B,2)
for (i in 1:B){
like.fn[i,] = c(beta.prior[i],exp(sum(dnorm(Y1,X1*beta.prior[i],1,log=T))))
}
like.fn[,2]=like.fn[,2]/max(like.fn[,2])     # Standardize the likelihood function so that it is not too small numerically. Likelihood can be scaled arbitrarily without any effect on the posterior. It is the shape and location that matters.

# Now use the JAGS program to get the posterior distribution.

reg.fn = function(){
	for (i in 1:N){
		Y1[i] ~ dnorm(mu[i],1)
		mu[i] <- X1[i]*beta0
	}
	beta0 ~ dunif(-2,3)
}

dat= list(X1=X1,Y1=Y1,N=N)

reg.fit = jags.fit(dat,"beta0",reg.fn)
summary(reg.fit)
Bayeshat = as.matrix(reg.fit)

# Figure 3: Prior, Likelihood and Posterior

plot(density(beta.prior),xlim=c(-3,3),ylim=c(0,2),col=4,lty=3,main="Prior, Likelihood and Posterior distribution",xlab="Slope parameter",ylab="")
par(new=T)
plot(like.fn,type="l",xlim=c(-3,3),ylim=c(0,2),main="",xlab="",ylab="",col=1,lty=1)
par(new=T)
plot(density(Bayeshat,from=-3,to=3),type="l",xlim=c(-3,3),ylim=c(0,1.6),main="",xlab="",ylab="",col=2,lty=2)
legend(-2,1.5,legend=c("Likelihood","Posterior","Prior"),lty=1:3,col=c(1,2,4),cex=0.7)


 


