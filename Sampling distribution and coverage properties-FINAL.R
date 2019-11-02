
# Coverage properties and contradiction proportions

# Date: November 1, 2019

# FIRST PART OF THIS CODE IS FOR SHOWING THE DIFFERENCE BETWEEN TRUE AND ESTIMATED SAMPLING DISTRIBUTION. IT WILL PLOT A FIGURE SIMILAR TO FIGURE 1 IN THE PAPER.

# SECOND PART OF THE CODE SHOWS THE REPRODUCIBILITY CRISIS THAT ARISES BECAUSE OF INCORRECT INTERPRETATION OF CONFIDENCE INTERVALS

# *************************************************************************************************************************

# Frequentist sampling distributions: True and Estimated

# ************************************************************************************************************************
N = 20   # Sample size
mu1 = 10
sigma1 = 1
B = 1000

# Sampling distribution (True)
out.true = rep(0,B)

for (i in 1:B){
Y = rnorm(N,mu1,sigma1)
out.true[i] = mean(Y)
}

 
true.CI = c(mu1 - 1.68*1/sqrt(N),mu1 + 1.68*1/sqrt(N))

# Parametric bootstrap sampling distribution

Y = rnorm(N,mu1,sigma1)
mu.hat = mean(Y)

# Bootstrap Sampling distribution (Parametric)
out.boot = rep(0,B)

for (i in 1:B){
Y.boot = rnorm(N,mu.hat,sigma1)
out.boot[i] = mean(Y.boot)
}

# Plot the true and bootstrap estimate of the sampling distribution


range.x = range(c(density(out.true)$x,density(out.boot)$x))
 range.y = range(c(density(out.true)$y,density(out.boot)$y))
plot(density(out.true),main="True and estimated sampling distribution",xlab="",ylab="Density",xlim=range.x,ylim=range.y,lty=1,col=1)
par(new=T)
plot(density(out.boot),main="",xlim=range.x,ylim=range.y,col=2,lty=2,ylab="",xlab="Sample means for different samples")
legend(8.5,1.5,legend=c("True","Estimated"),lty=1:2,col=1:2,cex=0.7)

boot.CI = c(mu.hat - 1.68*sd(out.boot),mu.hat + 1.68*sd(out.boot))    # This is slightly different than the true CI.

# ****************************************************************************************************************

# Coverage property (Betting on the CI): Suppose the experiment is replicated R times, how many times would the bootstrap confidence interval cover the true value? This should be approximately equal to the theoretical coverage. 

# *****************************************************************************************************************

R = 1000    # Number of different experimenters
B = 1000  # Number of bootstrap samples
coverage.R = 0
mu.hat = rep(0,R)
boot.CI = matrix(0,R,2)

for (j in 1:R){

# Generate a data set (replication of the experiment) from the true data generating mechanism
Y = rnorm(N,mu1,sigma1)
mu.hat[j] = mean(Y)

# Compute the parametric bootstrap interval using the estimated mean and see if it covers the true mean
for (i in 1:B){
Y.boot = rnorm(N,mu.hat[j],sigma1)
out[i] = mean(Y.boot)
}
boot.CI[j,] = quantile(out,c(0.05,0.95))    # This is a 90% percentile interval.
if (boot.CI[j,1] < mu1 & boot.CI[j,2] > mu1){coverage.R=coverage.R + 1}
}

coverage.R = coverage.R/R
coverage.R                   # This should be approximately 0.90 (the nominal coverage)

# *******************************************************************************************************

# Compute the rate of contradiction using all replicated experiments. (Replication crisis result)

# **********************************************************************************************************


contradict = rep(0,R)
for (j in 1:R){
	contradict[j] = sum(as.numeric(mu.hat[-j] < boot.CI[j,1])) + sum(as.numeric(mu.hat[-j] > boot.CI[j,2]))
}

Contradiction.rate = mean(contradict)/(R-1)
Contradiction.rate      # Notice that this is substantially higher than the nominal error rate of 0.1. 

 






