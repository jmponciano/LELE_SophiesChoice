# Frequentist sampling distributions: True and Estimated, Coverage properties and contradiction proportions


N = 20   # Sample size
mu1 = 10
sigma1 = 1
B = 1000

# Sampling distribution (True)
out = rep(0,B)

for (i in 1:B){
Y = rnorm(N,mu1,sigma1)
out[i] = mean(Y)
}

 range.x = range(density(out)$x)
 range.y = range(density(out)$y)
plot(density(out),main="True and Parametric Bootstrap sampling distribution",xlab="",ylab="Density",xlim=range.x,ylim=range.y)
true.CI = c(mu1 - 1.68*1/sqrt(N),mu1 + 1.68*1/sqrt(N))

# Parametric bootstrap sampling distribution

Y = rnorm(N,mu1,sigma1)
mu.hat = mean(Y)

# Bootstrap Sampling distribution (Parametric)
out = rep(0,B)

for (i in 1:B){
Y.boot = rnorm(N,mu.hat,sigma1)
out[i] = mean(Y.boot)
}

par(new=T)
plot(density(out),main="",xlim=range.x,ylim=range.y,col="red",ylab="",xlab="Sample means for different samples")

boot.CI = c(mu.hat - 1.68*sd(out),mu.hat + 1.68*sd(out))    # This is slightly different than the true CI.

# Coverage property (Betting on the CI): Suppose the experiment is replicated 100 times, how many times would the bootstrap confidence interval cover the true value? This should be approximately equal to the theoretical coverage. 

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
boot.CI[j,] = quantile(out,c(0.05,0.95))    # This is a percentile interval.
if (boot.CI[j,1] < mu1 & boot.CI[j,2] > mu1){coverage.R=coverage.R + 1}
}

coverage.R = coverage.R/R
coverage.R                   # This should be approximately 0.90 (the nominal coverage)



# Let us compute the rate of contradiction using all replicated experiments.

contradict = rep(0,R)
for (j in 1:R){
	contradict[j] = sum(as.numeric(mu.hat[-j] < boot.CI[j,1])) + sum(as.numeric(mu.hat[-j] > boot.CI[j,2]))
}

Contradiction.rate = mean(contradict)/(R-1)
Contradiction.rate      # Notice that this is substantially higher than the nominal error rate of 0.1. 

 






