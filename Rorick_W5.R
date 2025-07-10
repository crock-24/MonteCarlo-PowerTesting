#7.1 Estimate the MSE of the level k trimmed means for random samples of size 20 generated 
#from a standard Cauchy distribution. (The target parameter θ is the center or median; 
#the expected value does not exist.) Summarize the estimates of MSE in a table for 
#k = 1,2,...,9
set.seed(5)
n <- 20
K <- n/2 - 1
m <- 1000
mse.se <- matrix(0, K, 2)

trimmed.mse <- function(n, m, k)
{
    # MC est of mse for k-level trimmed mean
    tmean <- numeric(m)
    for (i in 1:m) 
    {
        x <- sort(rcauchy(n))
        tmean[i] <- sum(x[(k+1):(n-k)]) / (n-2*k)
    }
    mse.est <- mean((tmean-0)^2)
    se.mse <- sqrt(var((tmean-0)^2)/m)
    return(c(mse.est, se.mse))
}

for (k in 1:K) 
{
    mse.se[k, 1:2] <- trimmed.mse(n=n, m=m, k=k)
}

cbind(mse.se, seq(9))

#7.2 Plot the empirical power curve for the t-test in Example 7.9, changing the alternative hypothesis 
#to H1 : µ̸= 500, and keeping the significance level α= 0.05.
n <- 20
m <- 1000
mu0 <- 500
sigma <- 100
mu <- c(seq(450, 650, 10)) #alternatives
M <- length(mu)
power <- numeric(M)
for (i in 1:M) {
    mu1 <- mu[i]
    pvalues <- replicate(m, expr = {
        x <- rnorm(n, mean = mu1, sd = sigma)  #simulate under alternative mu1
        ttest <- t.test(x, alternative = "two.sided", mu = mu0)  
        ttest$p.value } )
    power[i] <- mean(pvalues <= .05)
}
se <- sqrt(power * (1-power) / m)

library(ggplot2)
df <- data.frame(mean=mu, power=power,
                 upper=power+2*se, lower=power-2*se)
ggplot(df, aes(x=mean, y=power)) +
    geom_line() +
    geom_vline(xintercept=500, lty=2) +
    geom_hline(yintercept=c(0,.05), lty=1:2) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2, lwd=1.5)


#7.3 Plot the power curves for the t-test in Example 7.9 for sample sizes 10, 20, 
#30, 40, and 50, but omit the standard error bars. Plot the curves on the same graph, 
#each in a diﬀerent color or diﬀerent line type, and include a legend. Comment on 
#the relation between power and sample size.
n <- c(10, 20, 30, 40, 50)
m <- 1000
mu0 <- 500
sigma <- 100
mu <- c(seq(450, 650, 10)) #alternatives
M <- length(mu)
power <- matrix(0, M, length(n))
for (j in 1:length(n)){
    for (i in 1:M) {
        mu1 <- mu[i]
        pvalues <- replicate(m, expr = {
            x <- rnorm(n[j], mean = mu1, sd = sigma)  #simulate under alternative mu1
            ttest <- t.test(x, alternative = "two.sided", mu = mu0)  
            ttest$p.value } )
        power[i, j] <- mean(pvalues <= .05)
    }
}

library(ggplot2)
df <- data.frame(mean=mu, power=power)
ggplot(df, aes(x=mean, y=power)) +
    geom_line(aes(y = power[,1], colour = "n = 10")) + 
    geom_line(aes(y = power[,2], colour = "n = 20")) +
    geom_line(aes(y = power[,3], colour = "n = 30")) +
    geom_line(aes(y = power[,4], colour = "n = 40")) +
    geom_line(aes(y = power[,5], colour = "n = 50")) +
    geom_vline(xintercept=500, lty=2) +
    geom_hline(yintercept=c(0,.05), lty=1:2)

#7.6 Suppose a 95% symmetric t-interval is applied to estimate a mean, but the sample 
#data are non-normal. Then the probability that the confidence interval covers the 
#mean is not necessarily equal to 0.95. Use a Monte Carlo experiment to estimate the 
#coverage probability of the t-interval for random samples of χ2(2) data with sample 
#size n = 20. Compare your t-interval results with the simulation results in Example 
#7.4. (The t-interval should be more robust to departures from normality than 
#the interval for variance.)
set.seed(5)
m = 1000
n = 20
pvalues <- replicate(m, expr = {
    x <- rchisq(n, df = 2)  #simulate 20 samples from χ2(2)
    ttest <- t.test(x, alternative = "two.sided", mu = 2)  #mean of χ2(2) is 2, see if t test captures
    ttest$p.value })
print(mean(pvalues >= .05))

#7.8 Estimate the power of the skewness test of normality against symmetric 
#Beta(α,α) distributions and comment on the results. Are the results diﬀerent for 
#heavy-tailed symmetric alternatives such as t(ν)?
a <- c(1, 2, 3, 4, 5, 10, 50, 100) #holds different values of alpha for beta distribution
pwr.beta <- numeric(length(a)) #holds the power of rejecting normality at that particular alpha
alpha <- 0.1 #confidence in normality
n <- 30
m <- 2500
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3)))) # critical value for the skewness test

sk <- function(x) #computes the sample skewness coeff.
{
    xbar <- mean(x)
    m3 <- mean((x - xbar)^3)
    m2 <- mean((x - xbar)^2)
    return( m3 / m2^1.5 )
}


for (j in 1:length(a))  #for each value of alpha
{
    sktests <- numeric(m)
    for (i in 1:m)        #for each replicate
    {
        x <- rt(n, df)
        sktests[i] <- as.integer(abs(sk(x)) >= cv)
    }
    pwr.beta[j] <- mean(sktests)
}
cbind(a, pwr.beta)

df <- c(1, 2, 3, 4, 5, 10, 50, 100) #holds different degrees of freedom for t distribution
pwr.t <- numeric(length(a)) #holds the power of rejecting normality at that particular df

for (j in 1:length(df))  #for each value of alpha
{
    sktests <- numeric(m)
    for (i in 1:m)        #for each replicate
    {
        x <- rt(n, df[j])
        sktests[i] <- as.integer(abs(sk(x)) >= cv)
    }
    pwr.t[j] <- mean(sktests)
}
cbind(df, pwr.t)

#7.9 Refer to Example 7.16. Repeat the simulation, but also compute the F test of 
#equal variance, at significance level α. = 0.055. Compare the power of the 
#Count Five test and F test for small, medium, and large sample sizes. 
#(Recall that the F test is not applicable for non-normal distributions.)
count5test <- function(x, y) {
    X <- x - mean(x)
    Y <- y - mean(y)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))
    # return 1 (reject) or 0 (do not reject H0)
    return(as.integer(max(c(outx, outy)) > 5))
}

n <- c(10, 20, 30, 40, 50, 100)
m <- 1000
sigma1 <- 1
sigma2 <- 1.5
count5 <- numeric(m)
f <- numeric(m)
power.count5 <- numeric(length(n))
power.f <- numeric(length(n))
for (j in 1:length(n))  #for each value of n
{
    for (i in 1:m)  #for each value of n
    {
        x <- rnorm(n[j], 0, sigma1)
        y <- rnorm(n[j], 0, sigma2)
        count5[i] <- count5test(x, y)
        f[i] <- as.integer(var.test(x, y, alternative = "two.sided")$p.value < 0.055)
    }
    power.count5[j] <- mean(count5)
    power.f[j] <- mean(f)
}

cbind(power.count5, n)
cbind(power.f, n)
