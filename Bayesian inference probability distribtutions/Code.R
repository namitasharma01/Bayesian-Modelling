########################################################################
# 1. Bernoulli
########################################################################

# a) Posterior distribution of theta
# Sample data
n      <- 20
s      <- 5
f      <- n-s

# Prior distribution parameters
alpha0 <- 2
beta0  <- 2

# Simulation from posterior of theta
N    <- 10000
mean <- numeric()
sd   <- numeric()
for (nDraws in 1:N) {
  theta        <- rbeta(nDraws, alpha0+s, beta0+f)
  mean[nDraws] <- mean(theta)
  sd[nDraws]   <- sd(theta)  
}

# Plot mean and SD of simulations of theta with increasing number of draws 
plot(mean, ylab="Posterior Mean", xlab="Simulation runs", col="blue", pch=16, cex=0.2)
plot(sd, ylab="Posterior SD", xlab="Simulation runs", col="blue", pch=16, cex=0.2)

# Analytical Mean and SD
alpha_new <- alpha0+s
beta_new  <- beta0+f
true_mean <- alpha_new / (alpha_new+beta_new)
true_SD   <- sqrt(alpha_new*beta_new / ((alpha_new+beta_new)^2 * (alpha_new+beta_new+1)))

cat("Analytical mean of posterior of theta = ", 
    true_mean,
    "\n\nAnalytical SD of posterior of theta = ",
    true_SD)

# b) Posterior probability
theta     <- rbeta(N, alpha_new, beta_new)
prob      <- sum(theta>0.3) / length(theta)
true_prob <- pbeta(q=0.3, alpha_new, beta_new, lower.tail=FALSE)

cat("Pr(theta > 0:3|y) using Simulation = ", 
    prob,
    "\n\nPr(theta > 0:3|y) using true distribution = ",
    true_prob)

# c) Posterior distribution of log-odds
theta <- rbeta(N, alpha_new, beta_new)
phi   <- log(theta / (1-theta))

# Distribution of phi posterior using hist() function
hist(phi, breaks=100, probability=TRUE)

# Density plot of phi posterior using density() function
phi_pdf <- density(phi)
plot(phi_pdf$x, phi_pdf$y, xlab="phi", ylab="Density", type="l")

########################################################################
# 2. Log-normal distribution and the Gini coefficient
########################################################################
library("invgamma")

# a) Simulate from posterior of Sigma square
obs <- c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
n   <- length(obs)

# LogNormal posterior from non-Informative prior
logNormal_nonInfoPrior <- function(nDraws, data, mu) {
  n      <- length(data)
  tao_sq <- sum((log(data)-mu)^2)/n
  
  # Draw from chisq(n) since we are not losing any df in calculating the mean, 
  # i.e. we are given the mean  
  X      <- rchisq(nDraws, n)     
  
  # Draw from posterior of sig_sq ~ Inv-chisq(n, tao_sq)
  sig_sq <- n*tao_sq /X   
  
  return(sig_sq)
}

# Empirical posterior distribution of sig_sq  
postdraws <- logNormal_nonInfoPrior(nDraws=10000, data=obs, mu=3.7)
hist(postdraws, breaks=100, probability=TRUE)

# Density plot of sig_sq posterior using density() function
sigsq_pdf <- density(postdraws)
plot(sigsq_pdf$x, sigsq_pdf$y, type="l")

# Theoretical posterior distribution of sig_sq
sigsq_true <- rinvchisq(10000, n)
hist(sigsq_true, breaks=100, probability=TRUE)

# b) Gini Coefficient 
gini_coeff <- function(sig_sq) {
  G <- 2 * pnorm(sqrt(sig_sq/2), 0, 1) - 1
  return(G)
}
G <- gini_coeff(postdraws)
hist(G, breaks=100, probability=TRUE)

# Density plot of G posterior using density() function
G_pdf <- density(G)
plot(G_pdf$x, G_pdf$y, type="l", xlab="Gini Coefficient", ylab="Density")

# c) 90% Credible Interval and Highest Posterior Density for G
credInterval <- function(credible=0.9) {
  eq_tail  <- (1-credible)/2                     # tail region
  tail     <- eq_tail*length(G)                  # tail % of the posterior samples
  CredI    <- G[order(G)][tail:(length(G)-tail)] # credible interval
  CredI_L  <- CredI[1]                           # lower limit  
  CredI_U  <- CredI[length(CredI)]               # upper limit              
  
  return(c(CredI_L, CredI_U))
}

CredI <- credInterval(credible=0.9)
cat("90% equal tail interval = [", CredI[1], ":", CredI[2], "]")

# c) Highest Posterior Density Interval for G
HPD <- function(emp_pdf, prob=0.9) {
  # Mathematical integration of the empirical density curve can be approximated using Riemann sum 
  x    <- emp_pdf$x                   # 512 points where the density is estimated
  y    <- emp_pdf$y                   # 512 density values estimated for x
  dx   <- emp_pdf$x[2]-emp_pdf$x[1]   # Spacing or bin size 
  C    <- sum(emp_pdf$y) * dx         # Normalizing constant C = Riemann sum approx. of area under the curve ~ very close to 1
  mode <- x[which.max(y)]             # Mode of the density curve
  
  # Area under the curve to the right of x=a and left of x=b or Pr(a<x>b)
  prob_ab <- function(a, b) {
    p.unscaled <- sum(y[x > a & x < b]) * dx
    p.scaled   <- p.unscaled / C
    
    return(p.scaled)
  }
  
  # Cost function to find the interval that contains prob probability 
  O <- function(d) {
    p.area <- prob_ab(a=mode-d, b=mode+d)
    cost <- (prob - p.area)^2 
    
    return(cost)
  }
  
  # Search over range x for the 90% HPD
  d_optim <- optimize(O, range(x))$minimum

  HPD_Int  <- c(mode-d_optim, mode+d_optim)       # Highest posterior density 
  HPD_prob <- prob_ab(a=HPD_Int[1], b=HPD_Int[2]) # Verify Interval contains 90% of the highest posterior density 
  cat(HPD_Int, "covers", 
      HPD_prob, "% of the posterior density curve")
  
  return(HPD_Int)
}

HPD_Int  <- HPD(emp_pdf=G_pdf, prob=0.9)
cat("90% Highest Posterior Density Interval = [", HPD_Int[1], ":", HPD_Int[2], "]")

# Graphical representation of CredI and HPD on G distribution curve
plot(G_pdf$x, G_pdf$y, type="l", xlab="Gini Coefficient", ylab="Density")
abline(v=CredI, col="blue", lty="dashed")
abline(v=HPD_Int, col="red", lty="dashed")
legend("topright", c("Credible Interval", "HPD Interval"), fill=c("blue", "red"))


########################################################################
# 3. Bayesian inference for the concentration parameter in the von 
#    Mises distribution
########################################################################

########################################################################
# 3. Bayesian inference for the concentration parameter in the von 
#    Mises distribution
########################################################################
y <- c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

prior <- function(k, rate=1) {
  return(dexp(k, rate))
}

likelihood <- function(y, mu=2.39, k) {
  n <- length(y)
  l <- exp(k*sum(cos(y-mu))) / (2*pi*besselI(k, nu=0))^n
  return(l)
}

k_grid      <- seq(0, 5, length.out=10000)
k_posterior <- likelihood(y, k=k_grid) * prior(k=k_grid, rate=1)
k_posterior <- k_posterior / sum(k_posterior)

# This plot of prior does not give much insight because this is the just pdf curve 
# of exponential distribution. Because for all x>0, the pdf function will evaluate 
# to a point on this curve.
k_grid      <- seq(0, 5, length.out=10)
plot(k_grid, prior(k=k_grid, rate=1), main="Prior distribution assumption of k",
     xlab="k", ylab="prior(k)", col="red")

k_grid      <- seq(0, 5, length.out=100)
plot(k_grid, prior(k=k_grid, rate=1), main="Prior distribution assumption of k",
     xlab="k", ylab="prior(k)", col="red")

k_grid      <- seq(0, 5, length.out=10000)
plot(k_grid, prior(k=k_grid, rate=1), main="Prior distribution assumption of k",
     xlab="k", ylab="prior(k)=exponential(k, rate=1)", col="red")


# This likelihood curve tells us which is the MLE estimate of k
plot(k_grid, likelihood(y, k=k_grid), main="Likelihood curve of concentration parameter k",
     xlab="k", ylab="Likelihood(y|k, mu)", col="red")

# This plot of posterior again follows the pdf of posterior
plot(k_grid, k_posterior, main="Posterior  curve of concentration parameter k",
     xlab="k", ylab="Posterior(k|y, mu)", col="red")


# Compare with the plot of prior(k) above. This histogram gives us the density  
# ditribution of the k values. K=0 is most dense, k=0.8 is less dense in the prior 
# distribution of k

# If the values of k follow an exponential distribution, then the pdf(k)=0 for 
# most values of k. pdf(k) is in the range [0.8, 0.1] for very few values of k 
# which can be seen in this density curve
hist(prior(k=k_grid, rate=1), breaks=100, probability=TRUE)

kprior_dens <- density(prior(k=k_grid, rate=1))
plot(kprior_dens$x, kprior_dens$y, col="red", type="l")

# This plot of likelihood does not make sense because we want the likelihood value 
# at each value of k to check the maximum likelihood. Not the ditribution of the 
# likelihood values themselves.
hist(likelihood(y, k=k_grid), breaks=100, probability=TRUE)

# Posterior density
hist(k_posterior, breaks=100, probability=TRUE)

kpost_dens <- density(k_posterior)
plot(kpost_dens$x, kpost_dens$y, col="red", type="l")

# It is bimodal

index  <- order(kpost_dens$y)[1]
k_mode <- k_grid[index]

k_posterior[index] # values of posterior pdf at k_mode 

dexp(k_mode)
