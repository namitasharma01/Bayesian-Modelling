##################################################################
# 1. Time series models in Stan
##################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#-------------------------------------------
# (a) Simulate data from AR process
simulate_AR <- function(mu, phi, sigma2, T) {
  x    <- numeric(T)
  x[1] <- mu
  
  for (t in 2:T) {
    epsilon_t <- rnorm(n=1, mean=0, sd=sqrt(sigma2))
    x[t]      <- mu + phi*(x[t-1]-mu) + epsilon_t
  }
  return(x)
}

PHI = seq(-0.5, 1, length=4)
i = 1
par(mfrow = c(2,2))
for (phi in PHI) {
  # Simulate from AR process
  x <- simulate_AR(mu=10, phi=phi, sigma2=2, T=200)
  
  # Plot the realizations
  plot(1:200, x, type="l", col=i, xlab="time", ylab="x", main="Effect of phi on x")
  i <- i+1
}

#-------------------------------------------
# (b) Estimate parameters using MCMC
x <- simulate_AR(mu=10, phi=0.3, sigma2=2, T=200)
y <- simulate_AR(mu=10, phi=0.95, sigma2=2, T=200)

# Stan code to sample from posterior mu, phi and sigma2
stanNormal = '
data {
  int<lower=0> T;
  real x[T];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
}
model {
  mu ~ normal(0, 100);
  sigma2 ~ scaled_inv_chi_square(0.01, 10);
  phi ~ normal(0, 1);
  
  for (t in 2:T)
    x[t] ~ normal(mu + phi*(x[t-1]-mu), sqrt(sigma2));
}'

# Data (synthetic)
N <- length(x)
dataX <- list(T=N, x=x)
dataY <- list(T=N, x=y)

# Fit stan model
burnin <- 1000
niter  <- 2000

# For AR(1)-process 1
fitX <- stan(model_code=stanNormal, data=dataX, warmup=1000, iter=2000)
print(fitX)

# Posterior of mu, phi and sigma2
postDrawsX <- extract(fitX)
# Trace plot of mu, phi and sigma2
traceplot(fitX)

# For AR(1)-process 2
fitY <- stan(model_code=stanNormal, data=dataY, warmup=1000, iter=2000)
print(fitY)

# Posterior of mu, phi and sigma2
postDrawsY <- extract(fitY)
# Trace plot of mu, phi and sigma2
traceplot(fitY)

# Summary Stats
summaryX = summary(fitX, 
                   pars=c("mu", "phi", "sigma2"), 
                   probs = c(0.025, 0.975))$summary

summaryY = summary(fitY, 
                   pars=c("mu", "phi", "sigma2"), 
                   probs = c(0.025, 0.975))$summary

knitr::kable(summaryX)
knitr::kable(summaryY)

# 95% credible interval for mu, phi and sigma2
quanX = round(summaryX[, c(4,5)], 2)
quanY = round(summaryY[, c(4,5)], 2)

paste("\nAR(1)-process for x")
paste("\n95% credible interval for mu = [", quanX[1,1], quanX[1,2], "]")
paste("95% credible interval for phi = [", quanX[2,1], quanX[2,2], "]")
paste("95% credible interval for sigma = [", quanX[3,1], quanX[3,2], "]")

paste("\nAR(1)-process for y")
paste("\n95% credible interval for mu = [", quanY[1,1], quanY[1,2], "]")
paste("95% credible interval for phi = [", quanY[2,1], quanY[2,2], "]")
paste("95% credible interval for sigma = [", quanY[3,1], quanY[3,2], "]")

# Plot of first chain of mu, phi and sigma2 with 95% credible interval
par(mfrow = c(3,2))

muX_chain1 = postDrawsX$mu[1:(niter-burnin)]
hist(muX_chain1, breaks=100, col="gray", main="95% CI for Mu AR(1) process X")
abline(v=quanX[1,1], lty="dashed", col="red")
abline(v=quanX[1,2], lty="dashed", col="red")

muY_chain1 = postDrawsY$mu[1:(niter-burnin)]
hist(muY_chain1, breaks=100, col="gray", main="95% CI for Mu AR(1) process Y")
abline(v=quanY[1,1], lty="dashed", col="red")
abline(v=quanY[1,2], lty="dashed", col="red")

phiX_chain1 = postDrawsX$phi[1:(niter-burnin)]
hist(phiX_chain1, breaks=100, col="gray", main="95% CI for Phi AR(1) process X")
abline(v=quanX[2,1], lty="dashed", col="red")
abline(v=quanX[2,2], lty="dashed", col="red")

phiY_chain1 = postDrawsY$phi[1:(niter-burnin)]
hist(phiY_chain1, breaks=100, col="gray", main="95% CI for Mu AR(1) process Y")
abline(v=quanY[2,1], lty="dashed", col="red")
abline(v=quanY[2,2], lty="dashed", col="red")

sigmaX_chain1 = postDrawsX$sigma[1:(niter-burnin)]
hist(sigmaX_chain1, breaks=100, col="gray", main="95% CI for Sigma2 AR(1) process X")
abline(v=quanX[3,1], lty="dashed", col="red")
abline(v=quanX[3,2], lty="dashed", col="red")

sigmaY_chain1 = postDrawsY$sigma[1:(niter-burnin)]
hist(sigmaY_chain1, breaks=100, col="gray", main="95% CI for Sigma2 AR(1) process X")
abline(v=quanY[3,1], lty="dashed", col="red")
abline(v=quanY[3,2], lty="dashed", col="red")

# Convergence of the samplers- Trace plot of mu, phi and sigma2
traceplot(fitX)
traceplot(fitY)

# Joint posterior of mu and phi
jointPostX = cbind(mu_x=postDrawsX$mu, phi_x=postDrawsX$phi)
jointPostY = cbind(mu_y=postDrawsY$mu, phi_y=postDrawsY$phi)

# Plot joint posterior
pairs(jointPostX, main="Joint posterior of Mu and Phi in AR process x")
pairs(jointPostY, main="Joint posterior of Mu and Phi in AR process y")

# Yes, we are able to estimate the true values of the parameters for both the AR processes 
# X and Y quite accurately. For AR process x, the 95% credible intervals are narrow indicating 
# less uncertainity in the estimated parameter values. The effective sample size is large 
# which means that the efficiency of MCMC sampling was high for process X. However, for AR 
# process y, the 95% credible intervals for the parameter mu is very large indicating high 
# uncertainity in the estimated value of mean whereas the other two parameters are tightly 
# bound in their intervals. The effective sample sizes are also quite small which tells us 
# that the MCMC sampling was not very efficient for process Y (due to high auto correlation).

# In the joint posterior plot of mu and sigma2 for AR process Y, it can be seen that the 
# uncertainity along the mu axis is quite high which makes the joint distribution skewed. 
# The joint posterior plot for AR process x shows less uncertainity for both parameters 
# and hence a better convergence.

#-------------------------------------------
# (c) Estimate Poisson model using Stan
campy <- read.table(
  file="C:/Users/namit/Downloads/Bayesian Learning/R files/Lab4/campy.dat",
  header=TRUE)

# Stan code to sample from posterior of x[t] 
stanPoisson = '
data {
  int<lower=0> T;
  int c[T];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
  real x[T];
}
model {
  mu ~ normal(0, 100);
  sigma2 ~ scaled_inv_chi_square(0.01, 10); 
  phi ~ normal(0, 1);
  
  for (t in 2:T){
    x[t] ~ normal(mu + phi*(x[t-1]-mu), sqrt(sigma2));
    c[t] ~ poisson(exp(x[t]));
  }
}'

# Data (campy file)
N <- length(campy$c)
data <- list(T=N, c=campy$c)

# Fit stan model to poisson process
burnin <- 1000
niter  <- 2000
fit    <- stan(model_code=stanPoisson, data=data, warmup=1000, iter=2000)
print(fit)

# Posteriors of mu, phi, sigma2, x
postDraws <- extract(fit)

# Posterior of Latent intensity theta
thetaPost <- exp(postDraws$x)

# Mean of posterior of theta
thetaPostMean = apply(thetaPost, 2, mean)

# 95% credible interval for theta
credI = apply(thetaPost, 2, quantile, probs=c(0.025, 0.975))

# Plot data and posterior of theta 
plot(campy$c, type="l", main="Data and Theta Posterior", 
     xlab="Time", ylab="infections")
lines(thetaPostMean, col="blue")
lines(credI[1, ], lty="dashed", col="red")
lines(credI[2, ], lty="dashed", col="red")
legend("topleft", c("Data", "Posterior Mean", "Credible Interval"), 
       col=c("black", "blue", "red"), 
       lty=c("solid", "solid", "dashed"))

#-------------------------------------------
# (d) Estimate Poisson model with informative prior using Stan

# Correction------------------------------------------
# Use a stronger prior with more df=100 (Seen more datapoints) and less sigma=1

# Stan code to sample from posterior of x[t] 
stanPoisson = '
data {
  int<lower=0> T;
  int c[T];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
  real x[T];
}
model {
  mu ~ normal(0, 100);
  sigma2 ~ scaled_inv_chi_square(100, 1);
  phi ~ normal(0, 1);
  
  for (t in 2:T){
    x[t] ~ normal(mu + phi*(x[t-1]-mu), sqrt(sigma2));
    c[t] ~ poisson(exp(x[t]));
  }
}'

# Data (campy file)
N <- length(campy$c)
data <- list(T=N, c=campy$c)

# Fit stan model to poisson process
burnin <- 1000
niter  <- 2000
fit    <- stan(model_code=stanPoisson, data=data, warmup=1000, iter=2000)
print(fit)

# Posteriors of mu, phi, sigma2, x
postDraws <- extract(fit)

# Posterior of Latent intensity theta
thetaPost <- exp(postDraws$x)

# Mean of posterior of theta
thetaPostMean = apply(thetaPost, 2, mean)

# 95% credible interval for theta
credI = apply(thetaPost, 2, quantile, probs=c(0.025, 0.975))

# Plot data and posterior of theta 
plot(campy$c, type="l", main="Data and Theta Posterior", 
     xlab="Time", ylab="infections")
lines(thetaPostMean, col="blue")
lines(credI[1, ], lty="dashed", col="red")
lines(credI[2, ], lty="dashed", col="red")
legend("topleft", c("Data", "Posterior Mean", "Credible Interval"), 
       col=c("black", "blue", "red"), 
       lty=c("solid", "solid", "dashed"))

