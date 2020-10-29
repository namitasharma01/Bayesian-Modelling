######################################################################
# 1. Normal model, mixture of normal model with semi-conjugate prior
######################################################################
library("geoR")

rainfall <- read.table(
  file="C:/Users/namit/Downloads/Bayesian Learning/R files/Lab3/rainfall.dat",
  header=FALSE)

# (a) Normal model

# (i) Gibbs sampler 
fullCondPostMu <- function(mu0, tausq0, sigsq, y) {
  n       <- length(y)
  ybar    <- mean(y)
  
  # Compute mu_n and tausq_n
  tausq_n <- sigsq*tausq0 / (n*tausq0 + sigsq)
  w       <- (n/sigsq) / (n/sigsq + 1/tausq0)
  mu_n    <- w*ybar + (1-w)*mu0
  
  # Sample from full conditional posterior of mu
  mu_post <- rnorm(1, mean=mu_n, sd=sqrt(tausq_n))
  
  return(mu_post)
}

fullCondPostSig <- function(nu0, sigsq0, mu_post, y) {
  n <- length(y)
  
  # Compute nu_n and sigsq_n
  nu_n       <- nu0 + n
  sigsq_n    <- (nu0*sigsq0 + sum((y-mu_post)**2)) / (n+nu0)
  
  # Sample from full conditional posterior of sigsq
  sigsq_post <- geoR::rinvchisq(1, df=nu_n, scale=sigsq_n)
  
  return(sigsq_post)
}

gibbsSampler <- function(iter=100, mu0, tausq0, nu0, sigsq0, sigsq_init, data=rainfall$V1) {
  mu_post       <- numeric(iter)
  sigsq_post    <- numeric(iter)
     
  # First iteration of Gibbs sampler - Start with a random value for mu or sigsq
  #sigsq_init    <- geoR::rinvchisq(1, df=nu0, scale=sigsq0)
  sigsq_post[1] <- sigsq_init
  mu_post[1]    <- fullCondPostMu(mu0=mu0, tausq0=tausq0, sigsq=sigsq_post[1], y=data)

  # Gibbs sampler for remaining iter-1 samples
  for (i in 2:iter) {
    mu_post[i]    <- fullCondPostMu(mu0=mu0, tausq0=tausq0, sigsq=sigsq_post[i-1], y=data)
    sigsq_post[i] <- fullCondPostSig(nu0=nu0, sigsq0=sigsq0, mu_post=mu_post[i], y=data)
  }
  
  return(list(mu=mu_post, sigsq=sigsq_post))
}

# (ii) Convergence of gibbs sampler

iter=1000
gibbSample <- gibbsSampler(iter=iter, mu0=10, tausq0=1, nu0=1, sigsq0=1, sigsq_init=1)

plot(1:iter, gibbSample$mu, type="l", col="blue", 
     xlab="Iterations", ylab="Conditional Posterior Mu")
plot(1:iter, gibbSample$sigsq, type="l", col="blue",
     xlab="Iterations", ylab="Conditional Posterior Sigma Square")

# (b) Mixture normal model


# (c) Graphical comparison






######################################################################
# 2. Metropolis Random Walk for Poisson regression
######################################################################

# (a)

# (b)

# (c)

# (d)




#-NOTES-----------------------------------------------------------
## Multinomial is a generalization of Binomial

# Multinomial Binomial
rmultinom(n=4, size=4, prob=c(0.2, 0.8))

# Binomial - 4 independent binomial trial outcomes i.e. draws from the 
# distribution (same as above but the outcomes for the other class is implicit)
rbinom(n=4, size=4, prob=c(0.2, 0.8))

# Multinomial Bernoulli - 4 independent bernoulli trial outcomes
# (together will be one draw from a binomial distribution)
rmultinom(n=4, size=1, prob=c(0.2, 0.8))

# Bernoulli - 4 independent bernoulli trial outcomes (draws from the distribution)
rbinom(n=4, size=1, prob=c(0.2, 0.8))

# NOTE: n in the distribution formula is size here.
#       n here is the number of times we repeat the experiment of repeated trials

# Individual component of multinomial random vector are binomial and have 
# a binomial distribution 
rmultinom(n=1, size=5, prob=c(0.2, 0.3, 0.5))

# Above is same as 
rbinom(1, size=5, prob=c(0.2))
rbinom(1, size=5, prob=c(0.3))
rbinom(1, size=5, prob=c(0.5))

# NOTE: Each trial of multinomial is a categorical distribution i.e. the outcome 
#       tells us a particular class an object belongs to 
rmultinom(1, size=1, prob=c(0.2, 0., 0.5))


## Dirichlet is a generalization of Beta

# Beta distribution is a distribution of things that could be probabilities. Each 
# random draw gives us a single probability.
# Generalize that to probabilities over k different outcomes. Dirichlet 
# distribution gives us k probabilities (vector over k things that sum to 1)

# Prior:Posterior
# Beta:Bernoulli :: Dirichlet:Multinomial


# Beta and Gamma function 
# Beta function defines binomial coefficients for continuous variables likewise Gamma 
# function defines factorial for continuous variables