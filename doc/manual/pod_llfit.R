library(stochbb)

#
# Processing model
#
d  <- 500
V  <- gamma(5,30)
S1 <- gamma(10,50)
S2 <- gamma(1,120)
M  <- gamma(1,150) 

# And for the experimental condition 
#   R = V + min(S1, S2) + M.
# This might be considered as a two parallel proceses S1 (identical to the S1 under control
# condition) and S2 which is delayed by d=300 ms but once started, much faster than S1.
# The resulting distribution is not analytic on (0, infinity)! It is identical to the distribution
# of R under control condition (Rc) on the interval (0, d] and diverges for T>d.
Re <- V %+% minimum(S1, affine(S2, 1, d)) %+% M

Nstep <- 1000;
Nsam <- 30000;
# Get some samples
sam <- new(ExactSampler, c(Re))
res <- array(0, c(Nsam,1))
sam$sample(res)

tmin <- 0; tmax <- 2000
t <- seq(tmin, tmax, length.out=Nstep)
tmin2 <- min(res[,1]); tmax2 <- max(res[,1])
t2 <- seq(tmin2, tmax2, length.out=Nstep)
cdf <- array(0, c(Nstep)); cdf2 <- array(0, c(Nstep))
Re$density$evalCDF(tmin, tmax, cdf)
Re$density$evalCDF(tmin2, tmax2, cdf2)

costFunc <- function(theta) {
  # update S2 stage with theta
  S2 <- gamma(1,theta[1]);
  # Re-assemble result
  R <- V %+% minimum(S1, affine(S2, 1, d)) %+% M
  # Evaluate logLikelihood of data given theta
  return( -logLikelihood(R, Nstep, res[,1]) )
}

# Find estimate
opt <- optim(c(100), costFunc, method="L-BFGS-B", lower=c(1), upper=c(Inf))
print(opt)

theta = seq(30, 200, length.out=200)
ll <- sapply(theta, costFunc)
ll <- ll-min(ll)

pdf("likelihoodProf.pdf", width=4.5, height=3.5)
plot(theta, exp(-ll), type="l", xlab="theta (S2)", ylab="")
abline(v=120)
abline(v=opt$par[[1]], lty=2)
dev.off()
