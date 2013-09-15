################################################################################
# Test Statistic Functions for Randomization Distributions
#
# A test statistic is a function of three arguments:
# - y: the data, after adjustment
# - z: a treatment indicator (should be numeric, but calling as.numeric(z) is safe)
#
# The function should return a scalar real value
#
# A test statistic may also be a AsymptoticTestStatistic, a class that descends from
# a function. This class includes a `asymptotic` slot which is a backend function for
# the randomizationDistributionEngine. For example, the harmonic.mean.difference test
# statistic can use xBalance for the asymptotic results.
################################################################################

############################### Mean Differences ###############################
fastmean<-function(x){
  ## Skip checks and dispatch etc.
  .Internal(mean(x))
}

mean.difference <- function(ys, z) {
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z)
  return(mean(ys[z]) - mean(ys[!z]))
}

mean.diff.lsfit<-function(ys,z){ 
  ##Try using something that calls compiled code
  ##Gives same answer as mean.difference for balanced blocks and should be like harmonic.mean.difference for unbalanced blocks.
  lsfit(x=model.matrix(ys~z),y=ys,intercept=FALSE)[["coefficients"]][["z"]]
}

mean.diff.vect<-function(ys,z){
 X<-model.matrix(ys~z)
 solve(qr(X, LAPACK=TRUE), ys)[2] ## qr.coef(qr(X,LAPACK=TRUE),ys) ## to handle near singular X
}

############################## Rank Based Functions ##############################

rank.sum <- function(ys, z) {
  z <- as.logical(z)
  ys.ranks <- rank(ys)
  return(sum(ys.ranks[z]))
}

paired.sgnrank.sum <- function(ys, z){
  stopifnot(length(unique(z)) == 2) ##require binary treatment for now

  # we assume data are in block order, eg. 11223344...
  blocks <- gl(length(ys)/2, 2)

  Y<-sapply(split(data.frame(r=ys,z=z),blocks),function(dat){with(dat,r[z==1]-r[z==0])})
  sgn<-as.numeric(Y>0)
  q<-rank(abs(Y))
  sum(sgn*q)
}

# mean.rank.diff.strata<-function(ys, z, blocks){
#   ##average difference in ranks between control and treated obs
#   ##currently weighted equally across strata (i.e. good for pairs not good for other designs)
#   mean(
#        mapply(FUN = function(ys, z) {
#          R<-rank(ys)
#          mean(R*z)-mean(R*(1-z))
#          ##((R*z)/sum(z))-((R*(1-z))/sum(1-z))
#        },
#               ys = split(ys, blocks),
#               z = split(z, blocks))
#        )
# }


############################## Misc Statistics ##############################

odds.ratio <- function(y, z) {
  c11 <- sum(y == 1 & z == 1)
  c10 <- sum(y == 1 & z == 0)
  c01 <- sum(y == 0 & z == 1)
  c00 <- sum(y == 0 & z == 0)

  return((c11 * c00) / (c10 * c01))
}



############################# d and d^2 stats ###############################

d.stat <- function(y, z) {
  tmp <- xBalance(z ~ y,
           data = data.frame(y, z),
           report = "adj.mean.diffs")

  return(tmp$results[, "adj.diff",])
}

###################### Asymptotics ######################

setClass("AsymptoticTestStatistic",
  representation = c(asymptotic = "function"),
  contains = "function")


### Harmonic Mean Difference -- uses xBalance for asymptotics ###

.xBalanceBackEnd <- function(
  adjusted.data,
  treatment) {

  report <- c("adj.mean.diffs", "p.values")
  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.data, 2, function(y) {
    df <- data.frame(z = treatment, y = y)

    # ignoring blocks for now
    res <- xBalance(z ~ y,
                  data = df,
                  report = report)
    return(res$results[,,])
  })

  tmp <- as.data.frame(t(tmp))
  tcnms <- colnames(tmp)

  positions <- match(c("adj.diff", "p"), tcnms)
  results <- tmp[,positions]

  colnames(results) <- c("statistic", 'p.value')

  return(new("RandomizationDistribution",
      results,
      test.statistic = xBalance,
      z = as.numeric(treatment)))
}


# a little helper function
.h.weights <- function(n, m) {
  (m * (n-m)) / n
}

harmonic.mean.difference <- new("AsymptoticTestStatistic",
# the main implementation
function(ys, z) {
  z <- as.numeric(z)
  stopifnot(length(unique(z)) == 2) # require binary treatment for now

  f <- function(r,z) {
    mean(r[z==1]) - mean(r[z==0])
  }

  return(f(ys, z)) # just one block for now
  # use this when blocking is included again.
  # h.b <- tapply(z, blocks, function(z) { .h.weights(length(z), sum(z))})
  # d.b <- mapply(f, split(ys,blocks), split(z,blocks))

  # notice this is the same as "adj.diff" from xBalance
  # tests/test.testStatistic.R a demonstration
  # return((1/sum(h.b)) * sum(h.b * d.b))
  # end of stuff to turn on when blocking enabled

}, asymptotic = .xBalanceBackEnd)

### Mann-Whitney U, with wilcox.test as a backend ###

.wilcoxBackEnd <- function(
  adjusted.data,
  treatment) {


  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.data, 2, function(y) {
    res <- wilcox.test(y[!!treatment], y[!treatment])
    return(res[c("statistic", "p.value")])
  })

  tmp <- as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp) <- c("statistic", "p.value")

  return(new("RandomizationDistribution",
      tmp, # RD inherits from data.frame
      test.statistic = wilcox.test,
      z = as.numeric(treatment)))
}

mann.whitney.u <-  new("AsymptoticTestStatistic",
# the basic implementation
function(ys, z) {
  x <- ys[!!z]
  y <- ys[!z]

  # next stolen from wilcox.test
  r <- rank(c(x, y))
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))
  return(sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
}, asymptotic = .wilcoxBackEnd)

### Quantile Differences
### Similar to the KS test -- but evaluates differences at discrete points of
### the eCDF

quantileAbsoluteDifference <- function(quantiles) {
  function(y, z) {
    z1 <- ecdf(y[z == 1])
    z0 <- ecdf(y[z == 0])

    max(abs(quantile(z1, quantiles) - quantile(z0, quantiles)))
  }
}


iqrDiff <- function(ys,z,q1=.25,q2=.75,type=7){
  ## Inter-quartile/quantile difference
  ## a test statistic focusing on differences in scale
  qZtrted<-abs(diff(quantile(as.numeric(ys[!!z]), c(q1, q2), na.rm = FALSE, names = FALSE,
		    type = type)))
  qZctrl<-abs(diff(quantile(as.numeric(ys[!z]), c(q1, q2), na.rm = FALSE, names = FALSE,
		    type = type)))
  return(qZtrted-qZctrl)
}

madDiff <- function(ys,z){
  ## Median absolute deviation (which is scaled to be asymp. same as sd of a Normal)
  ## see help("mad")
 mad(ys[!!z]) - mad(ys[!z])
}


quantileDifference <- function(ys,z,q=.5){
  ## difference at a quantile
  quantile(ys[!!z],q)-quantile(ys[!z],q)
} 

### The ksTestStatistic with ks.test as a potential backend

# this borrowed from ks.test
.ksBackEnd <- function(
  adjusted.y,
  z) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.y, 2, function(y) {
    res <- stats::ks.test(y[!!z], y[!z]) # small problems may still be figured exactly
    return(res[c("statistic", "p.value")])
  })

  tmp <- as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp) <- c("statistic", "p.value")

  return(new("RandomizationDistribution",
      tmp, # RD inherits from data.frame
      test.statistic = ks.test,
      z = as.numeric(z)))
}

ksTestStatistic <- new("AsymptoticTestStatistic",
function(y, z) {
  # next borrowed from KS test
  # first, set up using the KS test var names
  x <- y[z == 1]
  y <- y[z == 0]

  x <- x[!is.na(x)]
  n <- length(x)

  y <- y[!is.na(y)]
  n.x <- as.double(n)
  n.y <- length(y)

  n <- n.x * n.y/(n.x + n.y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  if (length(unique(w)) < (n.x + n.y)) {
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
  }

  return(max(abs(z)))
}, asymptotic = .ksBackEnd)

ksTestStatistic.ranked <- new("AsymptoticTestStatistic",
function(y, z) {
  # next borrowed from KS test
  # first, set up using the KS test var names
  y<-rank(y)

  x <- y[z == 1]
  y <- y[z == 0]

  x <- x[!is.na(x)]
  n <- length(x)

  y <- y[!is.na(y)]
  n.x <- as.double(n)
  n.y <- length(y)

  n <- n.x * n.y/(n.x + n.y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  if (length(unique(w)) < (n.x + n.y)) {
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
  }

  return(max(abs(z)))
}, asymptotic = .ksBackEnd)


### The adTestStatistic with ad.test as a potential backend
.adBackEnd <- function(
  adjusted.y,
  z) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.y, 2, function(y) {
    res <- adk.test(y[!!z], y[!z]) # small problems may still be figured exactly
    return(res$adk[2,c(1,2)])
  })

  tmp <- as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp) <- c("statistic", "p.value")

  return(new("RandomizationDistribution",
      tmp, # RD inherits from data.frame
      test.statistic = ad.test,
      z = as.numeric(z)))
}

adTestStatistic <- new("AsymptoticTestStatistic",
		       function(y, z) {
			 require(adk)
                         x <- y[z == 1]
                         y <- y[z == 0]

                         x <- x[!is.na(x)]

                         y <- y[!is.na(y)]

                         return(adk.test(x,y)$adk[2,1])
                       }, asymptotic = .adBackEnd)

adTestStatistic.ranked <- new("AsymptoticTestStatistic",
		       function(y, z) {
                         require(adk)

                         y<-rank(y)

                         x <- y[z == 1]
                         y <- y[z == 0]

                         x <- x[!is.na(x)]

                         y <- y[!is.na(y)]

                         return(adk.test(x,y)$adk[2,1])
                       }, asymptotic = .adBackEnd)
### The ssrTestStatistic with the F-test as a potential asymp backend
### (and using F as the test statistic)

.ssrBackEnd <- function( adjusted.y, z, d) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.y, 2, function(y) {
    tmp <- summary(lm(y~z+d))$fstatistic
    res <- c("statistic"=tmp["value"],
             "p.value"=1-pf(tmp["value"],tmp["numdf"],tmp["dendf"]))
    return(res[c("statistic", "p.value")])
  })

  tmp <- as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp) <- c("statistic", "p.value")

  return(new("RandomizationDistribution",
      tmp, # RD inherits from data.frame
      test.statistic = ks.test,
      z = as.numeric(z)))
}

ssrTestStatistic <- new("AsymptoticTestStatistic",
                        function(y, z, d) {
                        ##  return(summary(lm(y~z+d))$sigma)
                          return(summary(lm(y~z+d))$fstatistic["value"])
                        }, asymptotic = .ssrBackEnd)



### The Cramer-von Mises test statistic and asymp backend (using W2 version)
### The ksTestStatistic with ks.test as a potential backend

# this borrowed from cvm.test
.cvmBackEnd <- function(
  adjusted.y,
  z) {

  require(dgof)
  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.y, 2, function(y) {
    res <- dgof::cvm.test(y[!!z], y[!z]) # small problems may still be figured exactly
    return(res[c("statistic", "p.value")])
  })

  tmp <- as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp) <- c("statistic", "p.value")

  return(new("RandomizationDistribution",
      tmp, # RD inherits from data.frame
      test.statistic = cvm.test,
      z = as.numeric(z)))
}

cvmTestStatistic <- new("AsymptoticTestStatistic",
                        function(y, z) {
			  ## next borrowed from  cvm.stat.disc
			  ## first, set up using the cvm.stat.disc test var names
			  x <- y[z == 1]
			  y <- y[z == 0]

			  x <- x[!is.na(x)]
			  y <- y[!is.na(y)]

			  I <- knots(y)
			  N <- length(x)
			  e <- diff(c(0, N * y(I)))
			  obs <- rep(0, length(I))
			  for (j in 1:length(I)) {
			    obs[j] <- length(which(x == I[j]))
			  }
			  S <- cumsum(obs)
			  T <- cumsum(e)
			  H <- T/N
			  p <- e/N
			  t <- (p + p[c(2:length(p), 1)])/2
			  Z <- S - T
			  Zbar <- sum(Z * t)
			  S0 <- diag(p) - p %*% t(p)
			  A <- matrix(1, length(p), length(p))
			  A <- apply(row(A) >= col(A), 2, as.numeric)
                          E <- diag(t)
                          One <- rep(1, nrow(E))
                          K <- diag(0, length(H))
                          diag(K)[-length(H)] <- 1/(H[-length(H)] * (1 - H[-length(H)]))
                          Sy <- A %*% S0 %*% t(A)
                          M <- E
                          ##switch(type, W2 = E, U2 = (diag(1, nrow(E)) - E %*%
                          ##                    One %*% t(One)) %*% E %*% (diag(1, nrow(E)) - One %*%
                          ##                                               t(One) %*% E), A2 = E %*% K)
                          lambda <- eigen(M %*% Sy)$values
                          STAT <- sum(Z^2 * t)/N
                            ## switch(type, W2 = sum(Z^2 * t)/N,
                            ## U2 = sum((Z -  Zbar)^2 * t)/N, A2 = sum((Z^2 * t/(H * (1 - H)))[-length(I)])/N)
                          return(STAT)

                        }, asymptotic = .cvmBackEnd
                        )


### Subset/subgroup Analysis
### using a logical vector, limit the test statistic to a smaller group

subsetStatistic <- function(statistic, ...) {
  function(y, z) {
    statistic(subset(y, ...), subset(z, ...))
  }
}

