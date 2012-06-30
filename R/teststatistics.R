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

mean.difference <- function(ys, z) { 
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z) 
  return(mean(ys[z]) - mean(ys[!z]))
}

mean.diff.lsfit<-function(ys,z){ ##Try using something that calls compiled code
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

### The ksTestStatistic with ks.test as a potential backend

# this borrowed from ks.test
.ksBackEnd <- function(
  adjusted.y,
  z) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.y, 2, function(y) {
    res <- ks.test(y[!!z], y[!z]) # small problems may still be figured exactly
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


### Subset/subgroup Analysis
### using a logical vector, limit the test statistic to a smaller group

subsetStatistic <- function(statistic, ...) {
  function(y, z) {
    statistic(subset(y, ...), subset(z, ...))
  }
}

