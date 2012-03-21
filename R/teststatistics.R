################################################################################
# Test Statistic Functions for Randomization Distributions
#
# A test statistic is a function of three arguments:
# - y: the data, after adjustment
# - z: a treatment indicator (should be numeric, but calling as.numeric(z) is safe)
# - b: an indicator for the blocking factor, may contain only one level
#
# The function should return a scalar real value
#
# A test statistic may also be a AsymptoticTestStatistic, a class that descends from
# a function. This class includes a `asymptotic` slot which is a backend function for
# the randomizationDistributionEngine. For example, the harmonic.mean.difference test
# statistic can use xBalance for the asymptotic results.
################################################################################

############################### Mean Differences ###############################

mean.difference <- function(ys, z, blocks) { 
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z) 

  tmeans <- aggregate(ys[z], by = list(blocks[z, drop=F]), mean)$x 
  cmeans <- aggregate(ys[!z], by = list(blocks[!z, drop=F]), mean)$x
  tcdiff <- tmeans - cmeans
  block.weights <- table(blocks) / length(blocks)

  return(sum(block.weights * tcdiff))
}

mean.diff.lsfit<-function(ys,z,blocks){ ##Try using something that calls compiled code
  ##Gives same answer as mean.difference for balanced blocks and should be like harmonic.mean.difference for unbalanced blocks.
  lsfit(x=model.matrix(ys~z+blocks),y=ys,intercept=FALSE)[["coefficients"]][["z"]]
}

mean.diff.vect<-function(ys,z,blocks){
 X<-model.matrix(ys~z+blocks)
 solve(qr(X, LAPACK=TRUE), ys)[2] ## qr.coef(qr(X,LAPACK=TRUE),ys) ## to handle near singular X
}


mean.diff.noblocks<-function(ys, z, blocks) { 
  z <- as.numeric(z)
  stopifnot(all(sort(unique(z))==c(0,1))) # require binary treatment for now
  return((sum(ys*z)/sum(z))-(sum(ys*(1-z))/sum(1-z))) # I supect that arithmetic is faster than indexing.
}
############################## Rank Based Functions ##############################

rank.sum <- function(ys, z, blocks) {
  # ignore blocks for now, could probably use the "aggregate" function
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z)  
  ys.ranks <- rank(ys)
  return(sum(ys.ranks[z]))
}

paired.sgnrank.sum<-function(ys,z,blocks){ stopifnot(length(unique(z))==2) ##require binary treatment for now
  Y<-sapply(split(data.frame(r=ys,z=z),blocks),function(dat){with(dat,r[z==1]-r[z==0])})
  sgn<-as.numeric(Y>0) 
  q<-rank(abs(Y))
  sum(sgn*q)
}

mean.rank.diff.strata<-function(ys, z, blocks){
  ##average difference in ranks between control and treated obs
  ##currently weighted equally across strata (i.e. good for pairs not good for other designs)
  mean(
       mapply(FUN = function(ys, z) {
         R<-rank(ys)
         mean(R*z)-mean(R*(1-z))
         ##((R*z)/sum(z))-((R*(1-z))/sum(1-z))
       },
              ys = split(ys, blocks),
              z = split(z, blocks))
       )
}


############################## Misc Statistics ##############################

odds.ratio <- function(y, z, b) {
  c11 <- sum(y == 1 & z == 1)
  c10 <- sum(y == 1 & z == 0)
  c01 <- sum(y == 0 & z == 1)
  c00 <- sum(y == 0 & z == 0)
  
  return((c11 * c00) / (c10 * c01))
}     



############################# d and d^2 stats ###############################

d.stat <- function(y, z, b) {
  tmp <- xBalance(z ~ y,
           data = data.frame(y, z, b),
           strata = list(b = ~b),
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
  treatment, 
  blocks, 
  samples, p.value, # these will be ignored
  summaries, ...) {

  if (length(levels(blocks)) > 1) {
    strata <- blocks  
  } else {
    strata <- NULL  
  }

  report <- unlist(c("adj.mean.diffs", "p.values", summaries))
  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp <- apply(adjusted.data, 2, function(y) {
    df <- data.frame(z = treatment, y = y)

    # ignoring blocks for now
    res <- xBalance(z ~ y, 
                  data = df,
                  report = report,
                  strata = blocks)
    return(res$results[,,])
  })

  tmp <- as.data.frame(t(tmp))
  tcnms <- colnames(tmp)

  positions <- match(c("adj.diff", "p"), tcnms)
  important <- tmp[,positions]
  requested <- tmp[, -positions, drop = F]
  
  colnames(important) <- c("statistic", 'p.value')
  colnames(requested) <- tcnms[-positions]

  return(new("RandomizationDistribution", 
      cbind(important, requested), # RD inherits from data.frame
      test.statistic = xBalance,
      z = as.numeric(treatment),
      blocks = as.numeric(blocks)))
}


# a little helper function
.h.weights <- function(n, m) {
  (m * (n-m)) / n
}

harmonic.mean.difference <- new("AsymptoticTestStatistic",
# the main implementation
function(ys, z, blocks) {
  z <- as.numeric(z)
  stopifnot(length(unique(z)) == 2) # require binary treatment for now

  f <- function(r,z) {
    mean(r[z==1]) - mean(r[z==0])
  }

  h.b <- tapply(z, blocks, function(z) { .h.weights(length(z), sum(z))})
  d.b <- mapply(f, split(ys,blocks), split(z,blocks))  
  
  # notice this is the same as "adj.diff" from xBalance
  # tests/test.testStatistic.R a demonstration
  return((1/sum(h.b)) * sum(h.b * d.b)) 
}, asymptotic = .xBalanceBackEnd)

### Mann-Whitney U, with wilcox.test as a backend ###

.wilcoxBackEnd <- function(
  adjusted.data,
  treatment, 
  blocks, 
  samples, p.value, summaries, # these will be ignored
  ...) {
  
  # TODO check blocks so see if they are paired
  # if (length(levels(blocks)) > 1) {
  #   strata <- blocks  
  # } else {
  #   strata <- NULL  
  # }
  
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
      z = as.numeric(treatment),
      blocks = as.numeric(blocks)))
}

mann.whitney.u <-  new("AsymptoticTestStatistic",
# the basic implementation
function(ys, z, blocks) {
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
  function(y, z, b) {
    z1 <- ecdf(y[z == 1])
    z0 <- ecdf(y[z == 0])

    max(abs(quantile(z1, quantiles) - quantile(z0, quantiles)))    
  }
}


### Subset/subgroup Analysis
### using a logical vector, limit the test statistic to a smaller group

subsetStatistic <- function(statistic, ...) {
  function(y, z, b) {
    statistic(subset(y, ...), subset(z, ...), b)
  }
}
