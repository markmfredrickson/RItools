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
##' Mean functions
##'
##' Various functions to compute means
##' @param ys ys
##' @param z Vector that can be made logical.
##' @return mean
##' @name test.stat.means
fastmean <- function(ys){
  ## Skip checks and dispatch etc.
  .Internal(mean(ys))
}

##' @rdname test.stat.means
mean.difference  <-  function(ys, z) {
  # z is usually a vector of 1s and 0s. make it logical
  z  <-  as.logical(z)
  return(fastmean(ys[z]) - fastmean(ys[!z]))
}

##' @rdname test.stat.means
mean.diff.lsfit <- function(ys,z){
  ##Try using something that calls compiled code
  ##Gives same answer as mean.difference for balanced blocks and should be like harmonic.mean.difference for unbalanced blocks.
  lsfit(x=model.matrix(ys~z),y=ys,intercept=FALSE)[["coefficients"]][["z"]]
}

##' @rdname test.stat.means
mean.diff.vect <- function(ys,z){
  X <- model.matrix(ys~z)
  solve(qr(X, LAPACK=TRUE), ys)[2] ## qr.coef(qr(X,LAPACK=TRUE),ys) ## to handle near singular X
}

##' @rdname test.stat.means
t.mean.difference  <-  function(ys, z) {
  # z is usually a vector of 1s and 0s. make it logical
  z  <-  as.logical(z)
  return( {fastmean(ys[z]) - fastmean(ys[!z])} / sqrt( {fastvar(ys[z]) / sum(z)} + {fastvar(ys[!z]) / sum(!z)}) )
}

##' @rdname test.stat.means
fastvar <- function(ys){
  ## See var for this. Right now it removes missing values
  .Call(stats:::C_cov,ys,NULL,5,FALSE)
}

############################## Rank Based Functions ##############################

##' Rank based functions
##'
##' @param ys ys
##' @param z Vector that can be made logical
##' @return sum
##' @name rank.sum.functions
rank.sum  <-  function(ys, z) {
  z  <-  as.logical(z)
  ys.ranks  <-  rank(ys)
  return(sum(ys.ranks[z]))
}

##' @rdname rank.sum.functions
paired.sgnrank.sum  <-  function(ys, z){
  stopifnot(length(unique(z)) == 2) ##require binary treatment for now

  # we assume data are in block order, eg. 11223344...
  blocks  <-  gl(length(ys) / 2, 2)

  Y <- sapply(split(data.frame(r=ys,z=z),blocks),function(dat){with(dat,r[z == 1] - r[z == 0])})
  sgn <- as.numeric(Y > 0)
  q <- rank(abs(Y))
  sum(sgn * q)
}

# mean.rank.diff.strata <- function(ys, z, blocks){
#   ##average difference in ranks between control and treated obs
#   ##currently weighted equally across strata (i.e. good for pairs not good for other designs)
#   mean(
#        mapply(FUN = function(ys, z) {
#          R <- rank(ys)
#          mean(R * z)-mean(R * (1-z))
#          ##((R * z) / sum(z))-((R * (1-z)) / sum(1-z))
#        },
#               ys = split(ys, blocks),
#               z = split(z, blocks))
#        )
# }


############################## Misc Statistics ##############################

##' Odds ratio
##'
##' @param y y
##' @param z z
##' @return odds ratio
odds.ratio  <-  function(y, z) {
  c11  <-  sum(y == 1 & z == 1)
  c10  <-  sum(y == 1 & z == 0)
  c01  <-  sum(y == 0 & z == 1)
  c00  <-  sum(y == 0 & z == 0)

  return((c11 * c00) / (c10 * c01))
}



############################# d and d^2 stats ###############################

##' d stat
##'
##' @param y y
##' @param z z
##' @return d
d.stat  <-  function(y, z) {
  tmp  <-  xBalance(z ~ y,
                  data = data.frame(y, z),
                  report = "adj.mean.diffs")

  return(tmp$results[, "adj.diff",])
}

###################### Asymptotics ######################

##' S4 AsymptoticTestStatistic class
##'
##' @slot asymptotic function
##' @import kSamples
setClass("AsymptoticTestStatistic",
         representation = c(asymptotic = "function"),
         contains = "function")


### Harmonic Mean Difference -- uses xBalance for asymptotics ###

##' xBalanceBackEnd
##'
##' @param adjusted.data adjusted.data
##' @param treatment treatment
##' @return RandomizationDistribution
##' @name xBalanceBackEnd
.xBalanceBackEnd  <-  function(
                             adjusted.data,
                             treatment) {

  report  <-  c("adj.mean.diffs", "p.values")
  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp  <-  apply(adjusted.data, 2, function(y) {
               df  <-  data.frame(z = treatment, y = y)

               # ignoring blocks for now
               res  <-  xBalance(z ~ y,
                               data = df,
                               report = report)
               return(res$results[,,])
                             })

  tmp  <-  as.data.frame(t(tmp))
  tcnms  <-  colnames(tmp)

  positions  <-  match(c("adj.diff", "p"), tcnms)
  results  <-  tmp[,positions]

  colnames(results)  <-  c("statistic", 'p.value')

  return(new("RandomizationDistribution",
             results,
             test.statistic = xBalance))
}


##' Helper function
##'
##' @param n n
##' @param m m
##' @return value
##' @name h.weights
.h.weights  <-  function(n, m) {
  (m * (n-m)) / n
}


harmonic.mean.difference  <-  new("AsymptoticTestStatistic",
                                # the main implementation
                                function(ys, z) {
                                  z  <-  as.numeric(z)
                                  stopifnot(length(unique(z)) == 2) # require binary treatment for now

                                  f  <-  function(r,z) {
                                    fastmean(r[z == 1]) - fastmean(r[z == 0])
                                  }

                                  return(f(ys, z)) # just one block for now
                                  # use this when blocking is included again.
                                  # h.b  <-  tapply(z, blocks, function(z) { .h.weights(length(z), sum(z))})
                                  # d.b  <-  mapply(f, split(ys,blocks), split(z,blocks))

                                  # notice this is the same as "adj.diff" from xBalance
                                  # tests / test.testStatistic.R a demonstration
                                  # return((1 / sum(h.b)) * sum(h.b * d.b))
                                  # end of stuff to turn on when blocking enabled

                                }, asymptotic = .xBalanceBackEnd)

### Mann-Whitney U, with wilcox.test as a backend ###

##' Wilcox BackEnd
##'
##' @param adjusted.data adjusted.data
##' @param treatment treatment
##' @return RandomizationDistribution
##' @name wilcoxBackEnd
.wilcoxBackEnd  <-  function(
                           adjusted.data,
                           treatment) {


  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp  <-  apply(adjusted.data, 2, function(y) {
               res  <-  wilcox.test(y[!!treatment], y[!treatment])
               return(res[c("statistic", "p.value")])
                           })

  tmp  <-  as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp)  <-  c("statistic", "p.value")

  return(new("RandomizationDistribution",
             tmp, # RD inherits from data.frame
             test.statistic = wilcox.test))
}

mann.whitney.u  <-   new("AsymptoticTestStatistic",
                       # the basic implementation
                       function(ys, z) {
                         x  <-  ys[!!z]
                         y  <-  ys[!z]

                         # next stolen from wilcox.test
                         r  <-  rank(c(x, y))
                         n.x  <-  as.double(length(x))
                         n.y  <-  as.double(length(y))
                         return(sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2)
                       }, asymptotic = .wilcoxBackEnd)

###
### Similar to the KS test -- but evaluates differences at discrete points of
### the eCDF

##' Quantile Differences
##'
##' Similar to the KS test -- but evaluates differences at discrete
##' points of the eCDF
##' @param quantiles quantiles
##' @return value
quantileAbsoluteDifference  <-  function(quantiles) {
  function(y, z) {
    z1  <-  ecdf(y[z == 1])
    z0  <-  ecdf(y[z == 0])

    max(abs(quantile(z1, quantiles) - quantile(z0, quantiles)))
  }
}

##' IQR Difference
##'
##' @param q1 Lower quartile
##' @param q2 Upper quartile
##' @param type type
##' @return IQR
iqrDiff  <-  function(q1=.25,q2=.75,type=7){
  function(ys,z){
    ## Inter-quartile / quantile difference
    ## a test statistic focusing on differences in scale
    qZtrted <- abs(diff(quantile(as.numeric(ys[!!z]), c(q1, q2), na.rm = FALSE, names = FALSE,
                               type = type)))
    qZctrl <- abs(diff(quantile(as.numeric(ys[!z]), c(q1, q2), na.rm = FALSE, names = FALSE,
                              type = type)))
    return(qZtrted-qZctrl)
  }
}

##' MAD
##'
##' @param ys ys
##' @param z z
##' @return diff
madDiff  <-  function(ys,z){
  ## Median absolute deviation (which is scaled to be asymp. same as sd of a Normal)
  ## see help("mad")
  mad(ys[!!z]) - mad(ys[!z])
}

##' Quantile difference
##'
##' @param q quantiles
##' @return diff
quantileDifference  <-  function(q=.5){
  function(ys,z){
    ## difference at a quantile
    quantile(ys[!!z],q)-quantile(ys[!z],q)
  }
}

### The ksTestStatistic with ks.test as a potential backend

## KS background
##'
##' this borrowed from ks.test
##' @param adjusted.y y
##' @param z z
##' @return RandomizationDistribution
##' @name ksBackEnd
.ksBackEnd  <-  function( adjusted.y, z) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp  <-  apply(adjusted.y, 2, function(y) {
               res  <-  stats::ks.test(y[!!z], y[!z]) # small problems may still be figured exactly
               return(res[c("statistic", "p.value")])
                       })

  tmp  <-  as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp)  <-  c("statistic", "p.value")

  return(new("RandomizationDistribution",
             tmp, # RD inherits from data.frame
             test.statistic = ks.test))
}

ksTestStatistic  <-  new("AsymptoticTestStatistic",
                       function(y, z) {
                         # next borrowed from KS test
                         # first, set up using the KS test var names
                         x  <-  y[z == 1]
                         y  <-  y[z == 0]

                         x  <-  x[!is.na(x)]
                         n  <-  length(x)

                         y  <-  y[!is.na(y)]
                         n.x  <-  as.double(n)
                         n.y  <-  length(y)

                         n  <-  n.x * n.y / (n.x + n.y)
                         w  <-  c(x, y)
                         z  <-  cumsum(ifelse(order(w) <= n.x, 1 / n.x, -1 / n.y))
                         if (length(unique(w)) < (n.x + n.y)) {
                           z  <-  z[c(which(diff(sort(w)) != 0), n.x + n.y)]
                         }

                         return(max(abs(z)))
                       }, asymptotic = .ksBackEnd)

##' The adTestStatistic with ad.test as a potential backend
##'
##' @param adjusted.y y
##' @param z z
##' @return RandomizationDistribution
##' @name adBackEnd
.adBackEnd  <-  function( adjusted.y, z) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp  <-  apply(adjusted.y, 2, function(y) {
               res  <-  ad.test(y[!!z], y[!z])
               return(res$ad[1,c(1,3)])
                              })

  tmp  <-  as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp)  <-  c("statistic", "p.value")

  return(new("RandomizationDistribution",
             tmp, # RD inherits from data.frame
             test.statistic = ad.test))
}

adTestStatistic  <-  new("AsymptoticTestStatistic",
                       function(y, z) {
                         x  <-  y[z == 1]
                         y  <-  y[z == 0]

                         x  <-  x[!is.na(x)]

                         y  <-  y[!is.na(y)]

                         return(ad.test(x,y)$ad[1,1])
                       }, asymptotic = .adBackEnd)

##' ssrTestStatistic
##'
##' The ssrTestStatistic with the F-test as a potential asymp backend
##' (and using F as the test statistic)
##' @param S S
##' @return AsymptoticTestStatistic
ssrTSmaker <- function(S){
  force(S)
  return(new("AsymptoticTestStatistic",
             function(y,z){
               d <- as.vector(z %*% S)
               return(summary(lm(y~z+d))$fstatistic[["value"]])
             },
             function(adjusted.y,z){
               d <- as.vector(z %*% S)
               # adjusted.data should be a matrix, where each column is an adjusted data
               tmp  <-  apply(adjusted.y, 2, function(y) {
                            tmp  <-  summary(lm(y~z+d))$fstatistic
                            res  <-  c("statistic"=tmp[["value"]],
                                     "p.value"=1-pf(tmp[["value"]],tmp[["numdf"]],tmp[["dendf"]]))
                            return(res[c("statistic", "p.value")])
                              })

               tmp  <-  as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
               colnames(tmp)  <-  c("statistic", "p.value")

               return(new("RandomizationDistribution",
                          tmp, # RD inherits from data.frame
                          test.statistic = ssrTestStatistic))
             }
             ))
}

##' ssr backEnd
##'
##' @param adjusted.y y
##' @param z z
##' @param d d
##' @return RandomizationDistribution
##' @name ssrBackEnd
.ssrBackEnd  <-  function( adjusted.y, z, d) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp  <-  apply(adjusted.y, 2, function(y) {
               tmp  <-  summary(lm(y~z+d))$fstatistic
               res  <-  c("statistic"=tmp["value"],
                        "p.value"=1-pf(tmp["value"],tmp["numdf"],tmp["dendf"]))
               return(res[c("statistic", "p.value")])
                              })

  tmp  <-  as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp)  <-  c("statistic", "p.value")

  return(new("RandomizationDistribution",
             tmp, # RD inherits from data.frame
             test.statistic = ssrTestStatistic))
}

ssrTestStatistic  <-  new("AsymptoticTestStatistic",
                        function(y, z, d) {
                          ##  return(summary(lm(y~z+d))$sigma)
                          return(summary(lm(y~z+d))$fstatistic["value"])
                        }, asymptotic = .ssrBackEnd)



### The Cramer-von Mises test statistic and asymp backend (using W2 version)
### The ksTestStatistic with ks.test as a potential backend

#
##' cvm BackEnd
##'
##' this borrowed from cvm.test
##' @param adjusted.y y
##' @param z z
##' @return RandomizationDistribution
##' @import dgof
##' @name cvmBackEnd
.cvmBackEnd  <-  function(
                        adjusted.y,
                        z) {

  # adjusted.data should be a matrix, where each column is an adjusted data
  tmp  <-  apply(adjusted.y, 2, function(y) {
               res  <-  dgof::cvm.test(y[!!z], y[!z]) # small problems may still be figured exactly
               return(res[c("statistic", "p.value")])
                        })

  tmp  <-  as.data.frame(matrix(unlist(tmp), ncol = 2, byrow = T))
  colnames(tmp)  <-  c("statistic", "p.value")

  return(new("RandomizationDistribution",
             tmp, # RD inherits from data.frame
             test.statistic = cvm.test))
}

cvmTestStatistic  <-  new("AsymptoticTestStatistic",
                        function(y, z) {
                          ## next borrowed from  cvm.stat.disc
                          ## first, set up using the cvm.stat.disc test var names
                          x  <-  y[z == 1]
                          y  <-  y[z == 0]

                          x  <-  x[!is.na(x)]
                          y  <-  y[!is.na(y)]

                          I  <-  knots(y)
                          N  <-  length(x)
                          e  <-  diff(c(0, N * y(I)))
                          obs  <-  rep(0, length(I))
                          for (j in 1:length(I)) {
                            obs[j]  <-  length(which(x == I[j]))
                          }
                          S  <-  cumsum(obs)
                          T  <-  cumsum(e)
                          H  <-  T / N
                          p  <-  e / N
                          t  <-  (p + p[c(2:length(p), 1)]) / 2
                          Z  <-  S - T
                          Zbar  <-  sum(Z * t)
                          S0  <-  diag(p) - p %*% t(p)
                          A  <-  matrix(1, length(p), length(p))
                          A  <-  apply(row(A) >= col(A), 2, as.numeric)
                          E  <-  diag(t)
                          One  <-  rep(1, nrow(E))
                          K  <-  diag(0, length(H))
                          diag(K)[-length(H)]  <-  1 / (H[-length(H)] * (1 - H[-length(H)]))
                          Sy  <-  A %*% S0 %*% t(A)
                          M  <-  E
                          ##switch(type, W2 = E, U2 = (diag(1, nrow(E)) - E %*%
                          ##                    One %*% t(One)) %*% E %*% (diag(1, nrow(E)) - One %*%
                          ##                                               t(One) %*% E), A2 = E %*% K)
                          lambda  <-  eigen(M %*% Sy)$values
                          STAT  <-  sum(Z^2 * t) / N
                          ## switch(type, W2 = sum(Z^2 * t) / N,
                          ## U2 = sum((Z -  Zbar)^2 * t) / N, A2 = sum((Z^2 * t / (H * (1 - H)))[-length(I)]) / N)
                          return(STAT)

                        }, asymptotic = .cvmBackEnd
                        )


## Neymans Smooth Test

##' Quasi Relative Distance
##'
##' @param y y
##' @param yo y0
##' @param binn binn
##' @param location location
##' @return vector
quasireldist <- function(y,yo,binn=100,location="median"){
  ### This version follows the methods in 9.1.2.3 and 9.2 of Handcock and Morris, Relative Dist.
  ## first part is rmdata (from reldist.R)
  nas  <-  is.na(y)
  m  <-  length(y)
  ywgt <- rep(1,length=m) / m
  y  <-  y[!nas]

  nas  <-  is.na(yo)
  n  <-  length(yo)
  yowgt <- rep(1,length=n) / n
  yo  <-  yo[!nas]
  n  <-  length(yo)
  m  <-  length(y)

  #  sy is 1:n
  #  ys is the ordered y's
  #  ry is the ranks of ys in the joint vector of {yo, y}
  #  ry - sy is the number of yo's le ys
  #
  sy  <-  seq(along = y)
  ys  <-  sort.list(y)
  ywgt  <-  ywgt[ys]
  ys  <-  y[ys] ## ordered y's
  # ry  <-  sort.list(sort.list(c(yo, ys)))[sy + n]
  ry  <-  sort.list(sort.list(c(yo, ys)))
  ry  <-  ry[sy + n] ## take the last n of the ry's
  yos  <-  sort.list(yo) ## order of yo's
  yowgts  <-  yowgt[yos]
  #
  #   x is the new relative data, but note it has weight ywgt
  #
  x  <-  c(0,cumsum(yowgts))[ry - sy + 1]
  x[x > 1]  <-  1
  x[x < 0]  <-  0
  #  Use x[order(order(y))]
  #  to get the relative data back in original order
  #
  #  returned values are in the original order
  #
  rorder  <-  order(order(y)) ## integer rank of y basically rank(y) but no partials
  ##THis is what is returned from rmdata to rcdist
  ##invisible(list(x=x[rorder],wgt=ywgt[rorder],y=y,yo=yo,ywgt=ywgt[rorder],yowgt=yowgt,n=n,m=m))

  x <- x[rorder]
  wgt <- ywgt[rorder]
  ywgt <- ywgt[rorder]

  ## Next part is from rcdist
  r  <-  seq(0, 1, length = binn + 1)[-1] - 0.5 / binn

  ywgt  <-  ywgt / sum(ywgt)
  ## Break x into categories
  xxx  <-  cut(x, breaks = seq(0, 1, length = binn + 1),
             include.lowest=TRUE,labels=FALSE)

  ## Sum up density of y within bins
  xx  <-  rep(0,length=binn)
  for (i in unique(xxx)){
    xx[i]  <-  sum(ywgt[xxx == i],na.rm=TRUE)
  }
  return(as.vector(xx))
}

##' ddst test Statistic
##'
##' @param binn binn
##' @param location location
##' @return function
##' @import ddst
ddstTestStatistic <- function(binn=100,location="median"){
  function(ys,z){
    ## Neyman's Smooth test with data driven tuning parameter choice from Ledwina
    qdat <- quasireldist(y=ys[!!z],yo=ys[!z])
    ## theddst <- ddst.uniform.test(qdat)
    ## return(theddst$statistic)
    ## this next from ddst.uniform.test.R
    coord <- ddst:::ddst.uniform.Nk(qdat,ddst:::ddst.base.legendre,10)
    l <- ddst:::ddst.IIC(coord,length(qdat),2.4)
    return(coord[l])
  }
}

##' entropy relative distance
##'
##' @param smooth smooth
##' @param binn binn
##' @param location location
##' @return function
##' @import mgcv
entropyRelDistTestStatistic <- function(smooth=.01,binn=100,location="median"){
  function(ys,z){
  ## from Reldist rcdist
    ## n  <-  length(ys[!z])
    m  <-  length(ys[!!z])
    qdat <- quasireldist(y=ys[!!z],yo=ys[!z])
    r  <-  seq(0, 1, length = binn + 1)[-1] - 0.5 / binn
    is.wholenumber  <-  function(x, tol = .Machine$double.eps^0.5){
      abs(x - round(x)) < tol
    }
    family  <-  ifelse(is.wholenumber(m * qdat),"poisson","quasi")
    yl  <-  mgcv::gam(y ~ s(x, bs = "cr"), sp = smooth, family = family,
              data = data.frame(x = r, y = m * qdat))
    gpdf  <-  predict(yl, type = "response")
    scalef  <-  binn / sum(gpdf)
    gpdf  <-  gpdf * scalef

    egv  <-  gpdf[gpdf > 0]
    ##entropy  <-  sum(egv * log(egv)) / binn

    return(sum(egv * log(egv)) / binn)
  }
}








## Entropy (From reldist)

##' Subset / subgroup Analysis
##'
##' using a logical vector, limit the test statistic to a smaller group
##' @param statistic statistic
##' @param ... Add'l arguments to subset
##' @return function
subsetStatistic  <-  function(statistic, ...) {
  function(y, z) {
    statistic(subset(y, ...), subset(z, ...))
  }
}
