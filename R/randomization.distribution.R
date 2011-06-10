# require(gmp)

# an arbitrarily exact n choose k algorithm
# this could be replace by an approximation, as we only need to 
# know if nCk is larger than samples, where samples is probably@#
# a normal int or a long.
# MF: are you aware that R provides "lchoose()", which unlike choose()
# doesn't get choked up on big integers?  You might use lchoose()
# and compare to log(samples). But maybe you knew that and have 
# a reason not to. -BH
bigchoose <- function(n, k) {
  if (n < 1 || k < 1 || k > n) {
    return(as.bigz(0))  
  }

  if (n == k) {
    return(as.bigz(1))  
  }

  if( k > (n / 2)) {
    k <- n - k;
  }

  numer <- as.bigz(1)
  for (i in n:(n - k + 1)) {
    numer <- numer * i
  }

  denom <- as.bigz(1)
  for (i in 1:k) {
    denom <- denom * i  
  }

  return(numer / denom)
}
stopifnot(choose(10, 3) == bigchoose(10, 3))
stopifnot(bigchoose(3, 0) == 0)

###################################################
### Test statistics and helper functions
###################################################

mean.difference <- function(ys, z, blocks) { ##This does not use the harmonic mean. Perhaps include such a function
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z);stopifnot(length(unique(z))==2) ##require binary treatment for now
  tmeans <- aggregate(ys[z], by = list(blocks[z, drop=F]), mean)$x 
  cmeans <- aggregate(ys[!z], by = list(blocks[!z, drop=F]), mean)$x
  tcdiff <- tmeans - cmeans
  block.weights <- table(blocks) / length(blocks)

  return(sum(block.weights * tcdiff))
}

mean.difference.2 <- function(ys, z, blocks) { ##This does not use the harmonic mean. Perhaps include such a function
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z);stopifnot(length(unique(z))==2) ##require binary treatment for now
  tmeans <- aggregate(ys[z], by = list(blocks[z, drop=F]), mean)$x 
  cmeans <- aggregate(ys[!z], by = list(blocks[!z, drop=F]), mean)$x
  tcdiff <- tmeans - cmeans
  block.weights <- table(blocks) / length(blocks)

  return(sum(block.weights * tcdiff))
}


adj.mean.diff<-function(ys,z,blocks){
  if(is.logical(z)){z<-as.numeric(z)}
  stopifnot(length(unique(z))==2) ##require binary treatment for now
  h.fn<-function(n,m){(m*(n-m))/n} ##harmonic weighting function
  
  d.b.fn<-function(ys,z,blocks){ ##blockwise mean diffs function is this better or worse that above?
    mapply(function(r,z){mean(ys[z==1])-mean(ys[z==0])},split(ys,blocks),split(z,blocks))
  } 

  h.b<-tapply(z,blocks,function(z){h.fn(n=length(z),m=sum(z))})
  d.b<-d.b.fn(ys=ys,z=z,blocks=blocks)
  return((1/sum(h.b))*sum(h.b*d.b)) ##notice this is the same as "adj.diff" from xBalance, that should be the test
}

rank.sum <- function(ys, z, blocks) {
  # ignore blocks for now, could probably use the "aggregate" function
  # z is usually a vector of 1s and 0s. make it logical
  z <- as.logical(z);  stopifnot(length(unique(z))==2) ##require binary treatment for now
  ys.ranks <- rank(ys)
  return(sum(ys.ranks[z]))

}


paired.sgnrank.sum<-function(ys,z,blocks){ stopifnot(length(unique(z))==2) ##require binary treatment for now
  Y<-sapply(split(data.frame(r=ys,z=z),blocks),function(dat){with(dat,r[z==1]-r[z==0])})
  sgn<-as.numeric(Y>0) 
  q<-rank(abs(Y))
  sum(sgn*q)
}



get.cis<-function(distobj,thelevels,p.value.function){ ##get several confidence intervals from a distribution object
  thecis<-sapply(thelevels,function(l){
    confint(distobj,p.value.function=p.value.function,level=l)$ci
  })
  colnames(thecis)<-thelevels
  return(thecis)
}


###################################################
### Modes of effect/hypotheses, residuals
###################################################

# a helper function to create constant additive models of effect
constant.additive.hypothesis.factory <- function(hypothesized.value) {
  function(ys, z) { ys + (z * hypothesized.value) } ##This ia a problem. All of the math implies subtraction rather than addition. 
}

# a helper to create a list of constant hypos.
constant.hypotheses <- function(range,factory=constant.additive.hypothesis.factory) {
  n <- length(range)
  functions <- vector("list", n)
  for (i in 1:n) {
    functions[[i]] <- factory(range[i])

    # this is a bizarre bug: if we don't call the function on some
    # data, the list will overwrite each i with the last value of i
    # e.g. const.hyp(1:5) would equivalent to const.hyp(c(5,5,5,5,5))
    # It is hard to tell if this is a bug in the list code or in how
    # closures are generated (i.e. functions that save their environment,
    # specifically the hypothesized.value variable in c.a.h.f()
    functions[[i]](0,0) 
  }
  names(functions) <- range
  return(functions)
}

# the default model of effect for make.conf.interval
sharp.null.hypothesis <- constant.additive.hypothesis.factory(0)

# lm.residualizer take a formula of the form y ~ ..., with all values
# included in the ... included in the data argument (a data.frame).
# it returns a function that takes new y values and applies the model to them.
# the inner function returns the residuals of the model for later test stats.
lm.residualizer <- function(fmla, data){
  # it might be more efficient to build the model once, and then use
  # update(). I'm proceeding with the simplier implementation for now
  function(adjusted.outcomes) {
    model <- lm(fmla, cbind(y = adjusted.outcomes, data))
    residuals(model)
  }
}


###################################################
### Randomization distribution producing function
###################################################


setClass("MVRandomizationDistribution",
  representation(moe = "list",
                 observed.outcomes = "data.frame",
                 observed.treatment = "numeric",
                 test.statistic = "function",
                 secondary.statistic.generators = "list",
                 secondary.statistic.results = "data.frame",
                 blocks = "numeric",
                 residualizer = "function",
                 samples = "numeric"),
  contains = "namedList")

# the only initialization I care about is setting the length of the object
# itself to the length of the models of effects
# this is commented out because of cryptic errors about returning a list
# not a RandomizationDistribution
#setMethod("initialize", "RandomizationDistribution",
#  function(.Object, moe, ...) {
#    length(.Object) <- length(moe)
#    # class(.Object) <- "RandomizationDistribution"
#    callNextMethod(.Object, moe = moe, ...)
#    .Object
#  })

# note: it might be better to move the r.d() function below into the initializer
# above and then have r.d() just call new with some defaults as necessary
produceRandomizations <- function(observed.treatment, blocks, samples) {
  if(!is.factor(blocks)){blocks<-factor(blocks)}
  # block size and the number of treated per block
  block.size <- table(blocks)
  block.treated <- aggregate(observed.treatment, list(blocks), sum)$x
  # these offsets will be needed later to turn relative into absolute indices
  block.starts <- append(0, cumsum(block.size))
  block.starts <- block.starts[1:(length(block.starts) - 1)]

  # overall statistics
  total.treated <- sum(observed.treatment)
  total.randomizations <- 1 # thers is a bug in BigIntegers that requires a loop. ??Contact them? Fix it?
  for (i in 1:(nlevels(blocks))) {
    total.randomizations = 
      total.randomizations * bigchoose(block.size[i], block.treated[i])
    NULL
  }
  # randomizations is a matrix (often abbreviated omega) of 
  # possible randomziations, for now ignoring blocks
  # it is generated either by direct enumeration (if small enough)
  # or by drawing from the distribution of randomizations
  if (total.randomizations > samples) {
    randomizations <- matrix(nrow = total.treated, ncol = samples)

    for (i in 1:samples) {
      raw.draws <- mapply(sample.int, block.size, block.treated, SIMPLIFY = F)
      randomizations[, i] <- unlist(mapply(function(b,o) { b + o }, raw.draws, block.starts))
    NULL }
  } else { # small enough to figured exactly 
    # we first get all the block level samples (these vary in size
    block.combinations <- mapply(combn, block.size, block.treated, SIMPLIFY = F)
        
    block.combinations <- mapply(function(b,o) { b + o }, block.combinations, block.starts, SIMPLIFY = FALSE)
    # then all the indexes we need to get unique permutations of block.combs
    expansions <- 
      expand.grid(
        mapply(seq, 
               rep(1, length(block.size)), 
               sapply(block.combinations, function(b) { dim(b)[2] }),
               SIMPLIFY = F))
    
    # loop thru, creating a unique set of combinations
    exp.count <- dim(expansions)[1]
    stopifnot(exp.count == total.randomizations)

    randomizations <- matrix(nrow = total.treated, ncol = exp.count)
    for(i in 1:exp.count) {
      # from each block, grab the combination indexed by the row in expansions
      # combine all such items into a row in randomizations
      randomizations[,i] <- unlist(
        mapply(function(block, i) { block[, i] }, 
               block.combinations, 
               expansions[i,]))
      NULL
    }
  }

  return(randomizations)
}

# note: this could be more efficient to have the moe and residualizer be NULL
# by default. if not specified the calls to moe() and residualizer() which are
# effectively no-ops, could be avoided.

# the RandomizationDistribution object
# for a single test statistic
# results of a call to randomizationDistribution, along with some of the
# each distribution is with respect to a given test statistic
setClass("RandomizationDistribution",
  representation(test.statistic = "function",
                 observed.test.stat = "numeric", # applying the test stat to
                                                 # the data
                 models.of.effect = "list", # of functions
                 treatment = "numeric", # 1/0 vector 
                 blocks = "numeric", # indicator vector
                 samples = "numeric", # number of samples
                 distribution = "matrix"))

# are the observed treatment and blocks necessary? Should I include data as
# well? What about ... arguments to the rDE() function?

# the engine 
randomizationDistributionEngine <- function(
  data,
  treatment,
  models, # a list of list(testStatistic, moe1, moe2, ...) all functions
  blocks = NULL,
  samples = 5000,
  ...) {

  n <- length(treatment)
  if (is.null(blocks)) blocks <- rep(TRUE, n)
  blocks <- as.factor(blocks)
  stopifnot(n == length(treatment) && n == length(blocks))
  
  # note: it might be better to pass total pool, number of treatment, and a
  # vector of the size of the blocks. Rather than a vector indicating
  # treatment and a vector indicating block membership. keeping this signature
  # for now, but will think about it...
  randomizations <- produceRandomizations(treatment, blocks, samples)
  
  expand.z <- function(i) { a <- numeric(n); a[i] <- 1; return(a) }

  k <- length(models)
  distributions <- vector("list", k)
   
  for (i in 1:k) {
    this.model <- models[[i]]
    test.statistic <- this.model[[1]]
    moes <- this.model[-1]

    this.distrib <- apply(randomizations, 2, function(z) {
      z <- expand.z(z)
      sapply(moes, function(m) { test.statistic(
                                   m(data, z, blocks, ...),
                                   z, blocks, ...)})})

    distributions[[i]] <- new("RandomizationDistribution",
      test.statistic = test.statistic,
      observed.test.stat = test.statistic(data, treatment, blocks, ...),
      models.of.effect = moes,
      treatment = as.numeric(treatment),
      blocks = as.numeric(blocks),
      samples = samples,
      distribution = this.distrib)   
  }
  # temporary hack until I can get the intialize() method working
  names(distributions) <- names(models)

  return(distributions)
}

# front end as per Ben's recommendation. parameters are specified, not fns.
setClass("ParameterizedRandomizationDistribution",
  representation(
    call = "call",
    params = "data.frame"),
  contains = "RandomizationDistribution")

parameterizedRandomizationDistribution <- function(
  data,
  treatment,
  test.stat,
  moe, # single function with signature f(data, z, blocks, param1, param2, etc.)
  parameters, # list of name = values, name = values, ...
  blocks = NULL,
  samples = 5000) {
  
  stopifnot(inherits(moe, "function"))
  stopifnot(inherits(parameters, "list"))

  parameter.space <- do.call(expand.grid, parameters)

  functions <- apply(parameter.space, 1, function(params) {
    force(params)
    do.call(moe, as.list(params))})

  rds <- randomizationDistributionEngine(data, treatment, models=list(c(test.stat,
  functions)), blocks, samples)

  rd <- as(rds[[1]], "ParameterizedRandomizationDistribution")
  rd@params <- parameter.space
  rd@call <- match.call()

  return(rd)
}

setClass("ParameterizedRandomizationDistributionSummary",
  representation(
    randomizationDistribution = "ParameterizedRandomizationDistribution",
    point.estimate = "data.frame",
    showCall = "logical"))

setMethod("summary", "ParameterizedRandomizationDistribution", function(object, 
  p.value.function = general.two.sided.p.value, showCall = T, ...) {

  # point estimate(s)
  pvs <- p.values(object, p.value.function)
  maxp <- max(pvs$p)
  point.estimate <- pvs[pvs$p == maxp, ]
  rownames(point.estimate) <- NULL
  
  # observed test statistic is in the PRD object 

  return(new("ParameterizedRandomizationDistributionSummary", randomizationDistribution = object,
    point.estimate = point.estimate, showCall = showCall))
})

setMethod("show", "ParameterizedRandomizationDistributionSummary", function(object) {
  if (object@showCall) {
    dc <- deparse(object@randomizationDistribution@call)
    cat("Call: ", dc[1], "\n")
    cat(paste("      ", dc[-1] , "\n", sep = ""))
    cat("\n")
  }

  cat("Observed Test Statistic: ")
  cat(object@randomizationDistribution@observed.test.stat)
  cat("\n\n")

  # save and show the call to x@randomizationDistribution
  cat("Hodges-Lehmann Point Estimate(s):\n")
  print(object@point.estimate)


  invisible(object)
})

setMethod("show", "ParameterizedRandomizationDistribution", function(object) {
  show(summary(object))
})

#### clippings, to be put into a front end

#   # check inputs
#   # if moe is not a list, make it one so we sapply later
#   if (!inherits(moe, "list")) {
#     moe <- c(moe)
#   }
#   
#   # ensure moe is composed of functions, top to bottom
#   for (i in 1:(length(moe))) {
#     stopifnot(inherits(moe[[i]], "function"))  
#   }
#   stopifnot(inherits(test.statistic, "function"),
#             is.null(secondary.statistic.generators) || all(sapply(secondary.statistic.generators, function(x) inherits(x, "function")))
#             )
#   
#   if (is.data.frame(observed.outcomes)) {n <- nrow(observed.outcomes)} else {n <- length(observed.outcomes)}
# 

#   secondaries <- if (is.null(secondary.statistic.generators)) NULL else
#   {
#     as.data.frame(lapply(secondary.statistic.generators,
#                          function(FUN)
#                          {
#                            apply(randomizations, 2, function(z) {
#                              z <- expand.z(z)
#                              FUN(observed.outcomes, z, blocks)}
#                                  )
#                          }
#                          )
#                   )
#   }



###################################################
### Using distributions
###################################################
# computes a one sided p-value, to compute double sided p-value you must
# use twice, with lower.tail set to TRUE, and propbably manually
# transforming your value to the equivalent lower tail value
simple.p.value <- function(value, distribution, lower.tail = FALSE) {
  sapply(distribution, function(d) {
    if (lower.tail) {
      length(d[d <= value]) / length(d) ##?is this better than mean(d<=value)
    } else { 
      length(d[d >= value]) / length(d)
    }
  })
}

manual.symmetric.p.value <- function(value, distribution, pivot = 0) {
  lower.weight <- simple.p.value(value, distribution)
  upper.weight <- simple.p.value(pivot - value, distribution, lower.tail = F)
  return(lower.weight + upper.weight)
}

auto.symmetric.p.value <- function(value, distribution) {
  sapply(distribution, function(d) {
    pivot <- mean(d)
    if (pivot > value) {
      lowerv <- value
      upperv <- pivot + abs(pivot - value)
    } else {
      upperv <- value
      lowerv <- pivot - abs(pivot - value)
    }

    upper.weight <- length(d[d >= upperv]) / length(d)
    lower.weight <- length(d[d <= lowerv]) / length(d)
    return(upper.weight + lower.weight)
  }) 
}

##There is some debate about two-sided p-values. I've been going with Rosenbaum 2009, Chapter 2, Footnote 2

##"In general, if you want a two-sided P-value, compute both one-sided P-values, double the smaller one, and take the minimum of this value and 1. This approach views the two-sided P-value as a correction for testing twice [10]. Both the sample mean in the current section and Wilcoxon’s signed rank statistic in the next section have randomization distributions under the null hypothesis that are symmetric, and with a symmetric null distribution there is little ambiguity about the mean- ing of a two-sided P-value. When the null distribution is not symmetric, different definitions of a two-sided P-value can give slightly different answers. As discussed in [10], the view of a two- sided P-value as a correction for testing twice is one sensible approach in all cases. For related results, see [55]."

general.two.sided.mid.p.value<-function(value,distribution){
  sapply(distribution,function(dist){
    high.mid<-mean(dist>value)+mean(dist==value)/2
    low.mid<-mean(dist<value)+mean(dist==value)/2
    return(2*min(low.mid,high.mid)) ##don't worry about bounding by 1 for now
  })
}

general.two.sided.p.value<-function(value,distribution){
  upper.tail <- mean(distribution >= value)
  lower.tail <- mean(distribution <= value)
  2 * min(upper.tail, lower.tail)
}

setClass("ParameterPvals", contains = "data.frame")

p.values <- function(object, p.value.function = general.two.sided.p.value) {
  stopifnot(inherits(object, "ParameterizedRandomizationDistribution"))

  k <- dim(object@distribution)[1] # total number of distributions to check
  pvs <- vector("numeric", k)

  for (i in 1:k) {
    pvs[i] <- p.value.function(object@observed.test.stat,
      object@distribution[i,])
  }

  return(as(cbind(object@params, p = pvs), "ParameterPvals"))
}

plot.ParameterPvals <- function(object, ...) {
  library(lattice)
  # for each pair of parameters, plot somethign like this:
  width <- dim(object)[2]  
  params <- colnames(object)[1:(width - 1)] # TODO function to get
  # it from the PRD
  fmla <- as.formula(paste("p ~ ", params[1], "+", params[2]))
  levelplot(fmla, data = object, ...)
}


##The confidence interval function ought to refuse to provide levels that are not supported (see the behavior of wilcox.test when you ask for an exact CI but the discreteness of the randomization distribution does not allow it

setClass("ParameterizedRandomizationDistributionConfInt",
  representation(
    prd = "ParameterizedRandomizationDistribution",
    rejections = "data.frame"))

confint.ParameterizedRandomizationDistribution <- function(
  object,
  param,
  level = 0.95,
  p.value.function = general.two.sided.p.value) {
  
  pvs <- p.values(object, p.value.function) 
  
  reject <- (1 - pvs$p) > level
  results <- cbind(object@params, reject)
   
  return(new("ParameterizedRandomizationDistributionConfInt", prd = object,
    rejections = results))
}

summary.ParameterizedRandomizationDistributionConfInt <- function(object, ...)
{
  width <- dim(object@rejections)[2]
  params <- colnames(object@rejections)[1:(width - 1)]
  minmax <- t(sapply(params, function(p) {
    tmp1 <- object@rejections[, p]
    tmp2 <- object@rejections[!object@rejections$reject, p]
    c(min(tmp1), min(tmp2), max(tmp2), max(tmp1))
  }))
  colnames(minmax) <- c("Tested Min", "Int. Min", 
                        "Int. Max", "Tested Max")
  return(minmax)
}

plot.ParameterizedRandomizationDistributionConfInt <- function(object, ...) {
  library(lattice)
  # for each pair of parameters, plot somethign like this:
  width <- dim(object@rejections)[2]  
  params <- colnames(object@rejections)[1:(width - 1)] # TODO function to get
  # it from the PRD
  fmla <- as.formula(paste("reject ~ ", params[1], "+", params[2]))
  levelplot(fmla, data = object@rejections, colorkey = F)
}

## displaying the results
## TODO: implement a lattice::histogram method
rd.histograms <- function(distribs) {
  for (i in 1:(length(distribs))) {
    print("Press <return> for next histogram")
    readline()
    hist(distribs[[i]])
  }
}

### While the combinatic stuff is neat, it is not needed right now
### I'm keeping it here for lack of a better place (Mark 20100619)

#combinatic.decoder.factory <- function(n, k) {
#  # n, k are fixed at the start to define the space of combinations
#  # i is the combinatic index (0 < i < n) 

#  # precompute a few sequences we'll need
#  ks <- k:1 # the bottom of the choose tests
#  max.combinadic <- bigchoose(n, k)
#  function(i) {
#    i <- as.bigz(i)
#    if (i < 1 || `>.bigz`(i, max.combinadic)) {
#      stop(paste("Combinatic out of range for", n, "choose", k ))
#    }
    
#    # part of the frequent translation to R's 1... sequences
#    i <- i - 1
#    # initialize a vector to hold the values
    
#    remaining <- i
#    previous.candidate <- n
#    combination <- numeric(k)

#    # I wonder if this could be improved with a binary search?
#    for(j in ks) {
#      value <- remaining + 1
#      while (`>.bigz`(value, remaining)) {
#        current.candidate <- previous.candidate - 1
#        value <- bigchoose(current.candidate, j)
#        previous.candidate <- current.candidate
#      }
#      remaining <- sub.bigz(remaining, value)
#      combination[j] <- current.candidate
#    }

#    return(combination + 1) # translate to 1... counting
#  }  
#}

#combinatic.encoder.factory <- function(n, k) {
#  ks <- 1:k    
#  function(encoded) {
#    stopifnot(length(encoded) == k)
#    encoded <- encoded - 1 # translate from 1... counting
#    expanded <- as.bigz(0)
#    for (i in ks) {
#      expanded <- `+.bigz`(expanded, bigchoose(encoded[i], i))
#    }
#    return(expanded + 1)
#  }  
#}

#cdf <- combinatic.decoder.factory(5, 3)
#cef <- combinatic.encoder.factory(5, 3)
#ttt <- sapply(1:10, cdf)
#for (i in 1:10) {
#  stopifnot(identical(cef(ttt[, i]), as.bigz(i)))
#}

## large combinadic spaces
#bigdf <- combinatic.decoder.factory(100, 50)
#stopifnot(length(bigdf(as.bigz("10000"))) == 50)
#stopifnot(length(bigdf(bigchoose(100, 50))) == 50)
#stopifnot(!all(
#  bigdf(bigchoose(100, 50)) ==
#  bigdf(bigchoose(100, 50) - 1)))

#bigef <- combinatic.encoder.factory(100, 50)
#bigef(bigdf(bigchoose(100, 50)))
#stopifnot(class(bigdf(as.bigz("10000"))) == "numeric")

## very large uniform random distribution
## this is probably wrong, but I'm going to use it for now
## generating a bumch of answers, 1/2 are over the 1/2 point
## and 1/2 are under for many runs (with some random variation)
## at least in the large, this is a reasonable (if slow) algorithm
#bigrandom <- function(n, k) {
#  # better future implementation takes nb for urand.bigz as an argument
#  # will need to generate a bunch, discard the invalds, repeat for the
#  # the subset.
#  function(seed = 0) {
#    max.value <- bigchoose(n, k)
#    base2size <- sizeinbase(max.value, 2)
#    candidate <- 0 
#    while (candidate < 1 || `>.bigz`(candidate, max.value)) {
#      candidate <- urand.bigz(1, size = base2size, seed = seed) + 1  
#    }
#    return(candidate)
#  }
#}

#bg <- bigrandom(100, 50)
#stopifnot(class(bg) == "function")
#stopifnot(class(bg()) == "bigz")
#stopifnot(bg() !=  bg())

