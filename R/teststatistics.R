################################################################################
# Test Statistic Functions for Randomization Distributions
#
# A test statistic is a function of three arguments:
# - y: the data, after adjustment
# - z: a treatment indicator (should be numeric, but calling as.numeric(z) is safe)
# - b: an indicator for the blocking factor, may contain only one level
#
# The function should return a scalar real value
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


harmonic.mean.difference <- function(ys, z, blocks) {
  z<-as.numeric(z)
  stopifnot(length(unique(z))==2) # require binary treatment for now

  h.fn <- function(n,m){(m*(n-m))/n} # harmonic weighting function
  
  d.b.fn<-function(ys,z,blocks){ ##blockwise mean diffs function is this better or worse that above?
    mapply(function(r,z){mean(ys[z==1])-mean(ys[z==0])},split(ys,blocks),split(z,blocks))
  } 

  h.b<-tapply(z,blocks,function(z){h.fn(n=length(z),m=sum(z))})
  d.b<-d.b.fn(ys=ys,z=z,blocks=blocks)
  return((1/sum(h.b))*sum(h.b*d.b)) ##notice this is the same as "adj.diff" from xBalance, that should be the test
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

mann.whitney.u <- function(ys, z, blocks) {
  z <- as.logical(z)
  ys.ranks <- rank(ys)
  return(sum(ys.ranks[z]) - sum(ys.ranks[!z]))
}




