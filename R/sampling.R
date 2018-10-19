################################################################################
# Sampling: Drawing samples from sample spaces implied design choices (e.g.
# simple random sampling, srs with blocks, multinomial sampling, etc.)
################################################################################

##' Simple Random Sampler
##'
##' The simple random sampler creates a function that draws from the set of all 
##' randomizations with a fixed number of treated units within blocks
##' @param total Size of the experimental pool
##' @param treated Total number to be assigned treatment. One of \code{treated}
##' or \code{z} is required
##' @param z Treatment indicator vector, must be binary or logical
##' @param b Records the membership in blocks within which treatment is assigned
##' @return Function to use in \code{sampler} argument of \code{RItest}
##' @seealso \code{\link{RItest}}
##' @examples  
##' data("nuclearplants")
##' region.blocks <- simpleRandomSampler(total = 32, z = nuclearplants$pr, b = nuclearplants$ne)
##' @export
simpleRandomSampler <- function(total, treated, z, b) {
  naError <- function() { stop("NAs not allowed in arguments to simpleRandomSampler")}

  if (!missing(b)) {
    # since b is passed, it must be free of NA values
    if (any(is.na(b))) {
      naError()
    }

    total <- table(b)
    reordering <- as.vector(unlist(lapply(unique(b), function(nm) { which(nm == b) })))
  }

  if (!missing(z) && !missing(b)) {
    # b is checked above, check z for NAs
    if (any(is.na(z))) {
      naError()
    }

    treated <- aggregate(z, list(b), sum)$x
  }

  if (!missing(total) && !missing(z) && missing(treated)) {

    if (sum(total) != length(z)) {
      stop("Total units per blocks must agree with 'z' argument")
    }

    # tmp is the start of each block, tmp1, the end of each
    # from this, we can add up the 1's in each block
    # note: these could probably be unified with block.starts below.
    tmp <- c(0, cumsum(total)[1:(length(total) - 1)]) + 1
    tmp1 <- cumsum(total)
    treated <- mapply(function(a,b) { sum(z[a:b]) }, tmp, tmp1)
  }

  # both total and treated exist by now, so check them for NAs
  if (any(is.na(total)) || any(is.na(treated))) {
    naError()
  }

  if (length(total) != length(treated)) {
    stop("total and treated arguments must match lengths")
  }

  if (!all(total > treated)) {
    stop("Each block must contain strictly more observations than treated units")
  }

  if (!all(treated > 0)) {
    stop("Each block must have at least one treated unit")
  }

  if (length(treated) != length(total)) {
    stop("Arguments 'total' and 'treated' must be the same size")
  }

  # precompute some constants used in the process
  # the size of the omega sample space (log'ed)
  total.randomizations <- sum(lchoose(total, treated))

  # the indices of the block starting locations
  block.starts <- append(0, cumsum(total))
  block.starts <- block.starts[1:(length(block.starts) - 1)]

  # a useful little helper for turning indices into a 1/0 treatment vector
  n <- sum(total)
  expand.z <- function(i) { a <- numeric(n); a[i] <- 1; return(a) }

  # b wasn't passed, we'll need to create a reordering vector for later
  # but it will just be the "identity" reordering.
  if (missing(b)) {
    reordering <- 1:n
  }

  # total number of treated units
  total.treated <- sum(treated)

  function(samples) {

    if (total.randomizations > log(samples)) {
      randomizations <- matrix(nrow = total.treated, ncol = samples)

        for (i in 1:samples) {
          raw.draws <- mapply(sample.int, total, treated, SIMPLIFY = F)
          randomizations[, i] <- unlist(mapply(function(b,o) { b + o }, raw.draws, block.starts))
          NULL
        }
    } else { # small enough to figured exactly
      # we first get all the block level samples (these vary in size
      block.combinations <- mapply(combn, total, treated, SIMPLIFY = F)

      block.combinations <- mapply(function(b,o) { b + o }, block.combinations, block.starts, SIMPLIFY = FALSE)
      # then all the indexes we need to get unique permutations of block.combs
      expansions <-
      expand.grid(
        mapply(seq,
          rep(1, length(total)),
          sapply(block.combinations, function(b) { dim(b)[2] }),
          SIMPLIFY = F))

      # loop thru, creating a unique set of combinations
      exp.count <- dim(expansions)[1]

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

    # arg. I couldn't get any of the apply family to work with a matrix ->
    # matrix workflow (could do it with a cast to a data.frame, but trying to
    # be efficient
    actual.samples <- dim(randomizations)[2]
    tmp <- matrix(nrow = n, ncol = actual.samples)
    for (i in 1:actual.samples) {
      x <- expand.z(randomizations[,i])
      x[reordering] <- x
      tmp[, i] <- x
    }

    return(list(weight = 1, samples = tmp))
  } # end inner function
}


##' Independent Probability Sampler
##'
##' Draw from the 2^n possible randomizations for n units with
##' probability vector p (e.g. where each unit gets its own, possibly
##' biased, coin flip)
##' @param n Size of the experimental pool
##' @param p Vector of probabilities of treatment assignment for each unit. 
##' By default, every observation has the same probability of being assigned to
##' treatment and control
##' @return Function to use in \code{sampler} argument of \code{RItest}
##' @seealso \code{\link{RItest}}
##' @examples 
##' sameprob <- independentProbabilitySampler(n = 32)
##' @export
independentProbabilitySampler <- function(n, p = rep(0.5, n)) {

  function(samples) {

    if (log(samples, base = 2) >= n) {
      # enumerate
      zs <- matrix(0, nrow = n, ncol = 2^n)

      # this seems a little convoluted, but it does the job
      indexes <- lapply(1:n, function(k) { combn(n, k) })
      matrices <- lapply(indexes, function(idxgrp) {
        apply(idxgrp, 2, function(i) { a <- numeric(n); a[i] <- 1; a })
      })

      zs <- cbind(do.call(cbind, matrices), 0)

    } else {
      # generating a matrix of zero and ones is pretty easy. just compare to the
      # p matrix over and over.
      zs <- 0 + matrix(runif(n * samples) < p, nrow = n, ncol = samples)
    }

    # probablility of seeing a 1 is p, the probability of seeing a zero is 1 -
    # p, so the total probability of draw is:
    # (when p == 0.5, save some computation now and during p-value
    # computation)
    if (all(p == 0.5)) {
      weight <- 1
    } else {
      weight <- apply(p * zs + (1 - p) * (1 - zs), 2, prod)
    }

    return(list(weight = weight, samples = zs))
  }
}


##' Cluster Random Sampler
##' 
##' Creates a function that draws from the set of all randomizations, assuming
##' treatment assignment to all observations within a cluster
##' @param clusters Cluster membership vector
##' @param z Treatment assignment vector, must be binary or logical
##' @return Function to use in \code{sampler} argument of \code{RItest}
##' @seealso \code{\link{RItest}}
##' @examples 
##' clusters <-  sample(letters[1:3], 10, replace = TRUE)
##' treatment <- ifelse(clusters == sample(letters[1:3], 1), 1, 0)
##' clustered <- clusterRandomSampler(clusters = clusters, z = treatment)
##' @export
clusterRandomSampler <- function(clusters, z) {
  naError <- function() {stop("NAs not allowed in arguments to clusterRandomSampler")}
  
  # Check for NAs
  if (any(is.na(clusters)) || any(is.na(z))) {
    naError()
  }
  
  # Check length
  if (length(clusters) != length(z)){
    stop("'clusters' and 'z' arguments must match lengths")
  }
  
  # Treatment must be binary for now
  if(length(unique(z)) > 2){
    stop("treatment 'z' must be binary or logical")
  }
  
  # Number of observations by cluster
  total.cl <- table(clusters)

  # Reorder observations by cluster  
  reordering <- as.vector(unlist(lapply(sort(unique(clusters)), function(nm) {which(nm == clusters)})))
  
  # Determine which clusters get treatment
  treated.cl <- aggregate(z, list(clusters), sum)$x
  
  treated.cl <- ifelse(treated.cl > 0, 1, 0)
  
  # Useful quantities
  total.treated.cl <- sum(treated.cl)
  
  n.clusters <- length(unique(clusters))
  
  total <- length(z)
  
  # Check that treatment assignment and number of clusters are consistent
  if (total.treated.cl < 1) {
    stop("At least one cluster must be assigned to treatment")
  }
  
  if (n.clusters < 2) {
    stop("At least two unique 'clusters' must exist")
  }
  
  function(samples){# begin inner function
    # empty matrix of clusters x samples
    randomizations <- matrix(nrow = n.clusters, ncol = samples)
    
    # for each sample, randomly assign treatment to cluster
    for(i in 1:samples) {
      randomizations[,i] <- sample(treated.cl)
    }
    
    # each row is a cluster, turn into list
    list.cl <- rep(list(NULL), n.clusters)
    
    for(i in 1:n.clusters){
      list.cl[[i]] <- matrix(rep(randomizations[i,], total.cl[i], 
                                 ncol = samples, byrow = TRUE))
    }
    
    # Put them all together
    randomizations <- do.call(rbind, list.cl)
    
    # empty matrix total x samples
    tmp <- matrix(nrow = total, ncol = samples)
    
    # fill based on reordering
      for(i in 1:total){
      tmp[reordering[i],] <- randomizations[i,]
    }
    
    return(list(weight = 1, samples = tmp))
  } #end inner function
}
