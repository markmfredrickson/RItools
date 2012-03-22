################################################################################
# Sampling: Drawing samples from sample spaces implied design choices (e.g.
# simple random sampling, srs with blocks, multinomial sampling, etc.)
################################################################################


simpleRandomSampler <- function(total, treated, z, b) {
  
  if (!missing(b)) {
    total <- table(b)
    reordering <- as.vector(sapply(unique(b), function(nm) { which(nm == b) }))
  }

  if (!missing(z) && !missing(b)) {
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
