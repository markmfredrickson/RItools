################################################################################
# P-value Generating/Displaying Functions
################################################################################

############################## P-Value Generators ##############################

simple.p.value <- function(value, distribution, lower.tail = FALSE) {
  if (lower.tail) {
    mean(distribution <= value)
  } else { 
    mean(distribution >= value)
  }
}

upper.p.value <- function(value, distribution) { simple.p.value(value, distribution, FALSE) }
lower.p.value <- function(value, distribution) { simple.p.value(value, distribution, TRUE) }

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
  min(2 * min(upper.p.value(value, distribution), lower.p.value(value, distribution)))
}

############################## Using RDs ##############################

setClass("ParameterPvals", contains = "data.frame")

p.values <- function(object, p.value.function = general.two.sided.p.value) {
  stopifnot(inherits(object, "ParameterizedRandomizationDistribution"))

  k <- dim(object)[1] # total number of distributions to check
  pvs <- vector("numeric")

  # first row is the sharp null, that can be computed elsewhere
  for (i in 1:k) {
    pvs[i] <- p.value.function(object[i, 1],
      object[i, -1])
  }

  return(as(cbind(rbind(NA, object@params), p = pvs), "ParameterPvals"))
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
  results <- cbind(object@params, reject[-1]) # drop the null
   
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


