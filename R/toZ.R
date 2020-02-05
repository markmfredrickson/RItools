## Helper function to turn vectors into binary treatment indicators
##
## @param x The object to render as a treatment assignment vector of 1s and 0s
## @return A numeric vector of 1s and 0s
toZ <- function(x) { UseMethod("toZ") }

toZ.numeric <- function(x) { as.numeric(x > median(x)) }
toZ.logical <- function(x) { as.numeric(x) }
toZ.factor <- function(x) {
    fst <- levels(x)[1]
    as.numeric(fst != x)
}
