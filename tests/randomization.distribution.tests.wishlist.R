###########
### Accept multi-dimensional test statistic
###########

number.correct <- function(data, z, blocks) { sum(z == data$guess) / 2 }
temp.test <- function(data, z, blocks) { mean(data$temp[z == 1]) }

actual.cups <- c(1, 0, 0, 1, 1, 0, 1, 0)
ladydata <- data.frame(guess=c(1, 0, 0, 1, 1, 0, 0, 1), temp=rep(c(145,125,135,155),2))
stopifnot(all.equal(number.correct.and.mean.temp(ladydata, actual.cups, NULL), c(guess=3,temp=145)))

ldist2 <-  randomizationDistribution(ladydata, actual.cups,
test.statistic = number.correct, secondary.statistic.generators=list("temperature"=temp.test))

stopifnot(is.numeric(ldist2[[1]]), is.data.frame(ldist2@secondary.statistic.results) )
#stopifnot(is.matrix(ldist2[[1]]), all.equal(dimnames(ldist2[[1]])[[1]],c("guess", "temp")))
numcorr.table <- table(ldist2[[1]]) / length(ldist2[[1]])
stopifnot(all(numcorr.table == c(1/70, 16/70, 36/70, 16/70, 1/70)))

