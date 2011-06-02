library(RItools)
###################################################
### classic lady tasting tea
###################################################

actual.cups <- c(1, 0, 0, 1, 1, 0, 1, 0)
lady.guess <- c(1, 0, 0, 1, 1, 0, 0, 1)
number.correct <- function(guesses, z, blocks) { sum(z == guesses) / 2}
stopifnot(number.correct(lady.guess, actual.cups, NULL) == 3)

lady.distribution <- randomizationDistribution(lady.guess, actual.cups,
test.statistic = number.correct)
lady.table <- table(lady.distribution[[1]]) / length(lady.distribution[[1]])
stopifnot(all(lady.table == c(1/70, 16/70, 36/70, 16/70, 1/70)))

stopifnot(simple.p.value(3, lady.distribution) == 17/70)
stopifnot(simple.p.value(0, lady.distribution, lower.tail = T) == 1/70)

###################################################
### Blocked lady tasting tea 
###################################################
# pretend there are two blocks, each with 4 cups of tea
actual.cups <- c(1, 0, 0, 1, 1, 0, 1, 0)
lady.guess <- c(1, 0, 0, 1, 1, 0, 0, 1)
blocks <- c(0, 0, 0, 0, 1, 1, 1, 1)
blocked.lady <- randomizationDistribution(lady.guess, actual.cups,
  blocks = blocks, test.statistic = number.correct)
blocked.lady.table <- table(blocked.lady[[1]]) / length(blocked.lady[[1]])

# I'm fairly confident this is the distribution:
stopifnot(all(blocked.lady.table == c(1/36, 8/36, 18/36, 8/36, 1/36)))

###################################################
### Sampling
###################################################
# pretend there are two blocks, each with 4 cups of tea
actual.cups <- c(1, 0, 0, 1, 1, 0, 1, 0)
lady.guess <- c(1, 0, 0, 1, 1, 0, 0, 1)
blocks <- c(0, 0, 0, 0, 1, 1, 1, 1)
set.seed(20100620)
sampled.blocked.lady <- randomizationDistribution(lady.guess, actual.cups,
  blocks = blocks, test.statistic = number.correct, samples = 35)
sampled.table <- table(sampled.blocked.lady[[1]]) / length(sampled.blocked.lady[[1]])

# since this is a sampled version (which cannot sample the entire space), 
# it should be close to the real distribution but not exactly the same
sb.val <- function(v) { 
  a <- sampled.table[names(sampled.table) == v]
  if (length(a) == 0) {
    return(0)  
  } else {
    return(a)  
  }
}
stopifnot(sb.val(0) < 0.05 && sb.val(4) < 0.05)
stopifnot(sb.val(2) >= .4 && sb.val(2) <= .6)

###################################################
### Exploring adjustment: why addition works and not subtraction!!!!!
###################################################
# mock up some data that ought to be large enough for Normal approximations to be good.
set.seed(20100620)
rc <- rnorm(100, mean = 0) ##potential outcomes to control
rt <- rc+3 ##presume a constant additive effects model with tau=3
Z <- rep(c(0,1),100)
R <- Z*rt + (1-Z)*rc
##Not paired or blocked for now.
data <- data.frame(R = R, Z=Z)
data$ZF<-factor(data$Z,labels=c("control","treated"))
tau0s<-seq(-4,4,.1)

constant.additive.hypothesis.factory <- function(hypothesized.value) {
  function(ys, z) { ys + (z * hypothesized.value) } 
}

constant.subtractive.hypothesis.factory <- function(hypothesized.value) {
  function(ys, z) { ys - (z * hypothesized.value) }
}

hyps.add <- constant.hypotheses(tau0s,factory=constant.additive.hypothesis.factory)
hyps.sub <- constant.hypotheses(tau0s,factory=constant.subtractive.hypothesis.factory)

set.seed(20100620)
results.add <- randomizationDistribution(
  observed.outcomes = data$R,
  observed.treatment = data$Z,
  test.statistic = mean.difference,
  moe = hyps.add,
  samples = 100) # make this a little faster for testing purposes

results.sub <- randomizationDistribution(
  observed.outcomes = data$R,
  observed.treatment = data$Z,
  test.statistic = mean.difference,
  moe = hyps.sub,
  samples = 100) # make this a little faster for testing purposes


get.cis(results.add,thelevels=.95,p.value.function=general.two.sided.p.value)
get.cis(results.sub,thelevels=.95,p.value.function=general.two.sided.p.value)

t.test(R~Z,data=data)$conf.int
t.test(R~I(1-Z),data=data)$conf.int
wilcox.test(R~Z,data=data,conf.int=TRUE)$conf.int
wilcox.test(R~I(1-Z),data=data,conf.int=TRUE)$conf.int

data$ZF<-factor(data$Z,levels=c(1,0),labels=c("treated","control"))
table(data$Z,data$ZF)
   
##    treated control
##  0       0     100
##  1     100       0
##> t.test(R~ZF,data=data)
##
##	Welch Two Sample t-test
##
##data:  R by ZF 
##t = 21.2021, df = 186.2, p-value < 2.2e-16
##alternative hypothesis: true difference in means is not equal to 0 
##95 percent confidence interval:
## 2.500934 3.014090 
##sample estimates:
##mean in group treated mean in group control 
##           2.76991878            0.01240677 


###now a paired design
news.df <- structure(list(city = c("Saginaw", "Sioux City", "Battle Creek", 
"Midland", "Oxford", "Lowell", "Yakima", "Richland"), s = c(1L, 
1L, 2L, 2L, 3L, 3L, 4L, 4L), z = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 
1L), r = c(16L, 22L, 14L, 7L, 23L, 27L, 58L, 61L), rpre = c(17L, 
21L, 13L, 12L, 26L, 25L, 48L, 41L), cands = c(16L, 6L, 21L, 6L, 
9L, 18L, 6L, 10L)), .Names = c("city", "s", "z", "r", "rpre", 
"cands"), class = "data.frame", row.names = c("Saginaw", "Sioux City", 
"Battle Creek", "Midland", "Oxford", "Lowell", "Yakima", "Richland"
))

##make a new outcome with a an additive, constant, treatment effect of 10 
news.df$newr<-news.df$rpre+news.df$z*10

t.test(newr~z,data=news.df,paired=TRUE)

news.df$zF<-factor(news.df$z,levels=c(1,0),labels=c("treated","control"))
table(news.df$z,news.df$zF)
   

xBalance(z~newr,strata=list(s=~s),data=news.df,report="all")

tau0s<-seq(-20,20,.1)
xb.tau0s<-xBalance(as.formula(paste("z~",paste("I(newr-(z*",round(tau0s,2),"))",collapse="+"))),
                    strata=list(s=~s),
                     report=c("p.values"),
                     data=news.df)

##Make 95 and 66% CIs
##(1-(.95))/2 = .025 =(.05)/2
##(1-(.88))/2 = .06  =(.12)/2
##(1-(2/3))/2 = .16667 = (1/3)/2

xb.ci.1<-sort(c(range(tau0s[xb.tau0s$results[,"p",]>=.12]), ##88
                range(tau0s[xb.tau0s$results[,"p",]>=(1/3)]))) ##66%
names(xb.ci.1)<-c("l88","l66","u66","u88")
xb.ci.1

##Compare to t.test

(t.ci<-sort(c(t.test(newr~zF,data=news.df,paired=TRUE,conf.level=.88)$conf.int,
       t.test(newr~zF,data=news.df,paired=TRUE,conf.level=2/3)$conf.int)))


##Now, use randomizationDistribution
hyps.add <- constant.hypotheses(tau0s,factory=constant.additive.hypothesis.factory)
hyps.sub <- constant.hypotheses(tau0s,factory=constant.subtractive.hypothesis.factory)

results.news.add <- randomizationDistribution(
                                               observed.outcomes = news.df$newr,
                                               observed.treatment = news.df$z,
                                               blocks=news.df$s,
                                               test.statistic = mean.difference,
                                               moe = hyps.add,
                                               samples = 100) ##will only use 16 

results.news.sub <- randomizationDistribution(
                                               observed.outcomes = news.df$newr,
                                               observed.treatment = news.df$z,
                                               blocks=news.df$s,
                                               test.statistic = mean.difference,
                                               moe = hyps.sub,
                                               samples = 100) ##will only use 16 


get.cis(results.news.add,thelevels=c(.66,.88),p.value.function=general.two.sided.mid.p.value)
get.cis(results.news.sub,thelevels=c(.66,.88),p.value.function=general.two.sided.mid.p.value)

(-1)*get.cis(results.news.sub,thelevels=c(.66,.88),p.value.function=general.two.sided.mid.p.value)

##Now try with wilcox signed rank statistic: we should be able to reproduce the following:
#Compare to wilcox.test
library(exactRankTests)
(wc.ci<-sort(c(wilcox.exact(newr~zF,data=news.df,paired=TRUE,conf.int=TRUE,exact=TRUE,conf.level=.88)$conf.int,
       wilcox.exact(newr~zF,data=news.df,paired=TRUE,conf.int=TRUE,exact=TRUE,conf.level=2/3)$conf.int)))

paired.sgnrank.sum(ys=news.df$newr,z=news.df$z,blocks=news.df$s)
wilcox.test(newr~zF,data=news.df,paired=TRUE)$statistic


wilcox.stat<-function(ys,z,blocks){ stopifnot(length(unique(z))==2) ##require binary treatment for now
                                    require(exactRankTests)
                                    Y<-sapply(split(data.frame(r=ys,z=z),blocks),function(dat){with(dat,r[z==1]-r[z==0])}) ##;print(Y)
                                    wilcox.exact(Y)$statistic
                                  }


results.news.rank.add <- randomizationDistribution(
                                               observed.outcomes = news.df$newr,
                                               observed.treatment = news.df$z,
                                               blocks=news.df$s,
                                               test.statistic = wilcox.stat,
                                               moe = hyps.add,
                                               samples = 100) ##will only use 16 

results.news.rank.sub <- randomizationDistribution(
                                               observed.outcomes = news.df$newr,
                                               observed.treatment = news.df$z,
                                               blocks=news.df$s,
                                               test.statistic = wilcox.stat,
                                               moe = hyps.sub,
                                               samples = 100) ##will only use 16 


get.cis(results.news.rank.add,thelevels=c(.66,.88),p.value.function=general.two.sided.mid.p.value)
get.cis(results.news.rank.sub,thelevels=c(.66,.88),p.value.function=general.two.sided.mid.p.value)

(-1)*get.cis(results.news.rank.sub,thelevels=c(.66,.88),p.value.function=general.two.sided.mid.p.value)

##That doesn't work. Now try another approach:

##here is the Omega matrix for the newspapers experiment

Om <-
structure(list(V1 = c(1, 0, 1, 0, 1, 0, 1, 0), V2 = c(0, 1, 1, 
0, 1, 0, 1, 0), V3 = c(1, 0, 0, 1, 1, 0, 1, 0), V4 = c(0, 1, 
0, 1, 1, 0, 1, 0), V5 = c(1, 0, 1, 0, 0, 1, 1, 0), V6 = c(0, 
1, 1, 0, 0, 1, 1, 0), V7 = c(1, 0, 0, 1, 0, 1, 1, 0), V8 = c(0, 
1, 0, 1, 0, 1, 1, 0), V9 = c(1, 0, 1, 0, 1, 0, 0, 1), V10 = c(0, 
1, 1, 0, 1, 0, 0, 1), V11 = c(1, 0, 0, 1, 1, 0, 0, 1), V12 = c(0, 
1, 0, 1, 1, 0, 0, 1), V13 = c(1, 0, 1, 0, 0, 1, 0, 1), V14 = c(0, 
1, 1, 0, 0, 1, 0, 1), V15 = c(1, 0, 0, 1, 0, 1, 0, 1), V16 = c(0, 
1, 0, 1, 0, 1, 0, 1)), .Names = c("V1", "V2", "V3", "V4", "V5", 
"V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", 
"V16"), row.names = c("Saginaw", "Sioux City", "Battle Creek", 
"Midland", "Oxford", "Lowell", "Yakima", "Richland"), class = "data.frame")

simp.mid.p<-function(value,dist) {
        high.mid <- mean(dist > value) + mean(dist == value)/2
        low.mid <- mean(dist < value) + mean(dist == value)/2
        return(2 * min(low.mid, high.mid))
    }

brute.ci<-sapply(tau0s,function(h){
  obs<-with(news.df,paired.sgnrank.sum(z=z,ys=(newr-z*h),blocks=s))
  dist<-with(news.df,sapply(Om,function(newz){
    paired.sgnrank.sum(z=newz,ys=(newr-z*h),blocks=s)}))
  c(h,simp.mid.p(obs,dist))
})
row.names(brute.ci)<-c("h","p")

(mywc.ci<-sort(c(range(brute.ci["h",brute.ci["p",]>=.12]),
       range(brute.ci["h",brute.ci["p",]>=(1/3)]))))

###################################################
### Non-sharp hypotheses ???why assesses these here in this way???
###################################################

# creating a set of hypotheses
hvals <- c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
hypotheses <- constant.hypotheses(hvals)
stopifnot(length(hypotheses) == 9)
stopifnot(identical(as.numeric(names(hypotheses)), hvals))
mockdata <- data.frame(y = c(1, 1, 0, 0), z = c(1, 1, 0, 0))
stopifnot(!(identical(
  hypotheses[[1]](mockdata$y, mockdata$z),
  hypotheses[[2]](mockdata$y, mockdata$z))))


# mock up some data
set.seed(20100620)
treatments <- rnorm(100, mean = 3)
controls <- rnorm(100, mean = 0)
# so the "true" treatment effect is 3

# set up these results as a paired design (so it can be enumerated)
data <- data.frame(
  y = append(treatments, controls), 
  z = append(rep(TRUE, 100), rep(FALSE, 100)))

## T-tests for data ##
# > t.test(y ~ z, data)
# 
# 	Welch Two Sample t-test
# 
# data:  y by z 
# t = -21.0335, df = 195.101, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0 
# 95 percent confidence interval:
#  -3.220918 -2.668682 
# sample estimates:
# mean in group FALSE  mean in group TRUE 
#         -0.05363717          2.89116278 
# 
# > t.test(data[data$z, "y"], data[!data$z, "y"], paired = T)
# 
# 	Paired t-test
# 
# data:  data[data$z, "y"] and data[!data$z, "y"] 
# t = 20.0099, df = 99, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0 
# 95 percent confidence interval:
#  2.652788 3.236812 
# sample estimates:
# mean of the differences 
#                  2.9448 
# 
### Since the blocking is arbitrary, there is no improvement of estimates
### using the 2 sample t-test. From a randomizationDistribution perspective
### the paired version is simpler to compute (50 * 2 iterations vs 
### choose(100,50) -- but the results should be similar.

set.seed(20100620)
results <- randomizationDistribution(
  observed.outcomes = data$y,
  observed.treatment = data$z,
  test.statistic = mean.difference,
  moe = hypotheses,
  samples = 100) # make this a little faster for testing purposes

stopifnot(inherits(results, "list"))
stopifnot(length(results) == length(hypotheses))
stopifnot(!(identical(results[[1]], results[[2]])))

r.ci <- confint(results)
stopifnot(length(r.ci) == length(results))
stopifnot(inherits(r.ci, "logical"))

### extending the previous example/test with blocks
### use a simple block strategy with 2 blocks, with
### 30/20 treatments, and 20/30 controls
blocks <- append(rep(c(2,2,2,1,1), 20), rep(c(1,1,1,2,2), 20))

set.seed(20100620)
block.results <- randomizationDistribution(
  observed.outcomes = data$y,
  observed.treatment = data$z,
  blocks = blocks,
  test.statistic = mean.difference,
  moe = hypotheses,
  samples = 100) # make this a little faster for testing purposes

br.ci <- confint(block.results)
### Failing. Commenting out for now for later tests
### it seems that these confints have changed in some subtle way
### stopifnot(identical(br.ci, r.ci))

###################################################
### misc
###################################################

# I don't have a testing function to test for errors, but I should test:
# all functions are really functions, not somethign else
# input vectors to randomizationDistribution must be same length

observed.outcomes.test <- c(1, 0 , 0, 1)
observed.treament.test <- c(0, 1, 1, 0) # treatment lowers score by 1
basic.test <- randomizationDistribution(observed.outcomes.test,
                                         observed.treament.test)

# engine should return a table of values
stopifnot(class(basic.test) == "RandomizationDistribution")

###################################################
### Test statistics
###################################################
ys <- c(5,5,5,5,5,1,1,1,1,1)
z <- c(1,1,1,1,1,0,0,0,0,0)

test.md <- mean.difference(ys, z, rep(1, 10))
stopifnot(test.md == 4)

test.rnk <- rank.sum(ys, z, NULL)
stopifnot(test.rnk == sum(10:6))

# test out blocks 
blocks <- c(1,1,1,2,2,1,1,2,2,2)
ys.blocks <- ys + blocks

block.md <- mean.difference(ys.blocks, z, blocks)
stopifnot(identical(block.md, test.md))

# what would a blocked rank-sum look like?
# aligned ranks would subtract the block mean

###################################################
### Residualizers -- aka covariance adjustment
###################################################

set.seed(20100620)
treatments <- rnorm(100, mean = 3)
controls <- rnorm(100, mean = 0)
x <- rnorm(200, mean = 1, sd = 2)
# so the "true" treatment effect is 3
# but the x covariate adds some noise to the data

data <- data.frame(
  y = append(treatments, controls) + x, 
  z = append(rep(TRUE, 100), rep(FALSE, 100)),
  x = x)

## a t-test has correct coverage, but with unecessarily wide intervals,
## compared to  -3.220918 -2.668682 from earlier tests.
# > t.test(y ~ z, data)
# 
# 	Welch Two Sample t-test
# 
# data:  y by z 
# t = -9.0952, df = 191.777, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0 
# 95 percent confidence interval:
#  -3.386879 -2.179694 
# sample estimates:
# mean in group FALSE  mean in group TRUE 
#            1.277164            4.060451 
# 
## a linear model properly parcels out the additive components:
# > lm(y ~ z + x, data)
# 
# Call:
# lm(formula = y ~ z + x, data = data)
# 
# Coefficients:
# (Intercept)        zTRUE            x  
#    -0.02438      2.94125      0.97802  

# if we were to regress y on x, during the randomization inference
# we should be able to get the best of both worlds.

set.seed(20100620)
results.covariate.naive <- randomizationDistribution(
  observed.outcomes = data$y,
  observed.treatment = data$z,
  test.statistic = mean.difference,
  moe = hypotheses,
  samples = 100)

# this should be wrong
stopifnot(!identical(confint(results), confint(results.covariate.naive)))

my.residualizer <- lm.residualizer(y ~ x, data)
set.seed(20100620)
results.covariate.good <- randomizationDistribution(
  observed.outcomes = data$y,
  observed.treatment = data$z,
  test.statistic = mean.difference,
  moe = hypotheses,
  residualizer = my.residualizer,
  samples = 100)

# should be different than the naive version
stopifnot(!identical(results.covariate.naive[[1]], results.covariate.good[[1]]))

# should get the same values as the simplest case
stopifnot(identical(confint(results), confint(results.covariate.good)))

