require("RItools")

data(nuclearplants)

##################################################
### Basic uses
##################################################
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data=nuclearplants)

xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test'))

xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=list(pt=~pt), 
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test'))

(xb0 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=list(pt=~pt),
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test'))
)
##########################################################################################
### Oddness on LHS of formula
##########################################################################################Q
xBalance(I(pr==1) ~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=list(pt=~pt),
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test')
         )


#####################################################
######               xtable method                ###
#####################################################
if (require('xtable'))
  {
  xtablea <- xtable(xb0)
  xtableb <- xtable(xb0, caption="Caption!", label="thetable", digits=1,
       align=rep('l', prod(dim(xb0$result)[2:3])+1),
       display=c('s', rep(c(rep('fg', dim(xb0$result)[2]-1), 's'),
         dim(xb0$result)[3]) ) #,col.labels= do this one later
       )
  }


#####################################################
######               include.means                ###
#####################################################
###xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants,
###         covariate.scaling=1, include.means=TRUE)
###xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
###         ~factor(pt), nuclearplants, include.means=TRUE)
###
###
#####################################################
######  na.rm=FALSE with missing covariates       ###
#####################################################
### Should create a new variable (0=not missing,1=missing) and impute missing values with the mean (median is new default)

set.seed(123)
testdata<-nuclearplants
testdata$date[sample(1:32,10)]<-NA

xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, data = testdata,
    na.rm = FALSE,impfn=mean.default) ##first using the mean to match up with previous versions

#####################################################
######  handling factor with no obs for a level in a strata  ###
#####################################################

testdata<-nuclearplants
testdata$cum.nF<-factor(testdata$cum.n)

##Notice that for most levels of cum.n, there are no obs in one of the two strata
table(testdata$pt,testdata$cum.n)

##    1 2 3 5 6 7 8 11 12 14 15 16 17 18 19 20 21
##  0 4 3 4 2 1 1 1  0  2  1  1  1  1  1  1  1  1
##  1 0 0 0 0 0 1 2  3  0  0  0  0  0  0  0  0  0

##First no missing levels, same in both strata --- looks ok
xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.nF, strata= list(pt=~pt), data = testdata,impfn=mean.default)

##Second two missing levels, same in both strata
##This doesn't look as good --- we'd prefer to drop levels that don't exist.

testdata$cum.nF[testdata$cum.n>16]<-NA
testdata$cum.nF[testdata$cum.n==7]<-NA
table(testdata$pt,testdata$cum.nF,exclude=c()) ##Notice that the levels don't disappear by default.
xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.nF, strata=list(pt=~pt), data = testdata,na.rm=FALSE,impfn=mean.default)

##This isn't right either.
xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.nF, strata=list(pt=~pt), data = testdata,na.rm=TRUE,impfn=mean.default)


#####################################################
######  handling factor with no levels strata=argument  ###
#####################################################
testdata$badStrat <- rep(NA, dim(testdata)[1])
try(xBalance(pr ~ date, strata=list(~badStrat),
             data=testdata), FALSE)

#####################################################
######             WISHLIST                       ###
#####################################################
###
###
###
###
