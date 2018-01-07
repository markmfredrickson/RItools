set.seed(201801)
n <- 7L 
dat <- data.frame(y=rnorm(n), x=rnorm(n),
                  s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                  )
dat <- transform(dat, z=as.numeric( (x+rnorm(n))>0 ) )
## On-the-fly inverse probability of assignment weights
lm(y~z+x, weights=ipr(z,s), data=dat)

dat <- transform(dat, inv_pr_wt=ipr(z, s), inv_odds_wt=ipr(z, s, type = "odds vs 1"))
head(dat[3:6])

## Number of observations within a cluster
## doesn't affect inferred assignment probabilities  
dt <- transform(dat[3:6], clus=letters[1L:n])[rep(1L:n, rpois(n, 3)),]
with(dt, all.equal(inv_pr_wt, ipr(z, s, clus)))
dt <- data.frame(y=rnorm(nrow(dt)), x=rnorm(nrow(dt)),dt)
lm(y~z+x, weights=ipr(z,s,clus), data=dt) 
