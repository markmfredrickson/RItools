###These are tests of the basic function of xBalance under some
###restricted sitations and using functions that hew closely to the
###expressions in Hansen and Bowers 2008. For example, with only
###binary treatment. The idea here to show that the math in that
###article is equivalent to the output from xBalance.

require("RItools")

data(nuclearplants)

s2.fn<-function(x){sum((x-mean(x))^2)/(length(x)-1)} ##same as sd(x)
h.fn<-function(n,m){(m*(n-m))/n}

var1<-function(x,m){ ##var(d)
  h<-h.fn(n=length(x),m=m)
  (1/h)*s2.fn(x)
}

var2<-function(x,m){ ##var(Z'x) (i.e. var of the sum statistic)
  h<-h.fn(n=length(x),m=m)
  (h)*s2.fn(x)
}

#####First just looking at the unstratified calculations
xb1a<-xBalance(pr~ date+ t1 + t2 + cap + ne + ct + bw + cum.n,
               data=nuclearplants,
               report=c("adj.means","adj.mean.diffs","std.diffs","z.scores","p.values"))

##print(xb1a,digits=4)

testxb1a<-t(sapply(nuclearplants[,dimnames(xb1a$results)$vars],function(thevar){
  myssn<-with(nuclearplants,sum((pr-mean(pr))*thevar))
  myadjdiff<-myssn/h.fn(m=sum(nuclearplants$pr),n=length(nuclearplants$pr))
  mynullvar1<-var1(x=thevar,m=sum(nuclearplants$pr))
  myz<-myadjdiff/sqrt(mynullvar1)
  return(c(adj.diff=myadjdiff,adj.diff.null.sd=sqrt(mynullvar1),z=myz))
}))

##print(testxb1a,digits=4)

all.equal(xb1a$results[,c("adj.diff","z"),"Unstrat"],testxb1a,check.attributes = FALSE)


