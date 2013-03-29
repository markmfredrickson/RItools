xBalanceEngine <- function(ss,zz,mm,report, swt, s.p, normalize.weights,zzname)
{##ss is strata, zz is treatment, mm is the model matrix defined by the formula and data input to xBalance, swt is stratum weights, s.p. is the pooled sd, normalize.weights is logical (for creation of stratum weights)

	cnms <-
		c(
		  if ('adj.means'%in%report) c(paste(zzname,"0",sep="="),paste(zzname,"1",sep="=")) else character(0), ##c("Tx.eq.0","Tx.eq.1") else character(0),
		  if ('adj.mean.diffs'%in%report) 'adj.diff' else character(0),
		  if ('adj.mean.diffs.null.sd'%in%report) 'adj.diff.null.sd' else character(0),
		  if ('std.diffs'%in%report) 'std.diff' else character(0),
		  if ('z.scores'%in%report) 'z' else character(0),
		  if ('p.values'%in%report) 'p' else 'p'#character(0) turns out that it may be useful to have p-values in the object whether or not they are requested for printing
		  )


##Set up a data frame to contain the results    
ans <-
	if (length(setdiff(report,"chisquare.test")))
	{
		as.data.frame(matrix(0,dim(mm)[2], length(cnms),
				     dimnames= list(dimnames(mm)[[2]], cnms) )
		)
	} else
	{
		as.data.frame(matrix(0,0,0)
			      )
	}

	if (!length(zz))
		return(list(dfr=as.data.frame(matrix(0,0, length(cnms),
						     dimnames=list(character(0), cnms))),
			    chisq=c(pre.chisquare=0,pre.df=0,
				    post.chisquare=0,post.df=0)
			    )
	) 


	##Number of strata
	nlev <- nlevels(ss)

	### Calculate post.difference  
	ZtH <- unsplit(tapply(zz,ss,mean),ss) ##proportion treated (zz==1) in strata s.
	ssn <- drop(crossprod(zz-ZtH, mm*swt$wtratio)) ##weighted sum of mm in treated (using centered treatment)
	wtsum <- sum(unlist(tapply(zz,ss,function(x){var(x)*(length(x)-1)}))) ## h=(m_t*m_c)/m 
	post.diff <- ssn/(wtsum) ##diff of means
	if ('adj.mean.diffs'%in%report) ans[['adj.diff']] <- post.diff
	if ('std.diffs' %in% report) ans[['std.diff']] <- post.diff/s.p

	### Calculate post.Tx.eq.0, post.Tx.eq.1 --- now called "the treatment var"=0 and =1
	if ("adj.means"%in%report)
	{
		postwt0 <- unsplit(swt$sweights/tapply(zz<=0, ss, sum),
				   ss[zz<=0], drop=TRUE)
		ans[[paste(zzname,"0",sep="=")]] <- apply(mm[zz<=0,,drop=FALSE]*postwt0, 2,sum)
		ans[[paste(zzname,"1",sep="=")]] <- ans[[paste(zzname,"0",sep="=")]] + post.diff
	}


	msmn <- xBalance.make.stratum.mean.matrix(ss, (mm*swt$wtratio))

	tmat <- (mm*swt$wtratio - msmn)

	dv <- unsplit(tapply(zz,ss,var), ##sample variance of treatment
		      ss)
	ssvar <- apply(dv*tmat*tmat, 2, sum) ## for 1 column in  mm, sum(tmat*tmat)/(nrow(tmat)-1)==var(mm) and sum(dv*(mm-mean(mm))^2)=ssvar or wtsum*var(mm)

	##report (1/h)s^2. Since ssvar=(h)*s^2 multiply by (1/h)^2 to get (1/h)s^2.
	if ('adj.mean.diffs.null.sd' %in% report) ans[['adj.diff.null.sd']] <- sqrt(ssvar*(1/wtsum)^2) 

	if ('z.scores' %in% report) ans['z'] <- ifelse(ssvar<=.Machine$double.eps,0,
						       ssn/sqrt(ssvar))
	if (any(c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","std.diffs","z.scores","p.values") %in% report))  ##always produce a pvalue to use to create signif stars.
		ans['p'] <- ifelse(ssvar<=.Machine$double.eps,1,
				   2*pnorm(abs(ssn/sqrt(ssvar)),lower.tail=FALSE))

	if ("chisquare.test"%in%report)
	{
          require(svd)
          pst.svd <- try ( svd(tmat*sqrt(dv)) )
          if(inherits(pst.svd,'try-error')){
            pst.svd<-propack.svd(tmat*sqrt(dv))
          }
	##	pst.svd <- svd(tmat*sqrt(dv))
		Positive <- pst.svd$d > max(sqrt(.Machine$double.eps)*pst.svd$d[1], 0)
		Positive[is.na(Positive)]<-FALSE # JB Note: Can we imagine a situation in which we dont want to do this? 
		if (all(Positive)) ## is this faster? { ytl <- sweep(pst.svd$v,2,1/pst.svd$d,"*") }
		{ytl <- pst.svd$v *
		matrix(1/pst.svd$d, nrow=dim(mm)[2],ncol=length(pst.svd$d), byrow=T)
		} else{if (!any(Positive))
			ytl <- array(0, dim(mm)[2:1] )
		else ytl <- pst.svd$v[, Positive, drop = FALSE] * 
		matrix(1/pst.svd$d[Positive],ncol=sum(Positive),nrow=dim(mm)[2],byrow=TRUE)}

		mvz <- drop(crossprod(zz, tmat)%*%ytl)

		csq <- drop(crossprod(mvz))
		DF <- sum(Positive)

	} else csq <- DF <- 0

	list(dfr=ans,chisq=c('chisquare'=csq,'df'=DF))
}
