xBalance.makeMM <- function(tfm, dat)
{
attr(tfm, "intercept") <- 1
mf <- model.frame(tfm, dat, na.action=na.pass) ##na.pass leaves the NAs in dat in the model.frame
tlbl <- names(mf)
names(tlbl) <- as.character(tlbl)
clist <- lapply(mf,
                function(x) {
                    if (is.factor(x))
                        structure(diag(nlevels(x)),
                                  dimnames=list(levels(x), levels(x)))
                    else NULL })
clist <- clist[!sapply(clist, is.null)]
mm <- model.matrix(tfm,mf,contrasts.arg=clist)
matrix(mm[,dimnames(mm)[[2]]!="(Intercept)"],
       dim(mm)[1], dim(mm)[2]-1,
       dimnames=list(dimnames(mm)[[1]],dimnames(mm)[[2]][-1]) )

}
