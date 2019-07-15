## ----setup, include=FALSE------------------------------------------------
library(roxygen2)
library(haven)
library(tidyverse)
library(optmatch)
library(RItools)

## ------------------------------------------------------------------------
# Reading in aspirin data
aspirin <- read_xpt("./RXQASA_H.XPT")

eq1 <- function(x) {
  tmp <- x == 1 
  tmp[is.na(tmp)] <- FALSE
  return(tmp)
}

aspirin$taking_aspirin <- eq1(aspirin$RXQ515) | eq1(aspirin$RXQ520)

## Blood pressure exam
# https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BPX_H.htm 

bp <- read_xpt("./BPX_H.XPT")

bp$sys_mean <- rowMeans(bp[, c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")], na.rm = TRUE)
bp$dia_mean <- rowMeans(bp[, c("BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4")], na.rm = TRUE)

## now we need to match up the aspirin and bp tables


temp <- aspirin[aspirin$SEQN %in% bp$SEQN, ]

rownames(bp) <- bp$SEQN

combined <- cbind(temp, bp[as.character(temp$SEQN), ])

### TODO: load DEMO_H.XPT and merge it in
### https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Demographics&CycleBeginYear=2013


save.image(file = "analyzing_aspirin_data.rda")

demo <- read_xpt("./DEMO_H.XPT")

data <- aspirin %>% left_join(demo)


