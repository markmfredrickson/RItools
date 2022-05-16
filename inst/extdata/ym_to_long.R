## Create a version of the Yudkin and Moher ASSIST data that is disagregated

library(tidyverse)
library(randomizr)
library(RItools)

ym <- structure(list(Practice = 1:21, n = c(28L, 30L, 38L, 51L, 56L,
58L, 84L, 84L, 91L, 94L, 96L, 114L, 121L, 123L, 127L, 131L, 135L,
138L, 139L, 160L, 244L), assessed = c(14L, 23L, 16L, 35L, 23L,
33L, 20L, 19L, 25L, 18L, 22L, 40L, 12L, 18L, 46L, 22L, 44L, 49L,
32L, 21L, 38L), aspirin = c(79L, 73L, 79L, 96L, 84L, 66L, 90L,
77L, 66L, 74L, 72L, 75L, 77L, 72L, 81L, 76L, 85L, 77L, 81L, 67L,
74L), hypo = c(43L, 67L, 45L, 63L, 54L, 53L, 61L, 43L, 62L, 45L,
46L, 53L, 37L, 64L, 68L, 75L, 53L, 62L, 51L, 49L, 38L), lipid = c(50L,
33L, 16L, 31L, 21L, 28L, 26L, 14L, 24L, 21L, 19L, 31L, 20L, 16L,
24L, 16L, 29L, 41L, 24L, 22L, 26L)), class = "data.frame", row.names = c(NA,
-21L))

## Imagine that we are dividing the Practices into strata before randomizing
ym$assess_strata <- cut(ym$assess, c(0, 19.5, 34, 50))

## Assign a binary treatment rather than a three category treatment like the
## actual ASSIST Trial
set.seed(12345)
ym$trt <- block_ra(blocks=ym$assess_strata)

## Expand data to individual level for ease of analysis.
ym_long <- ym %>%
  mutate(ids = map(n, seq_len)) %>%
  unnest(cols = c(ids))

## Test to make sure that the expansion worked:
test1 <- ym_long %>% group_by(Practice) %>% summarize(n()==unique(n))
stopifnot(all(test1[,2]))

make_trt <- function(thevar){
    thelen <- length(thevar)
    prop_1 <- unique(thevar)
    num_1 <- round(thelen*prop_1/100)
    num_0 <- thelen - num_1
    rep(c(1,0),c(num_1,num_0))
}

ym_long <- ym_long %>% group_by(Practice) %>% mutate(assessed_bin=make_trt(assessed),
    aspirin_bin = make_trt(aspirin),
    hypo_bin = make_trt(hypo),
    lipid_bin = make_trt(lipid)
)

test2a <- ym_long %>% group_by(Practice) %>% summarize(prop_assessed=round(100*mean(assessed_bin)))
test2b <- left_join(ym,test2a) %>% select(Practice,assessed,prop_assessed)
stopifnot(all.equal(test2b$assessed,test2b$prop_assessed))

with(ym_long,table(Practice,trt,exclude=c()))

ym_long$assessed <- NULL
ym_long$aspirin <- NULL
ym_long$hypo <- NULL
ym_long$lipid <- NULL

ym_long <- rename(ym_long,n_practice = n, assessed = assessed_bin
,aspirin = aspirin_bin,hypo = hypo_bin,lipid = lipid_bin)

head(ym_long)

balanceTest(trt~n_practice+assessed+hypo+lipid+aspirin+strata(assess_strata)+cluster(Practice),
    data=ym_long)


save(ym_long,file="ym_long.rda")




