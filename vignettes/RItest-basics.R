## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, message = FALSE-----------------------------------------------
library(RItools)

## ----data----------------------------------------------------------------
data("nuclearplants")

head(nuclearplants)

## ----ttest---------------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$pr,
                     test.stat = t.mean.difference)

## ----qtest---------------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$pr, 
       test.stat = quantileDifference(.25))

## ----qtest2--------------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$pr, 
       test.stat = quantileDifference(.50))

## ----ks------------------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$pr, test.stat = ksTestStatistic)

## ----samples-------------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$pr, test.stat = t.mean.difference,
       samples = 100)

## ----setblock------------------------------------------------------------
region.blocks <- simpleRandomSampler(total = 32, z = nuclearplants$pr,
                                     b = nuclearplants$ne)

## ----testblock-----------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$pr, test.stat = t.mean.difference,
       sampler = region.blocks)

## ----fakedat-------------------------------------------------------------
set.seed(20180620)

nuclearplants$cl <- sample(letters[1:3], 32, replace = TRUE)

nuclearplants$trt <- with(nuclearplants,
                          ifelse(cl == sample(letters[1:3], 1), 1, 0))

## ----faketab-------------------------------------------------------------
with(nuclearplants, table(trt, cl))

## ----setcluster----------------------------------------------------------
sim.clusters <- clusterRandomSampler(clusters = nuclearplants$cl,
                                     z = nuclearplants$trt)

## ----cluster.test--------------------------------------------------------
RItest(nuclearplants$cost, nuclearplants$trt, test.stat = t.mean.difference,
       sampler = sim.clusters)

