## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load----------------------------------------------------------------
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

