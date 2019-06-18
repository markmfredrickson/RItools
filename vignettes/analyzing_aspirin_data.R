# https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/RXQASA_H.htm 
# https://wwwn.cdc.gov/Nchs/Nhanes/Search/DataPage.aspx?Component=Questionnaire&CycleBeginYear=2013 
 
library(haven)
library(tidyverse)
aspirin <- read_xpt("./RXQASA_H.XPT")

summary(aspirin)
View(aspirin)
### Questions:
# RXQ515 - Followed (doctor's) advice, took low-dose aspirin?
# RXQ520 - Taking low-dose aspirin on your own?  

eq1 <- function(x) {
  tmp <- x == 1 
  tmp[is.na(tmp)] <- FALSE
  return(tmp)
}

aspirin$taking_aspirin <- eq1(aspirin$RXQ515) | eq1(aspirin$RXQ520)

## Blood pressure exam
# https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BPX_H.htm 

bp <- read_xpt("./BPX_H.XPT")

## four systoltic and diastolic measurements
# BPXSY1 - Systolic: Blood pres (1st rdg) mm Hg
# BPXDI1 - Diastolic: Blood pres (1st rdg) mm Hg
# BPAEN1 - Enhancement used first reading
# BPXSY2 - Systolic: Blood pres (2nd rdg) mm Hg
# BPXDI2 - Diastolic: Blood pres (2nd rdg) mm Hg
# BPAEN2 - Enhancement used second reading
# BPXSY3 - Systolic: Blood pres (3rd rdg) mm Hg
# BPXDI3 - Diastolic: Blood pres (3rd rdg) mm Hg
# BPAEN3 - Enhancement used third reading
# BPXSY4 - Systolic: Blood pres (4th rdg) mm Hg
# BPXDI4 - Diastolic: Blood pres (4th rdg) mm Hg 

bp$sys_mean <- rowMeans(bp[, c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")], na.rm = TRUE)
bp$dia_mean <- rowMeans(bp[, c("BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4")], na.rm = TRUE)

## now we need to match up the aspirin and bp tables
sum(aspirin$SEQN %in% bp$SEQN)

temp <- aspirin[aspirin$SEQN %in% bp$SEQN, ]

rownames(bp) <- bp$SEQN

combined <- cbind(temp, bp[as.character(temp$SEQN), ])

### TODO: load DEMO_H.XPT and merge it in
### https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Demographics&CycleBeginYear=2013


save.image(file = "analyzing_aspirin_data.rda")

demo = read_xpt("./DEMO_H.XPT")

data = aspirin %>% left_join(demo) 
#TO DO combine demo file, create balance plots. 
plot(balanceTest(taking_aspirin ~ ., data))
plot(balanceTest(taking_aspirin ~ RXQ515 + SIAPROXY + WTINT2YR, data))
plot(balanceTest(taking_aspirin ~ RXQ515 + SIALANG + DMDEDUC2 + FIALANG, data))
plot(balanceTest(taking_aspirin ~ MIAINTRP + SDDSRVYR, data))


# Meeting
# Matching: find groups with at least one treatment and at least one control s.t. e(x) similar
# e(x) can be estimated by logistic regression (glm)
f = glm(taking_aspirin ~RXQ515 + SIAPROXY + WTINT2YR, data, family = "binomial" )
summary(f)
# strata is adding weights to balance plots
# use roxygen2 to make vignette on aspirin data
# optmatch
install.packages("optmatch")
library(optmatch)
weights = fullmatch(f,data = data)
plot(balanceTest(taking_aspirin ~ RXQ515 + SIAPROXY + WTINT2YR + strata(weights), cbind(data,weights = weights)))

# for vignette, in devtools, use document() to turn into vignette