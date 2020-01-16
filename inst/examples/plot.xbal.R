data(nuclearplants)
library(survival)

xb <- xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
               data = nuclearplants)

# Using the default grouping:
plot(xb, variable.labels = c(date = "Date",
             t1 = "Time 1",
             t2 = "Time 2",
             cap = "Capacity",
             ne = "In North East",
             ct = "Cooling Tower",
             bw = "Babcock-Wilcox",
             cum.n = "Total Plants Built"),
     strata.labels = c("--" = "Raw Data", "pt" = "Partial Turn-key"),
     absolute = TRUE)

# Using user supplied grouping
plot(xb, variable.labels = c(date = "Date",
             t1 = "Time 1",
             t2 = "Time 2",
             cap = "Capacity",
             ne = "In North East",
             ct = "Cooling Tower",
             bw = "Babcock-Wilcox",
             cum.n = "Total Plants Built"),
     strata.labels = c("--" = "Raw Data", "pt" = "Partial Turn-key"),
     absolute = TRUE,
     groups = c("Group A", "Group A", "Group A", "Group B",
                "Group B", "Group B", "Group A", "Group B"))
