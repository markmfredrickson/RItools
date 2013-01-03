data(nuclearplants)

xb <- xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
               data = nuclearplants,
               strata = list("none" = NULL,
                             "pt" = ~pt),
               groups = list("Temporal" = "t*"))

plot(xb, variable.labels = c(date = "Date",
                             t1 = "Time 1",
                             t2 = "Time 2",
                             cap = "Capacity",
                             ne = "In North East",
                             ct = "Cooling Tower",
                             bw = "Babcock-Wilcox",
                             cum.n = "Total Plants Built"),
    strata.labels = c("none" = "Raw Data", "pt" = "Partial Turn-key"),
    absolute = TRUE)

