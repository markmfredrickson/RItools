data(nuclearplants)

xb <- xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n, 
               data = nuclearplants,
               strata = list(unstratified = NULL, pt = ~ as.factor(pt)),
               groups = list("Timings" = list("t1", "t2")))

print(xb)

print(xb, show.pvals = TRUE)

print(xb, horizontal = FALSE)

# Print a subest of the variables, statistics, and strata
print(xb, which.vars = c("date", "t1"), 
          which.stats = c("pr=0", "pr=1", "z", "p"),
          which.strata = "pt",
          which.groups = NULL)

# Only look at the pt strata and the "Timings" group
print(xb, which.groups = "Timings")

