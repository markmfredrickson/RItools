set.seed(20121204)

# generate some balance data
vars.groups <- rep(c("Group A", "Group B", "Group C", "Group D"), each = 5)
names(vars.groups) <- paste("Variable", letters[1:20])
  
nvars <- length(vars.groups)

balance_data <- matrix(c(rnorm(n = nvars, mean = 1, sd = 0.5), 
                         rnorm(n = nvars, mean = 0, sd = 0.5)),
                       ncol = 2)

colnames(balance_data) <- c("Before Adjustment", "After Matching")

rownames(balance_data) <- names(vars.groups)

balanceplot(balance_data, vars.groups, xlab = "Balalnce Before/After Matching")

# base R graphics are allowed

abline(v = colMeans(balance_data), lty = 3, col = "grey")

