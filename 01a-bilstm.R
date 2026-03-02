# Clean R script with ReLU fix and debugging
library(reticulate)
#install.packages("epiworldRcalibrate")
library(epiworldRcalibrate)

#an example to see the package is working.
# Test prediction
incidence_vec <- c(
  103, 37, 60, 74, 108, 125, 138, 186, 215, 276, 318, 331, 414, 402, 446, 454, 405, 401, 373, 334, 285,
  241, 219, 156, 140, 108, 93, 82, 73, 78, 48, 38, 34, 22, 27, 22, 20, 11, 14, 8, 11, 14,
  8, 6, 6, 2, 0, 7, 2, 4, 4, 1, 1, 2, 2, 1, 0, 2, 1, 1, 0
)

n <- 7087
recov <- 0.203

result <- calibrate_sir(incidence_vec, n, recov)
#print(result)


