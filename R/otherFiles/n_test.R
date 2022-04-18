library (rapidmarkov)

states <- c("alive", "progressive", "dead")
startingStates <- c(1, 0, 0)
cycles <- 10
transitionsMatrix <- matrix(data = c(0.8, 0.1, 0.1,
                                      0.0, 0.7, 0.3,
                                     0.0, 0.0, 1.0),
                                     ncol = 3, byrow = TRUE)
runmarkov(transitionsMatrix, startingStates, cycles, states)
