
# test code
rm(list=ls(all=T))
#source("R/rrapidMarkovFunctions.R")
library(rrapidmarkov)
library(gtools)


# Example 4a
ex4a <- function(sims,cycles,cores) {
  start.time <- Sys.time()
  set.seed(34)
  probs = cbind(rdirichlet(sims,c(80, 10, 10)),
                rdirichlet(sims,c(0, 7, 3)),
                rdirichlet(sims,c(0, 0, 1)))
  colnames(probs) <- c("S1>S1", "S1>S2", "S1>S3",
                       "S2>S1", "S2>S2", "S2>S3",
                       "S3>S1", "S3>S2", "S3>S3")
  stateCosts <- cbind(rgamma(sims, 10, 1/50),rgamma(sims, 1, 1/1000),
                      rep(0, sims))
  stateUtilities <- cbind(rbeta(sims, 800, 200), rbeta(sims, 50, 50),
                          rep(0, sims))
  states <- c("alive", "progressive", "dead")
  startingStates <- c(1, 0, 0)
  transitionsMatrix <- converttomatrices(states, probs, cycles)
  output <- psamarkov(transitionsMatrix, startingStates, cycles, states, cores)
  costs <- psavalues(output, stateCosts)
  QALYs <- psavalues(output, stateUtilities, QALYs = TRUE)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

# Example 4b
ex4b <- function(sims,cycles) {
  start.time <- Sys.time()
  set.seed(100)
  states <- c("alive", "progressive", "dead")
  startingStates <- c(1, 0, 0)
  startAge <- 60
  lifeTable <- c(0.003548, 0.000224, 0.000127, 0.000098, 0.000073,
                 0.000081, 0.000075, 0.000060, 0.000060, 0.000062,
                 0.000059, 0.000076, 0.000069, 0.000078, 0.000101,
                 0.000119, 0.000153, 0.000152, 0.000218, 0.000196,
                 0.000197, 0.000224, 0.000219, 0.000220, 0.000226,
                 0.000260, 0.000252, 0.000286, 0.000330, 0.000314,
                 0.000374, 0.000394, 0.000482, 0.000500, 0.000545,
                 0.000586, 0.000654, 0.000738, 0.000720, 0.000846,
                 0.000882, 0.000993, 0.001051, 0.001183, 0.001328,
                 0.001436, 0.001540, 0.001700, 0.001823, 0.001935,
                 0.002136, 0.002363, 0.002581, 0.002756, 0.002952,
                 0.003265, 0.003621, 0.003896, 0.004324, 0.004730,
                 0.005104, 0.005600, 0.006303, 0.006832, 0.007346,
                 0.007995, 0.008860, 0.009523, 0.010379, 0.011361,
                 0.012601, 0.013781, 0.015916, 0.017545, 0.019298,
                 0.022011, 0.025052, 0.027787, 0.031365, 0.034408,
                 0.038922, 0.043938, 0.049786, 0.057500, 0.065049,
                 0.073792, 0.084250, 0.095303, 0.108358, 0.121616,
                 0.136979, 0.153256, 0.169425, 0.187195, 0.206281,
                 0.230368, 0.249110, 0.270828, 0.290344, 0.316440, 0.339686)
  WeibullSamples<- data.frame(shapes = rlnorm(sims, log(2), log(1.05)),
                              scales = rlnorm(sims, log(0.001), log(1.1)))
  p_AtoP <- t(apply(WeibullSamples, 1,
                      function(y) sapply((startAge + 1):(startAge + cycles - 1),
                                         function(x) {
                                           1 - exp(y[2] * ((x - 1)^y[1] - x^y[1]))
                                         }
                      )
                    )
              )

  colnames(p_AtoP) <- paste0("p_AtoP_Age", startAge:(startAge + cycles - 2))
  lifetableExtract <- lifeTable[(startAge+1):(startAge + cycles - 1)]
  p_AtoD <- matrix(rep(lifetableExtract, sims), nrow = sims, byrow = TRUE)
  colnames(p_AtoD) <- paste0("p_AtoD_Age", startAge:(startAge + cycles - 2))
  StayingAlive <- 1 - p_AtoP - p_AtoD
  colnames(StayingAlive) <- paste0("p_AtoA_Age", startAge:(startAge + cycles - 2))
  nStates <- length(states)
  probs <- matrix(data = NA, nrow = sims, ncol = nStates^2 * (cycles - 1))
  colnames(probs) <- unlist(lapply(1:(cycles - 1),
                                   function(i) c(paste0("A>>A", i),
                                                 paste0("A>>P", i),
                                                 paste0("A>>D", i),
                                                 paste0("P>>A", i),
                                                 paste0("P>>P", i),
                                                 paste0("P>>D", i),
                                                 paste0("D>>A", i),
                                                 paste0("D>>P", i),
                                                 paste0("D>>D", i))
                                   )
                            )
  probs[, seq(1, (cycles-1) * nStates^2, nStates^2)] <- StayingAlive
  probs[, seq(2, (cycles-1) * nStates^2, nStates^2)] <- p_AtoP
  probs[, seq(3, (cycles-1) * nStates^2, nStates^2)] <- p_AtoD
  FromProgressive <- rdirichlet(sims, c(0, 7, 3))
  for (i in 1:3) {
    probs[, seq(i + 3, (cycles - 1) * nStates^2, nStates^2)] <- FromProgressive[, i]
  }
  FromDead <- matrix(data = rep(c(0, 0, 1), sims), nrow = sims, byrow = TRUE)
  for (i in 1:3) {
    probs[, seq(i + 6, (cycles - 1) * nStates^2, nStates^2)] <- FromDead[, i]
  }
  transitionsMatrix <- converttomatrices(states, probs)
  output <- psamarkov(transitionsMatrix, startingStates, stateNames = states)
  stateCosts <- cbind(rgamma(sims, 10, 1/50), rgamma(sims, 1, 1/1000), rep(0, sims))
  stateUtilities <- cbind(rbeta(sims, 800, 200), rbeta(sims, 50, 50), rep(0, sims))
  costs <- psavalues(output, stateCosts)
  QALYs <- psavalues(output, stateUtilities, QALYs = TRUE)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

ex4a(10000,10,1)
ex4a(10000,10,8)
ex4a(100000,10,1)
ex4a(100000,10,8)
ex4a(1000000,10,1)
ex4a(1000000,10,8)

ex4b(10000,10)
ex4b(100000,10)
ex4b(1000000,10)
ex4b(10000,20)
ex4b(100000,20)
ex4b(1000000,20)

