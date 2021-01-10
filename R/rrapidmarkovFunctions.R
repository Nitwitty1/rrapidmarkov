#rrapidMarkov
# by Ed Wilson
# July 2020
#
# Quick and simple Markov modelling

#' Calculate Markov trace
#'
#' Returns Markov trace as matrix of proportion of cohort in each state at each
#'   time point.
#' @param transitionsMatrix 2d matrix or 3d array of transition probabilities.
#' @param startingStates Vector of starting states.
#' @param cycles Number of cycles to run (only need to specify if
#'   stationary matrix)
#' @param stateNames Vector of state names (optional)
#' @param quietly Boolean, if TRUE will print summary data about the model
#' @return n*s matrix containing proportion of cohort in each state each
#'   cycle, where n=number of cycles, s=number of states.
#' @export
#' @examples
#' # Three state stationary model (i.e. probabilities don't vary with time),
#' # ten cycles
#'
#' states <- c("alive", "progressive", "dead")
#' startingStates <- c(1, 0, 0)
#' cycles <- 10
#' transitionsMatrix <- matrix(data = c(0.8, 0.1, 0.1,
#'                                      0.0, 0.7, 0.3,
#'                                      0.0, 0.0, 1.0),
#'                                      ncol = 3, byrow = TRUE)
#' runmarkov(transitionsMatrix, startingStates, cycles, states)
#'
#'
#' # Three-state non-stationary model (probabilities vary with time), six
#' # cycles.
#'
#' states <- c("alive", "progressive", "dead")
#' startingStates <- c(1, 0, 0)
#' transitionsMatrix <- aperm(array(data = c(0.8, 0.1, 0.1,
#'                                           0.0, 0.7, 0.3,
#'                                           0.0, 0.0, 1.0,
#'
#'                                           0.7, 0.1, 0.2,
#'                                           0.0, 0.6, 0.4,
#'                                           0.0, 0.0, 1.0,
#'
#'                                           0.6, 0.1, 0.3,
#'                                           0.0, 0.5, 0.5,
#'                                           0.0, 0.0, 1.0,
#'
#'                                           0.5, 0.1, 0.4,
#'                                           0.0, 0.4, 0.6,
#'                                           0.0, 0.0, 1.0,
#'
#'                                           0.4, 0.1, 0.5,
#'                                           0.0, 0.3, 0.7,
#'                                           0.0, 0.0, 1.0),
#'                                  dim = c(3, 3, 5)),
#'                            c(2, 1, 3))
#' runmarkov(transitionsMatrix, startingStates, stateNames = states)
#'
#'
#' # Two state, four cycles, non-stationary model, state names not specified,
#' # 50/50 split between states at t=0
#'
#' startingStates <- c(0.5, 0.5)
#' transitionsMatrix <- aperm(array(data = c(0.8, 0.2,
#'                                           0.3, 0.7,
#'
#'                                           1, 0,
#'                                           0, 1,
#'
#'                                           0, 1,
#'                                           0.5, 0.5),
#'                                  dim = c(2, 2, 3)),
#'                             c(2, 1, 3))
#' runmarkov(transitionsMatrix, startingStates)
runmarkov <- function(transitionsMatrix, startingStates, cycles = NULL,
                      stateNames = NULL, quietly = FALSE) {

  #if stationary matrix, duplicate to number of cycles
  stationary <- FALSE
  if (length(dim(transitionsMatrix)) == 2) {
    stationary <- TRUE
    transitionsMatrix <- array(rep(transitionsMatrix, 2),
                               c(dim(transitionsMatrix)[1],
                                 dim(transitionsMatrix)[2], cycles))
  } else {
    cycles <- dim(transitionsMatrix)[3] + 1
  }

  #error traps
  if (dim(transitionsMatrix)[1] != dim(transitionsMatrix)[2]){
    stop("Not a square matrix - number of rows must equal number of columns")
  }
  if(abs(sum(apply(transitionsMatrix[,,], 1, sum))/
         (dim(transitionsMatrix)[2] * dim(transitionsMatrix)[3]) - 1)
            > 0.000001) {
    #Note this trap allows for a little bit of tolerance as rounding errors
    #occur when using gtools::rdirichlet()
    stop("Rows don't all sum to 1.\n
          Is your matrix the right orientation? (rows = from, columns = to)?\n
          Have you got a miscalculation in your probabilities?")
  }
  if(length(startingStates) != dim(transitionsMatrix)[1]) {
    stop("Starting states mis-specified.  Must have one value for each state")
  }
  if(sum(startingStates) != 1) {
    stop("Starting states do not sum to 1.")
  }

  nStates <- dim(transitionsMatrix)[1]
  if (is.null(stateNames)) {
    for (i in 1:nStates) {
      stateNames <- c(stateNames, paste0("state", i))
    }
  }

  if (quietly == FALSE) {
    cat("\nrrapidmarkov by Ed Wilson\nv0.1 July 2020\n")
    cat(paste("Number of states: ", nStates, "\n"))
    cat("State names: ")
    cat(stateNames)
    cat("\n")
    cat(paste("Cycles: ", cycles, "\n"))
    cat("State transitions matrix:\n")
    ifelse (stationary == TRUE, print(transitionsMatrix[ , , 1]),
                                print(transitionsMatrix))
    cat("\n")
  }
  results <- matrix(data = NA, ncol = nStates, nrow = cycles,
                    dimnames = list(paste0("cycle", 0:(cycles - 1)), stateNames))

  results[1, ] <- startingStates
   for (i in 1:(cycles - 1)) {
     results[i+1, ] <- (t(transitionsMatrix[,,i]) %*% results[i, ])
   }
  return(results)
}

#' Format state-transition matrices
#'
#' Helper function to format a 2d matrix of sampled transition probabilities
#'   into a list of 3d arrays suitable for input into runmarkov and psamarkov
#' @param states Either integer specifying number of health states or a vector
#'   of state names
#' @param probs Matrix of sampled probabilities.
#' @param cycles Number of cycles to run (only required if stationary matrix)
#' @return List of length sims, where $simNN = array of size nStates x nStates
#'   x cycles, containing state transition matrices, where nStates = number of
#'   health states.
#' @export
#' @examples
#' #Stationary model
#'
#' library(gtools)
#' sims <- 10
#' set.seed(34)
#' probs <- cbind(rdirichlet(sims,c(80, 10, 10)),
#'               rdirichlet(sims,c(0, 7, 3)),
#'               rdirichlet(sims,c(0, 0, 1)))
#' colnames(probs) <- c("S1>S1", "S1>S2", "S1>S3",
#'                      "S2>S1", "S2>S2", "S2>S3",
#'                      "S3>S1", "S3>S2", "S3>S3")
#' states <- c("alive", "progressive", "dead")
#' cycles <- 5
#' converttomatrices(states, probs, cycles)
#'
#' # without specifying state names:
#' states <- 3
#' converttomatrices(states, probs, cycles)
#'
#' #Non-stationary model
#'
#' library(gtools)
#' sims <- 10
#' set.seed(34)
#' probs = cbind(rdirichlet(sims,c(80, 10, 10)),
#'               rdirichlet(sims,c(0, 7, 3)),
#'               rdirichlet(sims,c(0, 0, 1)),
#'
#'               rdirichlet(sims,c(70, 20, 10)),
#'               rdirichlet(sims,c(0, 7, 3)),
#'               rdirichlet(sims,c(0, 0, 1)))
#' colnames(probs) <- c("c1_S1>S1", "c1_S1>S2", "c1_S1>S3",
#'                      "c1_S2>S1", "c1_S2>S2", "c1_S2>S3",
#'                      "c1_S3>S1", "c1_S3>S2", "c1_S3>S3",
#'
#'                      "c2_S1>S1", "c2_S1>S2", "c2_S1>S3",
#'                      "c2_S2>S1", "c2_S2>S2", "c2_S2>S3",
#'                      "c2_S3>S1", "c2_S3>S2", "c2_S3>S3")
#' states <- c("alive", "progressive", "dead")
#' converttomatrices(states, probs)
converttomatrices <- function(states, probs, cycles=NULL) {
  # work out state names and number of states
  if (typeof(states) != "character" && is.numeric(states)) {
    nStates <- states
    stateNames <- NULL
    for (i in 1:nStates) {
      stateNames <- c(stateNames, paste0("state", i))
    }
  } else {
    nStates <- length(states)
    stateNames <- states
  }

  # work out number of cycles
  if (is.null(cycles)) cycles <- ncol(probs) / (nStates^2) + 1

  #if a stationary model then expand out to the number of cycles
  if (ncol(probs) == nStates^2) {
    probs <- matrix(data = rep(probs, cycles - 1), nrow=nrow(probs))
  }

  output <- vector(mode = "list", length = nrow(probs))
  names(output) <- paste0("sim", 1:nrow(probs))
  pb <- utils::txtProgressBar(min = 0, max = nrow(probs), style = 3)
  for (i in 1:nrow(probs)) {
    utils::setTxtProgressBar(pb, i)
    output[[i]] <- aperm(array(data = probs[i,],
                               dim=c(nStates, nStates, cycles - 1),
                               dimnames = list(stateNames,
                                               stateNames,
                                               paste0("cycle", 1:(cycles - 1)))),
                         c(2, 1, 3))
      }
  close(pb)
  return(output)
}

#' Probabilistic analysis of Markov trace
#'
#' Calculates proportions of cohort in each health state at each time point for
#' multiple Monte Carlo simulations.
#' @param transitionsMatrix Output of converttomatrices()
#' @param startingStates Vector of starting states.  Must be of length n
#'   where n = number of states and sum(startingStates) = 1.
#' @param cycles Number of cycles to run (only required if stationary matrix)
#' @param stateNames Vector of state names (optional)
#' @param cores Number of processor cores to use
#' @return A list comprising:
#'
#'  $means - matrix of size cycles x states reporting mean proportions of
#'  cohort in each state, each cycle.
#'
#'  $LCL - matrix as above reporting Lower 95\% Credibility/Confidence Interval
#'
#'  $UCL - matrix as above reporting Upper 95\% Credibility/Confidence Interval
#'
#'  $raw - array of size cycles x states x sims, reporting raw data from
#'  each simulation
#' @export
#' @examples
#' #Stationary model
#' # Generate PSA samples
#' library(gtools)
#' sims <- 10
#' set.seed(34)
#' probs = cbind(rdirichlet(sims, c(80, 10, 10)),
#'               rdirichlet(sims, c(0, 7, 3)),
#'               rdirichlet(sims, c(0, 0, 1)))
#' colnames(probs) <- c("S1>S1", "S1>S2", "S1>S3",
#'                      "S2>S1", "S2>S2", "S2>S3",
#'                      "S3>S1", "S3>S2", "S3>S3")
#'
#' # set up Markov model
#' stateNames <- c("alive", "progressive", "dead")
#' startingStates <- c(1, 0, 0)
#' cycles <- 20
#' transitionsMatrix <- converttomatrices(stateNames, probs, cycles)
#' psamarkov(transitionsMatrix, startingStates, cycles, stateNames)
#'
#' #Non-stationary model
#' #Note: See userguide vignette for detailed explanation
#' #Assumptions:
#' #a) Probability of moving to progressive health state increases with age
#' #   following Weibull function.
#' #b) Probability of death from 'alive' state increases with age, using
#' #   life tables.
#' #c) Starting age of modelled cohort is 60.
#' #d) Probability of death from 'progressive' state is approximately 30% (time-
#' #   invariant).
#'
#' stateNames <- c("alive", "progressive", "dead")
#' cycles <- 5
#'
#' sims <- 10
#' set.seed(100)
#' WeibullSamples <- data.frame(shapes = rlnorm(sims, log(2), log(1.05)),
#'                              scales = rlnorm(sims, log(0.001), log(1.1)))
#' startAge <- 60
#' p_AtoP <- t(apply(WeibullSamples, 1,
#'                   function(y) sapply((startAge + 1):(startAge + cycles - 1),
#'                                      function(x) {
#'                                        1 - exp(y[2] * ((x - 1)^y[1] - x^y[1]))
#'                                      }
#'                                     )
#'                  )
#'            )
#' colnames(p_AtoP) <- paste0("p_AtoP_Age",
#'                            (startAge + 1):(startAge + cycles - 1))
#'
#' lifeTable <- c(0.005104, 0.005600, 0.006303, 0.006832)
#'
#' nStates <- length(stateNames)
#' probs <- matrix(data = NA, nrow = sims, ncol = nStates^2 * (cycles - 1))
#' colnames(probs) <- unlist(lapply(1:(cycles - 1),
#'                       function(i) c(paste0("A>>A", i), paste0("A>>P", i),
#'                                     paste0("A>>D", i), paste0("P>>A", i),
#'                                     paste0("P>>P", i), paste0("P>>D", i),
#'                                     paste0("D>>A", i), paste0("D>>P", i),
#'                                     paste0("D>>D", i))))
#'
#' StayingAlive <- 1 - p_AtoP - matrix(rep(lifeTable, sims),
#'                                     nrow = sims, byrow = TRUE)
#' colnames(StayingAlive) <- paste0("p_AtoA_Age",
#'                                  startAge:(startAge + cycles - 2))
#'
#' p_AtoD <- matrix(rep(lifeTable, sims), nrow = sims, byrow = TRUE)
#' colnames(p_AtoD) <- paste0("p_AtoD_Age", startAge:(startAge + cycles - 2))
#'
#' probs[,seq(1, (cycles - 1) * nStates^2, nStates^2)] <- StayingAlive
#' probs[,seq(2, (cycles - 1) * nStates^2, nStates^2)] <- p_AtoP
#' probs[,seq(3, (cycles - 1) * nStates^2, nStates^2)] <- p_AtoD
#'
#' library(gtools)
#' FromProgressive <- rdirichlet(sims,c(0, 7, 3))
#' for (i in 1:3) {
#'   probs[, seq(i + 3, (cycles - 1) * nStates^2, nStates^2)] <-
#'                                                      FromProgressive[, i]
#' }
#' FromDead <- matrix(data = rep(c(0, 0, 1), sims), nrow = sims, byrow = TRUE)
#' for (i in 1:3) {
#'   probs[, seq(i + 6, (cycles - 1) * nStates^2, nStates^2)] <- FromDead[,i]
#' }
#'
#' # convert to format for PSAmarkov()
#' transitionsMatrix <- converttomatrices(stateNames,probs,cycles)
#'
#' # define starting states and run the analysis
#' startingStates <- c(1,0,0)
#' psamarkov(transitionsMatrix,startingStates,stateNames = stateNames)
psamarkov <- function (transitionsMatrix, startingStates,
                       cycles = NULL, stateNames = NULL,
                       cores = 1){
  sims <- length(transitionsMatrix)

  #duplicate startingStates, cycles and states into list items
   x <- c()
   for (i in 1:sims) {
     x[[i]] <- list(transitionsMatrix = transitionsMatrix[[i]],
                    startingStates = startingStates,
                    cycles = cycles, stateNames = stateNames)
   }

  #run the thing
  #Only set up cluster if cores > 1 (avoids unnecessary overhead)
  if (cores > 1) {
    cores <- min(cores,parallel::detectCores())
    if(.Platform$OS.type == "windows") {
      type = "PSOCK"
    } else {
      type = "FORK"
    }
    cl <- parallel::makeCluster(parallel::detectCores(), type)
    parallel::clusterEvalQ(cl, library(rrapidmarkov))
    output <- parallel::parSapply(cl, x, function(x) {
      runmarkov(x[]$transitionsMatrix, x[]$startingStates,
                x[]$cycles, x[]$stateNames, quietly = TRUE)
    }, simplify = "array")
    parallel::stopCluster(cl)
  } else {
    pb <- utils::txtProgressBar(min = 0, max = sims, style = 3)
    counter <- 0
    output <- sapply(x, function(x) {
      counter <<- counter + 1
      utils::setTxtProgressBar(pb, counter)
      runmarkov(x[]$transitionsMatrix, x[]$startingStates,
                x[]$cycles, x[]$stateNames, quietly = TRUE)
    }, simplify = "array")
    close(pb)
  }
  outputRaw <- output
  outputMean <- apply(output, c(1, 2), mean)
  outputLCL  <- apply(output, c(1, 2), stats::quantile, 0.025)
  outputUCL  <- apply(output, c(1, 2), stats::quantile, 0.975)

  return(list(means = outputMean, LCL = outputLCL, UCL = outputUCL,
              raw = outputRaw))
}

#' Calculate values (eg costs and QALYs)
#'
#' Multiplies output of runmarkov() by state costs, utilities or other value
#' @param markovTrace Output of runmarkov()
#' @param stateValues Vector of state values (eg costs or health state
#'   utilities).  Must be of length = number of states.
#' @param transPeriod Model transition period.  Must be entered as a fraction
#'   of a year.  Eg day = 1/365, week = 1/52, calendar month = 1/12, lunar
#'   (4-week) month = 1/13.
#' @param discountRate Annual discount rate.
#' @param QALYs Set to true when calculating QALYs from health state utilities
#'   to adjust for the transition period (eg if transitionPeriod = 1/2, then 1
#'   cycle in full health = 0.5 QALYs)
#' @return A list comprising:
#'
#'   $undisc - matrix of size (cycles+1) x (states+1) reporting undiscounted
#'   values by state and cycle.  Final column = row totals, final row = column
#'   totals.
#'
#'   $disc - matrix as per $undisc, discounted at specified rate per annum.
#'
#'   $totalUndisc - total undiscounted value (= bottom right value of $undisc)
#'
#'   $totalDisc - total discounted value (= bottom right value of $disc)
#'
#'   $transPeriod - transition period used in model in years
#'
#'   $discountRate - annual discount rate used.
#' @export
#' @examples
#' states <- c("alive", "progressive", "dead")
#' startingStates <- c(1, 0, 0)
#' cycles <- 10
#' transitionsMatrix <- matrix(data = c(0.8, 0.1, 0.1,
#'                                      0.0, 0.7, 0.3,
#'                                      0.0, 0.0, 1.0),
#'                                      ncol = 3, byrow = TRUE)
#'
#' output <- runmarkov(transitionsMatrix, startingStates, cycles, states)
#' stateCosts <- c(100, 1000, 0)
#' stateUtilities <- c(0.8, 0.5, 0)
#'
#' calcvalues(output, stateCosts)
#' calcvalues(output, stateUtilities, QALYs = TRUE)
calcvalues <- function(markovTrace, stateValues,
                       transPeriod = 1, discountRate = 0.035,
                       QALYs = FALSE) {
  pb <- utils::txtProgressBar(min = 0, max = 4, style = 3)
  utils::setTxtProgressBar(pb, 0)

  cycles <- nrow(markovTrace)

  df <-  1 / (1 + discountRate)^(transPeriod * 0:(cycles - 1))
  df[1:(1 / transPeriod)] <- 1

  if (QALYs == FALSE) {
    undisc <- sweep(markovTrace, 2, stateValues, FUN="*")
  } else {
    undisc <- sweep(markovTrace, 2, stateValues, FUN="*") * transPeriod
  }
  utils::setTxtProgressBar(pb, 1)

  #column sums
  undisc <- cbind(undisc, apply(undisc, 1, sum))
  disc <- undisc * df
  utils::setTxtProgressBar(pb, 2)

  #row sums
  undisc <- rbind(undisc, apply(undisc, 2, sum))
  disc <- rbind(disc, apply(disc, 2, sum))
  utils::setTxtProgressBar(pb, 3)

  #colnames & rownames
  colnames(undisc) <- colnames(disc) <- c(colnames(markovTrace), "total")
  rownames(undisc) <- rownames(disc) <- c(paste0("cycle", 0:(cycles - 1)),
                                                 "total")

  utils::setTxtProgressBar(pb, 4)
  close(pb)

  return(list(undisc = undisc, disc = disc,
              totalUndisc = as.vector(undisc[nrow(undisc), ncol(undisc)]),
              totalDisc = as.vector(disc[nrow(disc), ncol(disc)]),
              transPeriod = transPeriod, discountRate = discountRate
             )
         )
}

#' Calculate expected values (eg costs and QALYs) from PSA output
#'
#' Multiplies output of PSAMarkov() by state costs, utilities or other value
#' @param markovTrace Output of psamarkov()
#' @param stateValues Matrix of state values (eg costs or health state
#'   utilities).  nCols must equal number of states, nrows must equal number of
#'   simulations
#' @param transPeriod Model transition period.  Must be entered as a fraction
#'   of a year.  Eg day = 1/365, week = 1/52, calendar month = 1/12, lunar
#'   (4-week) month = 1/13.
#' @param discountRate Annual discount rate.
#' @param QALYs Set to true when calculating QALYs from health state utilities
#'   to adjust for the transition period (eg if transitionPeriod = 1/2, then 1
#'   cycle in full health = 0.5QALYs)
#' @return A list comprising:
#'
#'   $undisc$mean - matrix of size (cycles+1) x (states+1) reporting
#'   undiscounted mean values by state and cycle.  Final column = row totals,
#'   final row = column totals.
#'
#'   $undisc$LCL - matrix as above reporting Lower 95\% Credibility/Confidence
#'     Interval
#'
#'   $undisc$UCL - matrix as above reporting Upper 95\% Credibility/Confidence
#'     Interval
#'
#'   $undisc$raw - array of size (cycles+1) x (states+1) x sims, reporting raw
#'     undiscounted values from each simulation
#'
#'   $undisc$total - matrix of size 1x3, reporting undiscounted mean, LCL and
#'     UCL values (=bottom right cell of $undisc$mean, $LCL and $UCL
#'     respectively)
#'
#'   $undisc$SAVI - vector of total values from each simulation (=bottom right
#'     cell of each slice of $undisc$raw), in a format which can easily be
#'     input into the
#'     \href{https://github.com/Sheffield-Accelerated-VoI/SAVI-package}{SAVI}
#'     package.
#'
#'   $disc - sub-lists with the same structure as $undisc, but discounted at
#'     discountRate per annum.
#' @export
#' @examples
#' library(gtools)
#' sims <- 10
#' set.seed(34)
#' probs = cbind(rdirichlet(sims, c(80, 10, 10)),
#'               rdirichlet(sims, c(0, 7, 3)),
#'               rdirichlet(sims, c(0, 0, 1)))
#' colnames(probs) <- c("S1>S1", "S1>S2", "S1>S3",
#'                      "S2>S1", "S2>S2", "S2>S3",
#'                      "S3>S1", "S3>S2", "S3>S3")
#' states <- c("alive", "progressive", "dead")
#' startingStates <- c(1, 0, 0)
#' cycles <- 10
#' transitionsMatrix <- converttomatrices(states, probs, cycles)
#' output <- psamarkov(transitionsMatrix, startingStates, cycles, states)
#'
#' stateCosts <- cbind(rgamma(sims, 10, 1/50),
#'                     rgamma(sims, 1, 1/1000),
#'                     rep(0, sims))
#' stateUtilities <- cbind(rbeta(sims, 800, 200),
#'                         rbeta(sims, 50, 50),
#'                         rep(0, sims))
#'
#' psavalues(output, stateCosts)
#' psavalues(output, stateUtilities, QALYs = TRUE)
psavalues <- function(markovTrace, stateValues,
                      transPeriod = 1, discountRate = 0.035,
                      QALYs = FALSE) {
  pb <- utils::txtProgressBar(min = 0, max = 6, style = 3)
  utils::setTxtProgressBar(pb, 0)

  cycles <- nrow(markovTrace$means)

  df <-  1 / (1 + discountRate)^(transPeriod * 0:(cycles - 1))
  df[1:(1 / transPeriod)] <- 1

  if (QALYs == F) {
    undisc <- sweep(markovTrace$raw, c(2, 3), t(stateValues), FUN = "*")
  } else {
    undisc <- sweep(markovTrace$raw, c(2, 3), t(stateValues), FUN = "*") *
              transPeriod
  }
  utils::setTxtProgressBar(pb, 1)

  #column & row sums
  # resizes undisc and applies column sums, creates disc then sums the
  # rows of disc and undisc
  undisc2 <- array(data = 0, dim = dim(undisc) + c(1, 1, 0))
  undisc2[1:(dim(undisc)[1]), 1:(dim(undisc)[2]), ] <- undisc
  undisc2[ ,dim(undisc2)[2], ] <- apply(undisc2[ ,1:(dim(undisc2)[2]), ],
                                        c(1, 3), sum)
  undisc <- undisc2
  rm(undisc2)
  utils::setTxtProgressBar(pb, 2)
  disc <- undisc * c(df, 0)
  undisc[dim(undisc)[1], , ] <- apply(undisc, c(2, 3), sum)
  utils::setTxtProgressBar(pb, 3)
  disc[dim(disc)[1], , ] <- apply(disc, c(2, 3), sum)
  utils::setTxtProgressBar(pb, 4)

  #colnames and rownames
  colnames(undisc) <- colnames(disc) <- c(colnames(markovTrace$raw), "total")
  rownames(undisc) <- rownames(disc) <- c(paste0("cycle", 0:(cycles - 1)),
                                          "total")

  undiscMean <- apply(undisc, c(1,2), mean)
  undiscLCL  <- apply(undisc, c(1,2), stats::quantile, 0.025)
  undiscUCL  <- apply(undisc, c(1,2), stats::quantile, 0.975)
  discMean <- apply(disc, c(1,2), mean)
  discLCL  <- apply(disc, c(1,2), stats::quantile, 0.025)
  discUCL  <- apply(disc, c(1,2), stats::quantile, 0.975)
  utils::setTxtProgressBar(pb, 5)
  #structure output
  undisc <- list(mean = undiscMean, LCL = undiscLCL, UCL = undiscUCL,
                 raw = undisc,
                 total = matrix(data = c(
                   undiscMean[nrow(undiscMean), ncol(undiscMean)],
                   undiscLCL[nrow(undiscLCL), ncol(undiscLCL)],
                   undiscUCL[nrow(undiscUCL), ncol(undiscUCL)]),
                   ncol = 3, dimnames = list(NULL, c("total", "LCL", "UCL"))),
                 SAVI = undisc[nrow(undisc), ncol(undisc), ])

  disc   <- list(mean = discMean, LCL = discLCL, UCL = discUCL,
                 raw = disc,
                 total = matrix(data = c(
                   discMean[nrow(discMean), ncol(discMean)],
                   discLCL[nrow(discLCL), ncol(discLCL)],
                   discUCL[nrow(discUCL), ncol(discUCL)]),
                   ncol = 3, dimnames = list(NULL, c("total", "LCL", "UCL"))),
                 SAVI = disc[nrow(disc), ncol(disc), ])
  utils::setTxtProgressBar(pb, 6)
  close(pb)

  return(list(undisc = undisc, disc = disc,
              transPeriod = transPeriod, discountRate = discountRate))
}
