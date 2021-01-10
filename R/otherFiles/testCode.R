# test code
rm(list=ls(all=T))

source("R/rrapidMarkovFunctions.R")

# general note: cycles always numbered from zero.  For non-stationary models,
# the number of sets of probabilities needed is cycles-1, as the state of the
# Markov chain in cycle 0 is defined by startingStates parameter.
#
# Why one row for every simulation then convertToMatrices just changes it all
# to a list with one item for every simulation set?  Because it's nice and neat
#  to have all the data for a complete simulation on one row of a table.
#  If you prefer to generate lists of the matrices yourself do go ahead!
#
#  Terminology: stationary = non-time varying - same probs each cycle
#  non-stationary = vary with time, eg increasing death rate with age,
#  changing hazard over time.
#
#  Markov models should be analysed probabilistically:
#  This example is a deterministic analysis. Deterministic analyses of
# Markov models should be used during model development only for debugging
# and QA purposes.  Markov models should always be analysed probabilistically
# (see Wilson [MS in development]).

# ------------------------ Example 1 ------------------------------#
# This example calculates the cost and QALYs from a three state, temporally
# stationary model (i.e. probabilities don't vary with time) over ten years
# (cycle length = 1 year).
#
# startingStates must be of length n where n = number of states and
# sum(startingStates) = 1
# transitionsMatrix is specified as a 2d matrix, hence cycles parameter must be
# specified.

rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")
states <- c("alive","progressive","dead")
startingStates <- c(1,0,0)
cycles <- 10
transitionsMatrix <- matrix(data = c(0.8, 0.1, 0.1,
                                     0.0, 0.7, 0.3,
                                     0.0, 0.0, 1.0), ncol=3, byrow=TRUE)

output <- runMarkov(transitionsMatrix, startingStates, cycles, states)
print(output)

# apply costs and outcomes
stateCosts <- c(100,1000,0)
stateUtilities <- c(0.8, 0.5, 0)
transPeriod <- 1
discountRate <- 0.035

costs <- calcvalues(output,stateCosts,transPeriod,discountRate)
QALYs <- calcvalues(output,stateUtilities,transPeriod,discountRate, QALYs=T)

print(costs)
print(QALYs)

# transPeriod must be entered as a fraction of a year.  Eg day = 1/365,
# week = 1/52, calendar month = 1/12, lunar (4-week) month = 1/13,
# quarter = 1/4 and so on.  Enter as a fraction rather than decimal.
#
# A note on discounting:
# By convention the first year is not discounted.  For transition periods of
# less than one year (i.e. where transPeriod<1), this means the discount factor
# for the first 1/transPeriod cycles will be 1.  Subsequent cycles are
# discounted by 1/((1+r)^(transPeriod*t)), so for a discountRate of 3.5% and a
# transPeriod of 1/2 (i.e. 6 months), cycles 0 and 1 will not be discounted,
# cycle 2 has a discount factor of 1/((1+0.035)^(0.5*2)) = 0.966, cycle 3 has a
# discount factor of 1/((1+0.035)^(0.5*3)) = 0.950 and so on.

# ------------------------ Example 2a -----------------------------#
# This example calculates the cost and QALYs from a three-state non-stationary
# model (probabilities vary with time), over six cycles (cycle length = 1 year)
#
# Note the 'dim' argument of array() for defining transitionsMatrix must equal
# (nStates, nStates, cycles-1).
#
# aperm() transposes the array to the correct structure.

rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")
states <- c("alive","progressive","dead")
startingStates <- c(1,0,0)
transitionsMatrix <- aperm(array(data = c(0.8, 0.1, 0.1,
                                          0.0, 0.7, 0.3,
                                          0.0, 0.0, 1.0,

                                          0.7, 0.1, 0.2,
                                          0.0, 0.6, 0.4,
                                          0.0, 0.0, 1.0,

                                          0.6, 0.1, 0.3,
                                          0.0, 0.5, 0.5,
                                          0.0, 0.0, 1.0,

                                          0.5, 0.1, 0.4,
                                          0.0, 0.4, 0.6,
                                          0.0, 0.0, 1.0,

                                          0.4, 0.1, 0.5,
                                          0.0, 0.3, 0.7,
                                          0.0, 0.0, 1.0),
                                 dim = c(3,3,5)), c(2,1,3))

output <- runMarkov(transitionsMatrix, startingStates, stateNames = states)
print(output)

# Costs and QALYs can be applied by running the 'calcvalues' function as
# described in Example 1:
stateCosts <- c(100,1000,0)
stateUtilities <- c(0.8, 0.5, 0)
costs <- calcvalues(output, stateCosts)
QALYs <- calcvalues(output,stateUtilities, QALYs=T)
print(costs)
print(QALYs)

# ------------------------ Example 2b -----------------------------#
# Two state, four cycles, non-stationary model, state names not specified,
# 50/50 split between states at t=0

rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")
startingStates <- c(0.5,0.5)
transitionsMatrix <- aperm(array(data = c(0.8, 0.2,
                                          0.3, 0.7,

                                          1, 0,
                                          0, 1,

                                          0, 1,
                                          0.5, 0.5),
                                 dim = c(2,2,3)), c(2,1,3))

output <- runMarkov(transitionsMatrix, startingStates)
print(output)

# Costs and QALYs can be applied by running the 'calcvalues' function as
# described in Example 1:
stateCosts <- c(100,1000)
stateUtilities <- c(0.8, 0)
costs <- calcvalues(output,stateCosts)
QALYs <- calcvalues(output,stateUtilities,QALYs=T)
print(costs)
print(QALYs)

# ------------------------ Example 3 -------------------------------#
# Using the convertToMatrices() helper function
#
# convertToMatrices() converts the matrix to a list where each list item
# contains the state-transitions matrix for each Monte Carlo simulation.  This
# can then be used as the input to PSAMarkov().
#
# ------------------------ Example 3a ------------------------------#
#Stationary models
# probs must contain a set of sampled transition probabilities as a
# matrix of size n*c, where n = number of Monte Carlo simulations and
# c = nStates^2 (where nStates = number of health states).
#
# The column order for, say, a 2 state model must be:
# FromState1toState1, FromState1toState2, FromState2toState1,
# FromState2toState2.

# Three state, five cycle example:
# Transition probabilities are sampled from Dirichlet distributions.

rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")
library(gtools)
sims <- 10
set.seed(34)
probs <- cbind(rdirichlet(sims,c(80,10,10)),
               rdirichlet(sims,c(0,7,3)),
               rdirichlet(sims,c(0,0,1)))
colnames(probs) <- c("S1>S1","S1>S2","S1>S3",
                    "S2>S1","S2>S2","S2>S3",
                    "S3>S1","S3>S2","S3>S3")
states <- c("alive","progressive","dead")
cycles <- 5
convertToMatrices(states,probs, cycles)

# ------------------------ Example 3b ------------------------------#

#Non-stationary models
# probs must contain a set of sampled transition probabilities as a
# matrix of size n*c, where n = number of Monte Carlo simulations and
# c = (nStates^2)*(cycles-1) columns, (where nStates = number of health
# states, cycles = number of cycles).
#
# The column order for, say, a 2 state model with 3 cycles must be:
# Cycle1FromState1toState1, Cycle1FromState1toState2,
# Cycle1FromState2toState1, Cycle1FromState2toState2,
# Cycle2FromState1toState1, Cycle2FromState1toState2,
# Cycle2FromState2toState1, Cycle2FromState2toState2.

# Three-state, two cycle example:
# Transition probabilities are sampled from Dirichlet distributions, with
# different functions each cycle.

rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")
library(gtools)
sims <- 10
set.seed(34)
probs = cbind(rdirichlet(sims,c(80,10,10)),
             rdirichlet(sims,c(0,7,3)),
             rdirichlet(sims,c(0,0,1)),

             rdirichlet(sims,c(70,20,10)),
             rdirichlet(sims,c(0,7,3)),
             rdirichlet(sims,c(0,0,1)))

colnames(probs) <- c("c1_S1>S1","c1_S1>S2","c1_S1>S3",
                    "c1_S2>S1","c1_S2>S2","c1_S2>S3",
                    "c1_S3>S1","c1_S3>S2","c1_S3>S3",

                    "c2_S1>S1","c2_S1>S2","c2_S1>S3",
                    "c2_S2>S1","c2_S2>S2","c2_S2>S3",
                    "c2_S3>S1","c2_S3>S2","c2_S3>S3")
states <- c("alive","progressive","dead")
convertToMatrices(states, probs)

# ------------------------ Example 4a ------------------------------#
# Probabilistic modelling, stationary model
rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")
library(gtools)

# Generate PSA samples
sims <- 10000
set.seed(34)
probs = cbind(rdirichlet(sims,c(80,10,10)),
              rdirichlet(sims,c(0,7,3)),
              rdirichlet(sims,c(0,0,1)))
colnames(probs) <- c("S1>S1","S1>S2","S1>S3",
                     "S2>S1","S2>S2","S2>S3",
                     "S3>S1","S3>S2","S3>S3")
stateCosts <- cbind(rgamma(sims,10,1/50),rgamma(sims,1,1/1000),rep(0,sims))
stateUtilities <- cbind(rbeta(sims,800,200), rbeta(sims,50,50), rep(0,sims))

# set up Markov model
states <- c("alive","progressive","dead")
startingStates <- c(1,0,0)
cycles <- 20
transitionsMatrix <- convertToMatrices(states,probs,cycles)
transPeriod <- 1
discountRate <- 0.035

# run the analysis
output <- PSAMarkov(transitionsMatrix,startingStates,cycles,states)

#view some output
print(output$means)
print(output$LCL)
print(output$UCL)


# Apply costs and outcomes
costs <- PSAvalues(output,stateCosts,transPeriod,discountRate)
QALYs <- PSAvalues(output,stateUtilities,transPeriod,discountRate, QALYs=T)

print(costs$undisc$mean)

print(costs$disc$total)
print(QALYs$disc$total)

# (discounted) costs and QALYs in format suitable for SAVI():
print(costs$disc$SAVI)
print(QALYs$disc$SAVI)

# See notes for Example 1 for how to define transPeriod (must be as a fraction
# of a year) and on how discounting is handled when transPeriod < 1.

# ------------------------ Example 4b ------------------------------#
# Probabilistic modelling, non-stationary model
# The idea here is to use your favourite data-generating process to sample the
# time-varying transition probabilities.  These are compiled into a wide matrix
# called probs, where each row represents one complete set of simulated values.  The order of
# the columns is cycle number, state 1 to every other state, state 2 to every
# other state and so on.  For example for a 2-state model with 3 cycles, the
# complete order of columns would be:
# C1_S1>S1, C1_S1>S2, C1_S2>S1, C1_S2>S2
# C2_S1>S1, C2_S1>S2, C2_S2>S1, C2_S2>S2
# (note cycles are numbered from 0, where the Markov chain is in its starting
# state, so you only need to specify cycles-1 sets of transitions).
#
# Example: suppose we're modelling the natural history of Disease X' with the
# following assumptions:
# The 'typical' patient is a 60 year old female.
# 1. Patients in the 'alive' state are at no increased risk of death compared
# with the general population.
# 2. The probability of a patient developing progressive disease follows a
# Weibull distribution.
# 3. Patients in the 'progressive' state have approximately 30% mortality risk
# per annum (this is time-invariant: a bit of an oversimplification to keep
# the example simple)
#
# Data sources for (1) are national life tables and for (2) some long term
# observational data on a cohort of patients with Disease X.

rm(list=ls(all=T))
source("R/rrapidMarkovFunctions.R")

sims <- 10000
set.seed(100)
states <- c("alive","progressive","dead")
startingStates <- c(1,0,0)
startAge <- 60
cycles <- 20
transPeriod <- 1
discountRate <- 0.035

# Background death rate:
# source: 2016-18 lifetables, female, age 60-79
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/
# lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables
# these are population level data so the standard error is small enough to
# ignore and enter these as constant, rather than resampling in the PSA
lifeTable <- c(0.005104,0.005600,0.006303,0.006832,0.007346,
               0.007995,0.008860,0.009523,0.010379,0.011361,
               0.012601,0.013781,0.015916,0.017545,0.019298,
               0.022011,0.025052,0.027787,0.031365)

# Assume that fitting a survival function to the KM plot of the cohort data
# found a Weibull with shape = 2, scale = 0.001 provided the best fit, with SE
# of ln(shape) = ln(1.05), and SE of ln(scale) = ln(1.1)

# sample set of shapes and scales:
WeibullSamples<- data.frame(shapes = rlnorm(sims, log(2), log(1.05)),
                      scales = rlnorm(sims, log(0.001), log(1.1)))

# convert to transition probabilities ("Alive" to "Progressive"), covering the
# time horizon of the model (note the nested apply functions take a few seconds
# to process):
p_AtoP <- t(apply(WeibullSamples,1,
                  function(y) sapply((startAge+1):(startAge+cycles-1),
                                     function(x) {
                                       1-exp(y[2]*((x-1)^y[1] -x^y[1]))
                                     }
                                    )
                  )
            )
colnames(p_AtoP) <- paste0("p_AtoP_Age",(startAge+1):(startAge+cycles-1))

# now generate the rest of the transition probabilities, and store in 'probs'
nStates <- length(states)
probs <- matrix(data=NA, nrow=sims, ncol=nStates^2 * (cycles-1))
colnames(probs) <- unlist(lapply(1:(cycles-1),
                                 function(i) c(paste0("A>>A",i),
                                               paste0("A>>P",i),
                                               paste0("A>>D",i),
                                               paste0("P>>A",i),
                                               paste0("P>>P",i),
                                               paste0("P>>D",i),
                                               paste0("D>>A",i),
                                               paste0("D>>P",i),
                                               paste0("D>>D",i))
                                 )
                          )

StayingAlive <- 1 - p_AtoP - matrix(rep(lifeTable,sims), nrow=sims, byrow=T)
colnames(StayingAlive) <- paste0("p_AtoA_Age",startAge:(startAge+cycles-2))

p_AtoD <- matrix(rep(lifeTable,sims),nrow=sims, byrow=T)
colnames(p_AtoD) <- paste0("p_AtoD_Age",startAge:(startAge+cycles-2))

probs[,seq(1,(cycles-1)*nStates^2,nStates^2)] <- StayingAlive
probs[,seq(2,(cycles-1)*nStates^2,nStates^2)] <- p_AtoP
probs[,seq(3,(cycles-1)*nStates^2,nStates^2)] <- p_AtoD

FromProgressive <- gtools::rdirichlet(sims,c(0,7,3))
for (i in 1:3) {
  probs[,seq(i+3,(cycles-1)*nStates^2,nStates^2)] <- FromProgressive[,i]
}

FromDead <- matrix(data=rep(c(0,0,1),sims), nrow=sims, byrow=T)
for (i in 1:3) probs[,seq(i+6,(cycles-1)*nStates^2,nStates^2)] <- FromDead[,i]


# convert to format for PSAmarkov()
transitionsMatrix <- convertToMatrices(states,probs,cycles)

# run the analysis
output <- PSAMarkov(transitionsMatrix,startingStates,states = states)
print(output$means)
print(output$LCL)
print(output$UCL)

#sample costs and health state utilities
stateCosts <- cbind(rgamma(sims,10,1/50),rgamma(sims,1,1/1000),rep(0,sims))
stateUtilities <- cbind(rbeta(sims,800,200), rbeta(sims,50,50), rep(0,sims))

# Apply costs and outcomes
costs <- PSAvalues(output,stateCosts,transPeriod,discountRate)
QALYs <- PSAvalues(output,stateUtilities,transPeriod,discountRate, QALYs=T)

print(costs$undisc$mean)

print(costs$disc$total)
print(QALYs$disc$total)

# (discounted) costs and QALYs in format suitable for SAVI():
print(costs$disc$SAVI)
print(QALYs$disc$SAVI)

# ----------- Quality / accuracy checks and Speed/performance tests -----------#
#
# Comparisons with long-hand Excel
# speed / performance tests
