library (readxl)
library (gtools)
library(rrapidmarkov)
library(tidyverse)
set.seed(2021)

"

The induvidual functions are called from a single sim function
which will bring the required values into the name space and update
those values for diiferent initial states, interventions, and testing
regimens. This relies on overloading some variable names which
may not constitute best practice, but is a pragmatic choice. Most functions
as written will work independently provided the required values are
in the namespace.

In general, I've used the global assignment parameter to
make data available outside of functions where required as
opposed to using return.

"

settings <- function (sims, cycles, age, prevalence){
  "Sets the base paramaters for the simulation as global"
  sims <<- sims
  cycles <<- cycles
  age <<- age
  prevalence <<- prevalence
}

init_state <- function(treatment){
  "
  Calculates the inital distribrution of patients in each NYHA class. However,
  startingStates does not accept a 2-d input. The only solution I could think of
  is to call calcvalues with a different starting state for each sim and append
  the results to an empty list. For now, i've run this function a few times and
  hard coded the mean which is what is returned. I've commented out the actual
  function below

  "

  "

  The if block pulls treats all patients receiving an implant indentically and
  does not discriminate between type of implant. This should be ammended. CRT-D
  is known to how more adverse events associated with implantation.

  "

  "

  In general, i've used if...elif blocks and without an escape else block to
  handle errorrs or exceptions. I should add these but I tested the functions as
  I wrote them and I can't imagine a scenario in which the programmed blocks receive
  anything but the expected input.

  "

  "

    if (treatment == 'icd' | treatment == 'crt-p' | treatment == 'crt-d'){
    #initial distribution of of nyha pataients
    dist_nyhaiii <-  .77
    dist_nyhaiv  <-  .23


    #implantation outcomes
    p1success <- rbeta(sims, 409, 60)
    p2success <- rbeta(sims, 62, 10 )
    p3success <- rbeta(sims, 8, 2   )
    pdeath    <- 0.002445

    #short-term tranisition between NYHA class
    x <- rdirichlet (sims, c(0.298 * dist_nyhaiii * 409, 0.459 * dist_nyhaiii * 409, 0.227 * dist_nyhaiii * 409,
                             0.016 * dist_nyhaiii * 409))

    y <- rdirichlet (sims, c(0.091 * dist_nyhaiv * 409, 0.455 * dist_nyhaiv * 409, 0.409 * dist_nyhaiv * 409,
                             0.045 * dist_nyhaiv * 409))

    stiiitoi   <-  x[, 1]
    stiiitoii  <-  x[, 2]
    stiiitoiii <-  x[, 3]
    stiiitoiv  <-  x[, 4]

    stivtoi   <- y[, 1]
    stivtoii  <- y[, 2]
    stivtoiii <- y[, 3]
    stivtoiv  <- y[, 4]

    firstsuccess  <- array(numeric(),c(sims))
    secondsuccess <- array(numeric(),c(sims))
    thirdsuccess  <- array(numeric(),c(sims))

    firstdeath  <- array(numeric(),c(sims))
    seconddeath <- array(numeric(),c(sims))
    thirddeath  <- array(numeric(),c(sims))

    nyhai   <- array(numeric(),c(sims))
    nyhaii  <- array(numeric(),c(sims))
    nyhaiii <- array(numeric(),c(sims))
    nyhaiv  <- array(numeric(),c(sims))

    for (i in 1:sims){
      firstsuccess[i] <- p1success[i]
      firstdeath[i] <- pdeath
      secondsuccess[i] <- (1 - (p1success[i] + pdeath)) * p2success[i]
      seconddeath[i] <- (1 - (p1success[i] + pdeath)) * pdeath
      thirdsuccess[i] <- (1 - (p1success[i] + pdeath)) * (1 - (p2success[i] + pdeath)) * p3success[i]
      thirddeath[i] <- (1 - (p1success[i] + pdeath)) * (1 - (p2success[i] + pdeath)) * pdeath

      nyhai[i] <- ((dist_nyhaiii * stiiitoi[i]) + (dist_nyhaiv * stivtoi[i])) * (
        firstsuccess[i] + secondsuccess[i] + thirdsuccess[i])
      nyhaii[i] <- ((dist_nyhaiii * stiiitoii[i]) + (dist_nyhaiv * stivtoii[i])) * (
        firstsuccess[i] + secondsuccess[i] + thirdsuccess[i])
      nyhaiii[i] <- ((dist_nyhaiii * stiiitoiii[i]) + (dist_nyhaiv * stivtoiii[i])) * (
        firstsuccess[i] + secondsuccess[i] + thirdsuccess[i])
      nyhaiv[i] <- ((dist_nyhaiii * stiiitoiv[i]) + (dist_nyhaiv * stivtoiv[i])) * (
        firstsuccess[i] + secondsuccess[i] + thirdsuccess[i])
    }
    return (list(mean(nyhai), mean(nyhaii) , mean(nyhaiii),mean(nyhaiv)))
  }
  else if (treatment == 'mt'){
    #initial distribution of of nyha pataients - mt
    dist_nyhaiii  <- .73
    dist_nyhaiv   <- .27

    #short-term tranisition between NYHA class - mt
    x <- rdirichlet(sims,c(0.103 * dist_nyhaiii * 404, 0.303 * dist_nyhaiii * 404, 0.528 * dist_nyhaiii * 404,
                           0.067 * dist_nyhaiii * 404))

    y <- rdirichlet(sims, c(0.0001 * dist_nyhaiv * 404, 0.200 * dist_nyhaiv * 404, 0.600 * dist_nyhaiv * 404,
                            0.200 * dist_nyhaiv * 404))

    stiiitoi   <-  x[, 1]
    stiiitoii  <-  x[, 2]
    stiiitoiii <-  x[, 3]
    stiiitoiv  <-  x[, 4]

    stivtoi   <- y[, 1]
    stivtoii  <- y[, 2]
    stivtoiii <- y[, 3]
    stivtoiv  <- y[, 4]

    nyhai   <- array(numeric(),c(sims))
    nyhaii  <- array(numeric(),c(sims))
    nyhaiii <- array(numeric(),c(sims))
    nyhaiv  <- array(numeric(),c(sims))

    for (i in 1:sims){
      nyhai[i]   <- ((dist_nyhaiii * stiiitoi[i]) + (dist_nyhaiv * stivtoi[i]))
      nyhaii[i]  <- ((dist_nyhaiii * stiiitoii[i]) + (dist_nyhaiv * stivtoii[i]))
      nyhaiii[i] <- ((dist_nyhaiii * stiiitoiii[i]) + (dist_nyhaiv * stivtoiii[i]))
      nyhaiv[i]  <- ((dist_nyhaiii * stiiitoiv[i]) + (dist_nyhaiv * stivtoiv[i]))
    }

    return (list(mean(nyhai), mean(nyhaii) , mean(nyhaiii),mean(nyhaiv)))
  }

  "

  if (treatment == 'icd' | treatment == 'crt-p' | treatment == 'crt-d'){
    return (list(0.248, 0.456,0.273,0.023))
  } else if (treatment == 'mt'){
    return (list(0.0757,0.275,0.548,0.1013))
  }

}

stable <- function(treatment){
  "
  Returns probabilities for tranistions between stable NYHA states (from nyha i-iv
  to nyha i-v). I'd assumed that ICD does not impact these transition probabilities
  in my diss and that is  reflected in the strucure of the if block. The variable names
  are lazy. p is from nyha i to nyha i-iv, q is from nyha ii to nyha i-iv and so on.

  "

  if (treatment == 'crt-d' | treatment == 'crt-p'){
    #initial distribution of of nyha pataients
    dist_nyhaiii <- 0.77
    dist_nyhaiv <- 0.23

    #long-term transition between NYHA class
    samplei <- (409 * dist_nyhaiii * 0.298) + (409 * dist_nyhaiv * 0.091)
    p <- rdirichlet (sims, c(0.906 * samplei, 0.075 * samplei, 0.016 * samplei, 0.003 * samplei))

    sampleii <- (409 * dist_nyhaiii * 0.459) + (409 * dist_nyhaiv * 0.455)
    q <- rdirichlet (sims, c(0.067 * sampleii, 0.896 * sampleii, 0.033 * sampleii, 0.004 * sampleii))

    sampleiii <- (409 * dist_nyhaiii * 0.227) + (409 * dist_nyhaiv * 0.409)
    r <- rdirichlet (sims, c(0.007 * sampleiii, 0.121 * sampleiii, 0.864 * sampleiii, 0.008 * sampleiii))

    sampleiv <- (409 * dist_nyhaiii * 0.016) + (409 * dist_nyhaiv * 0.045)
    s <- rdirichlet (sims, c(0.048 * sampleiv, 0.048 * sampleiv, 0.181 * sampleiv, 0.723 * sampleiv))
  } else{
    #initial distribution of of nyha pataients - mt
    dist_nyhaiii_m <- 0.73
    dist_nyhaiv_m <- 0.27

    #long-term transition between NYHA class - mt
    samplei <- (404 * dist_nyhaiii_m * 0.103) + (404 * dist_nyhaiv_m * 0.0001)
    p <- rdirichlet (sims, c(0.7956 * samplei, 0.1245 * samplei, 0.0738 * samplei, 0.0061 * samplei))

    sampleii <- (404 * dist_nyhaiii_m * 0.303) + (404 * dist_nyhaiv_m * 0.200)
    q <- rdirichlet (sims, c(0.0710 * sampleii, 0.8448 * sampleii, 0.0765 * sampleii, 0.0077 * sampleii))

    sampleiii <- (404 * dist_nyhaiii_m * 0.528) + (404 * dist_nyhaiv_m * 0.600)
    r <- rdirichlet (sims, c(0.0047 * sampleiii, 0.0893 * sampleiii, 0.8845 * sampleiii, 0.0216 * sampleiii))

    sampleiv <- (404 * dist_nyhaiii_m * 0.067) + (404 * dist_nyhaiv_m * 0.200)
    s <- rdirichlet (sims, c(0.00001 * sampleiv, 0.1064 * sampleiv, 0.1064 * sampleiv, 0.7872 * sampleiv))
  }
  return(list(p,q,r,s))
}

SCDWeib <- function(){
  "Weibull parameters and hazard ratios for SCD are set here. These are global."

  #sudden death
  scale_scd <<- exp(rnorm(sims, -5.361, 0.047))
  shape_scd <<- exp(rnorm(sims, -0.083, 0.009))

  hr_scd_icdvmt <<- exp(rnorm(sims, -1.109, 0.273))
  hr_scd_crtvmt <<- exp(rnorm(sims, -0.713, 0.253))

  hr_i <<- rep(1,sims)
  hr_scd_iivi <<- exp(rnorm(sims, -0.041, 0.329))
  hr_scd_iiivi <<- exp(rnorm(sims, -0.045, 0.341))
  hr_scd_ivvi <<- exp(rnorm(sims, -0.386, 0.65))
}

HospWeib <- function(){
  "Ditto for hospitalisations. There isn't a hazard ratio for ICD because I'd
  assumed these have no impact on hospitalisation. This is unlikely to be true."

  #hospitalisation
  scale_hosp <<- exp(rnorm(sims, -3.028, 0.072))
  shape_hosp <<- exp(rnorm(sims, -0.263, 0.044))

  hr_hosp_crtvmt <<- exp(rnorm(sims, -0.246, 0.13))

  hr_hosp_iivi <<- exp(rnorm(sims, 0.154, 0.189))
  hr_hosp_iiivi <<- exp(rnorm(sims, 0.154, 0.189))
  hr_hosp_ivvi <<- exp(rnorm(sims, 1.601, 0.264))

  p_hosp_die <<- rbeta(sims, 12, 162)
  p_hosp_die_mt <<- rbeta(sims, 39, 345)
}

Ht <- function(type = 'base'){
  "
  This function calcultes the cumulative hazard from the weibull parameters and
  then sets this as global. The optional argument scales the cumulative hazard to
  reflect a life time SCD prob. of 1. By default, the function DOES NOT scale the
  cumulative hazard.

  "
  "

  The way H(t) has been scaled means that the model probably should not be run
  for any less than a lifetime horizon as the value of H(t) at t will be 1 for
  any t.

  "

  #Cumulative hazard function
  #put the weibull parameters in memory, then calculations.
  SCDWeib()
  HospWeib()
  base_SCD <- array(numeric(),c(sims, cycles))
  scaled_SCD <- array(numeric(),c(sims, cycles))
  base_hosp <- array(numeric(),c(sims, cycles))
  for(n in 1:sims){
    for(c in 1:(cycles)){
      base_SCD[n,c] <- (scale_scd[n]*c)^shape_scd[n]
      base_hosp[n,c] <- (scale_hosp[n]*c)^shape_hosp[n]
    }
  }
  for(n in 1:sims){
    for(c in 1:(cycles)){
      scaled_SCD[n,c] <- base_SCD[n,c]/base_SCD[n,cycles]
    }
  }

  if (type == 'base'){
    Ht_SCD <<- base_SCD
    Ht_hosp <<- base_hosp
  } else if (type =='scaled'){
    Ht_SCD <<- scaled_SCD
    Ht_hosp <<- base_hosp
  }
}

failure <- function(hr, Ht){
  "A small convenience function to calculate number of failures from Ht"
  dead <- array (numeric(), c(sims, cycles) )
  for (n in 1:sims){
    for (c in 1:(cycles-1)){
      dead[n,c] <-  1-exp((hr[n]*Ht[n,c])-(hr[n]*Ht[n,(c+1)]))
    }
  }
  return (dead)
}

fail_SCD <- function (treatment, D = '++'){
  "
  Calculate the probability of fatal and non-fatal SCD events for each
  NYHA class. The 2 arguments, treatment and D(for dummy) specify the treatment
  the patient is on (mt, icd, crt-p or crt-d) and the classification of the patient
  as either true positive (++), false negative (+-), false positive (-+), and
  true negative (--).

  The function has no default value for treatment. The default for D is ++ as the
  logic for these patients (in code) is identical to that for calculating the
  same probabilities for patients when not considering the testing regime. The only
  differ in the Ht used.

  "
  #put all the NYHA hazard ratios togeather so I can use lapply
  hr_NYHAs <- list(hr_i, hr_scd_iivi, hr_scd_iiivi,  hr_scd_ivvi)

  if (D == '++'){

    #get the required hazard ratios for treatment together.
    if (treatment == 'crt-d'){
      hr_treatment <- hr_scd_icdvmt * hr_scd_crtvmt
    } else if (treatment == 'crt-p'){
      hr_treatment <- hr_scd_crtvmt
    } else if (treatment == 'icd'){
      hr_treatment <- hr_scd_icdvmt
    } else{
      hr_treatment <- hr_i#hr_i is an array of ones
    }

    "

    calculating the fatal events is straightforward (and fancy?). FSCD is an
    array of size (sims*cycles*4), i.e., sims*cycles results for each
    NYHA class.

    "
    FSCD <- lapply(hr_NYHAs, function(hr_NYHAs)failure(hr_NYHAs*hr_treatment, Ht_SCD))
    FSCD <- array(as.numeric(unlist(FSCD)), dim = c(sims,cycles,4))

    #Code looks very similar, but the failure function call has different args.
    NFSCD <- lapply(hr_NYHAs, function(hr_NYHAs)failure(hr_NYHAs, Ht_SCD))
    NFSCD <- (array(as.numeric(unlist(NFSCD)), dim = c(sims,cycles,4))) - FSCD

    "

    There are lots of these list to array conversions throughout. I don't know
    a better way. They are a bit wordy.

    "

  } else if (D == '+-'){

    #false negatives have received no treatment so their treatment hazard has been
    #set to 1. The rest is identical to ++.
    hr_treatment <- hr_i

    FSCD <- lapply(hr_NYHAs, function(hr_NYHAs)failure(hr_NYHAs*hr_treatment, Ht_SCD))
    FSCD <- array(as.numeric(unlist(FSCD)), dim = c(sims,cycles,4))

    NFSCD <- lapply(hr_NYHAs, function(hr_NYHAs)failure(hr_NYHAs, Ht_SCD))
    NFSCD <- (array(as.numeric(unlist(NFSCD)), dim = c(sims,cycles,4))) - FSCD

  } else if (D == '-+' | D == '--'){

    #Both false positives are true negatives have a 0 lifetime risk of SCD. This
    #has been set below. For sake of clarity, this should really be two blocks.

    FSCD <- array(as.numeric(0), dim = c(sims,cycles,4))
    NFSCD <- array(as.numeric(0), dim = c(sims,cycles,4))

  }

  return(list(NFSCD, FSCD))
}

fail_Hosp <- function (treatment){

  "
  Calculate the probability of fatal and non-fatal SCD events for each
  NYHA class. 1 arg to specify treatment allocation. No need to write seperate
  code blocks for ++. +-. -+, -- as the hospitalisation probabilities are
  identical, regardless of test classification

  "

  hr_NYHAh <- list(hr_i, hr_hosp_iivi, hr_hosp_iiivi, hr_hosp_ivvi)

  if (treatment == 'crt-d' | treatment == 'crt-p'){
    hr_treatment <- hr_hosp_crtvmt
    prob_of_death <- p_hosp_die
  } else{
    hr_treatment <- hr_i
    prob_of_death <- p_hosp_die_mt
  }

  FHosp <- lapply(hr_NYHAh, function(hr_NYHAh)(failure(hr_NYHAh*hr_treatment,
                                                       Ht_hosp)*(prob_of_death)))
  FHosp <- array(as.numeric(unlist(FHosp)), dim = c(sims,cycles,4))

  NFHosp <- lapply(hr_NYHAh, function(hr_NYHAh)(failure(hr_NYHAh*hr_treatment,
                                                        Ht_hosp)*(1-prob_of_death)))
  NFHosp <- array(as.numeric(unlist(NFHosp)), dim = c(sims,cycles,4))

  return (list(NFHosp, FHosp))
}

LT <- function(){
  "
  Fetch the values from the life table.  This function will not work untill the
  path has been changes (line 423). There should be a way to specify the path such
  that it is system agnostic but again, I've not been able to figure it out.
  Advice needed.

  "

  #Lifetable
  #Change path here
  all_cause_mortality <- read_excel ("C:\\Users\\Public\\nationallifetables3yearuk.xlsx", sheet="2017-2019")
  lx <- all_cause_mortality[['lx']]
  vec <- vector()
  for  (i in 1:length(lx)-1){
    some <- (lx[i] - lx[i + 1]) / lx[i]
    some <- 1 - (1-some)^(1/12)
    vec <- c(vec, some)
  }
  lifetable <- rep(vec[(age):length(vec)],each=12)
  return (lifetable)
}

make_probs <- function (treatment, D='++'){

  "Make the probs table. This function calls stable, fail_SCD, fail_Hosp and LT
  and puts it all togeather. "

  hold_stable <- array(as.numeric(unlist(stable(treatment))), dim=c(sims, 4, 4))
  SCD <- array(as.numeric(unlist(fail_SCD(treatment, D))), dim=c(sims, cycles, 4, 2))
  Hosp <- array(as.numeric(unlist(fail_Hosp(treatment))), dim=c(sims, cycles, 4, 2))

  p <- hold_stable[,,1]
  q <- hold_stable[,,2]
  r <- hold_stable[,,3]
  s <- hold_stable[,,4]
  NFSCD <- SCD[,,,1]
  FSCD <- SCD[,,,2]
  NFHosp <- Hosp[,,,1]
  FHosp <- Hosp[,,,2]
  lifetable <- LT()

  #using c bind in the loop has made it slow
  #would be better to assign values...
  #...to probs table using index numbers
  probs <- NULL
  for (time in 1:(cycles-1)){
    k1 <- as.vector(1 - (NFSCD[,time,1] + FSCD[,time,1] +  NFHosp[,time,1] +
                           FHosp[,time,1] + lifetable[time]))
    row1 <- cbind(k1*p, NFSCD[,time,1], FSCD[,time,1], NFHosp[,time,1],
                  FHosp[,time,1], lifetable[time])

    k2 <- as.vector(1 - (NFSCD[,time,2] + FSCD[,time,2] +  NFHosp[,time,2] +
                           FHosp[,time,2] + lifetable[time]))
    row2 <- cbind(k2*q, NFSCD[,time,2], FSCD[,time,2], NFHosp[,time,2],
                  FHosp[,time,2], lifetable[time])

    k3 <- as.vector(1 - (NFSCD[,time,3] + FSCD[,time,3] +  NFHosp[,time,3] +
                           FHosp[,time,3] + lifetable[time]))
    row3 <- cbind(k3*r, NFSCD[,time,3], FSCD[,time,3], NFHosp[,time,3],
                  FHosp[,time,3], lifetable[time])

    k4 <- as.vector(1 - (NFSCD[,time,4] + FSCD[,time,4] +  NFHosp[,time,4] +
                           FHosp[,time,4] + lifetable[time]))
    row4 <- cbind(k4*s, NFSCD[,time,4], FSCD[,time,4], NFHosp[,time,4],
                  FHosp[,time,4], lifetable[time])

    k5 <- as.vector(NFSCD[,time,1] + NFSCD[,time,2] + NFSCD[,time,3] +
                      NFSCD[,time,4])
    if (all(k5 == 0)){
      k5 <- c(0,0,0,0,1,0,0,0,0)
      row5 <- matrix(k5,nrow=sims,ncol=length(k5), byrow=TRUE)
    } else {
      row5 <- cbind(NFSCD[,time,1]/k5, NFSCD[,time,2]/k5, NFSCD[,time,3]/k5,
                    NFSCD[,time,4]/k5, 0, 0, 0, 0, 0)
    }

    k6 <- c(0,0,0,0,0,1,0,0,0)
    row6 <- matrix(k6,nrow=sims,ncol=length(k6), byrow=TRUE)

    k7 <- as.vector(NFHosp[,time,1] + NFHosp[,time,2] + NFHosp[,time,3] +
                      NFHosp[,time,4])
    if (all(k7 == 0)){
      k7 <- c(0,0,0,0,0,0,1,0,0)
      row7 <- matrix(k7,nrow=sims,ncol=length(k7), byrow=TRUE)
    } else {
      row7 <- cbind(NFHosp[,time,1]/k7, NFHosp[,time,2]/k7, NFHosp[,time,3]/k7,
                    NFHosp[,time,4]/k7, 0, 0, 0, 0, 0)
    }

    k8 <- c(0,0,0,0,0,0,0,1,0)
    row8 <- matrix(k8,nrow=sims,ncol=length(k8), byrow=TRUE)

    k9 <- c(0,0,0,0,0,0,0,0,1)
    row9 <- matrix(k9,nrow=sims,ncol=length(k9), byrow=TRUE)

    probs <- cbind(probs,row1, row2, row3, row4, row5, row6, row7, row8, row9)
  }
  return (probs)
}

utilities <- function (){
  "sets the utility values as global variables."

  #Utilities
  util_i <<- rbeta(sims, 395.792, 89.842)
  util_ii <<- rbeta(sims, 710.525, 276.315)
  util_iii <<- rbeta(sims, 359.881, 250.087)
  util_iv <<- rbeta(sims, 51.87, 50.236)
}

costs <- function(treatment){
  "
  I've put all the costs in here. Most of them are used to calculate the cost of
  hospitalisation by using the conditional probabilities of the possible procedures
  and their associated values. Treating costs and testing costs are also in here.
  These are set as global and used elsewhere."

  # days_ICU      <- np.random.normal(5.7, 0.009, size=width)
  # days_CCU      <- np.random.normal(6.8, 0.107, size=width)
  # days_hosp     <- np.random.normal(11.80, 0.510, size=width)
  cost_ICU        <- 1689.08  # CCU-01, XCO1Z-7Z
  cost_CCU        <- 1660.38  # CCU-06,XCO1Z-7Z
  cost_CABG       <- 12663.97  # ED26,27,28A-C
  cost_PTCA       <- 4414.38197  #EY40A-D, EY41A-D, EY44A-D
  cost_Transplant <- 79362.50  # ED04Z AND ED05Z weighted by FCE
  cost_h          <- 314.24 #VB01Z-VB09Z,VB11Z

  #cost_mt        <<- 0
  cost_ICD        <<- 3286.679171  #EY02A,EY02B
  cost_CRTP       <<- 7782.235254 #EY04A, EY04B
  cost_CRTD       <<- 7938.94858 #EY01A, EY01B

  cost_RF         <<- 191.27 + 130.26 #EY50Z, EY51Z
  cost_PEFA       <<- 2041.786715 + 350 #EY32A,EY32B

  if (treatment == 'crt-d' | treatment == 'crt-p'){
    #Conditioinal probability of procedure given hospitalisation
    n = rdirichlet(sims, c(0.0908 * 16.17, 0.3155 * 56.17, 0.0009 * 0.17,
                           0.0346 * 6.17, 0.0571 * 10.17, 0.5009 * 89.17))
    pICU        <- n[, 1]
    pCCU        <- n[, 2]
    pCABG       <- n[, 3]
    pPTCA       <- n[, 4]
    pTransplant <- n[, 5]
    pnone       <- n[, 6]

  } else {
    #Conditioinal probability of procedure given hospitalisation - mt
    n = rdirichlet(sims, c(0.0715 * 25.17, 0.2562 * 90.17, 0.0033 * 1.17,
                           0.0204 * 7.17, 0.0260 * 9.17, 0.6226 * 219.17))
    pICU        <- n[, 1]
    pCCU        <- n[, 2]
    pCABG       <- n[, 3]
    pPTCA       <- n[, 4]
    pTransplant <- n[, 5]
    pnone       <- n[, 5]
  }

  cost_hosp <-  array(numeric(),c(sims))
  for (i in 1:sims){
    cost_hosp[i] =  pICU[i] * cost_ICU + pCCU[i] * cost_CCU + pCABG[i] * cost_CABG + pPTCA[i] * cost_PTCA +
      pTransplant[i] * cost_Transplant + pnone[i]*cost_h
  }

  return (cost_hosp)
}

markov_chain <- function(treatment,  D = '++', type = 'base'){
  "
  This is the guts of the sim function. It sets the required values in memory and
  then calls rapidmarkov.

  "

  "

  I've commented out the call to settings as I call this function from sim() and
  have initialized things there. If you want to use this function to run your own
  sims, it will need to be brought back into the body of the code.

  "

  #settings(100,400,65,0.12)
  Ht(type)
  utilities()
  states <- c('NYHA I', 'NYHA II', 'NYHA III', 'NYHA IV',
              'NFSCD', 'FSCD', 'NFHosp', 'FHosp','Other')
  startingStates <- c(array(as.numeric(unlist(init_state(treatment)), dim = c(4))),
                      0, 0, 0, 0, 0)
  #startingStates <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)

  probs_table <- make_probs(treatment, D)

  transitionmatrix <- converttomatrices(states, probs_table)
  trace <- psamarkov(transitionmatrix, startingStates, cycles, states)
  stateUtilities <- cbind(util_i, util_ii, util_iii, util_iv, rep(0,sims),
                          rep(0,sims), rep(0,sims), rep(0,sims), rep(0,sims))
  stateCosts <- cbind(rep(0,sims), rep(0,sims), rep(0,sims), rep(0,sims),
                      rep(0,sims), rep(0,sims), costs(treatment),
                      costs(treatment), rep(0,sims))
  transPeriod <- (1/12)
  discount <- 0.035
  QALYs <- psavalues(trace, stateUtilities, transPeriod, QALYs = TRUE)
  Costs <- psavalues(trace, stateCosts, transPeriod)

  "

  this block adds the treatment cost to the simulated costs. I've done by taking
  the discounted SAVI output and just adding. I've also updated the mean and
  interval in the list of discounted costs.

  "
  if (treatment == 'mt'){
    Costs[[2]][[6]] <- Costs[[2]][[6]] + 0
    Costs[[2]][[5]] <- c(mean (Costs[[2]][[6]]), quantile(Costs[[2]][[6]], c(0.025,0.975)))
  }else if (treatment == 'icd'){
    Costs[[2]][[6]] <- Costs[[2]][[6]] + cost_ICD
    Costs[[2]][[5]] <- c(mean (Costs[[2]][[6]]), quantile(Costs[[2]][[6]], c(0.025,0.975)))
  }else if (treatment == 'crt-p'){
    Costs[[2]][[6]] <- Costs[[2]][[6]] + cost_CRTP
    Costs[[2]][[5]] <- c(mean (Costs[[2]][[6]]), quantile(Costs[[2]][[6]], c(0.025,0.975)))
  }else if (treatment == 'crt-d'){
    Costs[[2]][[6]] <- Costs[[2]][[6]] + cost_CRTD
    Costs[[2]][[5]] <- c(mean (Costs[[2]][[6]]), quantile(Costs[[2]][[6]], c(0.025,0.975)))
  }

  return (list(Costs = Costs, QALYs = QALYs))
}

ROC_specificity <- function(){
  "

  This function puts all the ROC curves together.

  "

  aa0  <- rep(1,sims)
  aa1  <- rep(1,sims)
  aa2  <- rep(1,sims)
  aa3  <- rep(1,sims)
  aa4  <- rep(1,sims)
  aa5  <- rep(1,sims)
  aa6  <- rbeta(sims,759.661008,3.817392)
  aa7  <- rbeta(sims,2235.344816,34.040784)
  aa8  <- rbeta(sims,3650.946	,93.614)
  aa9  <- rbeta(sims,580.524972	,43.695428)
  aa10 <- rep(0,sims)

  ab0  <- rep(1,sims)
  ab1  <- rep(1,sims)
  ab2  <- rep(1,sims)
  ab3  <- rep(1,sims)
  ab4  <- rbeta(sims,759.661008,	3.817392)
  ab5  <- rbeta(sims,2235.344816,	34.040784)
  ab6  <- rbeta(sims,736.914528,	15.039072)
  ab7  <- rbeta(sims,1415.207424,	58.966976)
  ab8  <- rbeta(sims,1094.731333,	88.762)
  ab9  <- rbeta(sims,481.050304,	91.62862933)
  ab10 <- rep(0,sims)

  ac0  <- rep(1,sims)
  ac1  <- rep(1,sims)
  ac2  <- rbeta(sims,759.661008,3.817392)
  ac3  <- rbeta(sims,2235.344816,	34.040784)
  ac4  <- rbeta(sims,3650.946,	93.614)
  ac5  <- rbeta(sims,1415.207424,	58.966976)
  ac6  <- rbeta(sims,2035.722656,	129.939744)
  ac7  <- rbeta(sims,777.024,	86.336)
  ac8  <- rbeta(sims,679.4186131,	124.6270829)
  ac9  <- rbeta(sims,448.213584,	161.600816)
  ac10 <- rep(0,sims)

  ad0  <- rep(1,sims)
  ad1  <- rbeta(sims,2235.344816,	34.040784)
  ad2  <- rbeta(sims,1415.207424,	58.966976)
  ad3  <- rbeta(sims,7969.2760196,	67.38282489)
  ad4  <- rbeta(sims,1327.56083,	139.3572142)
  ad5  <- rbeta(sims,993.576576,	161.745024)
  ad6  <- rbeta(sims,1196.415036,	280.640564)
  ad7  <- rbeta(sims,869.1860883,	297.5066477)
  ad8  <- rbeta(sims,463.918896,	233.703504)
  ad9  <- rbeta(sims,416.849744,	362.308656)
  ad10 <- rep(0,sims)

  ae0  <- rep(1,sims)
  ae1  <- rbeta(sims,6305.588432,	297.121968)
  ae2  <- rbeta(sims,3110.796,	345.644)
  ae3  <- rbeta(sims,1808.941968,	306.779632)
  ae4  <- rbeta(sims,1228.512,	307.128)
  ae5  <- rbeta(sims,869.1860883,	297.5066477)
  ae6  <- rbeta(sims,907.812679,	417.461305)
  ae7  <- rbeta(sims,618.8223893,	395.640544)
  ae8  <- rbeta(sims,553.4899413,	510.913792)
  ae9  <- rbeta(sims,230.096,	345.144)
  ae10 <- rep(0,sims)

  af0  <- rep(1,sims)
  af1  <- rbeta(sims,9859.982,	799.458)
  af2  <- rbeta(sims,17005.74533,	3119.397072)
  af3  <- rbeta(sims,5237.904672,	1564.568928)
  af4  <- rbeta(sims,2514.658761,	1103.555284)
  af5  <- rbeta(sims,2485.610896,	1556.032837)
  af6  <- rbeta(sims,1287.702144,	1096.931456)
  af7  <- rbeta(sims,710.5713658,	817.5390982)
  af8  <- rbeta(sims,381.7652942,	622.8802169)
  af9  <- rbeta(sims,240.665152,	618.853248)
  af10 <- rep(0,sims)

  ag0  <- rep(1,sims)
  ag1  <- rbeta(sims,1190.442816,	68.17513926)
  ag2  <- rbeta(sims,2407.440656,	357.0844993)
  ag3  <- rbeta(sims,1804.960417,	474.9895833)
  ag4  <- rbeta(sims,2011.833547,	795.3760533)
  ag5  <- rbeta(sims,1306.858733,	716.6644667)
  ag6  <- rbeta(sims,1224.6678,	952.5194)
  ag7  <- rbeta(sims,1037.82593,	1166.406134)
  ag8  <- rbeta(sims,390.3456933,	662.2719067)
  ag9  <- rbeta(sims,92.67048694,	260.3599395)
  ag10 <- rep(0,sims)

  PEFA1 <- cbind (aa0,aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9,aa10)
  PEFA2 <- cbind (ab0,ab1,ab2,ab3,ab4,ab5,ab6,ab7,ab8,ab9,ab10)
  PEFA3 <- cbind (ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10)
  PEFA4 <- cbind (ad0,ad1,ad2,ad3,ad4,ad5,ad6,ad7,ad8,ad9,ad10)
  PEFA5 <- cbind (ae0,ae1,ae2,ae3,ae4,ae5,ae6,ae7,ae8,ae9,ae10)
  PEFA6 <- cbind (af0,af1,af2,af3,af4,af5,af6,af7,af8,af9,af10)
  RF    <- cbind (ag0,ag1,ag2,ag3,ag4,ag5,ag6,ag7,ag8,ag9,ag10)

  rocs <- list(PEFA1, PEFA2, PEFA3, PEFA4, PEFA5, PEFA6, RF)
  names(rocs) <- c('PEFA1', 'PEFA2', 'PEFA3', 'PEFA4', 'PEFA5', 'PEFA6', 'RF')
  return(rocs)
}

weights <- function(itp, ifn, ifp, itn, specificity){
  "
  this function takes sim results for ++, +-, -+, -- and weights them using the
  ROC curve.

  "

  #calulcate weights
  specificity <- array(as.numeric(unlist(specificity)), dim=c(sims,11))
  sensitivity <- array(seq(0,1,length.out=11), dim=(11))
  tp <- array(rep(0,11), dim=c(11))
  fn <- array(rep(0,11), dim=c(11))
  fp <- array(rep(rep(0,11),sims), dim=c(sims,11))
  tn <- array(rep(rep(0,11),sims), dim=c(sims,11))
  for (j in 1:11){
    tp[j] = prevalence*sensitivity[j]
    fn[j] = prevalence*(1-sensitivity[j])
    for (i in 1:sims){
      fp[i,j] = (1-prevalence)*(1-specificity[i,j])
      tn[i,j] = (1-prevalence)*(specificity[i,j])
    }
  }
  #weights <- list(tp = tp,fn = fn,fp = fp,tn = tn)

  #get the costs and utilities from the input.
  c1 <- itp[[1]][[2]][[6]]
  c2 <- ifn[[1]][[2]][[6]]
  c3 <- ifp[[1]][[2]][[6]]
  c4 <- itn[[1]][[2]][[6]]

  u1 <- itp[[2]][[2]][[6]]
  u2 <- ifn[[2]][[2]][[6]]
  u3 <- ifp[[2]][[2]][[6]]
  u4 <- itn[[2]][[2]][[6]]

  #do the calculation
  c <- array(numeric(),c(sims, 11))
  u <- array(numeric(),c(sims, 11))
  for (j in 1:11){
    for (i in 1:sims){
      c[i,j] = tp[j]*c1[i] + fn[j]*c2[i] + fp[i,j]*c3[i] + tn[i,j]*c4[i]
      u[i,j] = tp[j]*u1[i] + fn[j]*u2[i] + fp[i,j]*u3[i] + tn[i,j]*u4[i]
    }
  }

  return (list(Costs = c, QALYs = u))
}

sim <- function(scenario ='treatment', sims = 100, cycles = 400, age = 65,
                prevalence = 0.12){

  "
  This is the main function. scenario 'treatment' returns the simulated outcomes
  for patients treated with each intervention ignoring testing.

  scenario 'typei' returns the simulated results when test positive patients are
  treated with icd and test negative patients with crt-p.

  scenario 'typeii' returns the simulated results when test positive patients are
  treated with crt-d and test negative patients with crt-p.

  "

  settings (sims,cycles,age,prevalence)
  if (scenario == 'treatment'){
    #run the simulation
    mt <- markov_chain('mt', '++', 'base')
    icd <- markov_chain('icd', '++', 'base')
    crtp <- markov_chain('crt-p', '++', 'base')
    crtd <- markov_chain('crt-d', '++', 'base')
    return (list(mt = mt, icd = icd, crtp = crtp, crtd = crtd))
  }
  else if (scenario == 'typei'){
    #run the simulation
    i_tp <-  markov_chain ('icd', '++' , 'scaled')
    i_fn <-  markov_chain ('mt' , '+-' , 'scaled')
    i_fp <-  markov_chain ('icd', '-+' , 'scaled')
    i_tn <-  markov_chain ('mt' , '--' , 'scaled')

    #weight the results using the ROCs
    spec <- ROC_specificity()
    empty <- list(1,2,3,4,5,6,7)
    result <- lapply(empty, function(empty)weights(i_tp, i_fn, i_fp, i_tn, spec[empty]))
    names(result) <- c('iPEFA1', 'iPEFA2', 'iPEFA3', 'iPEFA4', 'iPEFA5', 'iPEFA6', 'iRF')

    #add the test costs to the result. The solution is a bit btech.
    for (i in 1:7){
      if (i == 7){
        result[[i]][[1]] <- result[[i]][[1]] + cost_RF
        break
      }
      result[[i]][[1]] <- result[[i]][[1]] + cost_PEFA
    }
    return (result)
  }
  else if (scenario == 'typeii'){

    #same as for typei
    ii_tp <- markov_chain ('crt-d', '++', 'scaled')
    ii_fn <- markov_chain ('crt-p', '+-', 'scaled')
    ii_fp <- markov_chain ('crt-d', '-+', 'scaled')
    ii_tn <- markov_chain ('crt-p', '--', 'scaled')

    spec <- ROC_specificity()
    empty <- list(1,2,3,4,5,6,7)
    result <- lapply(empty, function(empty)weights(ii_tp, ii_fn, ii_fp, ii_tn, spec[empty]))
    names(result) <- c('iiPEFA1', 'iiPEFA2', 'iiPEFA3', 'iiPEFA4', 'iiPEFA5', 'iiPEFA6', 'iiRF')

    for (i in 1:7){
      if (i == 7){
        result[[i]][[1]] <- result[[i]][[1]] + cost_RF
        break
      }
      result[[i]][[1]] <- result[[i]][[1]] + cost_PEFA
    }
    return (result)
  }

}

pCE <- function (c2,c1,u2,u1){
  #wrote the funcion
  #haven't used it for anything
  #also haven't tested it
  #should work, no promises
  wtp <- seq(from=0,to=100000,by=100)
  nmbnew <- array(numeric(),c(length(wtp), sims))
  nmbold <- array(numeric(),c(length(wtp), sims))
  pCE <- array(numeric(),c(length(wtp), sims))
  for (i in 1:length(wtp)){
    for (j in 1:sims){
      nmbnew[i, j] = u2[j] * wtp[i] - c2[j]
      nmbold[i, j] = u1[j] * wtp[i] - c1[j]
      if (nmbnew[i, j] > nmbold[i, j]){
        pCE[i, j] = 1
      }else{
        pCE[i, j] = 0
      }
    }
  }
  return (pCE)
}

#treatment effect
treatment <- sim()
#typei results
testingi  <- sim('typei')
#typeii results
testingii <- sim('typeii')

#make the iso net benefit curve
#could use ggplot
#this worked and all I can see now is a blank plot.
plot.new()
for (i in 1:7){
  #change testingi to testingii in the next two lines for the other plot
  c2 <- array(as.numeric(unlist(testingii[[i]][[1]])), dim = c(sims,11))
  u2 <- array(as.numeric(unlist(testingii[[i]][[2]])), dim = c(sims,11))

  empty <- seq(from=1,to=11,by=1)
  cost <- lapply(empty,function(empty)mean(c2[,empty]))
  qalys <- lapply(empty,function(empty)mean(u2[,empty]))

  points (qalys,cost, type='b')
  lines (qalys,cost, type='b')
}













