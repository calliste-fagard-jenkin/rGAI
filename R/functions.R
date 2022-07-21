#' rGAI: Generalised Abundance Index fitting for seasonal invertibrates, with
#' covariates
#'
#' Fitting of the GAI from Dennis et al (2016), with normal mixture, stopover
#' and spline-based seasonal trend models, as well as poisson, negative binomial
#' and zero-inflated count models. Covariates can be fitted to mean brood
#' arrival time, mean brood dispersion and the proportion of the population
#' which finds itself in each brood. The package slightly extends the
#' supplementary materials of the paper to provide flexible model-fitting
#' options and covariate inclusion.
#' 
#' @details The package should be used through the fit_GAI function.
#' Options should be specified to give more detail of the model being fitted.
#' Bootstraps can be used to produce confidence interval estimates of the 
#' abundance index and parameter values.
#' @author C. Fagard-Jenkin, based on work by E. Dennis, et al.
#' @name rGAI
#' @useDynLib rGAI
#' @importFrom Rcpp sourceCpp
#' @importFrom splines bs
#' @importFrom mvtnorm rmvnorm
#' @import magrittr
#' @import parallel
#' @import doParallel
#' @docType package
NULL

#' Normal mixture sampling
#'
#' produces samples from a mixture of n normal distributions.
#'
#' @param n The number of samples to generate
#' @param par A list of parameter values. The first entry contains the
#' @param skeleton An object which helps to relist the vector of parameters into
#' the expected format.means of each component, the second the SDs, and the last
#' contains the weights. An exponential link is applied to the
#' SDs, and the probs_link function is used on the weights.
#' probs_link is simply a variant of the qlogis link
#'
#' @return A vector of n generated normal-mixture samples.
#' @export
sim_n_mixture <- function(n, par, skeleton = NULL){
  # purpose : produces samples from a mixture of n normal distributions.
  # inputs  : n   - The number of samples to generate
  #           par - A list of parameter values. The first entry contains the
  #                 means of each component, the second the SDs, and the last
  #                 contains the weights. An exponential link is applied to the
  #                 SDs, and the probs_link() function is used on the weights.
  #                 The probs_link() is simply a variant of the qlogis() link
  #                 which ensures the sum of all the probabilities is 1.
  # output  : A vector of n generated samples
  
  # Extract pars and apply links:
  if (!is.null(skeleton)){
    par %<>% relist(skeleton)
    means <- means_link(par$mu)
    sds <- exp(par$sigma)
    probs <- probs_link(par$w)
  }
  
  else{
    means <- means_link(par[[1]])
    sds <- exp(par[[2]])
    probs <- probs_link(par[[3]])
  }
  
  # Check dimensions are consistent:
  if (length(means)!=length(sds) | length(sds)!=length(probs))
    stop('inconsistent component numbers implied by parameter choices')
  
  # Generate the deviates:
  groups <- sample(length(means), n, replace = T, prob = probs)
  rnorm(n, means[groups], sds[groups]) %>% return
}

#' Normal mixture negative Log Likelihood
#'
#' Evaluates an n-mixture model negative log-likelihood given the parameters and
#' the observations
#'
#' @param par A list of parameters including mean, sd, and probs.
#' Alternatively, a vector of parameters, with a skeleton to relist the
#' parameters (parameter: skeleton).
#' @param obs A matrix of observed values.
#' @param skeleton The skeleton to \code{\link{relist}} the parameters
#' @param returnNLL If TRUE, returns the value of the NLL rather than the
#' likelihood for each point individually. If FALSE, also returns the output as
#' an nS*nT matrix.
#' @param DMs A named list of design matrices for covariate inclusion
#'
#' @return A matrix of evaluated densities for the count at each site and
#' location, or the single real number valued negative log likelihood
#'
#' @export
mixture_nll <- function(par, obs, skeleton, returnNLL = T, DMs = list()){
  # purpose : Evaluates an n-mixture model negative log-likelihood given the
  #           parameters and the observations
  # inputs  : par       - A list of parameters including mean, sd, and probs.
  #                       Alternatively, a vector of parameters, with a skeleton 
  #                       to relist the parameters (parameter: skeleton).
  #           obs       - A matrix of observed values.
  #           skeleton  - The skeleton to relist() the parameters 
  #           returnNLL - If TRUE, returns the value of the NLL rather than the
  #                       likelihood for each point individually. If FALSE, also
  #                       returns the output as an nS*nT matrix.
  #           DMs       - A named list of design matrices for covariate
  #                       inclusion
  # output  : A number on the real scale, the NLL of all the data points given
  #           the model parameters, for a normal mixture model
  
  # get parameter values with or without covariates, and apply the links:
  n <- length(obs)
  nS <- nrow(obs) ; nT <- ncol(obs)
  parameter_vals <- get_parameter_values(par, DMs, skeleton, n)
  means <- parameter_vals$means %>% t %>% as.vector
  sds <- parameter_vals$sds %>% t %>% as.vector
  probs <- parameter_vals$probs %>% t %>% as.vector
  
  B <- length(skeleton$mu)
  
  # if SDs for each brood are the same, replicate so dimensions match:
  if (length(sds) <= length(means)) sds <- rep(sds, length(means)/length(sds))
  
  # Now vectorise obs to only require one call to dnorm (far faster):
  obs <- rep(obs, rep(B, n))
  
  lfunc <- ifelse(isTRUE(returnNLL),
                  function(x) -sum(log(x), na.rm = T),
                  function (x) matrix(x, nrow = nS, ncol = nT))
  
  # Now calculate the likelihood:
  # note on the sequence of events. We: Calculate the likelihood for all data
  # points and all possible broods for each point %>% Split these probabilities
  # into a separate group for each data point, containing the likelihood value
  # for each brood it could belong to %>% add these values together to get the
  # unconditional likelihood of the point irrespective of brood %>% convert
  # these values from a list to a numeric vector %>% either -sum(log()) the
  # values to give the negative log likelihood, or prod the values to give the
  # likelihood %>% return the NLL or L
  dnorm(obs, means, sds) %>% `*`(probs) %>% split(rep(1:n, each = B)) %>%
    lapply(FUN = sum) %>% as.numeric %>% lfunc %>% return
}

#' Stopover model density evaluation.
#'
#' Evaluates the CDF of a mixture of normal distributions, then
#' returns the probability of being in a given bin.
#'
#' @param breaks The points to evaluate. The function return the probability of
#' being in each bin defined by the breaks.
#' @param par A vector of parameter values which should be
#' \code{\link{relist}}ed
#' @param skeleton The skeleton to \code{\link{relist}} the parameters
#' @param DMs A list of design matrices to be used when covaraites are included
#' in the model.
#' @param nS In the case where DMs is empty because we use no covariates, we
#' need to know the number of sites to determine the number of points.
#' @return A matrix of evaluated densities for the count at each site and
#' location, or the single real number valued negative log likelihood
#'
#' @export
evaluate_stopover <- function(breaks, par, skeleton, DMs = list(), nS){
  # purpose : Evaluates the CDF of a mixture of normal distributions, then 
  #           returns the probability of being in a given bin.
  # inputs  : breaks - The points to evaluate. The function return the
  #                    probability of being in each bin defined by the breaks.
  #           par    - A list of parameter values. The first entry contains the
  #                    means of each component, the second the SDs, and the last
  #                    contains the weights. An exponential link is applied to
  #                    the SDs, and the probs_link() function is used on the
  #                    weights. The probs_link() is simply a variant of the
  #                    qlogis() link which ensures the sum of all the
  #                    probabilities is 1.
  #           DMs    - A list of design matrices to be used when covaraites are
  #                    included in the model.
  #           nS     - In the case where DMs is empty because we use no
  #                    covariates, we need to know the number of sites to
  #                    determine the number of points.
  # output  : A vector of length breaks - 1 with the probability an individual
  #           enters population between the lower and upper limit of the 
  #           interval specified by the break points
  # Note    : The stopover model is incompatible with time-varying covariates
  #           for the means, SDs and weights of the arrival times, and so the 
  #           DMs should be checked for this before being passed ot this
  #           function. It is not checkled in this function since this would
  #           lead to unnecessary computation time during parameter estimation
  #           with optim.
  # get parameter values with or without covariates, and apply the links:
  n <- length(breaks)
  nT <- length(breaks) - 1
  parameter_vals <- get_parameter_values(par, DMs, skeleton, nS * nT)
  means <- parameter_vals$means
  sds <- parameter_vals$sds
  probs <- parameter_vals$probs
  phi <- par %>% relist(skeleton) %>% `$`(phi) %>% plogis
  B <- ncol(means)
  
  # if SDs for each brood are the same, replicate so dimensions match:
  if (ncol(sds) == 1) sds %<>% rep(B) %>% matrix(ncol = B)
  
  output <- matrix(NA, nrow = nS, ncol = nT)
  components <- rep(1:B, each = n)
  points <- rep(breaks, B)
  
  # Univoltine causes an edge case:
  if (length(probs) == 1) probs %<>% rep(nS * nT) %>% as.matrix
  
  for (i in 1:nS){
    # since we know the parameter values are the same for all time points for a 
    # given site, we simply loop through the sites and get the corresponding
    # row of "a_func":
    means_i <- means[i,] %>% `[`(components)
    sds_i <- sds[i, ] %>% `[`(components)
    probs_i <- probs[i, ]
    
    # Get CDF at each breakpoint for each component %>% break up the results
    # into a list where each entry is the calculation for all components for a
    # given data point %>% remove the unecessary points %>% multiply by the
    # probability of being in each component and add these probabilities together
    # %>% use C++ code to perform the for loop which multiplies by the retention
    # probabilities:
    output[i,] <- pnorm(points, means_i, sds_i) %>% diff %>%
      split(rep(1:n, B)[-n*B]) %>% `[`(-n) %>%
      lapply(function (x) sum(x*probs_i)) %>% unlist %>%
      stopover_for_loop(phi = phi)
  }
  
  return(output)
}

#' Spline model density evaluation.
#'
#' Calculates the seasonal flight curve using splines
#'
#' @param par The vector of estimated spline coefficients (and also potentially
#' the distributional parameters for ZIP and NB likelihoods)
#' @param skeleton The skeleton to \code{\link{relist}} the parameters
#' @param spline_specs The list containing the 'degree' of the splines and 
#' also df (their degrees of freedom)
#' @param nT The number of sampling occasions per site
#' @param nS The number of sampling locations during the survey
#' @return A vector of values of length nT, giving the relative density
#' of expected counts at each occasion
#' @export
evaluate_splines <- function(par, skeleton, spline_specs, nT, nS){
  # purpose : Calculates the seasonal flight curve using splines
  # inputs  : par          - The vector of estimated spline coefficients (and 
  #                          also potentially the distributional parameters for
  #                          ZIP and NB likelihoods)
  #           skeleton     - The skeleton which allows us to relist the vector
  #                          of parameters correctly
  #           spline_specs - The list containing the 'degree' of the splines and
  #                          also df (their degrees of freedom)
  #           nT           - The number of sampling occasions per site
  #           nS           - The number of sites in the study
  # outputs : A vector of values of length nT, giving the relative density
  #           of expected counts at each occasion
  
  # extract spline settings:
  degree <- spline_specs$degree
  df <- spline_specs$df
  
  # extract parameters
  par %<>% relist(skeleton) %>% `$`(spline_par)
  
  # calculate a_func:
  bsbasis <- bs(1:nT, df = df, degree = degree, intercept = T)
  a_func <- (bsbasis%*%par) %>% exp %>% sum_to_one %>% return
  return(a_func)
}

#' Generate random parameter values for non-spline models
#'
#' For ease of simulation studies or testing, generates some parameters for 
#' non-covariate including normal mixture and stopover models
#'
#' @param nT The number of sampling occasions in the experiment
#' @param skeleton The skeleton to \code{\link{relist}} the parameters
#' @param sd_max The maximum SD for broods
#' @param l The lower limit for the distributional parameter (before the link
#' function is applied)
#' @param u The upper limit for the distributional parameter (before the link
#' function is applied)
#' @return A vector of parameters for the model, such that it can be relisted
#' using skeleton
#' @export
generate_mixture_pars <- function(nT, skeleton, sd_max = 1, l = -5, u = 5){
  # purpose : Produces random parameter values for a mixture model as a means of
  #           introducing variation when testing and debugging
  # inputs  : nT       - The number of sampling occasions in the experiment
  #           skeleton - The skeleton for the model, which helps to indicate
  #                      the number of broods and SDs which are needed.
  #           sd_max   - The maximum SD for broods.
  #           l        - The lower limit for the distributional parameter
  #                      (before the link function is applied)
  #           u        - The upper limit for the distributional parameter
  #                      (before the link function is applied).
  # outputs : A vector of parameters for the model, such that it can be
  #           relisted using skeleton
  
  n_means <- length(skeleton$mu)
  n_sigmas <- length(skeleton$sigma)
  copy <- skeleton
  
  # generate the means, sds and weights for the mixture model:
  copy$mu <- runif(n_means, 0, nT) %>% sort %>% means_backtransform
  copy$sigma <- runif(n_sigmas, 0, sd_max) %>% log
  copy$w <- rbeta(n_means - 1, 1, (n_means - 1):1) %>% qlogis
  
  # generate a distributional parameter if it's required:
  if (!is.null(skeleton$dist.par)) copy$dist.par <- runif(1, l, u)
  
  # generate a phi for the stopover if required:
  if (!is.null(skeleton$phi)) copy$phi <- runif(1) %>% qlogis
  
  return(copy)
}
 
#' Generate random parameter values for spline models
#'
#' For ease of simulation studies or testing, generates some parameters for 
#' non-covariate including spline models
#'
#' @param skeleton The skeleton to \code{\link{relist}} the parameters
#' @param sigma The standard deviation of the normal distribution the parameters
#' @param l The lower limit for the dist.par parameter on the link scale
#' @param u The upper limit for the dist.par parameter on the link scale
#' are generated from
#' @return A vector of parameters for the model, such that it can be relisted
#' using skeleton
#' @export
generate_spline_pars <- function(skeleton, sigma = 2, l = -5, u = 5){
  # purpose : Creates a random list of parameter values for a spline model,
  #           including a distributional parameter if the model is ZIP or NB
  # inputs   : skeleton - A skeleton created by produce_skeleton()
  #            sigma    - The standard deviation of the normal distribution the
  #                       parameters are generated from
  #            l        - The lower limit for the dist.par parameter on the
  #                       link scale
  #            u        - The upper limit for the dist.par parameter on the link
  #                       scale
  # output  : A list of parameters which when passed to unlist() will be in the
  #           form fit_GAI() expects for an initial guess for parameter values
  output <- skeleton
  
  # Generate the spline parameters:
  if (is.null(skeleton$spline_par)) stop("spline skeleton must be provided")
  else output$spline_par <- rnorm(length(skeleton$spline_par), sd = sigma)
  
  # Generate the distributional parameter for NB or ZIP models:
  if (!is.null(skeleton$dist.par)) output$dist.par <- runif(1, l, u)
  return(output)
}

#' Simulate data from the GAI
#'
#' Produces a simulated data set for observed counts of individuals at each of
#' nS sites for nT occasions. Note, this function was implemented before 
#' covariate inclusion, and has no ability to deal with covariate-including
#' data simulation.
#'
#' @param par The vector of parameter values to be used
#' @param nT The integer number of sampling occasions
#' @param totals The integer site total for each site
#' @param a_type The character name of the seasonal distribution
#' @param dist_type The character name of the count distribution
#' @param options A list of options to specify behaviours of the seasonal
#' distribution. Inlcludes B, the number of broods. Currently covariates are
#' not supported.
#' @param thin The proportion of sampling opportunities which should be missed,
#' and replaced with an NA.
#' @param returnDF If TRUE, returns a data.frame with site, occasion, and count
#' columns, rather than a matrix of counts
#' @return A matrix of simulated counts for each site (rows) and each occasion
#' (columns).
#' @export
simulate_data <- function(par, nT = 25, totals = rep(500, 50),
                          a_type = "mixture", dist_type = "P",
                          options = list(), thin = 0.4, returnDF = F){
  # purpose : Produces a simulated data set for observed counts of individuals
  #           at each of nS sites for nT occasions
  # inputs  : par       - The vector of parameter values to be used
  #           nT        - The integer number of sampling occasions
  #           total     - The integer site total for each site
  #           a_type    - The character name of the seasonal distribution
  #           dist_type - The character name of the count distribution
  #           options   - A list of options to specify behaviours of the
  #                       seasonal distribution
  #           thin      - The proportion of sampling opportunities which should
  #                       be missed, and replaced with an NA.
  #           returnDF  - If TRUE, returns a data.frame with site, occasion, and
  #                       count columns, rather than a matrix of counts
  # output   : A matrix of simulated counts for each site (rows) and each
  #            occasion (columns).
  
  # extract distributional parameter if dist_type is not Poisson:
  skeleton <- produce_skeleton(a_type, dist_type, options)$skeleton
  dist_par <- relist(par, skeleton)$dist.par
  ns <- length(totals)
  
  if (is.null(dist_par) & dist_type != "P")
    stop("missing distributional parameter")
  
  # Aid the user in supplying the parameters correctly:
  if (length(par) != length(unlist(skeleton))){
    message <- paste(length(par), "parameters supplied when",
                     length(unlist(skeleton)), "are required.",
                     "The parameters required have been printed.")
    
    print(skeleton)
    
    stop(message)
  }
  
  nS <- length(totals)
  totals_matrix <- matrix(rep(totals, nT), nrow = nS)
  breaks <- c(-Inf, 1:(nT - 1), Inf)
  
  # determine the density of the flight pattern given the chosen parameters:
  flight_pattern <- switch(a_type,
                           mixture = mixture_nll(par, matrix(1:nT, nrow = 1),
                                                 skeleton, F),
                           stopover = evaluate_stopover(breaks, par, skeleton,
                                                        nS = nS),
                           splines = evaluate_splines(par, skeleton, options,
                                                      nT, nS) %>%
                             # easier to turn spline output into a matrix 
                             # instead of doing another call to switch for mu:
                             matrix(nrow = nS, ncol = nT, byrow = T))
  
  # scale the patterns by the site totals to get the expected count for 
  # each site/occasion combination, and then simulate the observed counts
  # before putting these back into a matrix for the user:
  mu <- flight_pattern %>% apply(2, function(x) x*totals) %>% unlist
  
  simulations <- switch(dist_type,
                        P = rpois(n = nS*nT, lambda = mu),
                        NB = rnbinom(n = nS*nT, size = exp(dist_par), mu = mu),
                        ZIP = rzip(n = nS*nT, mu, plogis(dist_par)))
  
  simulations %<>% matrix(nrow = nS) %>% thin_NA(alpha = thin)
  
  if (returnDF){
    data.frame(site = rep(1:nS, nT), occasion = rep(1:nT, each = nS),
               count = as.numeric(simulations)) %>% return}
  
  else return(simulations)
}

#' Zero-Inflated Poisson Density Evaluation
#'
#' Calculates the density of a zero-inflated poisson model at given data points
#' for a set of parameter values
#'
#' @param obs The observations. Can be a matrix or vector, and may contain NAs
#' @param lambda The vector of means for each observation. Must have the same 
#' length as obs
#' @param psi The probability of being from the poisson distribution as
#' opposed to being an 'inflated' 0
#' @param log A boolean, if TRUE the function returns the logs of the
#' densities instead
#' @return A vector (or matrix) of the same length (dimension) as obs, with
#' the density or log density of the zero-inflated poisson evaluated
#' at each point.
#' @export
dzip <- function(obs, lambda, psi, log = T){
  # purpose : Evaluates the density of a Zero=Infalted POisson distribution
  # inputs  : obs    - The observations. Can be a matrix or vector, and can
  #                    contain NAs
  #           lambda - The vector of means for each observation. Must be of the
  #                    same length as obs (can also be a matrix)
  #           psi    - The probability of being from the poisson distribution as
  #                    opposed to being an 'infalted' 0.
  #           log    - A boolean, if TRUE the function returns the logs of the
  #                    densities instead
  # outputs : A vector (or matrix) of the same length (dimesnion) as obs, with 
  #           the density or log density of the zero-inflated poisson evaluated
  #           at each point.
  
  # get the poisson likelihood of each point:
  output <- obs
  output <- try(dpois(obs, lambda), silent = T)
  zeros <- obs == 0 & !is.na(obs)
  others <- obs != 0 & !is.na(obs)
  
  output[zeros] <- (1 - psi) + psi*output[zeros]
  output[others] <- psi*output[others]
  
  if (isTRUE(log)) output %<>% log
  
  return(output)
}

#' Zero-Inflated Poisson Random Number Generation
#'
#' Simulates deviates from a zero-inflated poisson model given the model's
#' parameter values.
#'
#' @param n The integer number of deviates to generate
#' @param lambda The vector of means for each observation. Must have the same 
#' length n, or be of length 1 for constant means
#' @param psi The probability of being from the poisson distribution as
#' opposed to being an 'inflated' 0
#' @return A vector of n points from a zero inflated poisson distribution
#' @export
rzip <- function(n, lambda, psi){
  # purpose : Simulates values from a zero inflated poisson distribution
  # inputs  : n      - The number of elements to simulate
  #           lambda - The mean, or vector of means for each point
  #           psi    - The probability of being from the poisson distribution,
  #                    as opposed to being an 'inflated zero'.
  # outputs : A vector of n points from a zero inflated poisson distribution
  
  output <- rpois(n,lambda)
  output[weightedSelection(1:n, rep(1 - psi, n))] <- 0
  return(output)
}

nsolve_helper <- function(x, N, f, obs, a_func, dist_par, max_counter = 50){
  # purpose : Iteratives through uses of the uniroot() function to find the
  #           roots of the derivative of the log likelihood for the NB and ZIP
  #           cases.
  # inputs  : x           - An index for which site we are currently trying
  #                         to solve the problem for
  #           N           - The vector of current guesses for site totals
  #           f           - The function which evaluates the derivative of ll
  #           obs         - The matrix of observed counts
  #           a_func      - The matrix of arrival time densities
  #           dist_par    - The extra distributional parameter for the NB or ZIP
  #           max_counter - The maximum number of times to try using the 
  #                         uniroot function
  # output  : A single real number, the numerically calculated value of N which
  #           sets the derivative of the ll to 0.
  
  # fetch interval:
  upper <- ifelse(N[x] > 0, N[x]*2, max(N))
  lower <- ifelse(obs[x,] %>% na.omit %>% `>`(0) %>% any, 0.1, 0)
  if(class(lower)!="numeric"){print("lower error");print(lower)}
  
  # solve 1-D derivative:
  counter <- 0
  repeat{
    # attempt to solve:
    if (class(upper)!="numeric") print("upper error")
    output <- try(uniroot(f, interval = c(lower, upper), i = x,
                          obs = obs, a_func = a_func,
                          dp = dist_par)$root, silent = T)
    
    # increase upper limit if required:
    if (class(output) == "try-error" & counter < max_counter){
      counter %<>% `+`(1) ; upper %<>% `*`(2)}
    
    else {break}
  }
  
  return(output)
}

nsolve_u <- function(N, obs, a_func, dist_par, dist_choice){
  # purpose : Solves for N for the iterative MLE finder of the NB and ZIP
  #           likelihoods.
  # inputs  : N        - The vector of guesses for the site totals
  #           obs      - The matrix of obvservations
  #           dist_par - The r parameter of the NB, or the Psi for a ZIP
  #           a_choice - "NB" or "ZIP"
  # outputs : A vector of the same length as N, with the numerical solutions
  #           to the derivative of the LL given the MLEs\
  
  f_choice <- switch(dist_choice, NB = NB_mle_i, ZIP = ZIP_mle_i)
  
  new_n <- sapply(1:nrow(obs), nsolve_helper, f = f_choice, obs = obs, N = N, 
                  a_func = a_func, dist_par = dist_par)
  #cat('output N has length:', length(new_n))
  return(new_n)
}

NB_mle_i <- function(N, i, obs, a_func, dp){
  # purpose : Evaluates equation 7 of Dennis et al (2014), the partial deriv
  #           of the NB log likelihood with respect to the superpopulations at 
  #           each site
  # inputs  : N      - The vector of estimated site superpopulations
  #           i      - The occasion at which we solve for N
  #           obs    - The matrix of observed counts
  #           a_func - The estimated matrix of seasonal flight curve densities
  #           r      - The estimated number of events till failure for the NB
  # outputs : A vector of the same length as N with the numerically solved new
  #           estimates for N (at each site)
  # N %<>% `[`(i)
  
  # link function for r:
  dp %<>% exp
  
  # calculate and return the formula:
  LHS <- obs[i,]/N
  NUM <- (dp + obs[i,])*a_func[i,]
  DEN <- dp + N*a_func[i,]
  (LHS - NUM/DEN) %>% sum(na.rm = T) %>% return
}

ZIP_mle_i <- function(N, i, obs, a_func, dp){
  # purpose : Evaluates equation 7 of Dennis et al (2014), the partial deriv
  #           of the NB log likelihood with respect to the superpopulations at 
  #           each site
  # inputs  : N      - The vector of estimated site superpopulations
  #           obs    - The matrix of observed counts
  #           a_func - The estimated matrix of seasonal flight curve densities
  #           r      - The estimated number of events till failure for the NB
  # outputs : A vector of the same length as N with the numerically solved new
  #           estimates for N (at each site)
  # N %<>% `[`(i)
  
  # link function for phi:
  dp %<>% plogis
  
  # indicator for if observation is definitely not an 'inflated zero':
  delta <- (obs > 0) %>% `*`(1) %>% as.matrix
  
  # calculate the formula:
  exp_term <- exp(-N*a_func[i,])
  numer <- -dp*a_func[i,]*(1 - delta[i,])*exp_term
  denom <- 1 - dp + dp*exp_term
  
  result <- numer/denom - delta[i,]*a_func[i,] + (delta[i,]*obs[i,])/N
  result %>% sum(na.rm = T) %>% return
}

#' Obtain skeleton and Design Matrix
#'
#' Updates the skeleton and design matrix objects for a GAI model for a given 
#' input parameter type.
#'
#' @param base The character name of the type of parameter to be used. Can be
#' 'mu', 'sigma' or 'w'.
#' @param options A list of options given to the fit_GAI function to 
#' specify the model
#' @param DF The data frame with which covariate values are stored, for creating
#' design matrices
#' @param n The number of formulas to check (i.e. how many parameters of type
#' 'base' are we expecting?)
#' @param skeleton The object where skeleton items should be added
#' @param DMs The object where design matrix should be added
#' @return A named list containing the updated skeleton and DMs objects in 
#' entries 'skeleton' and 'DMs' respectively.
get_cov_skeleton_and_DM <- function(base, options, DF, n, skeleton, DMs){
  # purpose : Produces the correct entries in the skeleton and design matrix
  #           lists for stopover and mixture models
  # inputs  : base     - A character name of the type of parameter being dealt
  #                      with ("mu", "sigma" or "w")
  #           options  - A list of options given to the fit_GAI function to 
  #                      specify the model
  #           DF       - The data frame with which covariate values are stored,
  #                      for creating design matrices
  #           n        - The number of formulas to check
  #           skeleton - The object where skeleton items should be added
  #           DMs      - The object where design matrix should be added
  # outputs : A named list containing the updated skeleton and DMs objects
  # Extract the formulas from the options:
  formulas <- options[[paste(base, "_formula", sep = "")]]
  tvm <- "Covariate values must vary by site only, and not through time"
  
  # If the formula contains more than one thing, then it's a general formula 
  # shared among all of the 'base' type parameters:
  if (class(formulas) == "formula"){
    # Check formula doesn't involve NA covariate values or time-varying 
    # covariates:
    DMs[[base]] <- design_matrix(DF, formulas)
    if (is_time_varying(DMs[[base]], DF)) stop(tvm)
    skeleton[[paste(base, ".cov", sep = "")]] <- rep(NA, ncol(DMs[[base]]) - 1)
  }
  
  # Otherwise the formula is a list, containing formulas for each parameter
  # individually:
  else if (class(formulas) == "list"){
    # Only going up to n avoids odd behaviour where the formulas list contains
    # more entries than what we expect.
    for (i in 1:n){
      # Extract the formula, and do nothing if no covariate is specified:
      current_formula <- formulas[[i]]
      if (is.null(current_formula)) next
      
      # Produce the design matrix and update the output objects:
      DMs[[paste(base, i, sep = "")]] <- design_matrix(DF, current_formula)
      if (is_time_varying(DMs[[paste(base, i, sep = "")]], DF)) stop(tvm)
      n_par <- DMs[[paste(base, i, sep = "")]] %>% ncol %>% `-`(1)
      skeleton[[paste(base, i, ".cov", sep = '')]] <- rep(NA, n_par)
    }#for
  }#else if
  return(list(skeleton = skeleton, DMs = DMs))
}#get_cov_skeleton_and_DM

#' Produce a skeleton for GAI parameters
#' 
#' Produces a skeleton object which can be used by relist to restructure a 
#' vector of parameters into a named list of parameters for a GAI model
#' @param a_choice The character name for the seasonal trend function being
#' fitted. options include "splines", "stopover", and "mixture"
#' @param distribution The distribution of the counts, Poisson ("P"), Negative
#' Binomial ("NB") or Zero-Inflated Poisson ("ZIP")
#' @param options A list containing different specifications, which vary
#' depending on the model. For stopover and mixture models this contains B (the
#' number of broods), shared_sigma (boolean denoting if the SDs are the same for
#' each component), mu_formula (specifying a formula which describes a covariate
#' dependency for the mean arrival times for each brood), and sd_formula
#' (similar, for the SD of each brood).
#' @param DF If covariate relationships are specified, this is the data.frame
#' which contains these covariate values for each observation.
#' @return A named list containing 'skeleton', the skeleton of parameter values
#' used by the model to relist optim guesses, and 'DMs', the list of design
#' matrix objects which are required to obtain the correct parameter values 
#' in each case where the user has specified a formula
#' @export
produce_skeleton <- function(a_choice = "mixture", distribution = "P",
                             options = list(), DF = NULL){
  # purpose : Takes a variety of options for model specifications and determines
  #           the skeleton structure which is required for the model parameters.
  # inputs  : a_choice     - The character name for the seasonal trend function
  #                          being fitted. options include "splines",
  #                          "stopover", and "mixture".
  #           distribution - The distribution of the counts, Poisson ("P"),
  #                          Negative Binomial ("NB") or Zero-Inflated Poisson
  #                          ("ZIP")
  #           options      - A list containing different specifications,
  #                          which vary depending on the model. For stopover and
  #                          mixture models this contains B (the number of
  #                          broods), shared_sigma (boolean denoting if the SDs
  #                          are the same for each component), mu_formula
  #                          (specifying a formula which describes a covariate
  #                          dependency for the mean arrival times for each 
  #                          brood), and sd_formula (similar, for the SD of each
  #                          brood).
  #           DF           - If covariate relationships are specified, this is 
  #                          the data.frame which contains them.
  # output  : A named list of vectors containing NAs, which specificies the
  #           skeleton used by the nll function to relist() the vector of
  #           parameters given by optim().
  
  # input checks:
  if (! a_choice %in% c("splines", "stopover", "mixture"))
    stop("Invalid a_choice. Select 'splines', 'stopover' or 'mixture'.")
  
  if (! distribution %in% c("P", "ZIP", "NB"))
    stop("Invalid distribution choice. Select 'P', 'ZIP' or 'NB'.")
  
  skeleton <- list()
  DMs <- list()
  
  if (a_choice == "splines"){
    # determine the degrees of freedom of the splines (default of 15):
    df <- ifelse(is.null(options$df), 15, options$df)
    skeleton$spline_par <- rep(NA, df)
  }
  
  else if (a_choice %in% c("mixture", "stopover")){
    # determine number of broods, defaults to 1:
    B <- ifelse(is.null(options$B), 1, options$B)
    
    # determine the number of sigmas, defaults to 1:
    sigma_num <- ifelse(test = is.null(options$shared_sigma),
                        yes = 1,
                        no = B*!options$shared_sigma) %>% c(1) %>% max
    
    # Now make the skeleton:
    skeleton$mu <- rep(NA, B)
    skeleton$sigma <- rep(NA, sigma_num)
    if (B > 1) skeleton$w <- rep(NA, B - 1)
    
    # Update the design matrices and skeletons for the mean parameters:
    updates <- get_cov_skeleton_and_DM("mu", options, DF, B, skeleton, DMs)
    skeleton <- updates$skeleton ; DMs <- updates$DMs
    
    # Update the design matrices and skeletons for the sigma parameters:
    updates <- get_cov_skeleton_and_DM("sigma", options, DF,
                                       sigma_num, skeleton, DMs)
    skeleton <- updates$skeleton ; DMs <- updates$DMs
    
    # Update the design matrices and skeletons for the w parameters:
    if (B > 1){
      updates <- get_cov_skeleton_and_DM("w", options, DF, B - 1, skeleton, DMs)
      skeleton <- updates$skeleton ; DMs <- updates$DMs
    }
    
    # For the stopover model, add the retention prob:
    if (a_choice == "stopover") skeleton$phi <- NA 
  }
  
  # If the model is ZIP or NB, we need an extra distributional parameter:
  if (distribution != "P") skeleton$dist.par <- NA
  
  return(list(skeleton = skeleton, DMs = DMs))
}

#' Profile (concentrated) likelihood evaluation for GAI
#' 
#' A function which evaluations the profile negative log likelihood of our
#' observed count data, given the user model specifications of a GAI model.
#' @param par A vector of estimated parameters for the model
#' @param obs A matrix of observations. Rows for sites and columns for
#' occasions.
#' @param skeleton The skeleton which tells this function how to \code{relist}
#' the estimated parameters.
#' @param a_choice The character name of the method used to estimate the
#' seasonal trend. Either "mixture", "stopover" or "splines".
#' @param dist_choice The distribution for the observed counts, can be "P"
#' (poisson), "NB" (negative binomial), "ZIP" (zero inflated poisson)
#' @param spline_specs A list containing the df (degrees of freedom) and the
#' 'degree' for the splines.
#' @param N optionally, instead of using the profile likelihood approach to
#' estimate N (the superpopulation at each site), it can be supplied as a vector
#' or integer.
#' @param returnN if TRUE, returns the vector estimate of N and the matrix
#' a_func instead of the LL
#' @param DMs A list of design matrices if covariates have been included
#' @param returnA If TRUE, will also include the flight pattern densities
#' @export
profile_ll <- function(par, obs, skeleton, a_choice = "mixture", 
                       dist_choice = "P",
                       spline_specs = list(df = 10, degree = 3),
                       N = NULL, returnN = F, DMs = list(),
                       returnA = F){
  # purpose : Evaluates the likelihood of a poisson model, given a choice of
  #           seasonal trend function.
  # inputs  : par          - A vector of estimated parameters for the model
  #           obs          - A matrix of observations. Rows for sites and
  #                          columns for occasions.
  #           a_choice     - The character name of the method used to estimate
  #                          the seasonal trend. Either "mixture", "stopover" or 
  #                          "splines".
  #           dist_choice  - The distribution for the observed counts, can be
  #                          "P" (poisson), "NB" (negative binomial),
  #                          "ZIP" (zero inflated poisson)
  #           skeleton     - The skeleton which tells this function how to
  #                          relist() the estimated parameters.
  #           nS           - The number of sites in the study.
  #           spline_specs - A list containing the df (degrees of freedom) and
  #                          the degrees for the splines.
  #           N            - optionally, instead of using the profile likelihood
  #                          approach to estimate N (the superpopulation at each
  #                          site), it can be supplied as a vector or integer.
  #           returnN      - if TRUE, returns the vector estimate of N and the
  #                          matrix a_func instead of the LL
  #           DMs          - A list of design matrices if covariates have been
  #                          included
  #           returnA      - If TRUE, will return the A matrix as soon as is 
  #                         possible, without returning anything else
  # output  : A real number, the neg log-likelihood of all the data given the 
  #           model specification and parameters
  
  nS <- nrow(obs) ; nT <- ncol(obs)
  breaks <- c(-Inf, 1:(nT - 1), Inf)
  if (dist_choice %in% c("NB", "ZIP")) dist_p <- relist(par, skeleton)$dist.par
  
  # Get the seasonal flight pattern:
  a_func <- switch(a_choice,
                   splines = evaluate_splines(par,skeleton,spline_specs,nT,nS),
                   stopover = evaluate_stopover(breaks,par,skeleton,DMs,nS),
                   mixture = mixture_nll(par, matrix(1:nT,nrow = nS,ncol = nT,
                                                     byrow = T),skeleton,F,DMs))
  y_dot <- apply(obs, 1, sum, na.rm = T)
  
  # The splines function returns a vector, since no covariates are allowed,
  # so we convert the output to a matrix, to match "obs":
  if (a_choice == "splines") a_func %<>% matrix(nrow = nS, ncol = nT, byrow = T)
  a_funcNA <- a_func
  a_funcNA[is.na(obs)] <- NA
  
  # Get N from the specified argument or with a profile likelihood approach:
  if (is.null(N)) N <- y_dot / apply(a_funcNA, 1, sum, na.rm = T)
  else if (length(N) == 1) N %<>% rep(nS)
  lambda <- a_func * rep(N, nT)
  
  if (returnA) return(a_func)
  
  else{
  
    # Evaluate the negative log likelihood:
    lik <- switch(dist_choice,
                  P = dpois(obs, lambda, log = T),
                  NB = dnbinom(obs, mu = lambda, size = exp(dist_p), log = T),
                  ZIP = dzip(obs, lambda = lambda, psi = plogis(dist_p),
                             log = T))
    
    # Return the negative log lik for all points:
    if (returnN) return(list(N = N, a_func = a_func))
    else lik %>% sum(na.rm = T) %>% `*`(-1) %>% return
  }
}

check_GAI_inputs <- function(start, obs, a_choice, dist_choice, options,
                             tol, maxiter){
  # purpose : Performs sanity checks on the fit_GAI function inputs, to avoid
  #           cluttering the aforementioned function with if statements
  # inputs  : see fit_GAI
  # outputs : NULL, will cause an error if any checks fail
  if (a_choice == "splines" & (is.null(options$df) | is.null(options$degree)))
    stop("degree and df must be specified in options for a spline fit")
  
  if (! a_choice %in% c("mixture", "stopover", "splines"))
    stop("a_choice must be 'mixture', 'stopover', or 'splines'")
  
  if (! dist_choice %in% c("P", "NB", "ZIP"))
    stop("dist_choice must be 'P', 'NB' or 'ZIP'")
  
  if (class(maxiter) != "numeric" | maxiter <= 0)
    stop("maxiter must be a positive integer")
  
  if (class(options) != "list") stop("options must be a list")
  if (class(tol) != "numeric" | tol <= 0) stop("tol must be a positive number")
  if (class(obs) != "data.frame")
    stop("observations must be supplied as a data.frame")
  if (is.null(obs$site) | is.null(obs$occasion) | is.null(obs$count))
    stop("DF must contain columns 'site', 'occasion' and 'count'.")
}

#' GAI model fitting
#' 
#' Fits the GAI model to a set of data, with a choice of normal, stopover or
#' spline fits for the seasonal arrival times, as well as poisson, zero-inflated
#' poisson, or negative binomially distribted count observations.
#' @param start The vector of start guesses for the parameters.
#' @param DF The data.frame of observations. Should contain columns 'site', 
#' 'occasion' and 'count'. occasion and count should be integers, and missing 
#' entries in the count column are allowed to take the value NA. Any covariate
#' which is referred to in a formula specified in the options argument should
#' be present as a column in DF.
#' @param a_choice Choice of flight path function, which describes the seasonal
#' variation in counts over the monitoring period. Possible options are
#' "mixture": A normal mixture, with as many components as desired, will be used
#' to model the seasonal variation in counts.
#' "stopover": The Matechou et al. (2014) stopver model for seasonal variation 
#' will be used.
#' "splines": A general additive model, using splines, will model the seasonal 
#' variation. The degree and degrees of freedom of the splines can be specified
#' in the 'options' argument.
#' @param dist_choice The distribution for observed counts of individuals.
#' "P" indicates a poisson distribution, "NB" a, negative binomial, and "ZIP", 
#' the zero-inflated poisson distribution.
#' @param options A list containing different specifications, which vary
#' depending on the model. For stopover and mixture models this contains B (the
#' number of broods), shared_sigma (boolean denoting if the SDs are the same for
#' each component), mu_formula (specifying a formula which describes a covariate
#' dependency for the mean arrival times for each brood), and sd_formula
#' (similar, for the SD of each brood). When splines are being used, options 
#' 'degree' and 'DF' can be included to specify the degree of the polynomial used
#' by the splines, and the total degrees of freedom of the seasonal flight path
#' model, respectively.
#' @param tol The tolerance for the NB and ZIP iterative solvers. The solver
#' stops when the difference between the negative log likelihood from one
#' iteration to the next falls beneath this threshold.
#' @param maxiter The maximum number of iterations to perform when iterating for
#' the NB or ZIP likelihoods
#' @param verbose If TRUE, the function will print various updates while it
#' iterates for NB and for ZIP models.
#' @param hessian if TRUE, refits the model one last time with hessian = TRUE in
#' optim, usually so that this can be used to conduct a bootstrap
#' @param bootstrap Always NULL for the user.. Internally, bootstraps use
#' this argument to pass in sanitised precalculations to speed up the refitting
#' process of 'refit-the-model' bootstraps
#' @param method A character for the method which should be used to find MLEs.
#' See optim documentation for options and further detail. SANN (simulated 
#' annealing) can be a good method of fine-tuning starting points, if one is 
#' afraid their MLEs are caught in a local maximum.
#' @return The output of the `fit_GAI` function is a list with a few important
#' elements. The `par` element gives named estimates of the MLEs for the model
#' parameters, with `value` giving the value of the negative log likelihood
#' evaluated at the MLEs. `counts`, `convergence`, `message` and `hessian`
#' are all standard outputs given by the `optim` function with details on the
#' numerical process of estimating the MLEs. `spline_specs` contains the
#' user-specified options for fitting splines, and will be an empty list for
#' mixture and stopover models. `dist_choice` and `a_choice` contain the count
#' distribution and flight path distributions, respectively. `skeleton` contains
#'  the skeleton list of parameter values that the package uses to fit the
#'  model, with `options` being the list of options passed to `fit_GAI`.
#'  `maxiter` refers to the maximum number of iterations used to estimate the
#'  MLEs of ZIP and NB models. `A` and `N` contain the matrix of estimated
#'  seasonal densities at each site and occasion, and the vector of estimated
#'  site totals, respectively. `DMs` contains the list of design matrices used
#'  by the rGAI package to include the selected covariate formulas in the model.
#'  `obs` contains the count observations in matrix form, with sites as rows and
#'  occasions as columns. This is the same format as for `A`. `DF` contains the
#'  original `data.frame` supplied to `fit_GAI`, and finally, `tol`
#'  specififies the stopping condition used for the model (an epsilon such that 
#'  during an iterative process for fitting a ZIP or NB model, a difference of
#'  less than epsilon in the negative log likelihood between two iterations
#'  causes the process to terminate).
#' @examples
#' fit_GAI(rep(0, 3), example_data, a_choice = "stopover")
#' fit_GAI(rep(0, 6), example_data, a_choice = "mixture", options = list(B = 3))
#' fit_GAI(rep(0, 20), example_data, a_choice = "splines", options = list(df = 20, degree = 3))
#' fit_GAI(rep(0, 5), example_data, a_choice = "mixture", options = list(B = 2, shared_sigma = F))
#' @export
fit_GAI <- function(start, DF, a_choice = "mixture", dist_choice = "P",
                    options = list(), tol = 1e-3, maxiter = 1e3,
                    verbose = F, hessian = F, bootstrap = NULL,
                    method = "Nelder-Mead"){
  # purpose : Fits a GAI with either spline, mixture or stopover seasonal flight
  #           curves and either P, NB or ZIP distributed counts of individuals
  # inputs  : start       - The vector of start guesses for the parameters.
  #           a_choice    - "mixture", "stopover" or "splines".
  #           dist_choice - "P" (poisson), "NB" (negative binomial) or "ZIP"
  #                         (zero-inflated poisson).
  #           options     - Specifies some options for the model.
  #                         See produce_skeleton().
  #           tol         - The tolerance for the NB and ZIP iterative solvers.
  #                         The solver stops when the difference between the 
  #                         negative log likelihood from one iteration to the 
  #                         next falls beneath this threshold.
  #           maxiter     - The maximum number of iterations to perform when 
  #                         iterative for the NB or ZIP likelihoods.
  #           verbose     - If TRUE, the function will print various updates
  #                         while it iterates for NB and for ZIP models.
  #           hessian     - if TRUE, refits the model one last time with
  #                         hessian = TRUE in optim, usually so that this can
  #                         be used to conduct a bootstrap.
  #           skeleton    - Used by the bootstrap function to pass in 
  #                         arguments to avoid unnecessary calculations
  #           bootstrap   - If is not null, contains a set of pre-calculated
  #                         and sanitised objects.
  #           checkInpts  - If TRUE, will sanity check the input arguments
  bootFLAG <- !is.null(bootstrap)
  
  if (!bootFLAG){
    # We take a few shortcuts during a bootstrap to save on useless computation.
    check_GAI_inputs(start, DF, a_choice, dist_choice, options, tol, maxiter)
    extracted_counts <- extract_counts(DF)
    obs <- extracted_counts$matrix
    DF <- extracted_counts$DF

    # Produce a skeleton so the ll knows how to deal with the optim par vector:
    skeleton <- produce_skeleton(a_choice = a_choice,
                                 distribution = dist_choice,
                                 options = options, DF = DF)
    
    DMs <- skeleton$DMs ; skeleton <- skeleton$skeleton
  
    # Check the user's inputs match the expected number of parameters:
    if (length(start) != (skeleton %>% unlist %>% length)){
      paste("Incorrect number of parameters supplied.",
            "Please submit something similar to:") %>% print
      skeleton %>% unlist %>% print
      stop("Starting parameter value dimension mismatch")
    }
  
    # So that the output vector is more readable for the user:
    names(start) <- names(unlist(skeleton))
  
    spline_specs <- switch(a_choice, list(),
                           splines = with(options, list(degree = degree, 
                                                        df = df)))
  }#if
  
  else{
    # Just steal the precalculated values that are (mostly) constant for each
    # bootstrap iteration (other than the observation matrix).
    spline_specs <- bootstrap$spline_specs
    DMs <- bootstrap$DMs
    skeleton <- bootstrap$skeleton
    obs <- bootstrap$obs
  }
  
  # Get the first fit with optim:
  hessian_ff <- hessian & dist_choice == "P"
  fit <- optim(profile_ll, par = start, obs = obs, skeleton = skeleton,
               a_choice = a_choice, dist_choice = dist_choice,
               spline_specs = spline_specs, DMs = DMs, hessian = hessian_ff,
               method = method)
  
  # To ensure a bootstrap call has the estimated site totals:
  if (dist_choice == "P"){
    N_A <- profile_ll(fit$par, obs = obs, skeleton = skeleton,
                      a_choice = a_choice, dist_choice = dist_choice,
                      spline_specs = spline_specs, returnN = T, DMs = DMs)
    N <- N_A$N
  }
  
  if (dist_choice != "P") { # For non-P models, we solve iteratively
    # Get the initial profile likelihood estimate of N:
    N_A <- profile_ll(fit$par, obs = obs, skeleton = skeleton,
                      a_choice = a_choice, dist_choice = dist_choice,
                      spline_specs = spline_specs, returnN = T, DMs = DMs)
    
    N <- N_A$N
    if (length(N) == 1) N %<>% rep(nS)
    
    # The iterative process:
    fit <- iterate_GAI(N, fit, obs, skeleton, a_choice, dist_choice,
                       spline_specs, tol, maxiter, DMs, verbose, hessian,
                       method)
    
    # Add the GAI, seasonal component and estimated density:
    N_A <- profile_ll(fit$par, obs = obs, skeleton = skeleton,
                      a_choice = a_choice, dist_choice = dist_choice,
                      spline_specs = spline_specs, returnN = T, DMs = DMs)
  }
  
  if (bootFLAG) return(list(A = N_A$a_func, N = N, mle = fit$par))
  
  else{
    site_names <- DF$site %>% unique %>% names
    names(fit$N) <- site_names
    
    # Add the skeleton, options, tol, distribution choices, maxiter and 
    # data.frame for a clean bootstrap implementation:
    fit$spline_specs <- spline_specs
    fit$dist_choice <- dist_choice
    fit$a_choice <- a_choice
    fit$skeleton <- skeleton
    class(fit$skeleton) <- c(class(fit$skeleton), "GAIskel")
    fit$maxiter <- maxiter
    fit$options <- options
    fit$A <- N_A$a_func
    fit$N <- N_A$N
    fit$DMs <- DMs
    fit$obs <- obs
    fit$tol <- tol
    fit$DF <- DF
    
    # Make the output an S3 class, to allow for AIC and plotting general
    # methods to be implemented for user ease of use:
    class(fit) <- "GAI"
    return(fit)
  }#else
}#fit_gai

iterate_GAI <- function(N, fit, obs, skeleton, a_choice, dist_choice,
                        spline_specs, tol, maxiter, DMs, verbose, hessian,
                        method){
  # purpose : Performs the iterative process for the NB and ZIP models, given
  #           an initial fit
  # inputs  : see fit_gai
  # outputs : an optim output for the final fit
  
  val <- fit$value
  prev_val <- val + 2*tol
  iternum <- 1
  
  while (prev_val - val > tol & iternum < maxiter){
    
    if (verbose) cat("negative log likelihood for iteration", iternum, "is",
                     fit$value, "\n")
    
    # update prev_val and iternum:
    prev_val <- fit$value
    iternum <- iternum + 1
    
    # update distributional parameter and a_func:
    N_A <- profile_ll(fit$par, obs = obs, skeleton = skeleton,
                      a_choice = a_choice, dist_choice = dist_choice,
                      spline_specs = spline_specs, returnN = T, DMs = DMs)
    
    a_func <- N_A$a_func
    dist_par <- relist(fit$par, skeleton)$dist.par
    
    # update the site total estimates with a numerical solver:
    N <- nsolve_u(N, obs, a_func, dist_par, dist_choice)
    
    # refit the model using this version of N
    t.fit <- optim(fit$par, profile_ll, obs = obs, skeleton = skeleton,
                   a_choice = a_choice, dist_choice = dist_choice,
                   spline_specs = spline_specs, N = N, DMs = DMs,
                   method = method)
    
    if (t.fit$value > prev_val) break
    else fit <- t.fit
    
    # update value:
    val <- fit$value
  }
  
  if (verbose) cat("final negative log likelihood value is", fit$value, '\n')
  
  if (hessian){
    if (verbose) cat("performing final iteration to find hessian\n")
    fit <- optim(fit$par, profile_ll, obs = obs, skeleton = skeleton,
                 a_choice = a_choice, dist_choice = dist_choice,
                 spline_specs = spline_specs, N = N, DMs = DMs,
                 hessian = hessian, method = method)
  }
  
  return(fit)
}

#' Introduce NAs into matrices
#' 
#' Thin a matrix of count data by replacing each of its values by an NA
#' observation with a fixed probability
#' @param matvec The input matric or vector to be thinned
#' @param alpha The probability of a point being replaced
#' @return The object matvec, but with some elements replaced by NA values
#' @export
thin_NA <- function(matvec, alpha = 0.3){
  # purpose : Takes an input vector or matrix and replaces each value with NA
  #           with probability alpha
  # inputs  : matvec - The input matrix or vector to be thinned
  #           alpha  - The expected proportion of points to be removed
  # output  : matvec, with each point replaced by NA with probability NA
  
  to_replace <- rbinom(length(matvec), 1, alpha) %>% as.logical %>% which
  matvec[to_replace] <- NA
  return(matvec)
}

#' Backtransform for means
#' 
#' Performs the backtransform of the link function which is used to specify the
#' means of the broods
#' @param vec The vector of (positive) means which should be backtransformed to 
#' the link-scale space
#' @return A vector of the same length as vec, of backtransformed values
#' @export
means_backtransform <- function(vec){
  # purpose : Takes a vector of means and produces the vector of values which
  #           leads to those means when passed to means_link()
  # inputs  : vec - The vector of means to be backtransformed
  # outputs : A vector fo the same length as vec, of backtransformed values
  vec[-1] <- diff(vec)
  vec %>% log %>% return
}

#' Backtransform for weights
#' 
#' Performs the backtransform of the link function which is used to specify the
#' weights of the broods
#' @param vec The vector of weights between 0 and 1 which should be
#' backtransformed to the link-scale space
#' @return A vector of the same length as vec, of backtransformed values
#' @export
probs_backtransform <- function(vec, tol = 1e-10){
  # purpose : Takes a vector of weights and produces the vector of values which
  #           leads to those weights when passed to probs_link()
  # inputs  : vec - The vector of weights to be backtransformed
  #           tol - The tolerance for the maximum amount the sum of the input
  #                 weights is allowed to differ from 1.
  # outputs : A vector of the same length as vec, of backtransformed values
  n <- length(vec)
  output <- rep(NA, n)
  left <- 1
  
  for (i in 1:n){
    output[i] <- vec[i] / left
    left <- left - vec[i]
  }
  
  return(qlogis(output[1:(n - 1)]))
}


transform_values <- function(base, starting, skeleton){
  # purpose : Takes a given parameter type ("mu", "sigma" or "w") and performs
  #           the required checks to backtransform user supplied values onto the
  #           package's link scale
  # inputs  : base     - The character name of the parameter type "mu", "sigma"
  #                      or "w"
  #           starting - The list of starting values that the user gave
  #           skeleton - The skeleton for the specified model
  
  if (starting[[base]] %>% is.null %>% `!`){
    skel_vals <- skeleton[[base]]
    real_vals <- starting[[base]]
    lsv <- length(skel_vals)
    lrv <- length(real_vals)
    
    # For the weights, the user is expected to enter one more value than there
    # are final parameters:
    if (base == "w") lsv <- lsv + 1
    
    if (lsv != lrv){
      message <- paste(lsv, base, "values expected, but", lrv, "were given.",
                       "Therefore", base, "starting values are set to 0.")
      warning(message)
    }#if
    
    else{
      starting[[base]] <- switch(base,
                                 mu = means_backtransform(starting[[base]]),
                                 w = probs_backtransform(starting[[base]]),
                                 sigma = log(starting[[base]]))
                                 
    }#else
  }#if (starting_values[[base]])
  return(starting)
}#transform_values_one_param

#' Transform starting values to the link scale
#' 
#' Takes in starting guesses for parameter values, and maps them to the real 
#' valued link scale. Even parameters already on the real scale (such as
#' means) are transformed, to ensure they are uniquely determined. Starting
#' values for splines and covariates are always set to 0, as are starting
#' values for parameters that aren't specified.
#' @param starting_values A named list including starting values for "mu",
#' "sigma", "w" and "dist.par" (in the case of a non-poisson model)
#' @param a_choice The choice of flight path distribution
#' @param dist_choice The choice of observed count distribution
#' @param options The list of model options that specify the number of broods, 
#' spline parameters, covariate formulas, etc.
#' @param DF A data.frame which must be included when covariate formulas are
#' specified, so that the correct number of parameters required by each 
#' formula can be determined. Calling extract_counts on this data.frame first
#' to sanitize it is recommended
#' @return A named list of starting parameter values, all on the link scale.
#' @export
transform_starting_values <- function(starting_values, a_choice, dist_choice,
                                      options, DF = list()){
  # purpose : Takes a named list of inputs and outputs the vector of parameters
  #           expected by fit_GAI, on the link scale
  # inputs  : A named list of parameter values on the real scale
  # output  : A named vector of parameter values on the real scale.
  skeleton <- produce_skeleton(a_choice, dist_choice, options, DF)$skeleton
  
  dist_guess <- starting_values[["dist.par"]]
  sigma_guess <- starting_values[["sigma"]]
  phi_guess <- starting_values[["phi"]]
  mu_guess <- starting_values[["mu"]]
  w_guess <- starting_values[["w"]]
  
  # Add a default guess of 0 for each parameter value, and then relist the
  # the skeleton for convenience:
  output <- rep(0, skeleton %>% unlist %>% length) %>% relist(skeleton)
  names(output) <- names(skeleton)
  
  # Add in the user guesses for "mu", "sigma", "w" and "dist.par" if they exist,
  # but only if the model deems them necessary:
  if (!is.null(w_guess) & !is.null(skeleton[["w"]]))  output[["w"]] <- w_guess
  if (!is.null(sigma_guess))  output[["sigma"]] <- sigma_guess
  if (!is.null(mu_guess))  output[["mu"]] <- mu_guess
  
  if (!is.null(phi_guess) & !is.null(skeleton[["phi"]]))
    output[["phi"]] <- phi_guess
  
  if (!is.null(dist_guess) & !is.null(skeleton[["dist.par"]])) 
    output[["dist.par"]] <- dist_guess
  
  # A lazy guess of all 0s for splines:
  if (a_choice == "splines") return(output)
  
  else{
    # For mixtures and stopovers, go through means, sigmas and ws and transform 
    # the values given by the user with the link functions. Covariate parameters
    # will be given a default of 0, as well as anything that hasn't been
    # specified:
    output %<>% transform_values(base = "mu", skeleton = skeleton)
    output %<>% transform_values(base = "sigma", skeleton = skeleton)
    
    # In the case of a single brood, there won't be any weights, so we avoid
    # triggering any warning messages this way:
    if (skeleton[["w"]] %>% is.null %>% `!`) output %<>%
      transform_values(base = "w", skeleton = skeleton)
    
    # Add in phi if it's needed
    if (a_choice == "stopover" & phi_guess %>% is.null %>% `!`)
      output[["phi"]] <- qlogis(phi_guess)
    
    if (dist_choice != "P"){
      if (starting_values[["dist.par"]] %>% is.null %>% `!`){
        # Make sure the distributional parameter is on the right scale for a 
        # negative binomial (log link) and zero-inflated poisson (logistic link)
        linkfunc <- switch(dist_choice, NB = log, ZIP = qlogis)
        output[["dist.par"]] <- linkfunc(starting_values[["dist.par"]])
      }#if starting_values
    }#if dist_choice
  }#else
  output %>% unlist %>% return
}

#' Extract a count matrix from a data.frame
#' 
#' Loops through unique values of 'occasion' and 'site' in a data.frame in order
#' to extract the matrix of observed counts required for GAI model fitting
#' @param data_frame A data.frame with columns "site", "occasion" and "counts"
#' @param checks if TRUE, checks that all counts are specified
#' @param returnDF if TRUE returns instead a list with the matrix of counts as
#' well as a reordered version of the data.frame so that the rows match the
#' matrix output.
#' @return A matrix of counts, with rows as sites and columns as occasions.
#' @examples extract_counts(example_data)
#' @export
extract_counts <- function(data_frame, checks = T, returnDF = T){
  # purpose : Extracts from a data_frame the matrix of observed counts
  # inputs  : data_frame - A data.frame with columns "site", "occasion" and
  #                        "counts".
  #           checks     - if TRUE, checks that all counts are specified
  #           returnDF   - if TRUE returns instead a list with the matrix of 
  #                        counts as well as a reordered version of the
  #                        data.frame so that the rows match the matrix output.
  # outputs : A matrix of counts, with rows as sites and columns as occasions.
  if (checks){
    if (is.null(data_frame$site)) stop("'site' column must be present")
    if (is.null(data_frame$occasion)) stop("'occasion' column must be present")
    if (is.null(data_frame$count)) stop("'count' column must be present")
    
    if (! class(data_frame$count) %in% c("numeric", "integer"))
      stop('count must be numeric')
    if (! class(data_frame$occasion) %in% c("numeric", "integer", "factor"))
      stop('occasion must be numeric or factor')
    if (! class(data_frame$site) %in% c("numeric", "integer", "factor"))
      stop('site must be numeric or factor')
  }
  
  # Extract number of sites:
  unique_sites <- unique(data_frame$site)
  nS <- length(unique_sites)
  
  # Extract number of occasions:
  unique_occasions <- unique(data_frame$occasion)
  nT <- length(unique_occasions)
  
  # Get the names of extra columns:
  covariate_names <- data_frame %>% colnames %>%
    setdiff(c("count", "site", "occasion"))
  
  # To avoid the assumption that sites and occasions start at the number 1, or
  # that the counts column in the data_frame are well-ordered in terms of 
  # sites and occasions, we loop through to extract the correct count for each
  # site occasion pair:
  output <- matrix(NA, nrow = nS, ncol = nT)
  outputDF <- data.frame()
  over <- under <- t_counter <- s_counter <- 0
  
  for (t in unique_occasions){
    # Increment occasion index, and reset site index:
    t_counter <- t_counter + 1
    s_counter <- 0
    
    for (s in unique_sites){
      # Increment site index:
      s_counter <- s_counter + 1
      
      # Extract the relevant portion of the data.frame:
      count_st_sub <- subset(data_frame, data_frame$site == s &
                               data_frame$occasion == t)
      # Extract the count:
      count_st <- count_st_sub$count
      
      # Accept the count if it's the only one specified:
      if (length(count_st) == 1){output[s_counter, t_counter] <- count_st}
      
      # If multiple counts specified for the same point, take the first:
      else if (length(count_st) > 1){
        output[s_counter, t_counter] <- count_st[1]
        over <- over + 1
      }
      
     
      else {
        # If no count is specified, mark it as missing, and add back in NA
        # values for all covariates:
        output[s_counter, t_counter] <- NA; under %<>% `+`(1)
        count_st_sub <- list(site = s, occasion = t, count = NA)
        for (cov in covariate_names) count_st_sub[[cov]] <- NA
        count_st_sub %<>% as.data.frame
      }
      
      outputDF %<>% rbind(count_st_sub[1, ])
    }
  }
  
  # Let the user know we dealt with unexpected cases:
  if (over > 0 | under > 0){
    warning(paste(over, "counts were specified multiples times.", under,
                  "were not specified."))
    
    if (under > 0 & returnDF)
      warning("Covariates for missing count entries set to NA")
  }
  
  if (returnDF) return(list(matrix = output, DF = outputDF))
  else return(output)
}

design_matrix <- function(DF, covar_formula){
  # purpose : Produces a design matrix, given a data.frame and the specified
  #           relationship given by a formula object
  # inputs  : DF            - The data.frame containing the covariates
  #           covar_formula - The formula specifying the covariate relationship
  # output  : A design matrix for the specified parameter value given the
  #           covariate values.
  # We must remove the LHS of the formula, if it exists (it usually does):
  RHS <- try(covar_formula[[3]], silent = T)
  
  if (class(RHS) != 'try-error'){ # This means there is a LHS as well as a RHS
    covar_formula[[2]] <- NULL    # This deletes the LHS and fixes the indices
  }
  
  # Make the design matrix for the whole data set:
  current_action <- options("na.action")
  options(na.action = "na.pass")
  DM <- model.matrix(object = covar_formula, data = DF)
  if (DM %>% is.na %>% any) stop("Covariate values cannot include NA")
  options(na.action = current_action)
  return(DM)
}

lin_predictor <-  function(parameters, DM){
  # purpose : Produces the linear predictor for a parameter, given the covariate
  #           values it depends on, as described by a design matrix
  # inputs  : parameters - The parameters for the covariate relationship
  #           DM         - The design matrix
  # output  : A vector of values for the parameter, given the design matrix
  #           and the hyper-parameter values
  
  # Multiply DM by the column matrix of covariate parameters to obtain 
  # the linear predictor for the theta parameter of interest:
  parameters <- matrix(parameters, ncol = 1)
  return(DM %*% parameters)
}

is_time_varying <- function(DM, DF){
  # purpose : determines if a design matrix contains time-varying covariates
  # inputs  : DM - The design matrix for the chosen covariate.
  #           DF - The data.frame with the data and covaraite values.
  # output  : TRUE or FALSE
  output <- FALSE
  
  for (s in unique(DF$site)){
    site_rows <- DF$site == s # find out which rows refer to this site
    # find out if any covariate column has time-varying values:
    s_check <- DM[site_rows, ] %>% apply(2, unique) %>% lapply(length) %>%
      as.numeric %>% `>`(1) %>% any
    
    if (s_check){output <- TRUE; break}
  }
  
  return(output)
}


#' Parameter value finding helper for covariate inclusion
#'
#'  A helper function which obtains the linear predictor for a given set of
#'  parameters ('mu', 'sigma') etc, by scanning the list of design matrices and
#'  the model skeleton for the relevant information
#' @param brood the parameter number (integer) within the sequence. If being
#' used by the package internally, this is almost always the brood number
#' @param is_indiv A boolean indicating if the parameter for this given brood is
#'  specified with a brood specific formula
#' @param is_general A boolean indicating if the parameter for this given is the
#'  same as that of the other broods
#' @param par A vector of parameter estimates, typically supplied by optim
#' @param base The character name of the type of parameter. "mu", "sigma", or
#' "w"
#' @param DMs The list of design matrices which contains all of the design
#' matrices required to calculated the value of each parameter given its
#' covariate values
#' @return A matrix of parameter values, to be used by
#' get_parameter_values_all_rows 
lp_func <- function(brood, is_indiv, is_general, par, base, DMs){
  
  is_covariate <- is_indiv | is_general
  if (!is_covariate) par[[base]] %>% `[`(brood) %>% return
  
  else{
    # This section deals with the column as a covariate that's either using
    # the same formula as the other broods, or its own one:
    extension <- ifelse(is_indiv, brood, "")
    index <- paste(base, extension, ".cov", sep = "")
    # The most horrific piece of code in the whole package. This pulls the right
    # parameters from the design matrix depending on if the formula is general
    # to all broods, or specific to the individual brood, and then packages them
    # in a list to be able to use do.call to apply the linear predictor within
    # the same line of code:
    par[[base]] %>% `[`(brood) %>% c(par[[index]]) %>%
      list %>% c(DMs[[paste(base, extension, sep = "")]] %>% list) %>%
      do.call(what = lin_predictor) %>% return
    # No longer the most horrific piece of code, since it turns out this
    # approach was very versatile for performing smilar tasks at other points
    # in the package, too.
  }#else
}#lp_func


#' Parameter value finding helper for covariate inclusion
#'
#'  A helper function which uses lp_func to get the values of parameters with 
#'  a covariate formula, for 
#'  all sites.
#' @param base The character name of the type of parameter. "mu", "sigma", or
#' "w"
#' @param par A vector of parameter estimates, typically supplied by optiim
#' @param DMs The list of design matrices which contains all of the design
#' matrices required to calculated the value of each parameter given its
#' covariate values
#' @param B The number of parameters for which the base needs to be used.
#' Typically the same as the number of broods, i.e. one mean arrival time per
#' brood,
#' @return A matrix of parameter values to be used by get_parameter_values
get_parameter_values_all_rows <- function(base, par, DMs, B){
  # purpose : Helper which extracts the correct information from parameter 
  #           and design matrix object to calculate the correct final 
  #           parameter value for a given site
  # inputs  : base - The string indicating which type of parameter should be 
  #                  dealth with (mu, w, sigma, etc)
  #           par  - The relisted object containing the optim estimates for 
  #                  all parameters
  #           DMs  - The list of all design matrices for objects with covariates
  #           B    - The number of parameters for which the base needs to be
  #                  used. Typically the same as the number of broods, i.e.
  #                  one mean arrival time per brood,
  #
  # Note: If none of the broods have a covariate formula, this will cause an 
  #       edge case which isn't dealt with elegantly, so it's important that
  #       the rest of the package code performs as expected (i.e doesn't call
  #       this function when there are no covariates)
  
  # Seems convoluted, but in essence, because the user can specify either a
  # set of mu covariate parameters (using the mean brood arrival time as an
  # example), or a set of mu1, mu2, etc. covariates individually, we need to
  # check if the brood index in question has no mu covariate formula, a general
  # one (which is the same for all broods) or an individual one. If it's
  # individually specified, then there'll by a mu1.cov design matrix (for
  # brood 1, or a mu2.cov design matrix for brood 2, etc.). To simplify:
  # we're just checking to see if the DMs object has an entry with a name that
  # implies we're dealing with an individually specified covariate:
  is_general <- DMs[[base]] %>% is.null %>% `!`
  
  is_indiv_cov <- function(x){
    DMs[[paste(base, x, sep = "")]] %>% is.null %>% `!` %>% return()
  }
  
  # cbind can handle combining columns with 1 entry to columns with multiple.
  # However, it only likes doing so when all of these columns are specified in 
  # the same call to cbind, or it spits out a dimension error. So we need to
  # store all our columns as entries in a list, and call cbind on this entire
  # object at the end: brood, is_indiv, is_general, par, base
  output <- lp_func(1, is_indiv_cov(1), is_general, par, base, DMs) %>% list
  
  if (B > 1){
    for (i in 2:B){
      output[[i]] <- lp_func(i, is_indiv_cov(i), is_general, par, base, DMs)
    }
  }
  
  output %>% do.call(what = cbind) %>% return()
}#get_parameter_values_all_rows


#' Determine if a given parameter uses covariates
#'
#'  Parses through a GAI model's skeleton and design matrix list to determine
#'  which type of covariate specification has been used (No covariates, the
#'  same formula for each brood, or brood-specific formulas)
#' @param base The character name of the type of parameter. "mu", "sigma", or
#' "w"
#' @param par A vector of parameter estimates, typically supplied by optim
#' @param DMs The list of design matrices which contains all of the design
#' matrices required to calculated the value of each parameter given its
#' covariate values
#' @return A boolean. TRUE if the parameter "base" uses covariates, and 
#' FALSE otherwise.
contains_covariates <- function(base, par, DMs){
  # purpose : Scans the (unlisted) parameters and design matrices objects to 
  #           determine if a type of parameter (means, standard deviations, 
  #           weights, etc.) in the GAI model contains either a general or
  #           brood specific covariate specification
  #
  # Note : See get_parameter_values_all_rows for more detailed commenting on how
  #        this sequence of pipes is being used
  parnum <- par[[base]] %>% length
  is_general <- DMs[[base]] %>% is.null %>% `!`
  is_indiv_spec <- sapply(1:parnum,
                          function(i) DMs[[paste(base, i, sep = "")]] %>%
                            is.null) %>% sapply(`!`) %>% any
  
  return(is_general | is_indiv_spec)
}

#' Obtain parameter values and handle covariate formulas
#'
#'  Uses a suite of helper functions to calculate linear predictors for 
#'  parameter values, with and without covariate inclusion, for stopover and 
#'  mixture models
#' @param par A vector of parameter estimates, typically supplied by optim
#' @param DMs The list of design matrices which contains all of the design
#' matrices required to calculated the value of each parameter given its
#' covariate values
#' @param skeleton A named list of parameters which can be used to relist the 
#' vector of parameter guesses provided by optim into the structure expected by
#' GAI package functions
#' @param n The integer number of data points for which these parameter values
#' should be found
#' @return A name dlist containing the Mu, Sigma and W parameters for the 
#' data points
get_parameter_values <- function(par, DMs, skeleton, n){
  # purpose : Returns a list of parameter values given design matrices and
  #           a vector of hyper-parameters given by optim()
  # inputs  : par      - The vector of parameters and hyper-parameters given by
  #                      optim()
  #           DMs      - The list of design matrices
  #           skeleton - The skeleton with which we can relist the optim()
  #                      parameter values
  #           n        - The number of data points = nS*nT
  # outputs : A list with the mean, sd, and w parameters with links applied.
  #           Each component is a matrix. The rows are a site occasion pair, and
  #           the columns refer to the value of the parameter for each brood.
  #           The parameters are given for the first time step first, and then 
  #           the second, etc. i.e, a survey with 3 sites, 5 occasions and 2 
  #           broods would have a 'means' output which is a 15x2 matrix.
  par %<>% relist(skeleton)
  probs <- 1 # For single brood case
  
  single_brood_mod <- function(x){
    if (x %>% dim %>% `==`(1) %>% any) return(t(x))
    else return(x)
  }
  
  # Use the helper function to get the parameter values for the means, 
  # SDs, and weights one by one:
  if (contains_covariates("mu", par, DMs)) means <-
    get_parameter_values_all_rows("mu", par, DMs, length(par$mu)) %>%
    apply(1, means_link) %>% t %>% single_brood_mod
  
  # Without covariates, do it manually to save computation time:
  else means <- means_link(par$mu) %>%
      matrix(nrow = n, ncol = length(par$mu), byrow = T)
  
  # Standard Deviations:
  if (contains_covariates("sigma", par, DMs)) sds <-
    get_parameter_values_all_rows("sigma", par, DMs, length(par$sigma)) %>% exp
  
  else sds <- exp(par$sigma) %>% matrix(nrow = n, ncol = length(par$sigma),
                                        byrow = T)
  
  # Weights:
  if (contains_covariates("w", par, DMs)) probs <-
    get_parameter_values_all_rows("w", par, DMs, length(par$w)) %>%
    apply(1, probs_link) %>% t
  
  else if (!is.null(par$w)) probs <- probs_link(par$w) %>%
    matrix(nrow = n, ncol = length(par$w) + 1, byrow = T)
  
  else probs %<>% matrix(nrow = n, ncol = 1, byrow = T)
  
  return(list(means = means, sds = sds, probs = probs))
}

bootstrap_combine <- function(sapply_output){
  # purpose : Combines the result of two or more sets of parallel bootstrap
  #           calculations for the butterfly model
  # note : The function expects each output to be a list containing three items,
  #        the output is a list of three items with rbinded outputs from all
  #        the arguments.
  argg <- sapply_output
  n <- length(argg)
  output <- argg[[1]]

  for (i in 2:n){
    for (j in 1:length(argg[[1]])){
      output[[j]] = rbind(output[[j]], argg[[i]][[j]])
    }
  }

  return(output)
}

bootstrap_non_parallel <- function(R, MLE, sigma, obs, skeleton, a_choice,
                                   dist_choice, DMs, spline_specs, parameters,
                                   N_vals, A_vals, mean_vals){
  # purpose : Computes the for loop for a parameter resampling bootstrap on a
  #           single core
  # inputs  : R        -
  #           MLE      - The maximum likelihood estimate for the parameters
  #           sigma    - The variance-covariance matrix for the MLE normal distn
  #           obs      - The observations
  #           skeleton -
  #           
  # output : list of parameters, N_vals and mean vals for the GAI resulting
  #           from the bootstrap
  for (i in 1:R){
    # Resample the parameters and get new estimates of the GAI:
    iteration_par <- rmvnorm(1, MLE, sigma)
    iteration_N_A <- profile_ll(iteration_par, obs = obs, skeleton = skeleton,
                                a_choice = a_choice,
                                dist_choice = dist_choice, DMs = DMs,
                                spline_specs = spline_specs, returnN = T)
    
    iteration_N <- iteration_N_A$N
    iteration_A <- iteration_N_A$a_func %>% as.vector
    iteration_mean <- mean(iteration_N)
    
    # Store the iteration results:
    parameters[i,] <- iteration_par
    mean_vals[i,] <- iteration_mean
    N_vals[i,] <- iteration_N
    A_vals[i,] <- iteration_A
  }
  
  return(list(parameters = parameters, N_vals = N_vals, A_vals = A_vals,
              mean_vals = mean_vals))
}

reformDMs <- function(DMs, n_sites, sites){
  # purpose : Reorders the rows of the fitted design matrix, so that they are
  #           in the correct order for the resampled sites selected for a 
  #           bootstrap
  # inputs  : DMs      - The design matrices object 
  #           n_sites  - The number of sites in the study
  #           sites    - The resampled site indices
  # output  : The DMs object, but reordered.
  
  n_dms <- length(DMs)

  # Time-varying covariates aren't allowed, so we can take some shortcuts:
  if(length(DMs) > 0){
    n_occasions <- DMs[[1]] %>% dim %>% `[`(1) %>% `/`(n_sites)
    for (i in 1:n_dms) DMs[[i]] <- DMs[[i]][rep(sites, n_occasions),]
  }
  
  return(DMs)
}

resample_bootstrap_one_iter <- function(sites, DMs, iteration_object, obs,
                                        maxiter, tol, a_choice, dist_choice,
                                        start){
  # purpose : Performs a single iteration of the resample the data bootstrap, 
  #           for cleaner parallel implementations with sapply and par Sapply
  # inputs  :
  #           sites            -  The number of sites in the study
  #           DMs              -  A list of design matrices for covariate 
  #                               inclusion
  #           iteration-object - An object which stores pre-calculation which
  #                              limit the work needed by the fit_gai function
  #                              during reffiting
  #           obs              - A matrix of observations of counts (or NAs)
  #           maxiter          - The maximum number of iterations to iterate a 
  #                              non-Poisson model 
  #           tol              - The stopping condition for non-Poisson fitting
  #           a_choice         - The choice of seasonal trend
  #           dist_choice      - The choice of count data distribution
  #           start            - The maximum likelihood estimate for the
  #                              parameters, to be used as the starting guess
  #                              for optim
  indices <- sample(1:sites, sites, replace = T)
  iteration_object$DMs <- reformDMs(DMs, sites, indices)
  iteration_object$obs <- obs[indices, ]
  
  fitted <- fit_GAI(start = start, DF = NULL, a_choice = a_choice,
                    dist_choice = dist_choice, options = NULL,
                    verbose = F, hessian = F,  bootstrap = iteration_object,
                    tol = tol, maxiter = maxiter)
  
  return(list(iteration_par = fitted$mle,
              iteration_N = fitted$N,
              iteration_A = fitted$A %>% as.vector,
              iteration_mean = mean(fitted$N)))
}

resample_bootstrap <- function(cl, R, start, obs, skeleton, a_choice,
                               dist_choice, DMs, spline_specs, tol, maxiter,
                               parameters, N_vals, A_vals, mean_vals){
  # purpose : Computes the non-parametric bootstrap by refitting the model at
  #           each iteration.
  # inputs  : R            - The number of bootstrap iterations to be performed
  #           start        - The maximum likelihood estimate for the parameters,
  #                          to be used as the starting guess for optim
  #           obs          - The observations
  #           skeleton     - The skeleton object which allows parameter vectors
  #                          to be relisted correctly
  #           a_choice     - The choice of seasonal arrival pattern
  #           dist_choice  - The choice of count distribution for the model
  #           DMs          - A list of design matrices for covariate inclusion
  #           options      - A set of user-specified options that indicate 
  #                          things such as number of broods or spline settings
  #           tol          - The tolerance for the therative fitting process
  #           maxiter      - The maximum number of iterations for the fitting
  #                          process
  #           parameters   - Object in which to store the bootstrapped 
  #                          parameter values
  #           N_vals       - Object in which to store the bootstrapped site
  #                          totals
  #           mean_vals    - Object in which to store the bootstrapped mean 
  #                          estimated site total across sites
  #           
  # output  : list of parameters, N_vals and mean vals for the GAI resulting
  #           from the bootstrap
  
  sites <- nrow(obs)
  
  iteration_object <- list(obs = NULL, skeleton = skeleton, DMs = DMs,
                           spline_specs = spline_specs)
  if (!is.null(cl)){
    parSapply(cl, rep(sites, R), resample_bootstrap_one_iter, DMs = DMs,
              iteration_object = iteration_object, obs = obs, maxiter = maxiter,
              tol = tol, a_choice = a_choice, dist_choice = dist_choice,
              start = start, simplify = F) %>% bootstrap_combine %>% return
  }
  
  else{
    sapply(rep(sites, R), resample_bootstrap_one_iter, DMs = DMs,
           iteration_object = iteration_object, obs = obs, maxiter = maxiter,
           tol = tol, a_choice = a_choice, dist_choice = dist_choice,
           start = start, simplify = F) %>% bootstrap_combine %>% return
    }
}#resample_bootstrap

transform_params <- function(param_vector, GAI_fit, use_all){
  # purpose : Transforms the parameter values of fitted GAI model by using a 
  #           set of covariate values which has been resampled from the
  #           dataset
  GAI_fit$par <- param_vector
  if (use_all) DF <- GAI_fit$DF # If refitting, use the data covariate values
  else DF <- GAI_fit$DF[sample(1:nrow(GAI_fit$DF), 1),] # else sample one set
  transform_output(GAI_fit, DF) %>% apply(2, mean) %>% return
}

#' Transforms bootstrap confidence intervals for custom covariate values.
#' 
#' Transforms bootstrap confidence intervals of parameter values to the
#' parameter scale for a range of custom covariate values. Typically used to 
#' provide information to produce plots of parameter estimates against 
#' covariate value, with confidence intervals
#' @param bootstrap_param_output A table of parameter estimates obtained by
#' using rGAI's bootstrap function with transform = F. Custom values produced
#' manually can also be passed.
#' @param GAI_fit An object produced by using fit_GAI for model fitting
#' @param covariate_data_frame A data.frame containing the range of covariate
#' values for which parameter scale transformations are desired.
#' @export
transform_bootstrap_parameters <- function(bootstrap_param_output, GAI_fit,
                                           covariate_data_frame){
  
  N <- nrow(bootstrap_param_output)
  output <- list()
  
  for (i in 1:N){
    GAI_fit$par <- bootstrap_param_output[i, ]
    output[[i]] <- transform_output(GAI_fit, covariate_data_frame)
  }  
  
  names(output) <- rownames(bootstrap_param_output)
  return(output)
}

#' Bootstrapping for GAI models
#' 
#' Produces either a non-parameteric bootstrap by refitting the GAI at each
#' iteration, or produces a parametric resampling of the MLEs at each iteration
#' using an estimate of their asymptotic normal distribution.
#' @param GAI_fit An object produced by using fit_GAI for model fitting
#' @param R The number of resamples to produce for the bootstrap
#' @param refit If TRUE, resamples the observations from the sites and occasions
#' uniformly, and refits the model at each bootsrap iteration. If FALSE, the
#' bootstrap simply resamples the fitted parameter values from their asymptotic
#' normal distribution (therefore this option requires a Hessian to have been
#' produced during the model fitting stage).
#' @param alpha 1 - alpha gives the coverage the bootstrap confidence intervals
#' aim to produce.
#' @param parallel if TRUE, calculates the bootstraps in parallel, using the
#' maximum number of available cores.
#' @param cores If not NULL, this specifies the number of cores to use, if the
#' default of using all available cores is undesirable.
#' @param transform if TRUE, will return a bootstrap on the transformed
#' parameters, rather than on the link scale. It should be noted that when 
#' covariates are present in the model, transformed outputs are averaged across
#' all covariate values in the data, to avoid erroneous interpretation of 
#' covariate effects on the transformed scale.
#' @return A named list with entries "EC", "N", "A" and "par" giving the
#' \code{c(alpha / 2, 1 - alpha / 2)} confidence interval for the expected
#' observed count at each site on each occasion, the estimated 
#' site super-population, the seasonal component, and the estimated
#' parameter values, respectively.
#' @export
bootstrap <- function(GAI_fit, R = 100, refit = T, alpha = 0.05, parallel = T,
                      cores = NULL, transform = T){
  # note    : The resampling the perameters bootstrap is not performed in
  #           parallel, since the calculation is so quick that combining results
  #           from multiple cores is typically slower than running the single
  #           threaded version.
  if (class(alpha) != "numeric" | alpha < 0 | alpha > 1) stop("Invalid alpha")
  if (class(GAI_fit) != "GAI") stop("GAI_fit must be produced by fit_GAI()")
  if (class(parallel) != "logical") stop("Parallel must be TRUE or FALSE")
  if (class(refit) != "logical") stop("Refit must be TRUE or FALSE")
  if (switch (class(cores), NULL = F, integer = cores < 1 | cores %% 1 != 0,
      numeric = cores < 1 | cores %% 1 != 0, T)) stop("Invalid number of cores")
  if (! class(R) %in% c("numeric", "integer") | R < 1 | R %% 1 != 0)
    stop("Invalid bootstrap iteration number")
  
  require(mvtnorm)
  require(parallel)
  
  # Register the cluster for parallelising:
  if (parallel){
    if (is.null(cores)) cores <- max(1, detectCores() - 1)
    cl <- makeCluster(cores)
  }
  
  else {cl <- NULL}
  
  # Extract the original fit values:
  MLE <- GAI_fit$par
  
  # Create matrices and arrays to store outputs:
  dims <- dim(GAI_fit$obs)
  nS <- dims[1] ; nT <- dims[2]
  parameters <- matrix(NA, nrow = R, ncol = length(MLE))
  colnames(parameters) <- names(MLE)
  
  A_vals <- matrix(NA, nrow = R, ncol = nS * nT)
  mean_vals <- matrix(NA, nrow = R, ncol = 1)
  N_vals <- matrix(NA, nrow = R, ncol = nS)
  
  # Extract required objects:
  spline_specs <- GAI_fit$spline_specs
  dist_choice <- GAI_fit$dist_choice
  sigma <- solve(GAI_fit$hessian)
  skeleton <- GAI_fit$skeleton
  a_choice <- GAI_fit$a_choice
  maxiter <- GAI_fit$maxiter
  obs <- GAI_fit$obs
  tol <- GAI_fit$tol
  DMs <- GAI_fit$DMs
  
  if (refit){
    # Resample data bootstrap:
    result <- resample_bootstrap(cl, R, MLE, obs, skeleton, a_choice,
                                 dist_choice, DMs, spline_specs, tol,
                                 maxiter, parameters, N_vals, A_vals,
                                 mean_vals)
    }
  
  else {
    # Resample the parameters bootstrap
    if (GAI_fit$hessian %>% is.null) stop("Hessian required in fitted model")
    result <- bootstrap_non_parallel(R, MLE, sigma, obs, skeleton, a_choice,
                                     dist_choice, DMs, spline_specs,
                                     parameters, N_vals, A_vals, mean_vals)
  }
  
  parameters <- result[[1]]
  N_vals <- result[[2]]
  A_vals <- result[[3]]
  mean_vals <- result[[4]]
  
  if (transform){
    parameters %<>% apply(1, transform_params, GAI_fit = GAI_fit,
                          use_all = refit) %>% t
    
    # Remove w1 from parameters if not present in GAI_fit$par (i.e. single brood case)
    if(!"w1" %in% names(GAI_fit$par))
      parameters <- parameters[, !colnames(parameters) =="w1"]
    
    # Assume the transformation has removed parameters that are only used to 
    # calculate the effect of a covariate:
    colnames(parameters) <- names(GAI_fit$par) %>% {.[!grepl(".cov", .)]}
  }
  
  # Calculate the confidence intervals using the results:
  probs <- c(alpha / 2, 1 - alpha / 2)
  par_CI <- apply(parameters, 2, quantile, probs = probs)
  
  A_vals %<>% unlist %>% array(dim = c(R, nS, nT), 
                               dimnames = list(
                                 paste0("R_", 1:R),
                                 paste0("Site_", 1:nS),
                                 paste0("Occasion_", 1:nT))
                               )
  
  N_vals %<>% unlist %>% array(dim = c(R, nS, nT), 
                               dimnames = list(
                                 paste0("R_", 1:R),
                                 paste0("Site_", 1:nS),
                                 paste0("Occasion_", 1:nT))
  )
  
  expected_count <- A_vals * N_vals
  N_vals <- N_vals[1:R, 1:nS, 1]
  expected_count_CI <- apply(expected_count, c(2, 3), quantile, probs = probs)
  A_CI <- apply(A_vals, c(2, 3), quantile, probs = probs)
  N_CI <- apply(N_vals, 2, quantile, probs = probs)
  
  if (parallel) stopCluster(cl)
  
  return(list(
    EC = expected_count_CI,
    N = N_CI,
    A = A_CI,
    par = par_CI,
    EC_raw = expected_count,
    N_raw = N_vals,
    A_raw = matrix(A_vals, nrow = nS),
    par_raw = parameters))
}


#' Log likelihood extraction
#' 
#' Extracts the (postive) log likelihood value from an output from the fit_GAI
#' function, notably for use by the AIC function
#' @param GAIobj An object produced by using <- for model fitting
#' @return An object object of class "logLik" with an attribute "df" equal to 
#' the number of fitted parameters for the model at hand.
#' @export
logLik.GAI <- function(GAIobj){
  # purpose : A simple function which extracts the number of fitted parameters 
  #           and log likelihood of a fitted GAI model, so that the AIC generic
  #           function can be used
  # inputs  : GAIobj - A fitted GAI model output of class "GAI"
  # output  : An object of class "logLik"
  df <- length(GAIobj$par)
  ll <- -GAIobj$value
  attr(ll, "df") <- df 
  class(ll) <- "logLik"
  return(ll)
}

#' GAI model summary
#' 
#' Produces a set of model summaries for fitted GAI models. Includes the AIC, 
#' the MLEs and the average estmated site total across sites.
#' @param GAIobj An object produced by using fit_GAI for model fitting
#' @return An object object of class "summary.GAI"
#' @export
summary.GAI <- function(GAIobj){
  # purpose : Produces a table of summary outputs for objects fitted with the
  #           the fit_GAI function.
  # inputs  : GAIobj - A fitted GAI model output of class "GAI"
  # output  : A "summary.GAI" object
  output <- list()
  output$MLE <- GAIobj$par
  
  if (!is.null(GAIobj$hessian)){
    output$MLE.SE <- GAIobj$hessian %>% solve %>% diag %>% sqrt
  }
  
  output$AIC <- AIC(GAIobj)
  output$N <-GAIobj$N
  
  class(output) <- "summary.GAI"
  return(output)
}

#' GAI model summary printing
#' 
#' Prints the output of a GAI model summary, to display the information to the
#' end user.
#' @param obj An object produced by using the summary.GAI function an a model
#' output produce by fit_GAI.
#' @return An object object of class "summary.GAI"
#' @export
print.summary.GAI <- function(obj){
  # purpose : Prints the output generated from the summary.GAI function
  cat("Maximum Likelihood Estimates (MLEs):\n\n")
  obj$MLE %>% signif(4) %>% print
  
  if(obj$MLE.SE %>% is.null %>% `!`){
    cat("\n\n")
    cat("MLE Standard Errors:\n\n")
    obj$MLE.SE %>% signif(4) %>% print
  }
  
  cat("\n\nAkaike's Information Criterion (AIC):", obj$AIC, "\n\n")
  cat("Average estimated site super-population size:",
      obj$N %>% na.omit %>% mean,"\n")
}

#' GAI flight path plotting
#' 
#' Prints the flight path curves of a fitted GAI model with options to scale
#' by estimated site total, and avoid plotting all sites by smoothing through
#' results for each sites by taking a given quantile at each time point
#' @param GAIobj An object produced by using the summary.GAI function an a model
#' output produce by fit_GAI.
#' @param all_sites If TRUE, the curves for every single site are plotted rather
#'  than a smooth through the curves from each site. Defaults to FALSE.
#' @param quantiles For when all_sites = F, a vector of quantiles that will each
#' be used to produce a curve. Defaults to c(0.05, 0.5, 0.95)
#' @param scale_by_N If TRUE, then the seasonal flight component for each site
#' will be scaled by that site's estimated site total. Defaults to TRUE.
#' @param colours A list of colour values to be used for plotting different
#' curves. Defaults to the GGplot default 10 colour palette repeated 10 times to
#' accommodate up to 100 sites (even when ggplot is not used)
#' @param shift An optional parameter which can be used to shift the location
#' of the key when using base graphics, as this can often be misplaced.
#' @param useGGplot If TRUE, will use GGplot to produce the graphics instead of
#' the standard base R graphics
#' @return NULL
#' @examples 
#' # A plot with default settings:
#' fit_GAI(c(0, 0), example_data, "mixture") %>% plot
#' 
#' # A plot with custom quantiles:
#' fit_GAI(c(0, 0), example_data, "mixture") %>% plot(quantiles = c(0.01, 0.99))
#' 
#' # A plot with quantities not scaled by site totals:
#' fit_GAI(c(0, 0), example_data, "mixture") %>% plot(scale_by_N = F)
#' @export
plot.GAI <- function(GAIobj, all_sites = F, quantiles = c(0.05, 0.5, 0.95),
                     scale_by_N = T,
                     colours = hcl(h = seq(15, 375, length = 11),
                                   l = 65, c = 100)[1:10] %>% rep(10)){
  # purpose : Produces simple plots of the flight path distribution of a GAI
  # inputs  :GAIobj     - The fitted GAI model object from the fit_GAI function
  #          all_sites  - If TRUE, the curves for every single site are plotted
  #                       rather than a smooth through the curves from each site
  #          quantiles  - For when all_sites = F, a vector of quantiles that 
  #                       will each be used to produce a curve
  #          scale_by_N - If TRUE, then the seasonal flight component for each
  #                       site will be scaled by that site's estimated site 
  #                       total
  #          colours    - A list of colour values to be used for
  #                       plotting different curves. Defaults to the GGplot 
  #                       default 10 colour palette repeated 10 times to 
  #                       accommodate up to 100 sites (even when ggplot is not
  #                       used)
  #          useGGplot  - If TRUE, will use the GGplot package to produce the
  #                       plots instead
  # output  : NULL
  A <- GAIobj$A

  # Change the data based on user options:
  if (scale_by_N) A %<>% `*`(rep(GAIobj$N, ncol(A)))
  if (!all_sites) A %<>% apply(2, quantile, probs = quantiles, na.rm = T)
  
  # Determine the correct plot labels:
  plot_lines <- nrow(A)
  ylab <- "Seasonal flight path"
  if (scale_by_N) ylab %<>% paste("scaled by site total")
  
  require("ggplot2")
  require("reshape2")
    
  # Determine the correct label names:
  toplot <- A %>% t %>% as.data.frame
    
  if (all_sites){
    lab <- "Site"
    colnames(toplot) <- 1:ncol(toplot)
  }
    
  else lab <- "Quantile"
    
  # Reformat the data for ggplot:
  occasions <- 1:ncol(A)
  mes_vars <- colnames(toplot)
  toplot$occasions <- occasions
  toplot %<>% melt(id.vars = "occasions", measure.vars = mes_vars)
    
  # Produce the plot:
  ggplot(toplot, aes(x = occasions, y = value)) +
    geom_line(aes(color = variable)) +
    scale_color_manual(values = c(colours)) + 
    ylab(ylab) + xlab("Occasion") +
    labs(color = lab)
}

#' Parameter transformation
#' 
#' Transforms parameters from the real scale on which they are estimated to the
#' correct parameter-space scale using the appropriate link functions and 
#' covariate formulas
#' @param GAIoutput The output of the fit_GAI function, a fitted model
#' @param DF A data.frame containing rows of covariate values for which
#' transformed parameter values should be obtained. It is important that the 
#' names of columns in this data.frame are identical to those that were
#' specified in the options argument of the fit_GAI function. A blank data.frame 
#' is provided by default, which is assumed in a no covariate model.
#' @param provide_A If TRUE, will also return the A matrix for the supplied. It
#' should be noted that the A matrix this function returns is not scaled by site
#' population totals, since these are unavailable for new data.
#' points
#' @return A named vector of parameters, on the parameter-space scale
#' @examples 
#' fit_GAI(c(0, 0), example_data, "mixture") %>% transform_output(DF = example_data[1, ])
#' @export 
transform_output <- function(GAIoutput, DF = data.frame(), provide_A = F){
  # purpose : Produces a matrix of parameter outputs, transformed to the 
  #           parameter-space scale, for mixture and stopover models.
  # inputs  : GAIoutput - An object produced by the fit_GAI function, containing
  #                       a fitted model
  #           DF        - A data.frame that contains a named list of covariate
  #                       values. If null, covariate values are assumed to be
  #                       0, and a warning is produced.
  output <- list()
  
  # Sometimes users give a data.frame for transformation which is identitcal 
  # to the dataset (i.e with site, count and occasion columns, still). To avoid
  # the automatic checks for time-varying covariates by the DM producing
  # function, these columns are removed:
  for (col in c("site", "count", "occasion")) DF[col] <- NULL
  
  # Extract model settings:
  dist_choice <- GAIoutput$dist_choice
  a_choice <- GAIoutput$a_choice
  nT <- ncol(GAIoutput$obs)
  
  # For splines, there are no link functions:
  if (a_choice == "splines") return(GAIoutput$par)
  
  # It's simpler to produce the design matrices again from scratch, 
  # using the supplied covariate values. If missing covariates are present,
  # the produce skeleton function should produce an error:
  DMs <- produce_skeleton(a_choice, dist_choice, GAIoutput$options, DF)
  skeleton <- DMs$skeleton ; DMs <- DMs$DMs
  par <- GAIoutput$par
  
  # Get the means, sds, and weights using the same function as the rest of the 
  # package:
  vals <- get_parameter_values(par, DMs, skeleton, max(c(nrow(DF), 1)))
  
  # When only one row is given in the data.frame, the outputs can sometimes
  # erroneously be given as transposed:
  if ((dim(vals$means) %>% diff %>% `<`(0)) & nrow(DF) == 1) vals$means %<>% t
  if ((dim(vals$sds) %>% diff %>% `<`(0)) & nrow(DF) == 1) vals$sds %<>% t
  if ((dim(vals$probs) %>% diff %>% `<`(0)) & nrow(DF) == 1) vals$probs %<>% t
  
  # For the means, we get the number of means from the skeleton, and then 
  # add the appropriate columns to the output data.frame one by one:
  n_mu <- length(skeleton$mu)
  n_si <- length(skeleton$sigma)
  n_w <- length(skeleton$w)
  for (i in 1:n_mu) output[[paste("mu", i, sep = "")]] <- vals$means[,i]
  for (i in 1:n_si) output[[paste("sigma", i, sep = "")]] <- vals$sds[,i]
  for (i in 1:n_w) output[[paste("w", i, sep = "")]] <- vals$probs[,i]
  
  # If the model is a stopover, add in the fixed phi probability:
  par %<>% relist(skeleton = skeleton)
  if (a_choice == "stopover") output$phi <- plogis(par$phi) 
  
  # Add in the distributional parameter, if there is one:
  if (dist_choice == "ZIP") output$dist.par <- plogis(par$dist.par)
  else if (dist_choice == "NB") output$dist.par <- exp(par$dist.par)
  output %<>% do.call(what = cbind)
  
  if (provide_A){
    nT <- ncol(GAIoutput$A); nS <- nrow(DF)
    A <- profile_ll(GAIoutput$par, matrix(1, nrow = nS, ncol = nT), skeleton,
                    a_choice, dist_choice, list(), 1, F, DMs, T)
    output <- list(params = output, A = A)      
  }
  
  return(output)
}

#' Probability link function
#' 
#' Transforms a vector of numbers on the real line to a vector of weights that
#' sum to one
#' @param p A vector of real valued numerics
#' @return A vector weights that sum to one, distributed according to p
#' @examples probs_link(c(-1, 0, 1))
#' @export
probs_link <- function(p){
  p %<>% plogis %>% c(1)
  t <- 1
  n1 <- length(p)
  output <- rep(NA, n1)

  for (i in 1:n1){
    output[i] <- t * p[i]
    t <- t * (1 - p[i])
  }
  return(output)
}

#' Mean link function
#' 
#' Transforms a vector of numbers on the real line to a vector of means that
#' have ascending values
#' @param m A vector of real valued mean parameters
#' @return A vector ascending mean values, starting at m[1]
#' @examples means_link(c(-1, 0, 1))
#' @export
means_link <- function(m){
  m %>% exp %>% cumsum %>% return
}

#' Normalisation function
#' 
#' Scales an input vector to have a sum of 1
#' @param v A vector of real valued parameters
#' @return A vector of scaled values that add up to 1
#' @export
sum_to_one <- function(v){
  return(v / sum(v))
}
