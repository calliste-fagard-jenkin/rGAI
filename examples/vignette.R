## ----loading, results = 'hide', warning = F-----------------------------------------------------------------------------
#library(devtools, quietly = T)
#install_github("calliste-fagard-jenkin/GAI", quiet = F)
library(rGAI)


## ----data, fig.width = 7, fig.height = 5--------------------------------------------------------------------------------
# For the pipe operator:
library(magrittr)
data("example_data")
head(example_data)

# The GAI package comes with a handy function to convert data.frame data to 
# a matrix of counts, to facilitate plotting. This function also checks if 
# the survey has any missing occasions, so can be a useful way of sanitising
# a new dataset:
test_data <- extract_counts(example_data, returnDF = T)$matrix

# Plot the observed densities, averaged across all sites:
test_data %>% apply(2, mean, na.rm = T) %>% 
  plot(x = 1:ncol(test_data), type = 'l', col = rgb(1, 0, 0, 0.6), lty = 2, 
       xlab = 'Week', ylab = 'Observed count, averaged across all sites')


## ----fitting------------------------------------------------------------------------------------------------------------
# We can load the true parameter values of the simulated data set from the 
# GAI package, to more easily find starting parameter values for this example:
data("example_par")

# The 'options' argument of the fit_GAI function is how we specify the 
# specific model we'd like the package to fit. For mixture and stopover models
# we can specify if each brood should be modelled to have the same standard 
# deviation, and the number of broods our model should use:
my_options <- list(B = 3, shared_sigma = T)

# In the case of a spline-based model, the package needs to know the number of
# degrees of freedom and the degree of the splines we want to use:
options_for_splines <- list(df = 20, degree = 3)

# The fit_GAI function includes a 'verbose'  option, which prints the value of
# the negative log-likelihood at each iteration, and also the option to 
# estimate the hessian matrix of the parameters, to obtain their 
# variance covariance matrix:
my_mixture_GAI <- fit_GAI(start = example_par,
                          DF = example_data, a_choice = "mixture",
                          dist_choice = "ZIP", options = my_options,
                          verbose = T, hessian = T)

# Print the MLEs to show the model has fitted correctly:
my_mixture_GAI$par

# Not forgetting to add degrees of freedom and degree of splines information,
# and noting that specifying the number of broods will no longer have an 
# effect:
my_spline_GAI <- fit_GAI(start = rep(0, 20), DF = example_data,
                         a_choice = "splines", dist_choice = "P",
                         options = options_for_splines, verbose = T,
                         hessian = T)

# Print the MLEs to show the model has fitted correctly:
my_spline_GAI$par


## ----mistake------------------------------------------------------------------------------------------------------------
# We cut off some of our starting parameters on purpose, to cause the exception
# to be raised:
try({my_mixture_GAI <- fit_GAI(start = example_par[1:3], DF = example_data,
                               a_choice = "mixture", dist_choice = "ZIP",
                               options = my_options, verbose = T,
                               hessian = T)}, silent = T)


## ----covariates---------------------------------------------------------------------------------------------------------
# To specify a formula which will be identical for each brood, 
general_options <- list(B = 3, shared_sigma = T,
                        mu_formula = formula(~altitude))

# To be able to fit brood-specific covariate options for the standard deviations
# we require the standard deviation to be estimated individually for each brood,
# and so we set shared_sigma = F.
brood_specific_options <- list(B = 3, shared_sigma = F, 
                               sigma_formula =
                                 list(NULL,                     # for brood 1
                                      formula(~altitude),       # for brood 2
                                      formula(~I(altitude^2)))) # for brood 3

# The distributional parameter we had previously for the ZIP (Zero-Inflated 
# Poisson) distribution will disappear now that we're fitting a simple 
# Poisson. However, we are also adding a covariate formula which will require
# one additional parameter. We'll give this a starting value of 0, which 
# reflects an a priori belief that the covariate has no effect. To do this, 
# we simply replace the final value of example par with 0, since distributional
# parameters are placed at the end of the parameter vector, and covariate are 
# placed after these, when they are present:
general_fit_start <- example_par[1:(length(example_par) - 1)] %>% c(0)

# Having specified our covariate formulae, we can fit the model in exactly the
# same way as before:
general_fit <- fit_GAI(start = general_fit_start, DF = example_data,
                       a_choice = "mixture", dist_choice = "P",
                       options = general_options, hessian = T)

# Because these covariate data are dummies, we expect the fitted value to be
# very close to zero:
general_fit$par

# To make the starting values for the new model, we steal the known values from 
# the example_par vector (using the shared standard deviation across broods as
# the estimate for each brood individually). We then also add two 0s as our 
# starting values for the covariate parameters
brood_specific_start = c(example_par[c(1:3, rep(4, 3), 5:6)], rep(0, 2))
brood_specific_fit <- fit_GAI(start = brood_specific_start,
                              DF = example_data, a_choice = "mixture",
                              dist_choice = "P", hessian = T,
                              options = brood_specific_options)
# Checking the MLEs:
brood_specific_fit$par

# A quick example to show a stopover model with the same options:
general_fit_stopover <- fit_GAI(start = c(general_fit_start, 0),
                                DF = example_data, a_choice = "stopover",
                                dist_choice = "P", options = general_options)

general_fit_stopover$par


brood_specific_stopover <- fit_GAI(start = c(brood_specific_start, 0),
                                   DF = example_data, a_choice = "stopover",
                                   dist_choice = "P",
                                   options = brood_specific_options)

brood_specific_stopover$par

# And of course, we could also fit covariates to the simple univolitine case
# (although not for weights):
univ_options <- list(B = 1, sigma_formula = formula(~altitude))
univoltine_fit <- fit_GAI(start = c(2.2, 1, 0, 0), DF = example_data, 
                          a_choice =  "stopover", options = univ_options)

univoltine_fit$par


## ----covErrors----------------------------------------------------------------------------------------------------------
example_NA <- example_time_varying <- example_data

# Turn roughly 5% of our altitude data to NA values:
indices <- example_NA$altitude %>% length %>% rbinom(size = 1, prob = 0.05)
example_NA$altitude[indices %>% as.logical] <- NA

# Turn our altitude data into a time-varying covariate by making it a function 
# of the occasion:
example_time_varying$altitude <- example_data$altitude * example_data$occasion

# Trying to fit a covariate model with these new data will throw the relevant
# errors:
error_fit <- try(fit_GAI(start = general_fit_start, DF = example_NA,
                         a_choice = "mixture", dist_choice = "P",
                         options = general_options, hessian = T))

# Trying to fit with a time-varying formula: 
error_fit <- try(fit_GAI(start = general_fit_start, DF = example_time_varying,
                         a_choice = "mixture", dist_choice = "P",
                         options = general_options, hessian = T))



## ----transform1, fig.width = 7, fig.height = 5--------------------------------------------------------------------------
# Let's produce a slightly 'annotated' version of the plot of observed counts
# we had before. We'll add in some lines where we think the mean brood arrival
# times are, and add in horizontal lines to try to estimate thw standard 
# deviation of each brood's arrivals:
test_data %>% apply(2, mean, na.rm = T) %>% 
  plot(x = 1:ncol(test_data), type = 'l', col = rgb(1, 0, 0, 0.6), lty = 2, 
       xlab = 'Week', ylab = 'Observed count, averaged across all sites')

plot.col <- 'blue'
plot.base <- 1

# Mean brood arrival lines:
lines(rep(4.1, 2), c(plot.base, 4.6), col = plot.col)
lines(rep(13.0, 2), c(plot.base, 21.3), col = plot.col)
lines(rep(23.0, 2), c(plot.base, 8), col = plot.col) 

# Mean brood dispersion lines:
lines(c(1, 7.5), rep(plot.base, 2), col = plot.col)
lines(c(9, 18), rep(plot.base, 2), col = plot.col)
lines(c(19, 26), rep(plot.base, 2), col = plot.col)

# The length of the three horizontal bars that measure the width of our broods:
brood_widths <- c(6.5, 9, 7)

# Using a very rough rule of thumb that 95% of data in a normal distribution are
# found within 1.96 standard deviations from the means, we use the half width of
# our bars as an estimate of two sigma:
sigma_guesses <- brood_widths / 4

# Our guesses for the means can simply be read off of the lines of code that
# produced the plot:
mu_guesses <- c(4.1, 13, 23)

# To guess the weights, we make the crude assumption that the height of the peak
# of a brood relative to the others is roughly propoprtional to the percentage
# of the population that finds itself in this brood:
w_guesses <- c(4.6 - plot.base, 21.3 - plot.base, 8 - plot.base) %>%
  sum_to_one


## ----transform2---------------------------------------------------------------------------------------------------------
# We create a list of starting values for parameters with the same structure 
# as the options argument for fit_GAI:
my_starting_guesses <- list(mu = mu_guesses, sigma = sigma_guesses,
                            w = w_guesses)

# Note: When a ZIP or NB model is used, we can also include a guess for the
# distributional parameter by adding in a dist.par entry to this list, and 
# additionally a retention probability 'phi' can be added for stopover models 
new_brood_specific_start <-
  transform_starting_values(starting_values = my_starting_guesses,
                            a_choice = "mixture", dist_choice = "P",
                            options = brood_specific_options,
                            # When covariates are included in the model, we 
                            # need to include the data.frame of observations and
                            # covariate values, or the GAI package won't be able
                            # to figure out how many parameters the model 
                            # specification requires. This can be omitted for
                            # non-covariate-including models
                            DF = example_data)

# Let's print these new starting values, and refit the model:
print(new_brood_specific_start)
new_brood_specific_fit <- fit_GAI(start = new_brood_specific_start,
                                  DF = example_data, a_choice = "mixture",
                                  dist_choice = "P", hessian = T,
                                  options = brood_specific_options)
# Checking the MLEs:
new_brood_specific_fit$par


## ----bootstrap, warnings = F--------------------------------------------------------------------------------------------
# To produce a bootstrap, we must specify if we wish to refit the model
# (determining the type of bootstrap fitted), as well as the number of
# bootstrap resamples we would like to use (R):
general_fit_bootstrap <- bootstrap(general_fit, R = 500, refit = F,
                                   alpha = 0.01, transform = T)

# The more time-consuming bootstrap can be fitted in parallel, with a chosen
# number of cores. If this is not specified, the default number of cores 
# will be one less than the number available on the computer, so this must
# always be specified if running models on a large server:
refitting_bootstrap <- bootstrap(general_fit, R = 9, refit = T, parallel = T,
                                 cores = 3, alpha = 0.01, transform = T)
## ----intervals----------------------------------------------------------------------------------------------------------
# Taking a look at the results of the 'refit the model at each iteration'
# style of bootstrap:
refitting_bootstrap$par

## ----summary------------------------------------------------------------------------------------------------------------
# Get a basic summary of the model outputs:
summary(my_mixture_GAI)

# We can also obtain the AIC on its own using the standard AIC R function if 
# we aren't interested in producing the rest of the summary:
AIC(my_mixture_GAI)


## ----backtransform------------------------------------------------------------
# The easiest way to obtain transformed parameter values on the correct scale
# is to use the transform_output function. This function uses the skeleton and
# options present in the fitted model object to apply covariate formulas before
# applying link functions.

# We can create a data.frame with custom covariate values, or reuse values that
# we observed during the survey:
DF_to_transform <- data.frame(altitude = c(-100, 0, 100))
DF_to_transform <- general_fit$DF[1:3,]

# The transform_output function deals with all the covariate formulas and link
# functions by using the information contained in the fitted model object:
transform_output(general_fit, DF_to_transform)

# We can also use this function to get out the a_func matrix for a set of
# covariate values:
A <- transform_output(general_fit, DF_to_transform, provide_A = T)$A
matplot(t(A), type = 'l', col = c("blue", "darkblue", "darkgrey"),
        ylab = "Unscaled Flight Path Density", xlab = "Occasion")
    
# It's important to include all covariates exactly as they were named in the 
# call to fit_GAI, otherwise an error will be thrown, giving the name of the 
# missing covariate.
try(transform_output(brood_specific_fit, data.frame(altiitude = c(-10, 0, 10))))

# NA covariate values will throw an error, as always:
try(transform_output(brood_specific_fit, data.frame(altitude = c(NA, 0, 10))))


## ----plotting, fig.width = 7, fig.height = 5----------------------------------------------------------------------------
colours <- c("#33FFF9", "#33A8FF", "#4233FF")

# The default behaviour will use quantiles = c(0.05, 0.5, 0.95), and therefore 
# will not plot all sites individually.
plot(general_fit, scale_by_N = T, quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99))

# For the sake of demonstration, a plot with manual colours, with no scaling:
plot(brood_specific_fit, scale_by_N = F, quantiles = c(0.001, 0.5, 0.999),
     colours = colours)

