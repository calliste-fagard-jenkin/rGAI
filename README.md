# rGAI
Generalised Abundance Index for seasonal invertebrates

# Description
The rGAI package is an extension of functionality provided by Dennis et al (2016) https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12506. The package provides a user interface for calculating a Generalised Abundance Index (GAI) for seasonal invertebrates, with multiple options for covariate inclusion, seasonal flight patterns and observed count distributions, along with the ability to conduct bootstraps to obtain measures of uncertainty for estimated parameter and GAI values.

# Installation Instructions
rGAI can be installed with the `devtools` R library as follows:

```
library(devtools)
install_github("calliste-fagard-jenkin/rGAI", quiet = F)
```
# Bug Reports
Please report bugs via the GitHub interface, or by contacting the package's maintainer directly via email at cfj2@st-andrews.ac.uk. Please also contact the maintainer for requests to contribute. For questions related to the Generalised Abundance Index itself, please contact the lead author of the 2016 paper: E.B.Dennis@kent.ac.uk 

# Package Vignette
A full tutorial covering the package's range of features can be found in the /vignettes and /examples folders. For further descriptions and examples, please see the pre-print of the package's sister paper: https://www.authorea.com/doi/full/10.22541/au.163601791.14513063/v1

# Basic Usage

The GAI model has two major components. The first is a seasonal flight path model, which describes the seasonal variation in observed counts of individuals throughout a survey season. This is referred to as `a_choice` by the rGAI package, in reference to the use of the letter 'a' to denote this aspect of the model in its 2016 introductory paper (linked above). The options are the use of normal mixtures, a stopover model, or the use of splines. The second major component is the distribution of counts themselves. Either poisson, negative binomial, or zero-inflated poisson distributions can be selected. This aspect of the model is referred to as `dist_choice` by the rGAI package. The below set of simple examples are designed to illustrate using the `fit_GAI` function, to fit GAI models to your own data, as well as how to plot these outputs. Far more detail is available in the package vignette, including the use of covariates, and bootstraps, to produce confidence intervals on quantities of interest.

```
# The rGAI package comes with an example data set, which we will use to illustrate basic functionality:
data(example_data)

# Defining a_choice settings for mixture and stopover models. B refers to the number of broods
# in the model (i.e. the number of components in the mixture). shared_sigma = TRUE assumed 
# all mixture components have the same variance:
my_options <- list(B = 3, shared_sigma = T))

# Defining a_choice settings for spline based models:
options_for_splines <- list(df = 20, degree = 3)

# Now, fitting the model:
my_mixture_GAI <- fit_GAI(start = example_par,
                          DF = example_data, a_choice = "mixture",
                          dist_choice = "ZIP", options = my_options,
                          verbose = T, hessian = T)

# Print the MLEs to show the model has fitted correctly:
my_mixture_GAI$par

# We can plot the fitted model using R's inbuilt plot function. The scale_by_N option allows
# plotted points at each point in time to be scaled by estimated site population totals. 
# The quantiles option allows us to select which quantiles of the variation across sites
# we which to plot at each point in time:
plot(my_mixture_GAI, scale_by_N = T, quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99))

```
