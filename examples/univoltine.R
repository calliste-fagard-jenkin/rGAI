library(GAI)

# Load the data and look at the format:
data("example_data")
head(example_data)

# Simple options for a no-covariate univoltine model:
my_options <- list(B = 1)
options_for_splines <- list(df = 20, degree = 3)

# Mixture model:
my_mixture_GAI <- fit_GAI(start = c(2, 1, 1),
                          DF = example_data, a_choice = "mixture",
                          dist_choice = "ZIP", options = my_options,
                          verbose = T, hessian = T)
my_mixture_GAI$par

# Stopover model:
my_stopover_GAI <- fit_GAI(start = c(2, 1, 0),
                          DF = example_data, a_choice = "stopover",
                          dist_choice = "P", options = my_options,
                          verbose = T, hessian = T)
my_stopover_GAI$par