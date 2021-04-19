library(GAI)
library(magrittr)

data("example_data")
data("example_par")

test_data <- extract_counts(example_data, returnDF = T)$matrix

my_options <- list(B = 3, shared_sigma = T)
options_for_splines <- list(df = 20, degree = 3)

my_mixture_GAI <- fit_GAI(start = example_par,
                          DF = example_data, a_choice = "mixture",
                          dist_choice = "ZIP", options = my_options,
                          verbose = T, hessian = T)
my_mixture_GAI$par

my_spline_GAI <- fit_GAI(start = rep(0, 20), DF = example_data,
                         a_choice = "splines", dist_choice = "P",
                         options = options_for_splines, verbose = T,
                         hessian = T)

my_spline_GAI$par

general_options <- list(B = 3, shared_sigma = T,
                        mu_formula = formula(~altitude))

brood_specific_options <- list(B = 3, shared_sigma = F, 
                               sigma_formula =
                                 list(NULL,                     # for brood 1
                                      formula(~altitude),       # for brood 2
                                      formula(~I(altitude^2)))) # for brood 3

brood_specific_opt_mu <- list(B = 3, shared_sigma = F, 
                               sigma_formula =
                                 list(NULL,                     # for brood 1
                                      formula(~altitude),       # for brood 2
                                      formula(~I(altitude^2)))) # for brood 3

general_fit_start <- example_par[1:(length(example_par) - 1)] %>% c(0)

general_fit <- fit_GAI(start = general_fit_start, DF = example_data,
                       a_choice = "mixture", dist_choice = "P",
                       options = general_options, hessian = T)
general_fit$par

brood_specific_start = c(example_par[c(1:3, rep(4, 3), 5:6)], rep(0, 2))
brood_specific_fit <- fit_GAI(start = brood_specific_start,
                              DF = example_data, a_choice = "mixture",
                              dist_choice = "P", hessian = T,
                              options = brood_specific_options)

mu_specific_fit <- fit_GAI(start = brood_specific_start,
                              DF = example_data, a_choice = "mixture",
                              dist_choice = "P", hessian = T,
                              options = brood_specific_opt_mu)
brood_specific_fit$par

general_fit_stopover <- fit_GAI(start = c(general_fit_start, 0),
                                DF = example_data, a_choice = "stopover",
                                dist_choice = "P", options = general_options)

general_fit_stopover$par


brood_specific_stopover <- fit_GAI(start = c(brood_specific_start, 0),
                                   DF = example_data, a_choice = "stopover",
                                   dist_choice = "P",
                                   options = brood_specific_options)

brood_specific_stopover$par

univ_options <- list(B = 1, sigma_formula = formula(~altitude))
univoltine_fit <- fit_GAI(start = c(2.2, 1, 0, 0), DF = example_data, 
                          a_choice =  "stopover", options = univ_options)

univoltine_fit$par

plot.base <- 1
brood_widths <- c(6.5, 9, 7)
sigma_guesses <- brood_widths / 4
mu_guesses <- c(4.1, 13, 23)
w_guesses <- c(4.6 - plot.base, 21.3 - plot.base, 8 - plot.base) %>%
  sum_to_one

my_starting_guesses <- list(mu = mu_guesses, sigma = sigma_guesses,
                            w = w_guesses)
new_brood_specific_start <-
  transform_starting_values(starting_values = my_starting_guesses,
                            a_choice = "mixture", dist_choice = "P",
                            options = brood_specific_options,
                            DF = example_data)
print(new_brood_specific_start)
new_brood_specific_fit <- fit_GAI(start = new_brood_specific_start,
                                  DF = example_data, a_choice = "mixture",
                                  dist_choice = "P", hessian = T,
                                  options = brood_specific_options)
new_brood_specific_fit$par

# general_fit_bootstrap <- bootstrap(general_fit, R = 1000, refit = F,
#                                    alpha = 0.01)
# refitting_bootstrap <- bootstrap(general_fit, R = 10, refit = T, parallel = T,
#                                  cores = 3, alpha = 0.01)
# refitting_bootstrap$par
# 
# summary(my_mixture_GAI)
# AIC(my_mixture_GAI)


DF_to_transform <- data.frame(altitude = c(-10, 0, 10))
transform_output(brood_specific_fit, DF_to_transform)

A <- transform_output(general_fit, DF_to_transform, provide_A = T)$A
A <- transform_output(brood_specific_fit, DF_to_transform, provide_A = T)$A
matplot(t(A), type = 'l', col = c("blue", "darkblue", "darkgrey"),
        ylab = "Unscaled Flight Path", xlab = "Occasion")

# debug_options <- list(B = 3, mu_formula =
#                         list(~altitude, ~I(altitude^2) + altitude, ~altitude))
# 
# debug_fit <- fit_GAI(rep(0, 10), example_data, "mixture", "P", debug_options)
# transform_output(debug_fit, data.frame(altitude = c(-10, 0, 10)), provide_A = T)
# 
# colours <- c("#33FFF9", "#33A8FF", "#4233FF")
# 
# # The default behaviour will use quantiles = c(0.05, 0.5, 0.95), and therefore 
# # will not plot all sites individually.
# plot(general_fit, scale_by_N = T, quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99))
# plot(brood_specific_fit, scale_by_N = F, quantiles = c(0.001, 0.5, 0.999),
#      colours = colours)