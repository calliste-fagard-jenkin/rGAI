# Author E. B. Dennis

# Example 2 - modelling multiple broods
# This example demonstrates use of the rGAI package to a species with
# three broods per year by applications to Small Copper Lycaena phlaeas

# Load packages
library(rGAI)
library(ggplot2)
library(dplyr)

# Read in example data for Small Copper
data("smallcopper_data")
head(smallcopper_data)

# Summarise and plot the observed count data
count_summary <- smallcopper_data %>% group_by(occasion) %>%
                      summarise(mean = mean(count, na.rm=TRUE),
                                lower = quantile(count, 0.05, na.rm = TRUE),
                                upper = quantile(count, 0.95, na.rm = TRUE))
   
# Plot average counts across sites per week
ggplot(count_summary, aes(occasion, mean))+
                geom_point(size = 2)+
                geom_errorbar(aes(ymin = lower, ymax = upper))+
                geom_line()+
                ylab("Count")+
                xlab("Occasion (week)")+
                theme_classic()+
                theme(text = element_text(size = 18))+
                scale_x_continuous(breaks = seq(0,25,5)) 

# Fit the GAI model with several different options as shown in the paper

# 1. Fit for B = 1 brood, Poisson distribution
starting_guesses1 <- list(mu = 13, sigma = 2)
transformed_starting_guesses1 <- transform_starting_values(starting_values = starting_guesses1,
                                                          a_choice = "mixture", 
                                                          dist_choice = "P",
                                                          options = list(B = 1),
                                                          DF = smallcopper_data)
GAI_B1_P <- fit_GAI(start = transformed_starting_guesses1,
                          DF = smallcopper_data, 
                          a_choice = "mixture",
                          dist_choice = "P", 
                          options = list(B = 1),
                          verbose = T, hessian = T)

# 2. Fit for B = 2 broods with shared sigma, Poisson distribution
starting_guesses2 <- list(mu = c(7, 20), sigma = 2, w = c(.2, .8), shared_sigma = TRUE)
transformed_starting_guesses2 <- transform_starting_values(starting_values = starting_guesses2,
                                                          a_choice = "mixture", 
                                                          dist_choice = "P",
                                                          options = list(B = 2, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B2_P <- fit_GAI(start = transformed_starting_guesses2,
                          DF = smallcopper_data, 
                          a_choice = "mixture",
                          dist_choice = "P", 
                          options = list(B = 2, shared_sigma = TRUE),
                          verbose = T, hessian = T)


# 3. Fit for B = 3 broods with shared sigma, Poisson distribution
starting_guesses3 <- list(mu = c(8, 17, 25), sigma = 2, w = c(.3, .3, .4), shared_sigma = TRUE)
transformed_starting_guesses3 <- transform_starting_values(starting_values = starting_guesses3,
                                                          a_choice = "mixture", 
                                                          dist_choice = "P",
                                                          options = list(B = 3, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B3_P <- fit_GAI(start = transformed_starting_guesses3,
                           DF = smallcopper_data, 
                           a_choice = "mixture",
                           dist_choice = "P", options = list(B = 3, shared_sigma = TRUE),
                           verbose = T, hessian = T)



# 4. Fit for B = 3 broods with shared sigma, Zero-inflated Poisson distribution
starting_guesses4 <- list(mu = c(8, 17, 25), sigma = 2, w = c(.3, .3, .4), shared_sigma = TRUE)
transformed_starting_guesses4 <- transform_starting_values(starting_values = starting_guesses4,
                                                          a_choice = "mixture", 
                                                          dist_choice = "ZIP",
                                                          options = list(B = 3, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B3_ZIP <- fit_GAI(start = transformed_starting_guesses4,
                           DF = smallcopper_data, 
                           a_choice = "mixture",
                           dist_choice = "ZIP", options = list(B = 3, shared_sigma = TRUE),
                           verbose = T, hessian = T)


# 5. Fit for B = 3 broods with shared sigma, negative-binomial distribution
starting_guesses5 <- list(mu = c(8, 17, 25), sigma = 2, w = c(.3, .3, .4), shared_sigma = TRUE)
transformed_starting_guesses5 <- transform_starting_values(starting_values = starting_guesses5,
                                                          a_choice = "mixture", 
                                                          dist_choice = "NB",
                                                          options = list(B = 3, shared_sigma = TRUE),
                                                          DF = smallcopper_data)
GAI_B3_NB <- fit_GAI(start = transformed_starting_guesses5,
                               DF = smallcopper_data, 
                               a_choice = "mixture",
                               dist_choice = "NB", 
                               options = list(B = 3, shared_sigma = TRUE),
                               verbose = T, hessian = T)


# 6. Fit for B = 3 broods with separate sigma, negative-binomial distribution
starting_guesses6 <- list(mu = c(8, 17, 25), sigma = c(2, 2, 2), w= c(.3, .3, .4), shared_sigma = FALSE)
transformed_starting_guesses6 <- transform_starting_values(starting_values = starting_guesses6,
                                                          a_choice = "mixture", dist_choice = "NB",
                                                          options = list(B = 3, shared_sigma = FALSE),
                                                          DF = smallcopper_data)
GAI_B3_NB_sig <- fit_GAI(start = transformed_starting_guesses6,
                              DF = smallcopper_data, a_choice = "mixture",
                              dist_choice = "NB", options = list(B = 3, shared_sigma = FALSE),
                              verbose = T, hessian = T)

# Collate models for comparison by AIC (Table 2 in the paper)
modelcomparison <- data.frame(mod = 1:6,
                   npar = c(length(GAI_B1_P$par),
                            length(GAI_B2_P$par),
                            length(GAI_B3_P$par),
                            length(GAI_B3_ZIP$par),
                            length(GAI_B3_NB$par),
                            length(GAI_B3_NB_sig$par)),
                   AIC = c(AIC(GAI_B1_P),
                           AIC(GAI_B2_P),
                           AIC(GAI_B3_P),
                           AIC(GAI_B3_ZIP),
                           AIC(GAI_B3_NB),
                           AIC(GAI_B3_NB_sig)))
                   
modelcomparison$deltaAIC <- modelcomparison$AIC - min(modelcomparison$AIC)

# Now we focus on the "best" model
# Model summary 
summary(GAI_B3_NB_sig)

# Plot output using the in-built function
plot(GAI_B3_NB_sig)

# Get estimates on the parameter scale 
transform_output(GAI_B3_NB_sig)

# Parametric bootstrap to get confidence intervals for all parameters, including seasonal pattern
best_GAI_boot <- bootstrap(GAI_B3_NB_sig, R = 500, refit = FALSE,
                                   alpha = 0.05)
best_GAI_boot$par

# Produce Figure 3
# Get the predicted average count per occasion (week)
best_GAI_wcount <- data.frame(Week <- 1:26,
                     avcount = transform_output(GAI_B3_NB_sig, DF = smallcopper_data[1,], provide_A=TRUE)$A[1,]*mean(GAI_B3_NB_sig$N))

# Add 5% and 95% quantiles for the predicted counts per occasion (week)
best_GAI_wcount$predcount_upper <- apply(GAI_B3_NB_sig$N*GAI_B3_NB_sig$A, 2, quantile, 0.95)
best_GAI_wcount$predcount_lower <- apply(GAI_B3_NB_sig$N*GAI_B3_NB_sig$A, 2, quantile, 0.05)

ggplot(count_summary, aes(occasion, mean))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper))+
  geom_line(aes(Week, avcount), data = best_GAI_wcount, col = "blue", size = 1)+
  geom_ribbon(aes(Week, avcount, ymin = predcount_lower, ymax = predcount_upper), 
              alpha = .3, fill = "blue", data = best_GAI_wcount)+
  ylab("Count")+
  xlab("Occasion (week)")+
  theme_classic()+
  theme(text = element_text(size = 18))+
  scale_x_continuous(breaks = seq(0,25,5)) 


