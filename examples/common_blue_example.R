# Author E. B. Dennis

# Example 1 - incorporating covariates
# This example demonstrates use of the rGAI package to allow for seasonal
# variation to vary with space by application to data for Common Blue
# Polyommatus icarus.

# Load packages
library(rGAI)
library(ggplot2)
library(ggridges)
library(rgdal)
library(scales)
library(msm)
library(reshape2)

set.seed(123)

# Read in example data for Common Blue
data("commonblue_data")
head(commonblue_data)

# Scale covariates 
commonblue_data$east <- scale(commonblue_data$EAST)
commonblue_data$north <- scale(commonblue_data$NORTH)
commonblue_data$north2 <- commonblue_data$north^2

# Define model fitting choices 
stopover_options <- list(B = 2, 
                         shared_sigma = FALSE,
                         w_formula = list(~north + north2), 
                         mu_formula = list(~north, ~north))
# Provide guesses for starting values on the parameter scale
stopover_starting_guesses <- list(mu = c(10,18), sigma = rep(1.5,2), w = c(.7,.3))
# Produce starting values on the link scale
transformed_starting_guesses <- transform_starting_values(starting_values = stopover_starting_guesses,
                                                          a_choice = "stopover",
                                                          dist_choice = "P",
                                                          options = stopover_options,
                                                          DF = commonblue_data)

# Fit the GAI with stopover model and Poisson distribution, with the
# mean emergence time, mu, regressed linearly on northing for each brood and
# the weighting parameter, w, a quadratic function of northing
stopover_GAI <- fit_GAI(start = transformed_starting_guesses,
                          DF = commonblue_data, 
                          a_choice = "stopover",
                          dist_choice = "P", 
                          options = stopover_options,
                          verbose = T, hessian = T)

# Model summary 
summary(stopover_GAI)

# Plot output using the in-built function
plot(stopover_GAI)

# Example of using transform_output to get estimates on the parameter scale for custom covariate values
transform_output(stopover_GAI, DF = commonblue_data[c(1,100), c("north", "north2")])

# Get SEs for sigma and phi on the parameter scale using the delta method 
# (but see later code for using the in-built bootstrap functions in the package)
stopover_cov <- diag(solve(stopover_GAI$hessian))
deltamethod(~exp(x1), stopover_GAI$par["sigma1"], stopover_cov["sigma1"])
deltamethod(~exp(x1), stopover_GAI$par["sigma2"], stopover_cov["sigma2"])
deltamethod(~exp(x1)/(1+exp(x1)), stopover_GAI$par["phi"], stopover_cov["phi"])

# Figure 1 - estimated flight curve for different covariate values
uk <- readOGR("UK_outline.shp")

# Define a sequence of northing values for plotting and scale them
north_plot_vals <- seq(50000, 950000, 100000)
north_plot_vals_s <- (north_plot_vals - attr(commonblue_data$north, "scaled:center"))/attr(commonblue_data$north, "scaled:scale")
# Get transformed estimates of the seasonal pattern for these northing values
gai_output <- transform_output(stopover_GAI, 
                               DF = data.frame(north = north_plot_vals_s,
                                               north2 = north_plot_vals_s^2),
                               provide_A = T)
# Produce data frames for plotting
df_for_plot <- data.frame(NORTH = north_plot_vals,
                     occasion = rep(1:26, each = length(north_plot_vals_s)),
                     a = c(gai_output$A))
w_for_plot <- data.frame(NORTH = north_plot_vals,
                         w = gai_output$params[, "w1"],3)

ggplot(df_for_plot)+
  # Plot UK outline
  geom_polygon(data = uk, 
               aes(rescale(long, c(1, 27)), lat/1000, group = group), fill = "darkgrey")+
  # add the seasonal curves 
  geom_density_ridges(aes(x = occasion, y = NORTH/1000, 
                          group = factor(NORTH/1000),
                          height = a, fill=NORTH/1000), 
                      stat = "identity", alpha=.75)+
  # add the estimates for w
  geom_label(aes(x = 25, y = NORTH/1000+25, 
                label = format(round(w, 3), nsmall = 3), 
                color = NORTH/1000), 
             size = 6,
             data = w_for_plot) +
  scale_fill_gradient(low = "red", high = "blue")+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("Occasion (week number)")+
  ylab("Northing (km)")+
  theme(legend.position = "bottom")+
  theme_light()+
  coord_fixed(ratio = (27-1)/(bbox(uk)[1,2] - bbox(uk)[1,1])*1000, 
              ylim = c(0,1250),
              expand = FALSE)+
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 18),
        legend.position = "none")

# Figure 2 - transformed parameter estimates varying with the covariate northing

# Bootstrap on the untransformed scale
# note that we have set R=100 for demonstrate but this may need to be higher in applications
untransformed_bootstrap <- bootstrap(stopover_GAI, R = 100, refit = FALSE,
                                     alpha = 0.05, transform = FALSE)
# Get transformed bootstrap estimates
param_ci <- transform_bootstrap_parameters(untransformed_bootstrap$par, stopover_GAI, 
                               data.frame(north = north_plot_vals_s,
                                          north2 = north_plot_vals_s^2))
# Create a data frame with transformed parameter estimates, for defined covariate values
param_cov <- melt(gai_output$params, value.name = "estimate")
# 95% CI for transformed parameter estimates, for defined covariate values
param_cov_ci <- merge(melt(param_ci$`2.5%`, value.name = "lower"),
                      melt(param_ci$`97.5%`, value.name = "upper"))
param_cov <- merge(param_cov, param_cov_ci)
# Add northing
param_cov <- merge(param_cov, data.frame(Var1 = 1:10, NORTH = north_plot_vals))

ggplot(param_cov[param_cov$Var2 %in% c("mu1","mu2","w1"),], aes(NORTH/1000, estimate))+
  facet_wrap(~Var2, scale = "free_y")+
  geom_line(aes(color = NORTH), size = 1.2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3)+
  geom_line(aes(y = upper, color = NORTH))+
  geom_line(aes(y = lower, color = NORTH))+
  scale_color_gradient(low = "red", high = "blue")+
  scale_fill_gradient(low = "red", high = "blue")+
  xlab("Northing (km)")+
  ylab("Parameter estimate")+
  theme_bw()+
  theme(text = element_text(size = 20),
                      legend.position = "none")