% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{bootstrap}
\alias{bootstrap}
\title{Bootstrapping for GAI models}
\usage{
bootstrap(
  GAI_fit,
  R = 100,
  refit = T,
  alpha = 0.05,
  parallel = T,
  cores = NULL,
  transform = T
)
}
\arguments{
\item{GAI_fit}{An object produced by using fit_GAI for model fitting}

\item{R}{The number of resamples to produce for the bootstrap}

\item{refit}{If TRUE, resamples the observations from the sites and occasions
uniformly, and refits the model at each bootsrap iteration. If FALSE, the
bootstrap simply resamples the fitted parameter values from their asymptotic
normal distribution (therefore this option requires a Hessian to have been
produced during the model fitting stage).}

\item{alpha}{1 - alpha gives the coverage the bootstrap confidence intervals
aim to produce.}

\item{parallel}{if TRUE, calculates the bootstraps in parallel, using the
maximum number of available cores.}

\item{cores}{If not NULL, this specifies the number of cores to use, if the
default of using all available cores is undesirable.}

\item{transform}{if TRUE, will return a bootstrap on the transformed
parameters, rather than on the link scale. It should be noted that when 
covariates are present in the model, transformed outputs are averaged across
all covariate values in the data, to avoid erroneous interpretation of 
covariate effects on the transformed scale.}
}
\value{
A named list with entries "EC", "N", "A" and "par" giving the
\code{c(alpha / 2, 1 - alpha / 2)} confidence interval for the expected
observed count at each site on each occasion, the estimated 
site super-population, the seasonal component, and the estimated
parameter values, respectively.
}
\description{
Produces either a non-parameteric bootstrap by refitting the GAI at each
iteration, or produces a parametric resampling of the MLEs at each iteration
using an estimate of their asymptotic normal distribution.
}
