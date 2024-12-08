% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_mx.R
\name{glm_mx}
\alias{glm_mx}
\title{Maxent-like Generalized Linear Models (GLM)}
\usage{
glm_mx(formula, family = binomial(link = "cloglog"), data, weights = NULL, ...)
}
\arguments{
\item{formula}{A formula specifying the model to be fitted, in the format
used by \code{\link[stats]{glm}}.}

\item{family}{A description of the error distribution and link function to be
used in the model. Defaults to \code{binomial(link = "cloglog")},
which is commonly used for presence-background data.}

\item{data}{A \code{data.frame} containing the variables in the model. Must
include a column named \code{pr_bg} that indicates whether a record is
a presence (1) or background (0).}

\item{weights}{Optional. A numeric vector of weights for each observation. If
not provided, default weights of 1 for presences and 10000 for
background points will be used.}

\item{...}{Additional arguments to be passed to \code{\link[stats]{glm}}.}
}
\value{
A fitted \code{\link[stats]{glm}} object. The model object includes
the minimum and maximum values of the non-factor variables in the
dataset, stored as \code{model$varmin} and \code{model$varmax}.
}
\description{
This function fits a Maxent-like Generalized Linear Model (GLM) to binary
presence-background data. It allows for the specification of custom weights,
with a default weighting scheme where presences are weighted as 1 and
background points as 10000.
}
\details{
For more details, see Fithian and Hastie (2013) DOI: 10.1214/13-AOAS667.
}
