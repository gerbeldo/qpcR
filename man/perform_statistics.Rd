% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{perform_statistics}
\alias{perform_statistics}
\title{Perform Statistical Analysis}
\usage{
perform_statistics(wells, sample_type_levels = NULL, WT_string = "WT")
}
\arguments{
\item{wells}{A data frame containing wells information with `rel_exp_norm` computed.}

\item{sample_type_levels}{An optional character vector specifying the factor levels for `sample_type`. If `NULL`,
the wild-type level (matched by `WT_string`) is placed first.}

\item{WT_string}{A character string used to identify wild-type samples. Defaults to `"WT"`.}
}
\value{
A summary of the fitted linear mixed effects model.
}
\description{
Applies a linear mixed effects model (using lmer) to the normalized relative expression data,
analyzing the effect of sample type with biological replicates as random effects.
}
