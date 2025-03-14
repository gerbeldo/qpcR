% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{normalize_to_mean_WT}
\alias{normalize_to_mean_WT}
\title{Normalize Relative Expression to Mean Wild-Type}
\usage{
normalize_to_mean_WT(wells, target, WT_string = "WT")
}
\arguments{
\item{wells}{A data frame containing wells information with a `rel_exp` column.}

\item{target}{A character string specifying the target gene to use for normalization.}

\item{WT_string}{A character string used to identify wild-type samples. Defaults to `"WT"`.}
}
\value{
A data frame with an additional column `rel_exp_norm` representing the normalized relative expression.
}
\description{
Normalizes the relative expression values by dividing by the mean relative expression of wild-type samples.
}
