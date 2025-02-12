% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{cal_cov_inv_list_same}
\alias{cal_cov_inv_list_same}
\title{Calculate the inverse of covariance matrices for all tissue modules, same for all modules.}
\usage{
cal_cov_inv_list_same(r_const, num_comp, dist_mat, params)
}
\arguments{
\item{r_const}{constant r of the SE kernel for all tissue modules.}

\item{num_comp}{number of tissue modules.}

\item{dist_mat}{Euclidean distance matrix of the spots/cells.}

\item{params}{parameters of the estimation process.}
}
\value{
a list of inverse matrices of the covariance matrices of all tissue moduels.
}
\description{
Calculate the inverse of covariance matrices for all tissue modules, same for all modules.
}
