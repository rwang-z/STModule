% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{convert_res_to_matrix}
\alias{convert_res_to_matrix}
\title{Convert results of GPU versions to matrices.}
\usage{
convert_res_to_matrix(res, type = "all")
}
\arguments{
\item{res}{the results of 'run_STModule' or 'run_spatial_map_estimation' running on a tissue section.}

\item{type}{type of estimation.
\itemize{
\item 'all': for resutls of 'run_STModule'.
\item 'estimate': for results of 'run_spatial_map_estimation'.
}}
}
\value{
processed results.
}
\description{
Convert results of GPU versions to matrices.
}
