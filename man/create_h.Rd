% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_h.R
\name{create_h}
\alias{create_h}
\title{Create h matrix}
\usage{
create_h(U_and_W, data_obj, add_quad, nh, h_start)
}
\arguments{
\item{U_and_W}{Output from create_U_and_W}

\item{data_obj}{Output from create_embedded_data}

\item{add_quad}{Indicates whether to add a quadratic term to the linear ridge
regression}

\item{nh}{Number of hidden units}

\item{h_start}{Values to use for initializing h (vector of length nh)}
}
\description{
This function prepares the h matrices to use in the ridge regression.
}
