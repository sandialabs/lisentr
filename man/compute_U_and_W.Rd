% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_U_and_W.R
\name{compute_U_and_W}
\alias{compute_U_and_W}
\title{Create U and W matrices}
\usage{
compute_U_and_W(data_obj, nh, U_width, W_width, U_pi, W_pi, nu)
}
\arguments{
\item{data_obj}{Output from create_embedded_data}

\item{nh}{Number of hidden units}

\item{U_width}{Uniform distribution sampling range for U ("width" parameter)}

\item{W_width}{Uniform distribution sampling range for W ("width" parameter)}

\item{U_pi}{Probabilities of non-zeros in U (sparseness parameter)}

\item{W_pi}{Probabilities of non-zeros in W (sparseness parameter)}

\item{nu}{Spectral-radius parameter (scaling parameter for W-weight matrix)}
}
\description{
This function creates the U and W matrices for the model hidden state
}
