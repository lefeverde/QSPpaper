% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degs.R
\name{make_cont_fit}
\alias{make_cont_fit}
\title{Wrapper to create a fit object (see \code{\link[limma]{eBayes}}) using the contrast method}
\usage{
make_cont_fit(v, group_cont, mod)
}
\arguments{
\item{v}{a voom object}

\item{group_cont}{vector of contrasts}

\item{mod}{model matrix}
}
\value{
fit object
}
\description{
Wrapper to create a fit object (see \code{\link[limma]{eBayes}}) using the contrast method
}
