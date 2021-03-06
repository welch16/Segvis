% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-SegvisData.R
\docType{methods}
\name{plot_profile}
\alias{plot_profile}
\alias{plot_profile,SegvisData-method}
\alias{plot_profile}
\title{plot_profile method}
\usage{
plot_profile(object, ...)

\S4method{plot_profile}{SegvisData}(object, FUN = mean,
  nameFiles = basename(files(object)), type = "aggr", base = 1e+06,
  mc.cores = getOption("mc.cores", 2L), ...)
}
\arguments{
\item{object}{a \code{SegvisData} object.}

\item{...}{Any other additional parameters that \code{plot_profile} may need.}

\item{FUN}{a function to summarize the profiles by coordinate. The default
value is the \code{mean}.}

\item{nameFiles}{a character vector with the shortened names that are going
to be used in the plot.}

\item{type}{a character value indicating which strand to use. By default
uses the aggregate coverage between both strands.}

\item{base}{a numeric value indicating the number of aligned reads to which
the signal is going to be normalized. The default value is \code{1e6}.
\code{DT_profile} always normalizes the signal.}

\item{mc.cores}{a numeric value indicating the number of multi-cores to be
used.}
}
\value{
The \code{plot_profile} method returns a \code{ggplot2} object with
the summarized functional profiles.
}
\description{
plot_profile return a \code{ggplot2} object with the plot of the summarized
functional profiles (summarized by using the argument \code{FUN}). For each
coverage, file and distance from the anchor coordinate it applies FUN to
the tags vector for all regions.
}
\seealso{
\code{\link{SegvisData-class}}
}
