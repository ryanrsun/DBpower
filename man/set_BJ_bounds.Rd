% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_BJ_bounds.R
\name{set_BJ_bounds}
\alias{set_BJ_bounds}
\title{set_BJ_bounds.R}
\usage{
set_BJ_bounds(alpha, d)
}
\arguments{
\item{alpha}{Type I error of test.}

\item{d}{Number of elements in set.}
}
\value{
A d*1 vector of bounds on the magnitudes of the test statistics, where
the first element is the bound for |Z|_(1) and the last element is the bound for |Z|_(d).
}
\description{
Finds the boundary points of the rejection region for the BJ statistic when
all elements in a set are independent.
}
\examples{
set_BJ_bounds(alpha = 0.01, d=5)


}
