% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plsda.R
\name{gly_plsda}
\alias{gly_plsda}
\title{Partial Least Squares Discriminant Analysis (PLS-DA)}
\usage{
gly_plsda(
  exp,
  group_col = "group",
  ncomp = 2,
  scale = TRUE,
  add_info = TRUE,
  return_raw = FALSE,
  ...
)
}
\arguments{
\item{exp}{A \code{glyexp::experiment()} object containing expression matrix and sample information.}

\item{group_col}{A character string specifying the column name in sample information
that contains group labels. Default is "group".}

\item{ncomp}{An integer indicating the number of components to include. Default is 2.}

\item{scale}{A logical indicating whether to scale the data. Default is TRUE.}

\item{add_info}{A logical value. If TRUE (default), sample and variable information from the experiment
will be added to the result tibbles. If FALSE, only the PLS-DA results are returned.}

\item{return_raw}{A logical value. If FALSE (default), returns processed tibble results.
If TRUE, returns raw mixOmics plsda object.}

\item{...}{Additional arguments passed to \code{mixOmics::plsda()}.}
}
\value{
A list containing four tibbles (when return_raw = FALSE):
\itemize{
\item \code{samples}: PLS-DA scores for each sample with group information
\item \code{variables}: PLS-DA loadings for each variable
\item \code{variance}: PLS-DA explained variance information
\item \code{vip}: Variable Importance in Projection (VIP) scores for each variable
When return_raw = TRUE, returns the raw mixOmics plsda object.
}
}
\description{
Perform partial least squares discriminant analysis on the expression data.
The function uses \code{mixOmics::plsda()} to perform PLS-DA and returns tidy results.
}
\section{Required packages}{

This function requires the following packages to be installed:
\itemize{
\item \code{mixOmics} for PLS-DA analysis
}
}

\section{Sample size requirements}{

According to the Topliss ratio principle, the ratio of samples to variables (n/p)
should be at least 5 to avoid overfitting and ensure reliable results. This function
will throw an error if n/p < 5. For datasets with high dimensionality relative to
sample size, consider:
\itemize{
\item Feature selection before analysis
\item Collecting more samples
}
}

\seealso{
\code{\link[mixOmics:plsda]{mixOmics::plsda()}}
}
