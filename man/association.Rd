% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roxycode.R
\name{association}
\alias{association}
\title{association}
\usage{
association(
  object,
  what = c("p-value", "statistic", "z-score", "weighted"),
  cluster = "average",
  alpha = 0.05,
  sort = TRUE,
  zoom = FALSE,
  legend = TRUE,
  plot = TRUE,
  plot_name = "global_test",
  colors,
  alias,
  help.lines = FALSE,
  cex.labels = 0.6,
  ylim,
  pdf,
  trace,
  mirror
)
}
\arguments{
\item{object}{A gt.object, usually created by a call to gt function.}

\item{what}{Gives a choice between various presentations of the same plot.}

\item{cluster}{Type of the hierarchical clustering performed for the dendrogram.}

\item{alpha}{Parameter between 0 and 1. Set the level of the family-wise error control in the multiple testing procedure performed on the dendrogram.}

\item{sort}{If TRUE, the plot sorts the bars with the most significant covariates and subjects to the left, as far as is possible within the constraints of the dendrogram (if present).}

\item{zoom}{If TRUE, discards non-significant branches from the dendrogram with the corresponding covariates.}

\item{legend}{If TRUE, draws a legend in the plot.}

\item{plot}{If FALSE, suppress all plotting.}

\item{plot_name}{Title of the plot.}

\item{colors}{The colors to be used for the bars.}

\item{alias}{Optional alternative labels for the bars in the plots.}

\item{help.lines}{If TRUE, prints grey dotted lines that help connect the dendrogram to the bars.}

\item{cex.labels}{Magnification factor for the x-axis labels.}

\item{ylim}{Override for the y axis limits of the barplot.}

\item{pdf}{Optional filename (character) of the pdf file to which the plots are to be written.}

\item{trace}{If TRUE, prints progress information.}

\item{mirror}{If TRUE, plots the reverse of the scores for the subjects with negative residual response, so that "good" scores are positive for all subjects.}
}
\description{
This function creates visualization of the a global test results in terms of the contributions of the covariates.
}
\author{
Viktorian Miok
}
