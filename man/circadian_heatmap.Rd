% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roxycode.R
\name{circadian_heatmap}
\alias{circadian_heatmap}
\title{circadian_heatmap}
\usage{
circadian_heatmap(
  data,
  colinf,
  sumdat,
  title,
  filename,
  colors,
  distfun,
  hclustfun,
  legs,
  legc
)
}
\arguments{
\item{data}{numeric matrix of the values to be plotted.}

\item{colinf}{column information.}

\item{sumdat}{sum data information.}

\item{title}{main title of the heatmap.}

\item{filename}{name of the file.}

\item{colors}{colors used for the image. Defaults to heat colors (heat.colors).}

\item{distfun}{function used to compute the distance (dissimilarity) between both rows and columns.}

\item{hclustfun}{function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms.}

\item{legs}{information about the legends.}

\item{legc}{colors for the legends.}
}
\description{
This function performs circadian heatmap plot.
}
\author{
Viktorian Miok
}
