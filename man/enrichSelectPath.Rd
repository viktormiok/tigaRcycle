% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roxycode.R
\name{enrichSelectPath}
\alias{enrichSelectPath}
\title{enrichSelectPath}
\usage{
enrichSelectPath(genes, selectPathway, cut_off, name, orgdb, organism)
}
\arguments{
\item{genes}{genes for enrichment analysis.}

\item{selectPathway}{KEGG names, if NULL all pathway are included.}

\item{cut_off}{p-value cutoff on enrichment tests to report.}

\item{name}{name of the comparison.}

\item{orgdb}{organism data base.}

\item{organism}{supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'.}
}
\description{
This function performs KEGG enrichment analysis of selected pathways.
}
\author{
Viktorian Miok
}
