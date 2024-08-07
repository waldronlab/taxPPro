% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{getNsti}
\alias{getNsti}
\title{Get NSTI}
\usage{
getNsti(tree, annotated_tip_labels, nodes)
}
\arguments{
\item{tree}{A phylo object}

\item{annotated_tip_labels}{A character vector with the names of the tips}

\item{nodes}{A character vector with the names of nodes in the tree
that are annotated.}
}
\value{
A data.frame with two columns: tip_label and nsti.
}
\description{
\code{getNsti} gets the nsti described in picrust2. Code is based
on the source code of picrust2 using the castor package. NSTI
values are described in picrust and picrust2.
}
