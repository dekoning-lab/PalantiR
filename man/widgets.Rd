\name{widgets}
\alias{PhyloPlot}
\alias{PhyloPlotOutput}
\alias{PhyloPlotRender}
\alias{AlignmentPlot}
\alias{AlignmentPlotOutput}
\alias{AlignmentPlotRender}
\alias{plot.Alignment}
\alias{plot.Simulation}
\alias{plot.SubstitutionModel}
\title{Interactive phylogeny and alignment widgets}
\arguments{
\item{phylogeny}{A \code{Phylogeny}.}
\item{substitutions}{Optionally a simulation's substitution table, drawn on the
branches.}
\item{intervals}{Optionally the interval table, to shade mode changes.}
\item{alignment}{An \code{Alignment}.}
\item{...}{Passed to the widget.}
}
\description{
htmlwidgets for viewing a phylogeny with its substitution history, and for
viewing an alignment. The \code{*Output} and \code{*Render} functions are the
Shiny bindings.
}
\details{
These are reached through \code{plot} on a \code{Phylogeny}, \code{Simulation} or
\code{Alignment} object.
}
