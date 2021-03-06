% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/impute_expr.R
\name{tilling_lasso_preds_expr}
\alias{tilling_lasso_preds_expr}
\title{Tilling LASSO linear modal to impute expression of landmark genes}
\usage{
tilling_lasso_preds_expr(x, y, s.use = 20, y.name = NULL,
  do.print = FALSE, use.gram = FALSE, w = 0.8, loops = 10)
}
\arguments{
\item{x}{Matrix of predictors(columns), observed genes expression.}

\item{y}{Vector of response, observed expression of given landmarked gene.}

\item{s.use}{Numeric. See \link{predict.lars}'s parameter \code{s}.}

\item{y.name}{Chractor. Name of the given landmarked gene.}

\item{do.print}{Logic. Whether print \code{y.name} on screen.}

\item{use.gram}{Logic. See \link{lars}'s parameter \code{use.Gram}}

\item{w}{Window size to set the training data size relative to input
\code{x}. Default: 0.8.}

\item{loops}{Times to run tilling LASSO.}
}
\value{
Returns the imputed landmarked gene expession via LASSO modal.
}
\description{
Tilling LASSO linear modal to impute expression of landmark genes
}

