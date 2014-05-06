\name{org}
\alias{org}
\title{Convert object to ascii suitable for org-mode}
\usage{
  org(x, ..., ncol, include.rownames = TRUE,
    include.colnames = TRUE, header = TRUE,
    frame = "topbot", rownames = NULL, colnames = NULL,
    type = "org", tab = FALSE, margins = TRUE,
    print = TRUE, html, latex)
}
\arguments{
  \item{x}{R object}

  \item{...}{additional arguments to lower level functions}

  \item{ncol}{If \code{x} is a vector and \code{ncol} is
  given as argument, the resulting output will be a
  \code{matrix} with \code{ncol} columns}

  \item{include.rownames}{If \code{FALSE} row names are
  removed}

  \item{include.colnames}{If \code{FALSE} column names are
  removed}

  \item{header}{If TRUE the header is included}

  \item{frame}{Frame argument (see \code{ascii})}

  \item{rownames}{Optional vector of row names}

  \item{colnames}{Optional vector of column names}

  \item{type}{Type argument (see \code{ascii})}

  \item{tab}{Tabulate?}

  \item{margins}{Add margins to table?}

  \item{print}{print or return result}

  \item{html}{HTML prefix (added to ATTR_HTML)}

  \item{latex}{LaTeX prefix (added to ATTR_LaTeX)}
}
\description{
  Convert object to ascii suitable for org-mode
}
\author{
  Klaus K. Holst
}

