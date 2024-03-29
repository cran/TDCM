% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tdcm_emit}
\alias{tdcm_emit}
\alias{tdcm_warn}
\alias{tdcm_stop}
\title{Emit \code{tdcm}-related message}
\usage{
tdcm_emit(text, label = "INFO:", func = base::message, ...)

tdcm_warn(text, ...)

tdcm_stop(text, ...)
}
\arguments{
\item{text}{The message text to emit, an object that can be coerced to a
\code{"character"} via [\code{\link[base:character]{base::as.character()}}].}

\item{label}{The label that precedes the message. The default \code{label} is
\code{"INFO:"}.}

\item{func}{The function to use to emit the message. The default \code{func} is
[\code{\link[base:message]{base::message()}}].}
}
\description{
The \code{tdcm_emit} function is used internally by the \code{tdcm} package to ensure
that messages, warnings, and errors are emitted in a way that can be easily
suppressed by users of the package.
}
\section{Functions}{
\itemize{
\item \code{tdcm_warn()}: Emit \code{tdcm}-related warning message

\item \code{tdcm_stop()}: Emit \code{tdcm}-related stop message

}}
\keyword{internal}
