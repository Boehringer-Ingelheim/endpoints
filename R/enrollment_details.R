#' Specification format for \code{enrollment_details} used by \code{makeData()}
#' 
#' \code{enrollment_details} is a named list controlling administrative
#' censoring and stochastic enrollment in \code{\link{makeData}}.
#' 
#' If omitted, \code{makeData()} fills in defaults internally: \preformatted{
#' list( administrative_censoring = NULL, enrollment_distribution = "none",
#' enrollment_exponential_rate = NULL, piecewise_enrollment_cutpoints = NULL,
#' piecewise_enrollment_rates = NULL ) }
#' 
#' 
#' @section Purpose:
#' 
#' Administrative censoring and stochastic enrollment affect \emph{observed
#' follow-up}. This primarily affects time-to-event endpoints, but may also be
#' of use for possible interim analyses.
#' 
#' If no enrollment model is supplied (\code{enrollment_distribution =
#' "none"}), all subjects are treated as enrolling at time 0.
#' 
#' If enrollment is stochastic and administrative censoring is enabled, each
#' subject receives an \code{enrollTime}, and maximum available follow-up is
#' reduced from \eqn{\mathcal{A}} to \eqn{\mathcal{A} - T_E}, where
#' \eqn{\mathcal{A}} is the administrative censoring time and \eqn{T_E} is the
#' subject's enrollment time.
#' @section Arguments:
#' 
#' \describe{ \item{list("administrative_censoring")}{ Numeric scalar or
#' \code{NULL}. If non-\code{NULL} and positive, this defines the
#' administrative end of follow-up.  Any TTE outcome occurring after this limit
#' is administratively censored.  }
#' 
#' \item{list("enrollment_distribution")}{ Character string giving the
#' enrollment-time distribution. Must be one of: \describe{
#' \item{list("\"none\"")}{All subjects enroll at time 0.}
#' \item{list("\"uniform\"")}{Enrollment times are drawn from
#' \eqn{U(0,\mathcal{A})}.} \item{list("\"exponential\"")}{Enrollment times are
#' drawn from an exponential distribution, truncated at \eqn{\mathcal{A}}.
#' Requires \code{enrollment_exponential_rate}.}
#' \item{list("\"piecewise\"")}{Enrollment times are generated from a piecewise
#' exponential model over user-specified intervals. Requires
#' \code{piecewise_enrollment_cutpoints} and
#' \code{piecewise_enrollment_rates}.} } }
#' 
#' \item{list("enrollment_exponential_rate")}{ Numeric scalar (> 0). Used only
#' when \code{enrollment_distribution = "exponential"}.  Gives the exponential
#' rate governing enrollment times.  }
#' 
#' \item{list("piecewise_enrollment_cutpoints")}{ Numeric vector of strictly
#' increasing cutpoints defining the intervals for piecewise enrollment.  For
#' example, \code{c(0, 8, 16, 24)} defines three intervals: \code{[0,8)}, \code{[8,16)}, and
#' \code{[16,24]}.  }
#' 
#' \item{list("piecewise_enrollment_rates")}{ Numeric vector of positive
#' exponential rates, one for each interval defined by
#' \code{piecewise_enrollment_cutpoints}. Its length must be
#' \code{length(piecewise_enrollment_cutpoints) - 1}.  } }
#' @section Administrative censoring only:
#' 
#' If \code{administrative_censoring} is supplied and
#' \code{enrollment_distribution = "none"}, all subjects are assumed to enter
#' at time 0 and remain under observation until the administrative end of
#' follow-up (unless they experience an event or are independently censored
#' earlier).
#' @section Uniform enrollment:
#' 
#' Under uniform enrollment, subjects enter uniformly over the interval
#' \eqn{[0,\mathcal{A}]}. This induces heterogeneous follow-up, with later
#' enrollees having less potential follow-up time before administrative
#' censoring.
#' @section Exponential enrollment:
#' 
#' Under exponential enrollment, subject entry times are sampled from an
#' exponential distribution and truncated at the administrative censoring time.
#' This can be used to represent rapid early enrollment followed by a tapering
#' accrual pattern.
#' @section Piecewise exponential enrollment:
#' 
#' Under piecewise enrollment, the waiting time to enrollment is simulated
#' interval-by-interval. Within each interval, an exponential waiting time is
#' drawn using that interval's rate. If the draw exceeds the remaining interval
#' width, the process moves to the next interval.
#' 
#' This allows users to encode enrollment patterns such as slow-then-fast
#' accrual, fast-then-slow accrual, or more complex piecewise shapes. See the
#' \code{vignette("user_guide", package = "endpoints")} vignette for an example
#' of how to use this option to achieve certain enrollment proportions in
#' certain enrollment periods.
#' @name endpoint_details
#' @keywords documentation
#' @examples
#' 
#' ## Administrative censoring only
#' enrollment_details <- list(
#'   administrative_censoring = 24
#' )
#' 
#' ## Exponential enrollment
#' enrollment_details <- list(
#'   administrative_censoring    = 24,
#'   enrollment_distribution     = "exponential",
#'   enrollment_exponential_rate = 1 / 4
#' )
#' head(enrollment_details$data)
#' 
#' ## Piecewise exponential enrollment
#' enrollment_details <- list(
#'   administrative_censoring       = 24,
#'   enrollment_distribution        = "piecewise",
#'   piecewise_enrollment_cutpoints = c(0, 8, 16, 24),
#'   piecewise_enrollment_rates     = c(0.02, 0.06, 0.50)
#' )
#' 
NULL
