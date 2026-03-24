#' Specification format for \code{calibration_control} used by
#' \code{makeData()}
#' 
#' \code{calibration_control} is a named list of tuning parameters used by
#' \code{\link{makeData}} when \code{target_correlation = TRUE}.
#' 
#' These settings control the numerical search used to calibrate a latent
#' Gaussian correlation matrix so that, after transformation through the
#' endpoint-specific marginals, the observed Pearson correlation matrix is
#' approximately equal to the user-supplied \code{correlation_matrix}.
#' 
#' 
#' @section Default values: The default list is: \preformatted{ list( n_mc =
#' 10000, tol = 0.001, maxit = 100, rho_cap = 0.999, ensure_pd = TRUE,
#' conv_norm_type = "F" ) }
#' @section Arguments: \describe{ \item{list("n_mc")}{ Monte Carlo sample size
#' used internally when approximating the transformed correlation induced by a
#' candidate latent Gaussian correlation.  Larger values generally improve
#' numerical stability but increase computation time.  }
#' 
#' \item{list("tol")}{ Positive numeric tolerance passed to the one-dimensional
#' root-finding routine used in pairwise calibration.  Smaller values can
#' improve accuracy but may require more iterations and computation.  }
#' 
#' \item{list("maxit")}{ Positive integer giving the maximum number of
#' root-finding iterations.  }
#' 
#' \item{list("rho_cap")}{ Numeric scalar in (0,1) giving the maximum absolute
#' latent Gaussian correlation considered during calibration.  This avoids
#' numerical instability at the exact boundaries \eqn{(\pm 1}).  }
#' 
#' \item{list("ensure_pd")}{ Logical scalar. If \code{TRUE}, the calibrated
#' latent correlation matrix is repaired (if needed) to ensure it is a valid
#' positive-definite correlation matrix before simulation proceeds.  }
#' 
#' \item{list("conv_norm_type")}{ Character scalar passed through to
#' positive-definiteness repair routines (for example,
#' \code{Matrix::nearPD(..., conv.norm.type = ...)} when available).  The
#' default \code{"F"} (Frobenius)).  } }
#' @section When calibration is used:
#' 
#' These settings are used only when: \itemize{ \item \code{correlation_matrix}
#' is not \code{NULL}, and \item \code{target_correlation = TRUE}. }
#' 
#' If \code{target_correlation = FALSE}, the supplied \code{correlation_matrix}
#' is used directly as the latent Gaussian correlation matrix, and
#' \code{calibration_control} is ignored. In general, it is recommended that
#' \code{target_correlation = TRUE}.
#' @section Practical guidance: \itemize{ \item In many applications, the
#' defaults have been observed to be sufficient.  \item If the estimated
#' correlation matrices differ substantially from the target, the increasing
#' \code{n_mc} will usually have the most impact. Note that increasing this
#' value will increase computation time.  \item When all simulated endpoints
#' are continuous Gaussian, target_correlation = FALSE may be used as Pearson
#' correlation is preserved under the normal marginal transformation and no
#' calibration is needed. This may be especially important in simulation
#' studies as this will save computation time.  %\item For publication-quality
#' simulation studies where observed Pearson correlation matching is important,
#' % increasing \code{n_mc} may improve stability.  %\item If calibration
#' becomes slow, consider reducing \code{n_mc}, loosening \code{tol}, or using
#' % \code{target_correlation = FALSE} when the latent Gaussian correlation is
#' acceptable.  \item Even with calibration, finite-sample summaries from the
#' simulated data will not match the target matrix exactly.  \item Please note
#' that not all combinations of endpoint marginals and target Pearson
#' correlations are feasible, so some requested correlation structures may be
#' unattainable under the specified data-generating distributions. }
#' @seealso \code{\link{makeData}} for the main simulation routine.
#' 
#' \code{\link{endpoint_details}} for endpoint definitions.
#' 
#' \code{\link{enrollment_details}} for follow-up and enrollment settings.
#' @name calibration_control
#' @keywords documentation
#' @examples
#' 
#' ## Default settings
#' calibration_control <- list(
#'   n_mc = 10000,
#'   tol = 0.001,
#'   maxit = 100,
#'   rho_cap = 0.999,
#'   ensure_pd = TRUE,
#'   conv_norm_type = "F"
#' )
#' 
#' ## Higher-precision example
#' calibration_control <- list(
#'   n_mc = 50000,
#'   tol = 1e-4,
#'   maxit = 200,
#'   rho_cap = 0.999,
#'   ensure_pd = TRUE,
#'   conv_norm_type = "F"
#' )
#' 
NULL
