#' Specification format for \code{followup_details} used by \code{makeData()}
#'
#' \code{followup_details} is a named list controlling subject-level
#' administrative follow-up in \code{\link{makeData}}.
#'
#' If omitted, \code{makeData()} fills in defaults internally:
#' \preformatted{
#' list(
#'   min_followup = NULL,
#'   max_followup = NULL
#' )
#' }
#'
#' @section Purpose:
#' \code{followup_details} determines how much calendar follow-up each subject
#' can contribute once enrollment and the trial end time have been determined.
#'
#' It is distinct from:
#' \itemize{
#'   \item \code{enrollment_details}, which determines when subjects enter the
#'   trial calendar, and
#'   \item \code{trial_end_details}, which determines when the trial stops on
#'   the calendar scale.
#' }
#'
#' These three components define the realized trial calendar for any endpoint
#' type:
#' \enumerate{
#'   \item \code{enrollment_details} generates \code{enrollTime};
#'   \item \code{trial_end_details} determines the realized trial end time
#'   \eqn{\mathcal{T}};
#'   \item \code{followup_details} determines each subject's
#'   \code{availableFollowup}.
#' }
#'
#' For time-to-event endpoints, this same calendar structure also determines
#' administrative censoring of event times and statuses.
#'
#' A subject enrolled at calendar time \eqn{T_{E,i}} can contribute at most:
#' \deqn{
#'   \max(0, \mathcal{T} - T_{E,i})
#' }
#' units of calendar follow-up before the trial ends. If
#' \code{max_followup} is supplied, this is further capped at:
#' \deqn{
#'   \min\{M, \max(0, \mathcal{T} - T_{E,i})\},
#' }
#' where \eqn{M} denotes \code{max_followup}.
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{min_followup}}{
#'     Optional non-negative numeric scalar. This is the minimum follow-up
#'     duration that certain trial-end rules may require for the last
#'     randomized subject.\cr \cr
#'     By itself, \code{min_followup} does not censor subjects and does not
#'     force the trial to continue. It becomes operational when used together
#'     with a trial-end rule such as
#'     \code{trial_end_details$type = "last_patient_min_followup"} or
#'     \code{trial_end_details$require_min_followup = TRUE} in an event-driven
#'     trial.
#'   }
#'
#'   \item{\code{max_followup}}{
#'     Optional positive numeric scalar giving the maximum observable follow-up
#'     for any subject.\cr \cr
#'     This affects the realized trial calendar for any endpoint type by
#'     capping \code{availableFollowup}. For time-to-event endpoints, it also
#'     acts as a subject-level administrative censoring cap. If a subject would
#'     otherwise remain under observation longer than \code{max_followup},
#'     their observed TTE time is truncated at \code{max_followup} and the
#'     corresponding event indicator is administratively censored when
#'     appropriate. For non-TTE endpoints, \code{max_followup} changes the
#'     returned calendar metadata but does not modify the endpoint values
#'     themselves.
#'   }
#' }
#'
#' @section Common use cases:
#' \describe{
#'   \item{Only \code{max_followup}}{
#'     Use this when each subject can be followed for at most a fixed amount of
#'     time, regardless of how long the trial remains open.
#'   }
#'
#'   \item{Only \code{min_followup}}{
#'     Use this when the design requires the last randomized subject to have at
#'     least a certain amount of follow-up before final analysis, for example
#'     with \code{trial_end_details$type = "last_patient_min_followup"}.
#'   }
#'
#'   \item{Both \code{min_followup} and \code{max_followup}}{
#'     Use this when the design requires a minimum observation window for late
#'     enrollees, while also capping how long early enrollees can remain under
#'     observation.
#'   }
#' }
#'
#' @section Interaction with trial-end rules:
#' The effect of \code{followup_details} depends on the trial-ending rule.
#'
#' For a fixed-calendar trial:
#' \preformatted{
#' trial_end_details <- list(
#'   type = "fixed_calendar",
#'   trial_end_time = 36
#' )
#' }
#' subjects enrolled earlier will usually have more potential follow-up than
#' subjects enrolled later, but each subject is still capped at
#' \code{max_followup} if it is supplied.
#'
#' For a last-patient-minimum-follow-up design:
#' \preformatted{
#' followup_details <- list(
#'   min_followup = 24,
#'   max_followup = 36
#' )
#'
#' trial_end_details <- list(
#'   type = "last_patient_min_followup"
#' )
#' }
#' the realized trial end is:
#' \deqn{
#'   \mathcal{T} = \max_i(T_{E,i}) + 24.
#' }
#'
#' For non-TTE endpoints, these rules still determine \code{enrollTime},
#' \code{availableFollowup}, and the metadata in
#' \code{sim$meta$trial_calendar}, even though the endpoint values themselves
#' are not administratively censored.
#'
#' For an event-driven design, \code{max_followup} determines whether late TTE
#' events remain observable, and \code{min_followup} can optionally be enforced
#' through \code{trial_end_details$require_min_followup = TRUE}. Event-driven
#' trial end is specific to TTE endpoints.
#'
#' @section Returned data:
#' When a trial-calendar feature is active, \code{makeData()} returns an
#' \code{availableFollowup} column in the simulated dataset. This stores the
#' realized administrative follow-up available to each subject under the
#' combination of enrollment timing, trial end time, and any
#' \code{max_followup} cap.
#'
#' Subjects with \code{availableFollowup == 0} are subjects whose planned
#' enrollment occurs at or after the realized trial end time. They remain in
#' the simulated planned dataset but typically do not contribute to the
#' realized analysis population.
#'
#' @seealso
#' \code{\link{makeData}} for the main simulation function.
#'
#' \code{\link{enrollment_details}} for enrollment-time settings.
#'
#' \code{\link{trial_end_details}} for calendar stopping rules.
#'
#' @name followup_details
#' @keywords documentation
#'
#' @examples
#' ## Follow each subject for at most 36 time units
#' followup_details <- list(
#'   max_followup = 36
#' )
#'
#' ## Require at least 24 time units of follow-up for the
#' ## last randomized subject
#' followup_details <- list(
#'   min_followup = 24
#' )
#'
#' ## Combine a minimum and maximum follow-up rule
#' followup_details <- list(
#'   min_followup = 24,
#'   max_followup = 36
#' )
#'
#' ## Small makeData() example
#' ep_tte <- list(
#'   endpoint_type = "tte",
#'   baseline_rate = rate_from_prob(
#'     target_prob = 0.30,
#'     mode = "admin",
#'     admin_time = 12
#'   ),
#'   trt_effect = log(0.70),
#'   censoring_rate = 0.01
#' )
#'
#' sim <- makeData(
#'   correlation_matrix = NULL,
#'   sample_size_per_group = 10,
#'   SEED = 1,
#'   endpoint_details = list(ep_tte),
#'   enrollment_details = list(
#'     enrollment_distribution = "exponential",
#'     enrollment_exponential_rate = 4
#'   ),
#'   followup_details = list(
#'     min_followup = 12,
#'     max_followup = 18
#'   ),
#'   trial_end_details = list(
#'     type = "last_patient_min_followup"
#'   )
#' )
#'
#' dat <- as.data.frame(sim)
#' head(dat[, c("TTE_1", "Status_1", "enrollTime", "availableFollowup")])
NULL
