\name{enrollment_details}
\alias{enrollment_details}
\title{Specification format for \code{enrollment_details} used by \code{makeData()}}

\description{
\code{enrollment_details} is a named list controlling stochastic subject enrollment
in \code{\link{makeData}}.

Enrollment times are represented on the trial calendar scale. That is,
\code{enrollTime} is the calendar time at which a subject enters the trial,
measured relative to the first randomized subject.

If omitted, \code{makeData()} fills in defaults internally:
\preformatted{
list(
  enrollment_distribution        = "none",
  enrollment_exponential_rate    = NULL,
  piecewise_enrollment_cutpoints = NULL,
  piecewise_enrollment_rates     = NULL
)
}
}

\section{Purpose}{

Enrollment determines when subjects enter the trial calendar. It is distinct
from subject-level follow-up limits and trial-level stopping rules.

The trial calendar is generally constructed in three conceptual steps:
\enumerate{
  \item enrollment is simulated using \code{enrollment_details};
  \item the trial end time is determined using \code{trial_end_details};
  \item subject-level available follow-up is derived using \code{followup_details}.
}

For a subject \eqn{i}, \code{enrollTime_i} is measured on the trial calendar
scale. Time-to-event variables such as \code{TTE_1} are measured from that
subject's own enrollment/randomization time.

If the trial ends at calendar time \eqn{\mathcal{T}}, the subject's available
administrative follow-up is generally:
\deqn{
  \max(0, \mathcal{T} - T_{E,i}),
}
where \eqn{T_{E,i}} is the subject's enrollment time. If
\code{followup_details$max_followup} is supplied, available follow-up is further
capped at that value.

If no enrollment model is supplied
(\code{enrollment_distribution = "none"}), all subjects are treated as enrolling
at trial calendar time 0.
}

\section{Arguments}{

\describe{
  \item{\code{enrollment_distribution}}{
    Character string giving the enrollment process. Must be one of:
    \describe{
      \item{\code{"none"}}{
        All subjects enroll at time 0.
      }
      \item{\code{"exponential"}}{
        Homogeneous Poisson-process accrual. Inter-enrollment gaps are
        exponentially distributed with rate
        \code{enrollment_exponential_rate}. Requires
        \code{enrollment_exponential_rate}.
      }
      \item{\code{"piecewise"}}{
        Piecewise homogeneous Poisson-process accrual. Accrual rates vary by
        calendar interval. Requires \code{piecewise_enrollment_cutpoints} and
        \code{piecewise_enrollment_rates}.
      }
    }
  }

  \item{\code{enrollment_exponential_rate}}{
    Numeric scalar greater than 0. Used only when
    \code{enrollment_distribution = "exponential"}.

    This is the expected accrual rate per unit calendar time. For example,
    \code{enrollment_exponential_rate = 10} corresponds to an average of 10
    subjects enrolled per time unit. Inter-enrollment gaps are generated from
    an exponential distribution with this rate, and enrollment times are formed
    as cumulative sums of those gaps.
  }

  \item{\code{piecewise_enrollment_cutpoints}}{
    Numeric vector of strictly increasing cutpoints defining calendar intervals
    for piecewise accrual. The first cutpoint should be 0.

    For example, \code{c(0, 8, 16, 24)} defines three intervals:
    \eqn{[0,8)}, \eqn{[8,16)}, and \eqn{[16,24]}.
  }

  \item{\code{piecewise_enrollment_rates}}{
    Numeric vector of non-negative accrual rates, one for each interval defined
    by \code{piecewise_enrollment_cutpoints}. Its length must be
    \code{length(piecewise_enrollment_cutpoints) - 1}.

    At least one rate must be greater than 0. A rate of 0 can be used to encode
    an interval with no accrual. If the finite piecewise accrual window does not
    generate enough subjects, \code{makeData()} continues enrollment beyond the
    final cutpoint using the final interval's rate. Therefore, if the final rate
    is 0 and the target sample size has not been reached, simulation will fail.
  }
}
}

\section{No stochastic enrollment}{

When \code{enrollment_distribution = "none"}, all subjects are assigned:

\preformatted{
enrollTime = 0
}

In this case, if a trial end rule is used, all subjects share the same calendar
entry time. For example, if \code{trial_end_details$type = "fixed_calendar"}
and \code{trial_end_time = 80}, each subject's calendar-derived available
follow-up is 80, subject to any \code{max_followup} cap.
}

\section{Exponential enrollment}{

When \code{enrollment_distribution = "exponential"}, enrollment follows a
homogeneous Poisson-process style accrual model.

The inter-enrollment gaps are generated as:
\deqn{
  G_i \sim \mathrm{Exponential}(\lambda),
}
where \eqn{\lambda} is \code{enrollment_exponential_rate}. Enrollment times are
then generated as cumulative sums:
\deqn{
  T_{E,i} = \sum_{m=1}^{i} G_m.
}

The enrollment clock is shifted so the first enrolled subject has
\code{enrollTime = 0}. In \code{makeData()}, these calendar enrollment times
are randomly assigned across generated subject rows to avoid confounding
treatment arm with enrollment order.

This model is appropriate when the expected enrollment rate is approximately
constant over calendar time.
}

\section{Piecewise enrollment}{

When \code{enrollment_distribution = "piecewise"}, enrollment follows a
piecewise homogeneous Poisson-process style accrual model.

For each calendar interval, the accrual rate is assumed constant. For interval
\eqn{k}, with length \eqn{L_k} and accrual rate \eqn{\lambda_k}, the number of
enrollments in the interval follows:
\deqn{
  N_k \sim \mathrm{Poisson}(\lambda_k L_k).
}

Conditional on \eqn{N_k}, enrollment times within that interval are generated
as sorted uniform times over the interval. This construction supports
start-up, ramp-up, plateau, slowdown, and temporary enrollment pause patterns.

For example:
\preformatted{
enrollment_details <- list(
  enrollment_distribution        = "piecewise",
  piecewise_enrollment_cutpoints = c(0, 8, 24, 52),
  piecewise_enrollment_rates     = c(2, 8, 12)
)
}

This represents expected accrual of:
\itemize{
  \item 2 subjects per time unit from 0 to 8,
  \item 8 subjects per time unit from 8 to 24,
  \item 12 subjects per time unit from 24 to 52.
}

If the requested sample size is not reached by the final cutpoint, enrollment
continues beyond the last cutpoint using the final interval's rate.
}

\section{Interaction with follow-up and trial-end rules}{

\code{enrollment_details} only determines when subjects enter the trial. It
does not by itself censor time-to-event endpoints.

Administrative follow-up is determined by combining enrollment times with
\code{trial_end_details} and \code{followup_details}. For example:

\preformatted{
followup_details <- list(
  min_followup = 52,
  max_followup = 156
)

trial_end_details <- list(
  type = "last_patient_min_followup"
)
}

Under this design, the trial ends when the last randomized subject has at least
52 time units of follow-up:
\deqn{
  \mathcal{T} = \max_i(T_{E,i}) + 52.
}

Subject-level available follow-up is then:
\deqn{
  \min\{156, \max(0, \mathcal{T} - T_{E,i})\}.
}

For a fixed-calendar trial:
\preformatted{
trial_end_details <- list(
  type = "fixed_calendar",
  trial_end_time = 80
)
}

the trial ends at calendar time 80, and subject-level available follow-up is:


\deqn{
  \min\{M, \max(0, 80 - T_{E,i})\},
}


Subjects with \code{enrollTime} after the trial end time have
\code{availableFollowup = 0}. These rows may represent planned subjects who
would not contribute positive follow-up under the realized trial calendar.
}

\section{Event-driven trials}{

For event-driven trials, enrollment times are used to convert subject-level
event times to calendar event times:
\deqn{
  \text{calendar event time}_i = T_{E,i} + T_i,
}
where \eqn{T_i} is the subject-level time-to-event value.

For example:
\preformatted{
trial_end_details <- list(
  type = "event_driven",
  event_endpoint = "TTE_1",
  target_events = 150
)
}

The trial end time is the calendar time at which the target number of observed
events is reached, subject to any additional options in
\code{trial_end_details}, such as \code{require_min_followup} or
\code{max_trial_duration}.


In event-driven trials, \code{makeData()} returns the planned simulated sample.
If the event target is reached before all planned subjects have enrolled, some
rows may have \code{availableFollowup = 0}. These rows represent subjects who
would not contribute positive follow-up under the realized trial calendar. For
most analyses, users should restrict to subjects with positive follow-up:
\preformatted{
dat_analysis <- subset(dat, availableFollowup > 0) # or equivalently
dat_analysis <- subset(dat, enrollTime < sim$meta$trial_calendar$trial_end_time)
}


}

\section{Deprecated behavior}{

Earlier versions of \code{makeData()} allowed
\code{enrollment_details$administrative_censoring} and
\code{enrollment_distribution = "uniform"}.

These are no longer part of the recommended trial-calendar interface.
Administrative censoring is now represented through the combination of
\code{followup_details} and \code{trial_end_details}. Uniform enrollment has
been removed in favor of Poisson-process based accrual models.

If \code{administrative_censoring} is supplied in \code{enrollment_details}, it
is ignored and a warning may be issued.
}

\seealso{
\code{\link{makeData}} for the main simulation function.

\code{\link{endpoint_details}} for endpoint specifications.

\code{\link{calibration_control}} for correlation-calibration settings.
}

\examples{
## No stochastic enrollment
enrollment_details <- list(
  enrollment_distribution = "none"
)

## Homogeneous Poisson-process enrollment
## Average accrual rate of 10 subjects per time unit
enrollment_details <- list(
  enrollment_distribution     = "exponential",
  enrollment_exponential_rate = 10
)

## Piecewise Poisson-process enrollment
enrollment_details <- list(
  enrollment_distribution        = "piecewise",
  piecewise_enrollment_cutpoints = c(0, 8, 24, 52),
  piecewise_enrollment_rates     = c(2, 8, 12)
)

## Example with last-patient-minimum-follow-up trial end
followup_details <- list(
  min_followup = 52,
  max_followup = 156
)

trial_end_details <- list(
  type = "last_patient_min_followup"
)

## Example with fixed calendar trial end
followup_details <- list(
  max_followup = 104
)

trial_end_details <- list(
  type = "fixed_calendar",
  trial_end_time = 80
)

## Example with event-driven trial end
trial_end_details <- list(
  type = "event_driven",
  event_endpoint = "TTE_1",
  target_events = 150
)
}

\keyword{documentation}
