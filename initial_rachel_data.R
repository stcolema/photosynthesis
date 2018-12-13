
# === Libraries ================================================================

library(nlme)
library(lattice)
library(lme4)

library(tidyverse) # install.packages("tidyverse", dep = T)
library(magrittr)
library(attempt)

# === Functions ================================================================

# === Grouping assignments (ASIDE) =============================================
# https://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line
# Generic form
"%=%" <- function(l, r, ...) UseMethod("%=%")

# Binary Operator
"%=%.lbunch" <- function(l, r, ...) {
  Envir <- as.environment(-1)

  if (length(r) > length(l)) {
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  }

  if (length(l) > length(r)) {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }

  for (II in 1:length(l)) {
    do.call("<-", list(l[[II]], r[[II]]), envir = Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)

  # Assume that destin is a length when it is a single number and source is not
  if (d == 1 && s > 1 && !is.null(as.numeric(destin))) {
    d <- destin
  }

  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d / s))[1:d]
  }
  return(source)
}

# Grouping the left hand side
g <- function(...) {
  List <- as.list(substitute(list(...)))[-1L]
  class(List) <- "lbunch"
  return(List)
}

update_aci_data <- function(data,
                            T_var,
                            response_var,
                            PAR_var,
                            R = 0.008314472,
                            O_pres = 210,
                            Kelvin = FALSE) {
  T_leaf <- enquo(T_var)
  response <- enquo(response_var)
  PAR <- enquo(PAR_var)
  new_data <- data %>%
    dplyr::mutate(
      O = O_pres,
      R = R,
      Kc = exp(38.05 - 79.43 / (R * (!!T_leaf + 273.15 * (1 - Kelvin)))),
      Ko = exp(20.30 - 36.38 / (R * (!!T_leaf + 273.15 * (1 - Kelvin)))),
      Gstar = exp(19.02 - 38.83 / (R * (!!T_leaf + 273.15 * (1 - Kelvin)))),
      x = Ci - Gstar,
      Photosynthesis = !!response,
      y = Photosynthesis,
      PAR = !!PAR,
      c1 = Gstar + Kc * (1 + O / Ko),
      c2 = 3 * Gstar,
      x1 = c1 / x,
      x2 = c2 / x,
      I2 = PARi * (1 - f) / 2,
      Cc = Ci - !!response / gm_initial
    )
}

#' @title Rubisco limited curve
#'
#' @description Returns the expected CO2 assimilation for the curve restricted by Rubisco
#' assuming a saturating supply of RuBP. Limitation by Rubisco is associated
#' with low CO2 rather than V_max. This is part of FvCB model.
#' @param V_cmax (\eqn{V_{cmax}}) is the maximum velocity of Rubisco for
#' carboxylation.
#' @param C_c (\eqn{C_C}) is the \eqn{CO_2} partial pressure at Rubisco (within
#' chloroplast).
#' @param G_star (\eqn{\Gamma^*}) is the \eqn{CO_2} at which oxygenation proceeds at
#' twice the rate of carboxylaiton causing photosynthetic uptake of \eqn{CO_2} to be
#' exactly compensated by photorespiratory \eqn{CO_2} release.
#' @param K_c (\eqn{K_C}) is the Michaelis constant of Rubisco for carbon dioxide.
#' @param O is the partial pressure of oxygen at Rubisco.
#' @param K_o (\eqn{K_O}) is the inhibition constant of Rubisco for oxygen.
#' @param R_d (\eqn{R_d}) is respiratory \eqn{CO_2} release by methods other
#' than by photorespiration and is presumed to be primarily mitochondiral
#' respiraion.
#' @return Nothing at all.
#' @examples
#' Rubisco_limited(<examples>)
#' More examples and whatnot
Rubisco_limited <- function(V_cmax, C_c, G_star, K_c, O, K_o, R_d) {
  numerator <- V_cmax * (C_c - G_star)
  denom <- C_c + K_c * (1 + O / K_o)
  Rubisco_curve <- (numerator / denom) - R_d
}

#' @title RuBP limited curve
#'
#' @description Returns the expected CO2 assimilation for the curve restricted
#' by Ribulose 1,5-bisphosphate (RuBP) regeneration. This equation assumes
#' four electrons per carboxylation and oxygenation. Based on the number of
#' electrons required for NADP+ reduction the conservative values of 4 and 8
#' are used as defaults for the \eqn{C_C} and \eqn{\Gamma^*} multipliers in the
#' denominator, but 4.5 and 10.5 have been used.
#' @param C_c (\eqn{C_c}) is the \eqn{CO_2} partial pressure at Rubisco (within
#' chloroplast).
#' @param G_star (\eqn{\Gamma^*}) is the \eqn{CO_2} at which oxygenation proceeds at
#' twice the rate of carboxylaiton causing photosynthetic uptake of \eqn{CO_2} to be
#' exactly compensated by photorespiratory \eqn{CO_2} release.
#' @param J provided here is the rate of electron transport going to support
#' NADP+ reduction for RuBP regeneration at the measurement light intensity. J
#' is sometimes used to estimate a maximum rate that could be obtained at
#' saturating light, and this is called \eqn{J_{max}}.
#' @return Nothing at all.
#' @examples
#' RuBP_limited(<examples>)
#' More examples and whatnot
RuBP_limited <- function(J, C_c, G_star, R_d,
                         C_c_denom_multiplier = 4,
                         G_star_denom_multiplier = 8) {
  numerator <- C_c - G_star
  denom <- C_c_denom_multiplier * C_c + G_star_denom_multiplier * G_star
  return(J * numerator / denom - R_d)
}

# Need to relate this - possibly replaces RuBP_limiting effect?
#' @title Non rectangular hyperbola (from Gerrit)
#'
#' @description Returns the a fit for a non-rectangular hyperbola curve
#' @param PAR I think this is the assimilation rate?
#' @param Jmax Something to do with electrons
#' @param alpha Just a parameter
#' @param theta Ditto
non_rectangular_hyperbola <- function(PAR, J_max, alpha, theta) {
  # Solve for root of equation:
  # theta * x2 - (Jmax + alpha * PAR) * x + alpha * PAR * Jmax = 0
  # model approximates theta and alpha I think? Possibly not
  a <- theta
  b <- -(J_max + alpha * PAR)
  c <- alpha * PAR * J_max
  determ <- b^2 - 4 * a * c
  determ <- pmax(determ, 0)
  x <- (-b - sqrt(determ)) / (2 * a)

  # below is what Gerrit compared to Rubisco curve for deciding current limiting
  # factor - possily not part of the NRH function, but a represenation for RuBP
  # (x / 4) * (C_c - G_star) / (C_c  + 2 * G_star)
}

#' @title Gerrit FvCB model
#' @description FvCB model used in Didy's paper
#'
FvCB <- function(C_c, PAR, G_star, K_c, K_o, O, V_cmax, J_max, alpha, theta, R_d) {
  p1 <- Rubisco_limited(V_cmax, C_c, G_star, K_c, O, K_o, R_d) + R_d
  J <- non_rectangular_hyperbola(PAR, J_max, alpha, theta)
  # p2 <- J * (C_c - G_star) / (4.5 * C_c  + 10.5 * G_star)
  p2 <- RuBP_limited(J, C_c, G_star, R_d)
  ifelse(p1 < p2, p1, p2) - R_d
}


# === Data =====================================================================

# Read in data
my_wd <- "C:/Users/steph/Desktop/Bioinformatics/MAT80436 - Thesis/Data - Rachel"
setwd(my_wd)

LRC1 <- na.omit(read.table("LRC1.txt", sep = "\t", header = TRUE, na.strings = ""))
Dark <- read.table("Dark-F.txt", sep = "\t", header = TRUE, na.strings = "")
LRC2_F <- na.omit(read.table("LRC2-F.txt", sep = "\t", header = TRUE, na.strings = ""))

# Delete empty lines and any DIV/o errors - also empty column in LRC2
ACI1 <- na.omit(read.table("A-Ci1_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI2 <- na.omit(read.table("A-Ci2_clean.csv", sep = ",", header = TRUE, na.strings = ""))
LRC2 <- na.omit(read.table("LRC2_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI1_F <- na.omit(read.table("A-Ci1-F_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI2_F <- na.omit(read.table("A-Ci2-F_clean.csv", sep = ",", header = TRUE, na.strings = ""))

# === EDA ======================================================================

# LRC1 @ 21% O2
# LRC2 @ 2% (linear part of curve, R_d and J_max)
# ACI1 @ 21%
# ACI2 @ 2%
# Flourescence factor used to calculate the efficiency of Photosystem II
# (\phi_{PSII} or (1 - f) in my report)
# one g_m for each experimental unit (plant ID) - have flouresence measurements
# to help with this
# PAR is the incident light (Rachel has calculated the aborped PAR new column)
# DARK used Dark adapted flourescence used to calculate R_d (different every
# plant)
# Estimate quantum efficieny of PSII from strictly limited light level

# Most important R_d, V_cmax, J_max, \theta, g_m

# Light from above (AD), below (AB) and both (ADAB)
# Treatments A1 (from above) and A2 (from above and below)
# Therefore 6 treatments (3 within A1, 3 within A2)

#

# Constants
R <- 0.008314472 # R is the concentration of free (unbound) RuP 2
O_21 <- 210
O_2 <- 20
Tleaf <- quo(Tleaf)
response <- quo(Photo_corr_diff)
PAR <- quo(Parin.total.1)
factors <- c("TREATMENT", "SIDE")

# Random effects
random_effects <- rlang::syms(c(paste(factors, collapse = "."), factors[1]))
n_random_effects <- length(random_effects)

# random = pdDiag(list(Vcmax ~ leaf, Jmax ~ leaf, theta ~ leaf)),

# Initial parameter estimates based on my pdf "ACI curves"
Vcmax_initial <- 80
Jmax_initial <- 1.6 * Vcmax_initial
theta_initial <- 0.7
f <- 0.15
alpha_initial <- 0.85 * (1 - f) / 2 # not used as Rachel has done this bit - using value of 1
Rd_initial <- 0.015 * Vcmax_initial
gm_initial <- 0.0045 * Vcmax_initial


# Create indices for subgroups - these correspond to the normal (light from
# above) vs treated (50% above, 50% below) plants
normal_ind <- LRC1$TREATMENT == "A1"
treated_ind <- LRC1$TREATMENT == "A2"


ACI2 %>%
  select(contains("PAR")) %>%
  glimpse()


ACI2 %>%
  select(contains("photo")) %>%
  glimpse()

# Update the datasets to contain all the variables of interest and to have more 
# generic names
LRC1 <- update_aci_data(LRC1, Tleaf, Photo, Parin.total.1, O_pres = O_21)
ACI1 <- update_aci_data(ACI1, Tleaf, Photo_corr_diff, PARabs, O_pres = O_21)

# The non-constant values
LRC1 <- LRC1 %>%
  dplyr::mutate(
    O = O_21,
    R = R,
    Kc = exp(38.05 - 79.43 / (R * (!!Tleaf + 273.15))),
    Ko = exp(20.30 - 36.38 / (R * (!!Tleaf + 273.15))),
    Gstar = exp(19.02 - 38.83 / (R * (!!Tleaf + 273.15))),
    x = Ci - Gstar,
    Photosynthesis = !!response,
    y = Photosynthesis,
    PAR = !!PAR,
    c1 = Gstar + Kc * (1 + O / Ko),
    c2 = 3 * Gstar,
    x1 = c1 / x,
    x2 = c2 / x,
    I2 = PARi * alpha_initial,
    Cc = Ci - !!response / gm_initial
  )

# Indices for subgroups
A1_ind <- LRC1$TREATMENT == "A1"
A2_ind <- LRC1$TREATMENT == "A2"

A1_AD_ind <- as.logical((LRC1$TREATMENT == "A1") * (LRC1$SIDE == "AD"))
A1_AB_ind <- as.logical((LRC1$TREATMENT == "A1") * (LRC1$SIDE == "AB"))
A1_ADAB_ind <- as.logical((LRC1$TREATMENT == "A1") * (LRC1$SIDE == "ADAB"))

A2_AD_ind <- as.logical((LRC1$TREATMENT == "A2") * (LRC1$SIDE == "AD"))
A2_AB_ind <- as.logical((LRC1$TREATMENT == "A2") * (LRC1$SIDE == "AB"))
A2_ADAB_ind <- as.logical((LRC1$TREATMENT == "A2") * (LRC1$SIDE == "ADAB"))

# --- NLS ----------------------------------------------------------------------

nls(y ~ FvCB(Cc, PAR, Gstar, Kc, Ko, O, Vcmax, Jmax, 1.0, theta, Rd),
  data = ACI1,
  subset = A1_ind,
  start = list(
    Vcmax = 76.8,
    Jmax = 100,
    Rd = 1.4,
    theta = theta_0 # ,
    # alpha = alpha_initial
  ),
  control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-10)
)

# Grid search (kind of) for NLS in ACI1 data for each treatment effect
nls_aci1 <- NULL
for (V_cmax_0 in c(10, 50, 70, 80, 90, 100, 120, 150, 200, 500, 1000)) {
  for (J_max_0 in c(10, 50, 100, 120, 200, 300, 500, 800, 1000, 5000)) {
    for (R_d_0 in c(1:5)) {
      for (theta_0 in c(0.5, 1.0, 1.5, 2.0)) {
        nls_aci1 <- attempt::try_catch(
          expr = nls(y ~
          FvCB(
            Cc,
            PAR,
            Gstar,
            Kc,
            Ko,
            O,
            Vcmax,
            Jmax,
            1.0,
            theta,
            Rd
          ),
          data = ACI1,
          subset = A1_ind,
          start = list(
            Vcmax = V_cmax_0,
            Jmax = J_max_0,
            Rd = R_d_0,
            theta = theta_0 # ,
            # alpha = alpha_initial
          ),
          control = list(maxiter = 250, minFactor = 1e-10, printEval = F, tol = 1e-10)
          ),
          .e = function(a) return(NULL),
          .w = function(a) return(NULL)
        )
        if (!is.null(nls_aci1)) {
          print(paste(
            "SUCCESS!!! Vcmax:", V_cmax_0,
            "Jmax:", J_max_0,
            "Rd:", R_d_0,
            "theta:", theta_0
          ))
          Vc_work <- V_cmax_0
          J_work <- J_max_0
          break
        } else {
          print(paste(
            "Failed. Vcmax:", V_cmax_0,
            "Jmax:", J_max_0,
            "Rd:", R_d_0,
            "theta:", theta_0
          ))
        }
      }
    }
  }
}

nls_aci2 <- NULL
for (V_cmax_0 in c(10, 50, 70, 80, 90, 100, 120, 150, 200, 500, 1000)) {
  for (J_max_0 in c(10, 50, 100, 120, 200, 300, 500, 800, 1000, 5000)) {
    for (R_d_0 in c(1:5)) {
      for (theta_0 in c(0.5, 1.0, 1.5, 2.0)) {
        nls_aci2 <- attempt::try_catch(
          expr = nls(y ~
          FvCB(
            Cc,
            PAR,
            Gstar,
            Kc,
            Ko,
            O,
            Vcmax,
            Jmax,
            1.0,
            theta,
            Rd
          ),
          data = ACI1,
          subset = A2_ind,
          start = list(
            Vcmax = V_cmax_0,
            Jmax = J_max_0,
            Rd = R_d_0,
            theta = theta_0 # ,
            # alpha = alpha_initial
          ),
          control = list(maxiter = 250, minFactor = 1e-10, printEval = F, tol = 1e-10)
          ),
          .e = function(a) return(NULL),
          .w = function(a) return(NULL)
        )
        if (!is.null(nls_aci2)) {
          print(paste(
            "SUCCESS!!! Vcmax:", V_cmax_0,
            "Jmax:", J_max_0,
            "Rd:", R_d_0,
            "theta:", theta_0
          ))
          Vc_work <- V_cmax_0
          J_work <- J_max_0
          break
        } else {
          print(paste(
            "Failed. Vcmax:", V_cmax_0,
            "Jmax:", J_max_0,
            "Rd:", R_d_0,
            "theta:", theta_0
          ))
        }
      }
    }
  }
}

nls_aci1



# --- NLME ---------------------------------------------------------------------

CO2_grouped <- groupedData(y ~ 1 | ID,
  outer = ~TREAT, # !!random_effects[[n_random_effects]]
  inner = ~SIDE,
  data = ACI1
)

# Grid search for NLME
output <- list()
num_successes <- 0
nlme_aci1 <- NULL
for (V_cmax_0 in c(70, 90, 120, 150, 250)) {
  for (J_max_0 in c(50, 100, 120, 250, 500)) {
    for (R_d_0 in c(1.0, 2.0)) {
      for (V_cmax_ratio in c(0, 0.25, 0.50, 0.75, -0.25, -0.5, -0.75)) {
        for (J_max_ratio in c(0, 0.25, 0.50, 0.75, -0.25, -0.5, -0.75)) {
          nlme_aci1 <- attempt::try_catch(
            expr = nlme(y ~ FvCB(
              Cc,
              PAR,
              Gstar,
              Kc,
              Ko,
              O,
              Vcmax,
              Jmax,
              1.0,
              0.7,
              Rd
            ),
            data = CO2_grouped,
            fixed = list(Vcmax ~ TREAT, Jmax ~ TREAT, Rd ~ TREAT),
            random = pdSymm(list(Vcmax ~ SIDE, Jmax ~ 1, Rd ~ 1)),
            start = c(
              Vcmax = c(V_cmax_0, V_cmax_0 * V_cmax_ratio),
              Jmax = c(J_max_0, J_max_0 * J_max_ratio),
              Rd = c(R_d_0, 0)
            ),
            control = list(maxIter = 250, msVerbose = F, tolerance = 1e-4)
            ),
            .e = function(a) return(NULL),
            .w = function(a) return(NULL)
          )
          if (!is.null(nls_aci1)) {
            print(paste(
              "SUCCESS!!! Vcmax:", V_cmax_0,
              "Jmax:", J_max_0,
              "Rd:", R_d_0
            ))
            num_successes <- num_successes + 1
            level_name <- paste0("model_", num_successes)
            output[[level_name]] <- list()
            output[[level_name]]$model <- nlme_aci1
            output[[level_name]]$V_cmax <- V_cmax_0
            output[[level_name]]$J_max <- J_max_0
            output[[level_name]]$Rd <- R_d_0
            output[[level_name]]$V_ratio <- V_cmax_ratio
            output[[level_name]]$J_ratio <- J_max_ratio
            # Vc_work <- V_cmax_0
            # J_work <- J_max_0
            # break
          } else {
            print(paste(
              "Failed. Vcmax:", V_cmax_0,
              "Jmax:", J_max_0,
              "Rd:", R_d_0
            ))
          }
        }
      }
    }
  }
}


ggplot(data = Dark[-20, ], aes(x = Ci, y = Photo)) +
  # facet_wrap(vars(ID)) +
  geom_point() +
  geom_smooth()

ACI2 %>%
  dplyr::select(Ci, PAR, Gstar, Kc, Ko, O) %>%
  summary()

curve <- photosynthesisR::FvCB(
  LRC1$Ci,
  LRC1$PAR,
  LRC1$Gstar,
  LRC1$Kc,
  LRC1$Ko,
  LRC1$O,
  Vcmax_initial,
  Jmax_initial,
  alpha_initial,
  theta_initial,
  Rd_initial
)

LRC_test <- LRC1
LRC_test$FvCB <- curve
LRC_test$Cc

ggplot(data = LRC_test[LRC_test$ID == 1, ], aes(x = Ci, y = Photosynthesis)) +
  geom_point() +
  geom_point(aes(x = Cc, y = FvCB), colour = "blue")

CO2light



nls.all <- nls(y ~ FvCB2(x1, x2, PAR, Vcmax, Jmax, alpha, theta, Rd),
  data = LRC1,
  start = list(
    Vcmax = 76.8, Jmax = 100, alpha = 0.5, theta = -0.8, Rd = 1.4
    # Vm.d = 20, Jm.d = 300, al.d = 0, th.d = -1, Rd.d = -0.3
  ),
  control = list(maxiter = 250, minFactor = 0.000001, printEval = T, tol = 1e-04)
)
