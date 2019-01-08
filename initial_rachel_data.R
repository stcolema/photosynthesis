
# === Libraries ================================================================

library(nlme)
library(lattice)
library(lme4)

library(tidyverse) # install.packages("tidyverse", dep = T)
library(rlang)
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

# Add a new level to a factor
add_level <- function(x, new_level) {
  if (is.factor(x)) return(factor(x, levels = c(levels(x), new_level)))
  return(x)
}


update_aci_data <- function(data,
                            T_var,
                            response_var,
                            PAR_var,
                            Ci_var,
                            R = 0.008314472,
                            O_pres = 210,
                            Kelvin = FALSE,
                            reduce = FALSE,
                            other_vars_to_keep = NULL) {
  vars_of_interest <- enquos(T_var, response_var, PAR_var, Ci_var)

  kept_vars <- enquo(other_vars_to_keep)
  if (!quo_is_null(kept_vars)) {
    vars_of_interest <- c(vars_of_interest, kept_vars)
  }

  if (reduce) {
    data <- data %>%
      dplyr::select(!!!vars_of_interest)
  }

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
      # x = Ci - Gstar,
      Photosynthesis = !!response,
      # y = Photosynthesis,
      PAR = !!PAR # ,
      # c1 = Gstar + Kc * (1 + O / Ko),
      # c2 = 3 * Gstar,
      # x1 = c1 / x,
      # x2 = c2 / x,
      # I2 = PAR * (1 - f) / 2,
      # Cc = Ci - !!response / gm_initial
    )

  if (reduce) {
    new_data <- new_data %>%
      dplyr::select(-!!PAR) %>%
      dplyr::select(-!!response)
  }
  new_data
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

# === Method ===================================================================

# Estimate quantum efficieny of PSII from strictly limited light level
# (flourescence measurements)
# Flourescence factor used to calculate the efficiency of Photosystem II
# (\phi_{PSII} or (1 - f) in my report) - this should be plant specific

# Ideally calculate dark respiration for correction but does not tie in to model

# Calculate J by modelling Photosynthesis to I_inc * phi_psII * 0.25 (model for
# next line)
# use J = I_inc * phi_psII * slope of the above model
# J can be caluclated using the non-rectangular hyperbola too

# 1. LRC2 @ 2% (linear part of curve, R_d and J_max)
#      Lower intensities and O2 used to estimate J
#
# 2. ACI2 (@ 2% O2)
#      find compensation point (Gstar) per plant (Rd is not present as O2 is low)
#
# 3. LRC1 (@ ambient O2)
#      Calculate Cc using variable J
#

# Read Yin (2009) variable J method to get from J to C_c to g_m
# Yin adapts the variable J method in another paper from 2009 to estimate
# ambient CO2 levels

# Calculate gm using ambient CO2 ACI curve at low oxygen (ACI2) (at constant
# intensity) using variable J method

# one g_m for each experimental unit (plant ID) - have flouresence measurements
# to help with this

# Photosynthesic Absorbed Radiation (PAR) is the incident light
# (Rachel has calculated the aborped PAR new column)

# DARK used Dark adapted flourescence used to calculate R_dark (different to R_d)
# (different every plant)

# Most important R_d, V_cmax, J_max, \theta, g_m

# Combine ACI1 and LRC1 treating parameters as constant (as per Didy dataset)

# === Data =====================================================================

# Read in data
# set working directory based on computer
name <- Sys.info()["nodename"]
if (name == "D0132908") {
  my_wd <- "M:/My Documents/MAT80436 - Thesis/NLS_and_NLME_on_Rachel%27s_data"
} else {
  my_wd <- "C:/Users/steph/Desktop/Bioinformatics/MAT80436 - Thesis/Data - Rachel"
}
setwd(my_wd)

# LRC1 @ 21% O2
# LRC2 @ 2% O2
# ACI1 @ 21%
# ACI2 @ 2%

LRC1 <- na.omit(read.table("LRC1.txt", sep = "\t", header = TRUE, na.strings = ""))
Dark <- read.table("Dark-F.txt", sep = "\t", header = TRUE, na.strings = "")

# Delete empty lines and any DIV/o errors - also empty column in LRC2
ACI1 <- na.omit(read.table("A-Ci1_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI2 <- na.omit(read.table("A-Ci2_clean.csv", sep = ",", header = TRUE, na.strings = ""))
LRC2 <- na.omit(read.table("LRC2_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI1_F <- na.omit(read.table("A-Ci1-F_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI2_F <- na.omit(read.table("A-Ci2-F_clean.csv", sep = ",", header = TRUE, na.strings = ""))

# Light from above (AD), below (AB) and both (ADAB (not present in ACI))
# Treatments A1 (from above) and A2 (from above and below)
# Therefore 6 treatments (3 within A1, 3 within A2)


# Method of recording measurements (SIDE) is on the same individuals
# So plant 1 will have readings for AD, AB and ADAB for the same thing

# === EDA ======================================================================

# Constants
R <- 0.008314472 # R is the concentration of free (unbound) RuP 2
O_21 <- 210
O_2 <- 20
# Tleaf <- quo(Tleaf)
# response <- quo(Photo_corr_diff)
# PAR <- quo(Parin.total.1)
# factors <- c("TREATMENT", "SIDE")

# Random effects
# random_effects <- rlang::syms(c(paste(factors, collapse = "."), factors[1]))
# n_random_effects <- length(random_effects)

# random = pdDiag(list(Vcmax ~ leaf, Jmax ~ leaf, theta ~ leaf)),

# Initial parameter estimates based on my pdf "ACI curves"
Vcmax_initial <- 80
Jmax_initial <- 1.6 * Vcmax_initial
theta_initial <- 0.7
f <- 0.15
abs <- 0.85 # Rachel has corrected for absoprtance already so set this to 1
alpha_initial <- 1 * (1 - f) / 2 # normally by abs rather than 1
Rd_initial <- 0.015 * Vcmax_initial
gm_initial <- 0.0045 * Vcmax_initial


# Create indices for subgroups - these correspond to the normal (light from
# above) vs treated (50% above, 50% below) plants
normal_ind <- LRC1$TREATMENT == "A1"
treated_ind <- LRC1$TREATMENT == "A2"


LRC1 %>%
  select(contains("PAR")) %>%
  glimpse()


ACI1 %>%
  select(contains("photo")) %>%
  glimpse()

common_names_to_keep <- quos(MEAS, PLOT, ID, SIDE, TREATMENT, Photo, Ci, Tleaf)
LRC1_common <- LRC1 %>%
  dplyr::select(!!!common_names_to_keep)

# "TREAT"
# "PARabs"
# "Photo_corr_diff"
# "Parin.total"

# Update the datasets to contain all the variables of interest and to have more
# generic names

# === Data transformation ======================================================

LRC1_new <- update_aci_data(LRC1, Tleaf, Photo, Parin.total.1, Ci,
  O_pres = O_21,
  reduce = T,
  Kelvin = F,
  other_vars_to_keep = c(SIDE, TREATMENT, ID)
)

ACI1_new <- update_aci_data(ACI1, Tleaf, Photo_corr_diff, PARabs, Ci,
  O_pres = O_21,
  reduce = T,
  Kelvin = F,
  other_vars_to_keep = c(SIDE, TREAT, ID)
)
names(ACI1_new)
names(LRC1_new)

ACI1_new <- ACI1_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("CO2", labels = c("CO2")))

ACI1_new$SIDE <- add_level(ACI1_new$SIDE, "ADAB")
ACI1_new$Response <- add_level(ACI1_new$Response, "light")

LRC1_new <- LRC1_new %>%
  dplyr::mutate(Response = factor("light", labels = c("light")))
LRC1_new$Response <- add_level(LRC1_new$Response, "CO2")

combined_data <- dplyr::bind_rows(LRC1_new, ACI1_new)

# Indices for subgroups
select_A1 <- combined_data$TREATMENT == "A1"
select_A2 <- combined_data$TREATMENT == "A2"

select_A1_AD <- as.logical((combined_data$TREATMENT == "A1") * (combined_data$SIDE == "AD"))
select_A1_AB <- as.logical((combined_data$TREATMENT == "A1") * (combined_data$SIDE == "AB"))
select_A1_ADAB <- as.logical((combined_data$TREATMENT == "A1") * (combined_data$SIDE == "ADAB"))

select_A2_AD <- as.logical((combined_data$TREATMENT == "A2") * (combined_data$SIDE == "AD"))
select_A2_AB <- as.logical((combined_data$TREATMENT == "A2") * (combined_data$SIDE == "AB"))
select_A2_ADAB <- as.logical((combined_data$TREATMENT == "A2") * (combined_data$SIDE == "ADAB"))

select_CO2 <- combined_data$Response == "CO2"
select_light <- combined_data$Response == "light"

CO2_grouped <- groupedData(Photosynthesis ~ 1 | ID,
  outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
  inner = ~Ci + PAR + SIDE,
  data = combined_data
)

# === Transform 2 ==============================================================

LRC2 %>%
  select(contains("PAR")) %>%
  glimpse()

LRC2 %>%
  select(contains("photo")) %>%
  glimpse()

ACI2 %>%
  select(contains("PAR")) %>%
  glimpse()

ACI2 %>%
  select(contains("photo")) %>%
  glimpse()

LRC2_new <- update_aci_data(LRC2, Tleaf, Photo, PARabs, Ci,
                            O_pres = O_2,
                            reduce = T,
                            Kelvin = F,
                            other_vars_to_keep = c(SIDE, TREAT, ID)
)

ACI2_new <- update_aci_data(ACI2, Tleaf, Photo_corr_diff, PARabs, Ci,
                            O_pres = O_2,
                            reduce = T,
                            Kelvin = F,
                            other_vars_to_keep = c(SIDE, TREAT, ID)
)
names(ACI2_new)
names(LRC2_new)

summary(ACI2_new)
summary(LRC2_new)

ACI2_new <- ACI2_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("CO2", labels = c("CO2")))

ACI2_new$SIDE <- add_level(ACI2_new$SIDE, "ADAB")
ACI2_new$Response <- add_level(ACI2_new$Response, "light")

LRC2_new <- LRC2_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("light", labels = c("light")))

LRC2_new$SIDE <- add_level(LRC2_new$SIDE, "ADAB")
LRC2_new$Response <- add_level(LRC2_new$Response, "CO2")

combined_data_2 <- dplyr::bind_rows(LRC2_new, ACI2_new)

total_data <- dplyr::bind_rows(combined_data, combined_data_2)
# total_data$O <- factor(total_data$O)
summary(total_data)

CO2_grouped_full <- groupedData(Photosynthesis ~ 1 | ID,
                                outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
                                inner = ~Ci + PAR + SIDE,
                                data = total_data
)


# Indices for subgroups
select_A1_full <- total_data$TREATMENT == "A1"
select_A2_full <- total_data$TREATMENT == "A2"

select_A1_AD_full <- as.logical((total_data$TREATMENT == "A1") * (total_data$SIDE == "AD"))
select_A1_AB_full <- as.logical((total_data$TREATMENT == "A1") * (total_data$SIDE == "AB"))
select_A1_ADAB_full <- as.logical((total_data$TREATMENT == "A1") * (total_data$SIDE == "ADAB"))

select_A2_AD_full <- as.logical((total_data$TREATMENT == "A2") * (total_data$SIDE == "AD"))
select_A2_AB_full <- as.logical((total_data$TREATMENT == "A2") * (total_data$SIDE == "AB"))
select_A2_ADAB_full <- as.logical((total_data$TREATMENT == "A2") * (total_data$SIDE == "ADAB"))

select_CO2_full <- total_data$Response == "CO2"
select_light_full <- total_data$Response == "light"

select_O_2 <- total_data$O == 20
select_O_21 <- total_data$O == 210

# --- NLS ----------------------------------------------------------------------
names(combined_data)

nls_aci1 <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A1,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)

nls_aci2 <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A2,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)


model1_coef <- summary(nls_aci1)$coefficients
model2_coef <- summary(nls_aci2)$coefficients

comparison_nls <- bind_rows(
  as.data.frame(model1_coef),
  as.data.frame(model2_coef)
)

comparison_nls$Feature <- c(rownames(model1_coef), rownames(model2_coef))
comparison_nls$Treatment <- factor(c(
  rep("A1", nrow(model1_coef)),
  rep("A2", nrow(model2_coef))
),
level = c("A1", "A2")
)

comparison_nls_plot <- comparison_nls %>%
  mutate(Std_error = `Std. Error`) %>%
  select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nls_plot, aes(x = Feature, color = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std_error, ymin = Estimate - Std_error),
    position = "dodge"
  ) +
  labs(
    title = "NLS model estimate",
    y = "Value", x = "Parameters"
  )


model1_coef %<>%
  t %>%
  as.data.frame()

model2_coef %<>%
  t %>%
  as.data.frame()

g(Vcmax1, Jmax1, alpha1, theta1, Rd1) %=% model1_coef[1, ]
g(Vcmax2, Jmax2, alpha2, theta2, Rd2) %=% model2_coef[1, ]



#
# # Grid search (kind of) for NLS in ACI1 data for each treatment effect
# nls_aci1 <- NULL
# for (V_cmax_0 in c(50, 70, 90, 120, 200, 500, 1000)) {
#   for (J_max_0 in c(50, 100, 120, 200, 500, 800, 1000)) {
#     for (R_d_0 in c(1:5)) {
#       for (theta_0 in c(0.5, 1.0, 1.5, 2.0)) {
#         nls_aci1 <- attempt::try_catch(
#           expr = nls(y ~
#           FvCB(
#             Cc,
#             PAR,
#             Gstar,
#             Kc,
#             Ko,
#             O,
#             Vcmax,
#             Jmax,
#             1.0,
#             theta,
#             Rd
#           ),
#           data = combined_data,
#           subset = select_A1,
#           start = list(
#             Vcmax = V_cmax_0,
#             Jmax = J_max_0,
#             Rd = R_d_0,
#             theta = theta_0 # ,
#             # alpha = alpha_initial
#           ),
#           control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-10)
#           ),
#           .e = function(a) return(NULL),
#           .w = function(a) return(NULL)
#         )
#         if (!is.null(nls_aci1)) {
#           print(paste(
#             "SUCCESS!!! Vcmax:", V_cmax_0,
#             "Jmax:", J_max_0,
#             "Rd:", R_d_0,
#             "theta:", theta_0
#           ))
#           Vc_work <- V_cmax_0
#           J_work <- J_max_0
#           break
#         } else {
#           print(paste(
#             "Failed. Vcmax:", V_cmax_0,
#             "Jmax:", J_max_0,
#             "Rd:", R_d_0,
#             "theta:", theta_0
#           ))
#         }
#       }
#     }
#   }
# }

# --- NLS SIDE & TREATMENT -----------------------------------------------------

# Does not run
nls_A1_AB <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A1_AB,
start = list(
  Vcmax = Vcmax1, # with initial estimates still hit singular gradient
  Jmax = Jmax1,
  Rd = Rd1,
  theta = theta1,
  alpha = alpha1
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)

# Runs
nls_A1_AD <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A1_AD,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)

# Does not run (A1 ADAB does not exist)
nls_A1_ADAB <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A1_ADAB,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)

# Runs
nls_A2_AB <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A2_AB,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)

# Runs
nls_A2_AD <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A2_AD,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)

# Does not run
nls_A2_ADAB <- nls(Photosynthesis ~
FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = combined_data,
subset = select_A2_ADAB,
start = list(
  Vcmax = Vcmax_initial,
  Jmax = Jmax_initial,
  Rd = Rd_initial,
  theta = theta_initial,
  alpha = alpha_initial
),
control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)



nls_A1_AD
nls_A2_AB
nls_A2_AD

nls_A1_AD_coef <- summary(nls_A1_AD)$coefficients
nls_A2_AB_coef <- summary(nls_A2_AB)$coefficients
nls_A2_AD_coef <- summary(nls_A2_AD)$coefficients

comparison_nls_side <- bind_rows(
  as.data.frame(nls_A1_AD_coef),
  as.data.frame(nls_A2_AB_coef),
  as.data.frame(nls_A2_AD_coef)
)

comparison_nls_side$Feature <- c(
  rownames(nls_A1_AD_coef),
  rownames(nls_A2_AB_coef),
  rownames(nls_A2_AD_coef)
)

comparison_nls_side$Treatment <- factor(c(
  rep("A1", nrow(nls_A1_AD_coef)),
  rep("A2", nrow(nls_A2_AB_coef)),
  rep("A2", nrow(nls_A2_AD_coef))
),
level = c("A1", "A2")
)

comparison_nls_side$Side <- factor(c(
  rep("AD", nrow(nls_A1_AD_coef)),
  rep("AB", nrow(nls_A2_AB_coef)),
  rep("AD", nrow(nls_A2_AD_coef))
),
level = c("AB", "AD")
)

comparison_nls_side_plot <- comparison_nls_side %>%
  mutate(Std_error = `Std. Error`) %>%
  select(Estimate, Std_error, Feature, Treatment, Side)

ggplot(comparison_nls_side_plot, aes(x = Feature, color = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std_error, ymin = Estimate - Std_error),
    position = "dodge"
  ) +
  labs(
    title = "NLS model estimate",
    y = "Value", x = "Parameters"
  )

more_comparison <- comparison_nls_plot
more_comparison$Side <- "All"

more_comparison$Side <- add_level(more_comparison$Side, "AD")
more_comparison$Side <- add_level(more_comparison$Side, "AB")

more_comparison <- bind_rows(more_comparison, comparison_nls_side_plot)

ggplot(more_comparison, aes(x = Feature, color = Side)) +
  geom_errorbar(aes(ymax = Estimate + Std_error, ymin = Estimate - Std_error),
    position = "dodge"
  ) +
  facet_wrap(~Treatment, ncol = 2) +
  labs(
    title = "NLS model estimate",
    y = "Value", x = "Parameters"
  )

# --- NLME ---------------------------------------------------------------------

model_1 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ TREATMENT, alpha ~ TREATMENT, theta ~ TREATMENT),
random = pdDiag(list(Vcmax ~ SIDE, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 25
)
)

model_2 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT,
  Jmax ~ TREATMENT,
  Rd ~ TREATMENT,
  alpha ~ TREATMENT,
  theta ~ TREATMENT
),

random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

summary(model_2)

model_2a <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT,
  Jmax ~ TREATMENT,
  Rd ~ 1,
  alpha ~ TREATMENT,
  theta ~ TREATMENT
),

random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# === Model 2 investigation ====================================================

model_2_sum <- summary(model_2)


nlme_2_fixed <- model_2_sum$tTable
nlme_2_random <- summary(model_2)$coefficients$random

# Compare estimates for Treatments
fixed_comparison <- as.data.frame(nlme_2_fixed)
fixed_comparison$Feature <- rownames(fixed_comparison)
fixed_comparison$Feature <- c(
  rep("Vcmax", 2),
  rep("Jmax", 2),
  rep("Rd", 2),
  rep("alpha", 2),
  rep("theta", 2)
)

fixed_comparison$Treatment <- rep(c("A1", "A2"), nrow(fixed_comparison) / 2)


comparison_nlme_plot <- fixed_comparison %>%
  mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nlme_plot, aes(x = Feature, colour = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
    position = "dodge"
  ) +
  labs(
    title = "NLME model estimate",
    y = "Value", x = "Parameters"
  )

# more precision
model_2_precision <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT,
  Jmax ~ TREATMENT,
  Rd ~ TREATMENT,
  alpha ~ TREATMENT,
  theta ~ TREATMENT
),

random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = F, tolerance = 1e-3, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# === Back to NLME =============================================================


model_2b <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ 1,
  Jmax ~ 1,
  Rd ~ 1,
  alpha ~ 1,
  theta ~ 1
),
random = pdDiag(list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
)),
start = c(
  Vcmax = c(Vcmax1),
  Jmax = c(Jmax1),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

model_3 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ TREATMENT, alpha ~ TREATMENT, theta ~ TREATMENT),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = F, tolerance = 1e-1, msMaxIter = 100,
  pnlsTol = 1e-10, pnlsMaxIter = 25
)
)

# Try using output from model 2 with Rd independent of treatment
model_2_fixed_out <- fixed.effects(model_2)

model_4 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ TREATMENT, theta ~ TREATMENT),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# === Try SIDE rather than TREATMENT ===========================================

factorC <- with(CO2_grouped, interaction(TREATMENT, SIDE))
unique(factorC)

model_5 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, 0, 0, -30, -30, -30),
  Jmax = c(Jmax1, 0, 0, -10, -10, -10),
  Rd = c(Rd1, 0, 0),
  alpha = c(alpha1, 0, 0, 0, 0, 0),
  theta = c(theta1, 0, 0, 0, 0, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

model_6 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
),
random = pdDiag(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, 0, 0, -30, -30, -30),
  Jmax = c(Jmax1, 0, 0, -10, -10, -10),
  Rd = c(Rd1, 0, 0),
  alpha = c(alpha1, 0, 0, 0, 0, 0),
  theta = c(theta1, 0, 0, 0, 0, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

model_7 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1, 0, 0, -30, -30, -30),
  Jmax = c(Jmax1, 0, 0, -10, -10, -10),
  Rd = c(Rd1, 0, 0),
  alpha = c(alpha1, 0, 0, 0, 0, 0),
  theta = c(theta1, 0, 0, 0, 0, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)

model_8 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT,
  theta ~ TREATMENT
),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1, 0, 0, -30, -30, -30),
  Jmax = c(Jmax1, 0, 0, -10, -10, -10),
  Rd = c(Rd1, 0, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)

model_9 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ 1,
  alpha ~ TREATMENT,
  theta ~ TREATMENT
),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1, 0, 0, -30, -30, -30),
  Jmax = c(Jmax1, 0, 0, -10, -10, -10),
  Rd = c(Rd1),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)

model_10 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT,
  Rd ~ SIDE,
  alpha ~ TREATMENT,
  theta ~ TREATMENT
),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1, 0, 0, -30, -30, -30),
  Jmax = c(Jmax1, 20),
  Rd = c(Rd1, 0, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)
# === Try subsetting ===========================================================


# Try subsetting

model_A1_AD <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1),
  Jmax = c(Jmax1),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
subset = select_A1_AD,
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-5, msMaxIter = 100,
  pnlsTol = 1e-10, pnlsMaxIter = 25
)
)

# Try constant theta and alpha

# Try fixed ~ 1
model_const_ath <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  0.9654839,
  -0.213539,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1),
random = Vcmax + Jmax + Rd ~ 1,
start = c(
  Vcmax = c(Vcmax1),
  Jmax = c(Jmax1),
  Rd = c(Rd1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-5, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 25
)
)

summary(model_const_ath)

model_const_ath_TREAT <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  0.9654839,
  -0.213539,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1),
random = Vcmax + Jmax + Rd ~ 1,
start = c(
  Vcmax = c(51.42175, 0),
  Jmax = c(96.39980, 0),
  Rd = c(2.80716)
),
control = list(
  maxIter = 1000, msVerbose = T, tolerance = 1e-3, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 25
)
)

const_ath_Diag <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  0.9654839,
  -0.213539,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1),
random = pdDiag(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1)),
start = c(
  Vcmax = 51.42175,
  Jmax = 96.39980,
  Rd = 2.80716
),
# subset = select_A1,
control = list(
  maxIter = 1000, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-200, pnlsMaxIter = 50
)
)

# # Causes crashing if subsetting used
const_ath_symm <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  0.9654839,
  -0.213539,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1)),
start = c(
  Vcmax = 51.42175,
  Jmax = 96.39980,
  Rd = 2.80716
),
# subset = select_A1,
control = list(
  maxIter = 1000, msVerbose = F, tolerance = 1e-2, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)

# Causes crashing on my laptop
const_ath_symm <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  0.9654839,
  -0.213539,
  Rd
),
data = CO2_grouped,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1)),
start = c(
  Vcmax = c(51.421725, 10),
  Jmax = c(96.399792, 10),
  Rd = 2.807155
),
control = list(
  maxIter = 1000, msVerbose = T, tolerance = 1e-2, msMaxIter = 250,
  pnlsTol = 1e-20, pnlsMaxIter = 150
)
)


# # Grid search for NLME
# output <- list()
# num_successes <- 0
# nlme_aci1 <- NULL
# for (V_cmax_0 in c(70, 90, 120, 150, 250)) {
#   for (J_max_0 in c(50, 100, 120, 250, 500)) {
#     for (R_d_0 in c(1.0, 2.0)) {
#       for (V_cmax_ratio in c(0, 0.25, 0.50, 0.75, -0.25, -0.5, -0.75)) {
#         for (J_max_ratio in c(0, 0.25, 0.50, 0.75, -0.25, -0.5, -0.75)) {
#           nlme_aci1 <- attempt::try_catch(
#             expr = nlme(Photosynthesis ~ FvCB(
#               Ci,
#               PAR,
#               Gstar,
#               Kc,
#               Ko,
#               O,
#               Vcmax,
#               Jmax,
#               1.0,
#               0.7,
#               Rd
#             ),
#             data = CO2_grouped,
#             fixed = list(Vcmax ~ TREAT, Jmax ~ TREAT, Rd ~ TREAT),
#             random = pdSymm(list(Vcmax ~ SIDE, Jmax ~ 1, Rd ~ 1)),
#             start = c(
#               Vcmax = c(V_cmax_0, V_cmax_0 * V_cmax_ratio),
#               Jmax = c(J_max_0, J_max_0 * J_max_ratio),
#               Rd = c(R_d_0, 0)
#             ),
#             control = list(maxIter = 250, msVerbose = F, tolerance = 1e-4)
#             ),
#             .e = function(a) return(NULL),
#             .w = function(a) return(NULL)
#           )
#           if (!is.null(nlme_aci1)) {
#             print(paste(
#               "SUCCESS!!! Vcmax:", V_cmax_0,
#               "Jmax:", J_max_0,
#               "Rd:", R_d_0
#             ))
#             num_successes <- num_successes + 1
#             level_name <- paste0("model_", num_successes)
#             output[[level_name]] <- list()
#             output[[level_name]]$model <- nlme_aci1
#             output[[level_name]]$V_cmax <- V_cmax_0
#             output[[level_name]]$J_max <- J_max_0
#             output[[level_name]]$Rd <- R_d_0
#             output[[level_name]]$V_ratio <- V_cmax_ratio
#             output[[level_name]]$J_ratio <- J_max_ratio
#             # Vc_work <- V_cmax_0
#             # J_work <- J_max_0
#             # break
#           } else {
#             print(paste(
#               "Failed. Vcmax:", V_cmax_0,
#               "Jmax:", J_max_0,
#               "Rd:", R_d_0
#             ))
#           }
#         }
#       }
#     }
#   }
# }
#
# ggplot(data = LRC_test[LRC_test$ID == 1, ], aes(x = Ci, y = Photosynthesis)) +
#   geom_point() +
#   geom_point(aes(x = Cc, y = FvCB), colour = "blue")
#
# CO2light
#
#
#
# nls.all <- nls(y ~ FvCB2(x1, x2, PAR, Vcmax, Jmax, alpha, theta, Rd),
#   data = LRC1,
#   start = list(
#     Vcmax = 76.8, Jmax = 100, alpha = 0.5, theta = -0.8, Rd = 1.4
#     # Vm.d = 20, Jm.d = 300, al.d = 0, th.d = -1, Rd.d = -0.3
#   ),
#   control = list(maxiter = 250, minFactor = 0.000001, printEval = T, tol = 1e-04)
# )

# === NLME with full data ======================================================

full_model_1 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ TREATMENT, alpha ~ TREATMENT, theta ~ TREATMENT),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# === Model 1 investigation ====================================================

full_model_1_sum <- summary(full_model_1)


nlme_full_1_fixed <- full_model_1_sum$tTable
nlme_full_1_random <- summary(full_model_1)$coefficients$random

# Compare estimates for Treatments
fixed_comparison <- as.data.frame(nlme_full_1_fixed)
fixed_comparison$Feature <- rownames(fixed_comparison)
fixed_comparison$Feature <- c(
  rep("Vcmax", 2),
  rep("Jmax", 2),
  rep("Rd", 2),
  rep("alpha", 2),
  rep("theta", 2)
)

fixed_comparison$Treatment <- rep(c("A1", "A2"), nrow(fixed_comparison) / 2)


comparison_nlme_plot <- fixed_comparison %>%
  mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nlme_plot, aes(x = Feature, colour = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
    position = "dodge"
  ) +
  labs(
    title = "NLME model estimate",
    y = "Value", x = "Parameters"
  )

# === Back to NLME =============================================================

# Drop ADAB Side as not enough info
CO2_grouped_no_ADAB <- CO2_grouped_full[CO2_grouped_full$SIDE != "ADAB", ]

full_model_2 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1, 0, -(Vcmax1 - Vcmax2), -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, 0, -(Jmax1 - Jmax2), -(Jmax1 - Jmax2)),
  Rd = c(Rd1, 0),
  alpha = c(0.21989754, 0, 0, 0),
  theta = c(-1.17110061, 0, 0, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_3 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
),
random = pdDiag(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, 0, -(Vcmax1 - Vcmax2), -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, 0, -(Jmax1 - Jmax2), -(Jmax1 - Jmax2)),
  Rd = c(Rd1, 0),
  alpha = c(0.21989754, 0, 0, 0),
  theta = c(-1.17110061, 0, 0, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_4 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(
  Vcmax ~ TREATMENT / SIDE,
  Jmax ~ TREATMENT / SIDE,
  Rd ~ SIDE,
  alpha ~ TREATMENT / SIDE,
  theta ~ TREATMENT / SIDE
),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, 0, -(Vcmax1 - Vcmax2), -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, 0, -(Jmax1 - Jmax2), -(Jmax1 - Jmax2)),
  Rd = c(Rd1, 0),
  alpha = c(0.21989754, 0, 0, 0),
  theta = c(-1.17110061, 0, 0, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# Stop using SIDE for now

full_model_5 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ TREATMENT, alpha ~ TREATMENT, theta ~ TREATMENT),
random = pdDiag(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = c(Rd1, 0),
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# === Model 5 investigation ====================================================

full_model_5_sum <- summary(full_model_5)


nlme_full_5_fixed <- full_model_5_sum$tTable
nlme_full_5_random <- summary(full_model_5)$coefficients$random

# Compare estimates for Treatments
fixed_comparison <- as.data.frame(nlme_full_5_fixed)
fixed_comparison$Feature <- rownames(fixed_comparison)
fixed_comparison$Feature <- c(
  rep("Vcmax", 2),
  rep("Jmax", 2),
  rep("Rd", 2),
  rep("alpha", 2),
  rep("theta", 2)
)

fixed_comparison$Treatment <- rep(c("A1", "A2"), nrow(fixed_comparison) / 2)


comparison_nlme_plot <- fixed_comparison %>%
  mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nlme_plot, aes(x = Feature, colour = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
    position = "dodge"
  ) +
  labs(
    title = "NLME model estimate",
    y = "Value", x = "Parameters"
  )

anova(full_model_5, full_model_1)
summary(full_model_1)

# === Back to analysis =========================================================

full_model_6 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ TREATMENT, theta ~ TREATMENT),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1, Jmax1 - Jmax2),
  Rd = Rd1,
  alpha = c(alpha1, 0),
  theta = c(theta1, 0)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_7 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

anova(full_model_1, full_model_7)

# === Model 7 investigation ====================================================

full_model_7_sum <- summary(full_model_7)


nlme_full_7_fixed <- full_model_7_sum$tTable
nlme_full_7_random <- summary(full_model_7)$coefficients$random

# Compare estimates for Treatments
fixed_comparison <- as.data.frame(nlme_full_7_fixed)
fixed_comparison$Feature <- rownames(fixed_comparison)
fixed_comparison$Feature <- c(
  rep("Vcmax", 2),
  rep("Jmax", 1),
  rep("Rd", 1),
  rep("alpha", 1),
  rep("theta", 1)
)

fixed_comparison$Treatment <- c(c("A1", "A2"), rep("All", nrow(fixed_comparison) - 2))


comparison_nlme_plot <- fixed_comparison %>%
  mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nlme_plot, aes(x = Feature, colour = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
    position = "dodge"
  ) +
  labs(
    title = "NLME model estimate",
    y = "Value", x = "Parameters"
  )

anova(full_model_7, full_model_1)
summary(full_model_7)

# Rd is mad

# === Back to analysis =========================================================

# This model converges
full_model_8 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

summary(full_model_8)
anova(full_model_8, full_model_7)

full_model_7_no_ADAB <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_7_no_ADAB <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
  Jmax = c(Jmax1),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)

full_model_7_no_ADAB <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = Vcmax + Jmax + Rd + alpha + theta ~ 1,
start = c(
  Vcmax = c(Vcmax1),
  Jmax = c(Jmax1),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 150
)
)


full_model_9 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_no_ADAB,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ SIDE, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, 0),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)


full_model_10 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ SIDE, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, 0, 0),
  Rd = Rd1,
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_11 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ 1, Rd ~ SIDE, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1),
  Rd = c(Rd1, 0, 0),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

summary(full_model_11)
anova(full_model_11, full_model_7)

full_model_12 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ SIDE, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_13 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ SIDE, Jmax ~ SIDE, Rd ~ SIDE, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# === Model 7 investigation ====================================================

full_model_13_sum <- summary(full_model_13)

nlme_full_13_fixed <- full_model_13_sum$tTable
nlme_full_13_random <- summary(full_model_13)$coefficients$random

# Compare estimates for Treatments
fixed_comparison <- as.data.frame(nlme_full_13_fixed)
fixed_comparison$Feature <- rownames(fixed_comparison)
fixed_comparison$Feature <- c(
  rep("Vcmax", 2),
  rep("Jmax", 1),
  rep("Rd", 1),
  rep("alpha", 1),
  rep("theta", 1)
)

fixed_comparison$Treatment <- c(c("A1", "A2"), rep("All", nrow(fixed_comparison) - 2))


comparison_nlme_plot <- fixed_comparison %>%
  mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nlme_plot, aes(x = Feature, colour = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
    position = "dodge"
  ) +
  labs(
    title = "NLME model estimate",
    y = "Value", x = "Parameters"
  )

anova(full_model_13, full_model_1)
summary(full_model_13)
# Rd is mad

# === Back to analysis =========================================================

full_model_14 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ SIDE, Jmax ~ SIDE, Rd ~ SIDE, alpha ~ O, theta ~ O)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_15 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ SIDE, Jmax ~ SIDE, Rd ~ SIDE, alpha ~ SIDE, theta ~ SIDE)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_16 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ SIDE, Jmax ~ SIDE, Rd ~ SIDE, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_17 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ SIDE, Jmax ~ SIDE, Rd ~ SIDE, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_18 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdSymm(list(Vcmax ~ 1, Jmax ~ SIDE, Rd ~ SIDE, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

full_model_19 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  alpha,
  theta,
  Rd
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, Rd ~ 1, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ SIDE %in% TREATMENT, Jmax ~ SIDE %in% TREATMENT, Rd ~ SIDE, alpha ~ 1, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  Rd = c(Rd1),
  alpha = c(alpha1),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)


full_model_20 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  1,
  theta,
  1.5
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, theta ~ 1),
random = pdDiag(list(Vcmax ~ SIDE, Jmax ~ SIDE, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)

# This might work - R keeps aborting for me
full_model_21 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  1,
  theta,
  1.5
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, theta ~ 1),
random = pdSymm(list(Vcmax ~ SIDE, Jmax ~ SIDE, theta ~ 1)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  theta = c(theta1)
),
control = list(
  maxIter = 250, msVerbose = T, tolerance = 1e-1, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)


full_model_22 <- nlme(Photosynthesis ~ FvCB(
  Ci,
  PAR,
  Gstar,
  Kc,
  Ko,
  O,
  Vcmax,
  Jmax,
  1.2,
  2.0,
  1.5
),
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT),
random = pdSymm(list(Vcmax ~ SIDE, Jmax ~ SIDE)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2))
),
control = list(
  maxIter = 250, msVerbose = F, tolerance = 1e-3, msMaxIter = 250,
  pnlsTol = 1e-10, pnlsMaxIter = 50
)
)
