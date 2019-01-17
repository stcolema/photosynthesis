
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

LRC1 <- read.table("LRC1.txt", sep = "\t", header = TRUE, na.strings = "")
Dark <- read.table("Dark-F.txt", sep = "\t", header = TRUE, na.strings = "")
LRC2_F <- read.table("LRC2-F.txt", sep = "\t", header = TRUE, na.strings = "")
# ACI1 <- na.omit(read.table("A-Ci1.txt", sep = "\t", header = TRUE, na.strings = ""))
ACI2 <- read.table("A-Ci2.txt", sep = "\t", header = TRUE, na.strings = "", fileEncoding="UCS-2LE")
# LRC2 <- read.table("LRC1.txt", sep = "\t", header = TRUE, na.strings = "")


# Delete empty lines and any DIV/o errors - also empty column in LRC2
# ACI1 <- na.omit(read.table("A-Ci1_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI1 <- read.table("A-Ci1_clean.csv", sep = ",", header = TRUE, na.strings = "")
LRC2 <- read.table("LRC2_clean.csv", sep = ",", header = TRUE, na.strings = "")
# ACI2 <- read.table("A-Ci2_clean.csv", sep = ",", header = TRUE, na.strings = "")
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

ACI1_new <- update_aci_data(ACI1, Tleaf, Photo, PARabs, Ci,
                            O_pres = O_21,
                            reduce = T,
                            Kelvin = F,
                            other_vars_to_keep = c(SIDE, TREAT, ID, PhiPS2)
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

select_AD <- combined_data$SIDE == "AD"
select_AB <- combined_data$SIDE == "AB"
select_ADAB <- combined_data$SIDE == "ADAB"

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

unique(CO2_grouped$ID)
nrow(CO2_grouped[CO2_grouped$ID == 1 & CO2_grouped$SIDE == "AD", ])

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

LRC2_F_new <- update_aci_data(LRC2_F, Tleaf, Photo, PARabs, Ci,
                              O_pres = O_2,
                              reduce = T,
                              Kelvin = F,
                              other_vars_to_keep = c(SIDE, TREAT, ID, PhiPS2)
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

LRC2_F_new <- LRC2_F_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("light", labels = c("light")))

LRC2_F_new$SIDE <- add_level(LRC2_F_new$SIDE, "ADAB")
LRC2_F_new$Response <- add_level(LRC2_F_new$Response, "CO2")


combined_data_2 <- dplyr::bind_rows(LRC2_new, ACI2_new)

total_data <- dplyr::bind_rows(combined_data, combined_data_2)
# total_data$O <- factor(total_data$O)
summary(total_data)

CO2_grouped_full <- groupedData(Photosynthesis ~ 1 | ID,
                                outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
                                inner = ~Ci + PAR + SIDE,
                                data = total_data
)


# # Indices for subgroups
# select_A1_full <- total_data$TREATMENT == "A1"
# select_A2_full <- total_data$TREATMENT == "A2"
# 
# select_A1_AD_full <- as.logical((total_data$TREATMENT == "A1") * (total_data$SIDE == "AD"))
# select_A1_AB_full <- as.logical((total_data$TREATMENT == "A1") * (total_data$SIDE == "AB"))
# select_A1_ADAB_full <- as.logical((total_data$TREATMENT == "A1") * (total_data$SIDE == "ADAB"))
# 
# select_A2_AD_full <- as.logical((total_data$TREATMENT == "A2") * (total_data$SIDE == "AD"))
# select_A2_AB_full <- as.logical((total_data$TREATMENT == "A2") * (total_data$SIDE == "AB"))
# select_A2_ADAB_full <- as.logical((total_data$TREATMENT == "A2") * (total_data$SIDE == "ADAB"))
# 
# select_CO2_full <- total_data$Response == "CO2"
# select_light_full <- total_data$Response == "light"
# 
# select_O_2 <- total_data$O == 20
# select_O_21 <- total_data$O == 210

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

# === NLME example =============================================================

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


# === Yin step out model =======================================================
create_LRC2_combined <- function(LRC2, LRC2_F, A_param, phi_param, i_inc_param, id_param) {
  # enquose the variable names
  # A <- rlang::enquo(A_param)
  # phi <- rlang::enquo(phi_param)
  # I <- rlang::enquo(i_inc_param)
  # params_of_interest_LRC2 <- rlang::enquos(A_param, i_inc_param)
  params_of_interest_LRC2 <- c(A_param, i_inc_param, id_param)
  params_of_interest_LRC2_F <- c(phi_param, id_param)
  
  # select the relevant subsets of the datasets
  LRC2_rel <- LRC2 %>%
    dplyr::select(!!!params_of_interest_LRC2)
  
  # print("woof")
  LRC2_F_rel <- LRC2_F %>%
    dplyr::select(!!!params_of_interest_LRC2_F)
  
  # print("whos a good boy")
  # combine the dataframes and add the variable that the regression is on
  comb_data <- dplyr::bind_cols(LRC2_rel, LRC2_F_rel) %>%
    dplyr::mutate(New_var = 0.25 * !!phi_param * !!i_inc_param)
  
  comb_data
}

calc_s_rd <- function(LRC2, LRC2_F,
                      A_param = A,
                      phi_param = PHI,
                      i_inc_param = I,
                      id_var = ID,
                      ...) {
  # Enclose the ID variable
  A <- rlang::enquo(A_param)
  phi <- rlang::enquo(phi_param)
  I <- rlang::enquo(i_inc_param)
  ID <- rlang::enquo(id_var)
  
  # print("giragge")
  # Create the reduced, combined dataframe
  comb_data <- create_LRC2_combined(LRC2, LRC2_F, A, phi, I, ID)
  # print("hihihi")
  # return(comb_data)
  # Carry out regression using a mixed effects model
  lme_s_rd <- nlme::lme(formula(substitute(A_param ~ New_var)),
                        data = comb_data,
                        random = formula(substitute(~New_var | id_var)),
                        ...
  )
  
  # Return the summary of the above model
  sum_lm_s_rd <- summary(lme_s_rd)
  
  fix <- sum_lm_s_rd$coefficients$fixed
  rand <- sum_lm_s_rd$coefficients$random[[1]]
  
  Rd_ind <- rand[, 1] + fix[1]
  S_ind <- rand[, 2] + fix[2]
  
  ind_coef <- data.frame(ID = 1:length(Rd_ind), Rd = Rd_ind, S = S_ind)
  
  out <- list(
    data = comb_data,
    model = lme_s_rd,
    summary = sum_lm_s_rd,
    coefficients = ind_coef
  )
}



calc_J <- function(data, lump_var = S, i_inc_var = I, phi_psII_var = phi) {
  S <- rlang::enquo(lump_var)
  I <- rlang::enquo(i_inc_var)
  phi <- rlang::enquo(phi_psII_var)
  new_data <- data %>%
    dplyr::mutate(J = !!S * !!I * !!phi)
  new_data
}


# Delete empty lines and any DIV/o errors - also empty column in LRC2
LRC2 <- read.table("LRC2_clean.csv", sep = ",", header = TRUE, na.strings = "")
LRC2_F <- na.omit(read.table("LRC2-F.txt", sep = "\t", header = TRUE, na.strings = ""))

LRC2_new <- update_aci_data(LRC2, Tleaf, Photo, PARabs, Ci,
                            O_pres = O_2,
                            reduce = T,
                            Kelvin = F,
                            other_vars_to_keep = c(SIDE, TREAT, ID)
)

LRC2_F_new <- update_aci_data(LRC2_F, Tleaf, Photo, PARabs, Ci,
                              O_pres = O_2,
                              reduce = T,
                              Kelvin = F,
                              other_vars_to_keep = c(SIDE, TREAT, ID, PhiPS2)
)

LRC2_new <- LRC2_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("light", labels = c("light")))

LRC2_new$SIDE <- add_level(LRC2_new$SIDE, "ADAB")
LRC2_new$Response <- add_level(LRC2_new$Response, "CO2")

LRC2_F_new <- LRC2_F_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("light", labels = c("light")))

LRC2_F_new$SIDE <- add_level(LRC2_F_new$SIDE, "ADAB")
LRC2_F_new$Response <- add_level(LRC2_F_new$Response, "CO2")


out <- calc_s_rd(LRC2_new, LRC2_F_new,
                 A_param = Photosynthesis,
                 phi_param = PhiPS2,
                 i_inc_param = PAR,
                 id_var = ID,
                 control = list(
                   maxIter = 250,
                   msVerbose = F,
                   tolerance = 1e-2,
                   msMaxIter = 250,
                   pnlsTol = 1e-10,
                   pnlsMaxIter = 250,
                   niterEM = 1000
                 )
)


xyplot(Photosynthesis ~ New_var | as.factor(ID),
       group = as.factor(ID),
       data = out$data
)


names(LRC2_F_new)

LRC2_F_newer <- LRC2_F_new %>%
  left_join(out$coefficients, by = c("ID"))

ggplot2::ggplot(data = out$data, aes(x = New_var, y = Photosynthesis)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Regression for Rd, S",
    x = bquote(~frac(1, 4) ~ "I" ~ Phi[PSII]),
    y = "A"
  )


J_data <- calc_J(LRC2_F_newer, lump_var = S, i_inc_var = PAR, phi_psII_var = PhiPS2)

head(J_data)

total_data <- total_data %>%
  left_join(out$coefficients, by = c("ID"))

CO2_grouped_full <- groupedData(Photosynthesis ~ 1 | ID,
                                outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
                                inner = ~Ci + PAR + SIDE,
                                data = total_data
)

# === NLME again ===============================================================

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
data = CO2_grouped_full,
fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, alpha ~ 1, theta ~ 1),
random = pdDiag(list(Vcmax ~ ID, Jmax ~ ID)),
start = c(
  Vcmax = c(Vcmax1, -(Vcmax1 - Vcmax2)),
  Jmax = c(Jmax1, -(Jmax1 - Jmax2)),
  alpha = 0.8,
  theta = 1.5
),
subset = select_AD,
control = nlmeControl(
  opt = "nlminb",
  upper = c(
    Vcmax = c(200, 200),
    Jmax = c(300, 300),
    alpha = 2,
    theta = 2,
    Rd = 3
  ),
  lower = c(
    Vcmax = c(0, 0),
    Jmax = c(0, 0),
    alpha = 0,
    theta = 0,
    Rd = -2
  ),
  maxIter = 250,
  msVerbose = F,
  tolerance = 1e-2,
  msMaxIter = 250,
  pnlsTol = 1e-10,
  pnlsMaxIter = 250
)
)

# === Sco ======================================================================

calc_Sco <- function(ACI1, ACI2, Ci_threshold = 100, id_var = ID, A_var = A, ...) {
  # Enclose the ID variable
  ID <- rlang::enquo(id_var)
  A <- rlang::enquo(A_var)
  
  # Filter down to the linear part of the ACI curve in the high oxygen data
  ACI1_lin <- ACI1 %>%
    dplyr::filter(Ci < Ci_threshold)
  
  # Check what individuals are present
  ind_present_h <- ACI1_lin %>%
    dplyr::select(!!ID) %>%
    unique()
  
  ind_present_l <- ACI2 %>%
    dplyr::select(!!ID) %>%
    unique()
  
  # Return an error if any individuals are not present in both datasets
  if (any(unlist(ind_present_h) != unlist(ind_present_l))) {
    stop("Individuals present not matching.")
  }
  
  # Model for the bh and bl variables
  lme_h <- nlme::lme(formula(substitute(A_var ~ Ci)),
                     data = ACI1_lin,
                     random = formula(substitute(~Ci | id_var)),
                     ...
  )
  
  lme_l <- nlme::lme(formula(substitute(A_var ~ Ci)),
                     data = ACI2,
                     random = formula(substitute(~Ci | id_var)),
                     ...
  )
  
  # Get the slopes for each individual under high and low oxygen
  bh <- lme_h$coefficients$fixed[[2]]
  bl <- lme_l$coefficients$fixed[[2]]
  
  bh_vec <- lme_h$coefficients$random$ID[, 2] + bh
  bl_vec <- lme_l$coefficients$random$ID[, 2] + bl
  
  # Collect these in a dataframe with the corresponding ID
  id_str <- quo_name(ID)
  coef_df <- data.frame(
    ID = unlist(ind_present_h),
    bh = bh_vec,
    bl = bl_vec
  )
  names(coef_df)[1] <- id_str
  
  # Want the same number of observations in both high and low
  n_h <- ACI1_lin %>%
    dplyr::group_by(!!ID) %>%
    dplyr::count()
  
  n_l <- ACI2 %>%
    dplyr::group_by(!!ID) %>%
    dplyr::count()
  
  # Create reduced dataframes with renamed variables
  common_ci_data_h <- ACI1_lin %>%
    dplyr::mutate(Oh = O, Ah = !!A, Cih = Ci) %>%
    dplyr::select(-O, -!!A, -Ci) %>%
    dplyr::select(Oh, Ah, Rd, Cih, !!ID)
  
  common_ci_data_l <- ACI2 %>%
    dplyr::mutate(Ol = O, Al = !!A, Cil = Ci) %>%
    dplyr::select(-O, -!!A, -Ci) %>%
    dplyr::select(Ol, Al, Cil, !!ID)
  
  # Create empty dataframe which we will store the data in
  comb_data <- as.data.frame(matrix(nrow = 0, ncol = 8))
  names(comb_data) <- unique(c(names(common_ci_data_h), names(common_ci_data_l)))
  
  # Make sure we have the same number of observations for each ID
  for (i in unlist(ind_present_l)) {
    # Find the number of observations for the individual in both datasets
    n_pres_h <- n_h %>% dplyr::filter(!!ID == i) %>% .$n
    n_pres_l <- n_l %>% dplyr::filter(!!ID == i) %>% .$n
    
    # Use the lesser
    n_pres <- min(n_pres_h, n_pres_l)
    
    # Create subsets of the individual's data with n observations
    ind_subset_h <- common_ci_data_h %>%
      dplyr::filter(!!ID == i) %>%
      dplyr::top_n(n_pres, Cih)
    
    ind_subset_l <- common_ci_data_l %>%
      dplyr::filter(!!ID == i) %>%
      dplyr::top_n(n_pres, Cil) %>%
      dplyr::select(-!!ID)
    
    # Bind these together by column (hence the common number of observations)
    new_subset <- dplyr::bind_cols(ind_subset_h, ind_subset_l)
    
    # Store this in the full dataset
    comb_data <- dplyr::bind_rows(comb_data, new_subset)
  }
  
  # Add the bh and bl variables and calculate Sco (as per Yin et al 2009)
  full_data <- comb_data %>%
    dplyr::left_join(coef_df, by = quo_name(ID)) %>%
    dplyr::mutate(S_co = bh * bl * (Oh - Ol) /
                    (2 * bl * (Ah + Rd) - 2 * bh * (Al - Rd) + bh * bl * (Cih - Cil)))
}


# ACI1 <- read.table("A-Ci1_2.csv", sep = ",", header = TRUE, na.strings = "")
# data_h <- ACI1 %>% dplyr::filter(Ci < 200)

ggplot2::ggplot(data = data_h, aes(y = Photo, x = Ci)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Ci per individual",
    x = bquote(~C[i]),
    y = "A"
  )
ID_quo <- quo(ID)
ACI1_RD <- ACI1_new %>%
  left_join(out$coefficients, by = quo_name(ID_quo)) %>%
  filter(ID != 1)

names(ACI1_RD)
unique(ACI1_RD$ID)

unique(ACI2_new$ID)
ACI2_2 <-ACI2_new %>% dplyr::filter(ID != 1) %>% na.omit()
unique(ACI2_2$ID)

new_data_sco <- calc_Sco(ACI1_RD, ACI2_2,
                         Ci_threshold = 200,
                         id_var = ID,
                         A_var = Photosynthesis,
                         control = list(
                           maxIter = 250,
                           msVerbose = F,
                           tolerance = 1e-2,
                           msMaxIter = 250,
                           pnlsTol = 1e-10,
                           pnlsMaxIter = 250,
                           niterEM = 1000
                         )
)

ggplot2::ggplot(data = new_data_sco, aes(x = S_co)) +
  geom_density() +
  facet_wrap(~ID)

sco_data <- new_data_sco %>% 
  group_by(ID) %>%
  filter(!(abs(S_co - median(S_co)) > 2*sd(S_co)))

ggplot2::ggplot(data = sco_data, aes(x = S_co)) +
  geom_density() +
  facet_wrap(~ID)

ggplot2::ggplot(data = sco_data, aes(y = S_co, group = ID)) +
  geom_boxplot()

sco_data_2 <- sco_data %>% 
  group_by(ID) %>%
  filter(!(abs(S_co - median(S_co)) > 2*sd(S_co)))

ggplot2::ggplot(data = sco_data_2, aes(x = S_co)) +
  geom_density() +
  facet_wrap(~ID)

ggplot2::ggplot(data = sco_data_2, aes(y = S_co, group = ID)) +
  geom_boxplot()


new_data_sco$S_co

ACI1_RD <- dplyr::left_join(ACI1_new, out$coefficients, by = "ID")
ACI1_with_j <- calc_J(ACI1_RD, lump_var = S, i_inc_var = PAR, phi_psII_var = PhiPS2)


ACI1_final <- na.omit(ACI1_with_j)
ACI1_final %>% 
  dplyr::group_by(ID) %>% 
  dplyr::count()
names(ACI1_final)



ACI1_final <- ACI1_final %>% filter(ID != 1)
# Indices for subgroups
select_A1 <- ACI1_final$TREATMENT == "A1"
select_A2 <- ACI1_final$TREATMENT == "A2"

select_AD <- ACI1_final$SIDE == "AD"
select_AB <- ACI1_final$SIDE == "AB"
select_ADAB <- ACI1_final$SIDE == "ADAB"

select_A1_AD <- as.logical((ACI1_final$TREATMENT == "A1") * (ACI1_final$SIDE == "AD"))
select_A1_AB <- as.logical((ACI1_final$TREATMENT == "A1") * (ACI1_final$SIDE == "AB"))
select_A1_ADAB <- as.logical((ACI1_final$TREATMENT == "A1") * (ACI1_final$SIDE == "ADAB"))

select_A2_AD <- as.logical((ACI1_final$TREATMENT == "A2") * (ACI1_final$SIDE == "AD"))
select_A2_AB <- as.logical((ACI1_final$TREATMENT == "A2") * (ACI1_final$SIDE == "AB"))
select_A2_ADAB <- as.logical((ACI1_final$TREATMENT == "A2") * (ACI1_final$SIDE == "ADAB"))



grouped_ACI <- groupedData(Photosynthesis ~ 1 | ID,
                           outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
                           inner = ~SIDE,
                           data = ACI1_final
)

non_rectangular_hyperbola <- function(a, b, c){
  (-b - (b ** 2 - 4 * a * c) ** 0.5) / 2 * a
}

yin_a <- function(x_2, C_i, G_star, delta){
  x_2 + G_star + delta * (C_i + x_2)
}

yin_b <- function(x_1, x_2, C_i, G_star, R_d, g_mo, delta){
  -((x_2 + G_star) * (x_1 - R_d) + (C_i + x_2)
    * (g_mo * (x_2 + G_star) + delta * (x_1 - R_d))
    + delta * (x_1 * (C_i - G_star) - R_d * (C_i +x_2)))
}

yin_c <- function(x_1, x_2, C_i, G_star, R_d, g_mo, delta){
  c <- ((g_mo * (x_2 + G_star) + delta * (x_1 - R_d))
       * (x_1 * (C_i - G_star) - R_d * (C_i + x_2)))
  c
}

yin_FvCB_curves <- function(C_i, G_star, R_d, g_mo, delta, 
                            A_c = TRUE,
                            V_cmax = NULL, 
                            K_c = NULL, 
                            K_o = NULL,
                            O = NULL,
                            J = NULL){
  if(A_c) {
    if(any(sapply(list(V_cmax, K_c, K_o, O), is.null))) {
      num_null <- sum(sapply(list(V_cmax, K_c, K_o, O), is.null))
      stop(paste("Please give values for inputs.", 
                 num_null, 
                 " of V_cmax, K_c, K_o, and O are undeclared"
                 )
      )
    }
    x_1 <- V_cmax
    x_2 <- K_c * (1 + O / K_o)
  } else{
    if(is.null(J)){
      stop("Please declare J.")
    }
    x_1 <- J / 4
    x_2 <- 2 * G_star
  }
  a <- yin_a(x_2, C_i, G_star, delta)
  b <- yin_b(x_1, x_2, C_i, G_star, R_d, g_mo, delta)
  c <- yin_c(x_1, x_2, C_i, G_star, R_d, g_mo, delta)
  
  non_rectangular_hyperbola(a, b, c)
}

FvCB <- function(C_i, G_star, R_d, g_mo, delta, K_c, K_o, O, V_cmax, J) {
  p1 <- yin_FvCB_curves(C_i, G_star, R_d, g_mo, delta, 
                                    A_c = TRUE,
                                    V_cmax = V_cmax, 
                                    K_c = K_c, 
                                    K_o = K_o,
                                    O = O)
  p2 <- yin_FvCB_curves(C_i, G_star, R_d, g_mo, delta, 
                        A_c = FALSE,
                        J = J)
  ifelse(p1 < p2, p1, p2)
}

nls_aci1 <- nls(Photosynthesis ~
                  FvCB(Ci, Gstar, Rd, 0, 1.8, Kc, Ko, O, Vcmax, J),
                data = ACI1_final,
                subset = select_A1,
                start = list(
                  Vcmax = -1000
                ),
                control = list(maxiter = 250, minFactor = 1e-10, printEval = T, tol = 1e-6)
)



model_1 <- nlme(Photosynthesis ~ FvCB(Ci, Gstar, Rd, 0, 1.4, Kc, Ko, O, Vcmax, J),
data = grouped_ACI,
fixed = Vcmax ~ TREATMENT,
random = pdDiag(Vcmax ~ ID), # add %in% TREATMENT?
start = c(
  Vcmax = c(55, 0)
),
subset = select_AD,
control = nlmeControl(
  opt = "nlminb",
  upper = c(
    Vcmax = c(300, 100)
  ),
  lower = c(
    Vcmax = c(0, 0)
  ),
  maxIter = 250,
  msVerbose = T,
  tolerance = 1e-2,
  msMaxIter = 250,
  pnlsTol = 1e-10,
  pnlsMaxIter = 250,
  niterEM = 1000
)
)

ggplot(data = ACI1_final[select_AD,], aes(x = Ci, y = Photosynthesis, colour = TREATMENT, shape = SIDE)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ID)

