
# === Libraries ================================================================

library(nlme)
library(lattice)
library(lme4)

library(tidyverse) # install.packages("tidyverse", dep = T)
library(rlang)
library(magrittr)
library(attempt)

# === Functions ================================================================

# Add a new level to a factor
add_level <- function(x, new_level) {
  if (is.factor(x)) return(factor(x, levels = c(levels(x), new_level)))
  return(x)
}

# Function to strip datasets of unwanted variables and rename variables for ease 
# of use.
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
      # Gstar = exp(19.02 - 38.83 / (R * (!!T_leaf + 273.15 * (1 - Kelvin)))),
      Photosynthesis = !!response,
      PAR = !!PAR
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
  A <- J * (numerator / denom) - R_d
  A
}



FvCB <- function(C_c, G_star, K_c, K_o, O, V_cmax, J, R_d) {
  p1 <- Rubisco_limited(V_cmax, C_c, G_star, K_c, O, K_o, R_d) + R_d
  # J <- non_rectangular_hyperbola(PAR, J_max, alpha, theta)
  # p2 <- J * (C_c - G_star) / (4.5 * C_c  + 10.5 * G_star)
  p2 <- RuBP_limited(J, C_c, G_star, R_d)
  ifelse(p1 < p2, p1, p2) - R_d
}


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


FvCB_jmax <- function(C_c, G_star, K_c, K_o, O, V_cmax, J_max, R_d, PAR, alpha, theta) {
  p1 <- Rubisco_limited(V_cmax, C_c, G_star, K_c, O, K_o, R_d) + R_d
  J <- non_rectangular_hyperbola(PAR, J_max, alpha, theta)
  # p2 <- J * (C_c - G_star) / (4.5 * C_c  + 10.5 * G_star)
  p2 <- RuBP_limited(J, C_c, G_star, R_d)
  ifelse(p1 < p2, p1, p2) - R_d
}



# Calculate daytime respiration (Rd) and lump parameter (S)
calc_s_rd <- function(LRC2, LRC2_F,
                      A_param = A,
                      phi_param = PHI,
                      i_inc_param = I,
                      id_var = ID,
                      other_var_to_keep = NULL,
                      i_lower_bound = 0,
                      i_upper_bound = Inf,
                      ...) {
  # Enclose the variables used to avoid use of strings and enable use of dplyr
  A <- rlang::enquo(A_param)
  phi <- rlang::enquo(phi_param)
  I <- rlang::enquo(i_inc_param)
  ID <- rlang::enquo(id_var)
  other_var_kept <- rlang::enquos(other_var_to_keep)

  # select the relevant subsets of the datasets
  LRC2_rel <- LRC2 %>%
    dplyr::select(!!!c(A, I, ID, other_var_kept)) %>% 
    dplyr::filter(!! I > i_lower_bound, !! I < i_upper_bound)

  LRC2_F_rel <- LRC2_F %>%
    dplyr::select(!!!c(phi, I, ID, other_var_kept)) %>% 
    dplyr::filter(!! I > i_lower_bound, !! I < i_upper_bound) %>% 
    dplyr::select(- !! I)

  # combine the dataframes and add the variable that the regression is on
  comb_data <- dplyr::bind_cols(LRC2_rel, LRC2_F_rel)

  # Create the formula for LME
  formula_to_use <- formula(substitute(~New_var | id_var))

  # Create the variable to perform regression upon
  # we add an intercept to access it for boundary purposes
  # (can access the variable Intercept using lmeControl)
  comb_data <- comb_data %>%
    dplyr::mutate(New_var = 0.25 * !!phi * !!I) 

  # Carry out regression using a mixed effects model
  # Limit the intercept to being positive
  lme_s_rd <- nlme::lme(formula(substitute(A_param ~ New_var)),
    data = comb_data,
    random = formula_to_use,
    ...
  )

  # Return the summary of the above model
  sum_lm_s_rd <- summary(lme_s_rd)

  # Get the point estimates of the parameters for each individual
  fix <- sum_lm_s_rd$coefficients$fixed
  rand <- sum_lm_s_rd$coefficients$random[[1]]

  # Add the fixed component to the individuals random component to acquire the
  # parameters of interest
  Rd_ind <- rand[, 1] + fix[1]
  S_ind <- rand[, 2] + fix[2]

  # Put these parameters with their associated IDs in a dataframe for output
  ind_coef <- data.frame(
    ID = as.numeric(row.names(rand)),
    Rd = Rd_ind,
    S = S_ind
  )

  # Include them in the dataset they were calculated in, assigining the values
  # to the associated individuals
  comb_data <- comb_data %>%
    dplyr::left_join(ind_coef, by = c("ID"))

  # Return the data used, the model, its summary and the coefficient dataframe
  out <- list(
    data = comb_data,
    model = lme_s_rd,
    summary = sum_lm_s_rd,
    coefficients = ind_coef
  )
}


# Calculate the eltron transport rate (J) in accordance with Yin & Struik
calc_J <- function(data, lump_var = S, i_inc_var = I, phi_psII_var = phi) {
  # Enclose variables
  S <- rlang::enquo(lump_var)
  I <- rlang::enquo(i_inc_var)
  phi <- rlang::enquo(phi_psII_var)
  
  # Add new variable to data
  new_data <- data %>%
    dplyr::mutate(J = !!S * !!I * !!phi)
  new_data
}


# Calculate the relative CO2 / O2 specifcity factor for Rubisco
calc_Sco <- function(ACI1, ACI2,
                     Ci_threshold = 100,
                     id_var = ID,
                     A_var = A,
                     other_var = NULL,
                     ...) {
  # Enclose the ID variable
  ID <- rlang::enquo(id_var)
  A <- rlang::enquo(A_var)
  other_var <- rlang::enquo(other_var)

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

  if (!quo_is_null(other_var)) {
    key <- ACI1 %>%
      dplyr::group_by(!!ID) %>%
      dplyr::distinct(!!other_var)

    coef_df <- coef_df %>%
      left_join(key, by = quo_name(ID))
  }

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
    dplyr::select(Oh, Ah, Rd, Cih, !!ID, !!other_var)

  common_ci_data_l <- ACI2 %>%
    dplyr::mutate(Ol = O, Al = !!A, Cil = Ci) %>%
    dplyr::select(-O, -!!A, -Ci) %>%
    dplyr::select(Ol, Al, Cil, !!ID, !!other_var)

  # Create empty dataframe which we will store the data in
  n_col <- length(unique(c(names(common_ci_data_h), names(common_ci_data_l))))
  comb_data <- as.data.frame(matrix(nrow = 0, ncol = n_col))
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
    dplyr::left_join(coef_df, by = c(quo_name(ID))) %>% # , quo_name(other_var))) %>%
    dplyr::mutate(S_co = (Oh - Ol) /
      (2 * ((Al + Rd) / bl - (Ah + Rd) / bh - (Cil - Cih))))
}



# Calculate the Cc-based CO2 compensation point in the absence of Rd (G_star)
calc_G_star <- function(data, O_par = O, S_co_value = 3.022) {
  O <- rlang::enquo(O_par)
  data %>% dplyr::mutate(G_star = 0.5 * !!O / S_co_value)
}

# Caclulate msophyll diffusion conductance by Harley et al., 1992 Variable J 
# method
calc_gm <- function(data,
                    A_par = A,
                    Ci_par = C_i,
                    G_star_par = G_star,
                    J_par = J,
                    Rd_par = R_d,
                    num_const = 4,
                    denom_const = 8) {
  A <- enquo(A_par)
  C_i <- enquo(Ci_par)
  G_star <- enquo(G_star_par)
  J <- enquo(J_par)
  R_d <- enquo(Rd_par)
  
  # Calculate gm
  new_data <- data %>%
    dplyr::mutate(gm = !!A
    / (!!C_i -
        (!!G_star * (!!J + num_const * (!!A + !!R_d))
          / (!!J - denom_const * (!!A + !!R_d)))
      ))
  new_data
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
ACI2 <- read.table("A-Ci2.txt", sep = "\t", header = TRUE, na.strings = "", fileEncoding = "UCS-2LE")
# LRC2 <- read.table("LRC1.txt", sep = "\t", header = TRUE, na.strings = "")


# Delete empty lines and any DIV/o errors - also empty column in LRC2
# ACI1 <- na.omit(read.table("A-Ci1_clean.csv", sep = ",", header = TRUE, na.strings = ""))
ACI1 <- read.table("A-Ci1_clean.csv", sep = ",", header = TRUE, na.strings = "")
LRC2 <- read.table("LRC2_clean.csv", sep = ",", header = TRUE, na.strings = "")
# ACI2 <- read.table("A-Ci2_clean.csv", sep = ",", header = TRUE, na.strings = "")
# ACI1_F <- na.omit(read.table("A-Ci1-F_clean.csv", sep = ",", header = TRUE, na.strings = ""))
# ACI2_F <- na.omit(read.table("A-Ci2-F_clean.csv", sep = ",", header = TRUE, na.strings = ""))

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

# find the relevant parameters
LRC1 %>%
  select(contains("PAR")) %>%
  glimpse()


ACI1 %>%
  select(contains("photo")) %>%
  glimpse()

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

# Use the same variable names in both datasets
ACI1_new <- ACI1_new %>%
  dplyr::mutate(TREATMENT = TREAT) %>%
  dplyr::select(-TREAT) %>%
  dplyr::mutate(Response = factor("CO2", labels = c("CO2"))) # irrelevant

ACI1_new$SIDE <- add_level(ACI1_new$SIDE, "ADAB")
ACI1_new$Response <- add_level(ACI1_new$Response, "light") # irrelevant

LRC1_new <- LRC1_new %>%
  dplyr::mutate(Response = factor("light", labels = c("light")))
LRC1_new$Response <- add_level(LRC1_new$Response, "CO2")
# 
# combined_data <- dplyr::bind_rows(LRC1_new, ACI1_new)
# 
# # Indices for subgroups
# select_A1 <- combined_data$TREATMENT == "A1"
# select_A2 <- combined_data$TREATMENT == "A2"
# 
# select_AD <- combined_data$SIDE == "AD"
# select_AB <- combined_data$SIDE == "AB"
# select_ADAB <- combined_data$SIDE == "ADAB"
# 
# select_A1_AD <- as.logical((combined_data$TREATMENT == "A1") * (combined_data$SIDE == "AD"))
# select_A1_AB <- as.logical((combined_data$TREATMENT == "A1") * (combined_data$SIDE == "AB"))
# select_A1_ADAB <- as.logical((combined_data$TREATMENT == "A1") * (combined_data$SIDE == "ADAB"))
# 
# select_A2_AD <- as.logical((combined_data$TREATMENT == "A2") * (combined_data$SIDE == "AD"))
# select_A2_AB <- as.logical((combined_data$TREATMENT == "A2") * (combined_data$SIDE == "AB"))
# select_A2_ADAB <- as.logical((combined_data$TREATMENT == "A2") * (combined_data$SIDE == "ADAB"))
# 
# select_CO2 <- combined_data$Response == "CO2"
# select_light <- combined_data$Response == "light"
# 
# CO2_grouped <- groupedData(Photosynthesis ~ 1 | ID,
#   outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
#   inner = ~Ci + PAR + SIDE,
#   data = combined_data
# )
# 
# unique(CO2_grouped$ID)
# nrow(CO2_grouped[CO2_grouped$ID == 1 & CO2_grouped$SIDE == "AD", ])

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


# combined_data_2 <- dplyr::bind_rows(LRC2_new, ACI2_new)
# 
# total_data <- dplyr::bind_rows(combined_data, combined_data_2)
# # total_data$O <- factor(total_data$O)
# summary(total_data)
# 
# CO2_grouped_full <- groupedData(Photosynthesis ~ 1 | ID,
#   outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
#   inner = ~Ci + PAR + SIDE,
#   data = total_data
# )
# 
# 
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

# === Model 2 investigation ====================================================

# model_2_sum <- summary(model_2)
#
#
# nlme_2_fixed <- model_2_sum$tTable
# nlme_2_random <- summary(model_2)$coefficients$random
#
# # Compare estimates for Treatments
# fixed_comparison <- as.data.frame(nlme_2_fixed)
# fixed_comparison$Feature <- rownames(fixed_comparison)
# fixed_comparison$Feature <- c(
#   rep("Vcmax", 2),
#   rep("Jmax", 2),
#   rep("Rd", 2),
#   rep("alpha", 2),
#   rep("theta", 2)
# )
#
# fixed_comparison$Treatment <- rep(c("A1", "A2"), nrow(fixed_comparison) / 2)
#
#
# comparison_nlme_plot <- fixed_comparison %>%
#   mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# # select(Estimate, Std_error, Feature, Treatment)
#
# ggplot(comparison_nlme_plot, aes(x = Feature, colour = Treatment)) +
#   geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
#     position = "dodge"
#   ) +
#   labs(
#     title = "NLME model estimate",
#     y = "Value", x = "Parameters"
#   )
#
# # more precision
# model_2_precision <- nlme(Photosynthesis ~ FvCB(
#   Ci,
#   PAR,
#   Gstar,
#   Kc,
#   Ko,
#   O,
#   Vcmax,
#   Jmax,
#   alpha,
#   theta,
#   Rd
# ),
# data = CO2_grouped,
# fixed = list(
#   Vcmax ~ TREATMENT,
#   Jmax ~ TREATMENT,
#   Rd ~ TREATMENT,
#   alpha ~ TREATMENT,
#   theta ~ TREATMENT
# ),
#
# random = pdSymm(list(Vcmax ~ 1, Jmax ~ 1, Rd ~ 1, alpha ~ 1, theta ~ 1)),
# start = c(
#   Vcmax = c(Vcmax1, Vcmax1 - Vcmax2),
#   Jmax = c(Jmax1, Jmax1 - Jmax2),
#   Rd = c(Rd1, 0),
#   alpha = c(alpha1, 0),
#   theta = c(theta1, 0)
# ),
# control = list(
#   maxIter = 250, msVerbose = F, tolerance = 1e-3, msMaxIter = 250,
#   pnlsTol = 1e-10, pnlsMaxIter = 50
# )
# )
#

# === Yin step out model =======================================================

group_by_vars <- c("ID", "SIDE", "TREATMENT")
group_by_var <- "ID"

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

# === Calculate Rd, S ==========================================================

LRC2_side_AD <- LRC2_new$SIDE == "AD"
LRC2_AD <- LRC2_new %>%
  dplyr::filter(SIDE == "AD", TREATMENT == "A1")
# dplyr::filter(ID != 1)

LRC2_F_AD <- LRC2_F_new %>%
  dplyr::filter(SIDE == "AD", TREATMENT == "A1") # %>%
# dplyr::filter(ID != 1)

ggplot(data = LRC2_new, aes(x = PAR, y = Photosynthesis, colour = SIDE)) +
  geom_point() +
  facet_wrap(~ID)  +
  labs(
    title = "LCR curves",
    subtitle = "Coloured by SIDE", # google units cran package
    x = "Light intensity",
    y = "Net photosynthesis rate"
  )

LRC2_new %>%
  ggplot(aes(x = PAR, y = Photosynthesis, colour = SIDE)) +
  geom_point() +
  facet_wrap(~ID)  +
  labs(
    title = "LCR curves",
    subtitle = "All measurements",
    x = "Light intensity",
    y = "Net photosynthesis rate"
  )


LRC2_new %>%
  filter(PAR < 200) %>% 
  ggplot(aes(x = PAR, y = Photosynthesis, colour = SIDE)) +
    geom_point() +
    facet_wrap(~ID) +
  labs(
    title = "LCR curves",
    subtitle = "Intensity filtered below 200",
    x = "Light intensity",
    y = "Net photosynthesis rate"
  )


LRC2_new %>%
  filter(PAR < 100) %>% 
  ggplot(aes(x = PAR, y = Photosynthesis, colour = SIDE)) +
  geom_point() +
  facet_wrap(~ID)  +
  labs(
    title = "LCR curves",
    subtitle = "Coloured by SIDE, intensity filtered below 100", # google units cran package
    x = "Light intensity",
    y = "Net photosynthesis rate"
  )


LRC2_new %>%
  filter(PAR > 30, PAR < 100) %>% 
  ggplot(aes(x = PAR, y = Photosynthesis, colour = SIDE)) +
  geom_point() +
  facet_wrap(~ID) +
  labs(
    title = "LCR curves",
    subtitle = "Coloured by SIDE, intensity filtered below 100, above 30", # google units cran package
    x = "Light intensity",
    y = "Net photosynthesis rate"
  )

# Based on plots filter below 100 I

start <- list(Intercept = c(1), New_var = c(1))
lower <- list(Intercept = c(0), New_var = c(0))
upper <- list(Intercept = c(2), New_var = c(2))

i_lower_bound <- 30
i_upper_bound <- 100

out_AD <- calc_s_rd(LRC2_AD, LRC2_F_AD,
  A_param = Photosynthesis,
  phi_param = PhiPS2,
  i_inc_param = PAR,
  id_var = ID,
  other_var_to_keep = TREATMENT,
  i_lower_bound = i_lower_bound,
  i_upper_bound = i_upper_bound,
  control = lmeControl(
    opt = "nlminb",
    maxIter = 1000,
    msMaxIter = 1000,
    msVerbose = F,
    tolerance = 1e-4,
    pnlsTol = 1e-10,
    pnlsMaxIter = 500,
    niterEM = 1000,
    msMaxEval = 1000,
    sing.tol = 1e-20,
    start = start,
    lower = lower,
    upper = upper
  )
)


# small_model <- nlme::lme(Photosynthesis ~ New_var,
#   data = out_AD$data,
#   random = ~New_var | TREATMENT,
#   control = lmeControl(
#     opt = "nlminb",
#     maxIter = 1000,
#     msMaxIter = 1000,
#     msVerbose = F,
#     tolerance = 1e-4,
#     pnlsTol = 1e-10,
#     pnlsMaxIter = 500,
#     niterEM = 1000,
#     msMaxEval = 1000,
#     sing.tol = 1e-20,
#     lower = list(Intercept = 0)
#   )
# )

# small_model$coefficients$fixed
out_AD$model$coefficients$fixed

LRC2_AB <- LRC2_new %>%
  dplyr::filter(SIDE == "AB", TREATMENT == "A1")

LRC2_F_AB <- LRC2_F_new %>%
  dplyr::filter(SIDE == "AB", TREATMENT == "A1")


out_AB <- calc_s_rd(LRC2_AB, LRC2_F_AB,
  A_param = Photosynthesis,
  phi_param = PhiPS2,
  i_inc_param = PAR,
  id_var = ID,
  other_var_to_keep = TREATMENT,
  i_lower_bound = i_lower_bound,
  i_upper_bound = i_upper_bound,
  control = lmeControl(
    opt = "nlminb",
    maxIter = 1000,
    msMaxIter = 1000,
    msVerbose = F,
    tolerance = 1e-4,
    pnlsTol = 1e-10,
    pnlsMaxIter = 500,
    niterEM = 1000,
    msMaxEval = 1000,
    sing.tol = 1e-20,
    lower = list(Intercept = 0)
  )
)


LRC2_ADA2 <- LRC2_new %>%
  dplyr::filter(SIDE == "AD", TREATMENT == "A2")
# # dplyr::filter(ID != 1)
#
LRC2_F_ADA2 <- LRC2_F_new %>%
  dplyr::filter(SIDE == "AD", TREATMENT == "A2") # %>%
# dplyr::filter(ID != 1)
#
out_ADA2 <- calc_s_rd(LRC2_ADA2, LRC2_F_ADA2,
  A_param = Photosynthesis,
  phi_param = PhiPS2,
  i_inc_param = PAR,
  id_var = ID,
  other_var_to_keep = TREATMENT,
  i_lower_bound = i_lower_bound,
  i_upper_bound = i_upper_bound,
  control = list(
    maxIter = 1000,
    msMaxIter = 1000,
    msVerbose = F,
    tolerance = 1e-4,
    pnlsTol = 1e-10,
    pnlsMaxIter = 500,
    niterEM = 1000,
    msMaxEval = 1000,
    sing.tol = 1e-20,
    lower = list(Intercept = 0)
  )
)
#
LRC2_ABA2 <- LRC2_new %>%
  dplyr::filter(SIDE == "AB", TREATMENT == "A2")

LRC2_F_ABA2 <- LRC2_F_new %>%
  dplyr::filter(SIDE == "AB", TREATMENT == "A2")


out_ABA2 <- calc_s_rd(LRC2_ABA2, LRC2_F_ABA2,
  A_param = Photosynthesis,
  phi_param = PhiPS2,
  i_inc_param = PAR,
  id_var = ID,
  other_var_to_keep = TREATMENT,
  i_lower_bound = i_lower_bound,
  i_upper_bound = i_upper_bound,
  control = list(
    maxIter = 1000,
    msMaxIter = 1000,
    msVerbose = F,
    tolerance = 1e-4,
    pnlsTol = 1e-10,
    pnlsMaxIter = 500,
    niterEM = 1000,
    msMaxEval = 1000,
    sing.tol = 1e-20,
    lower = list(Intercept = 0)
  )
)


# head(out_AD$data)
# head(out_AB$data)
out_AB$data$SIDE <- "AB"
out_AD$data$SIDE <- "AD"

out_ABA2$data$SIDE <- "AB"
out_ADA2$data$SIDE <- "AD"

out <- dplyr::bind_rows(out_AD$data, out_AB$data)
out_A2 <- dplyr::bind_rows(out_ADA2$data, out_ABA2$data)
out$SIDE <- as.factor(out$SIDE)
out_A2$SIDE <- as.factor(out_A2$SIDE)
summary(out)
summary(out_A2)
out <- out %>% dplyr::select(-TREATMENT1)
out_A2 <- out_A2 %>% dplyr::select(-TREATMENT1)

out_full <- dplyr::bind_rows(out, out_A2)

ggplot(data = out_full, aes(y = S, x = SIDE, group = SIDE, colour = TREATMENT)) +
  facet_wrap(~ID) +
  geom_boxplot()

# ggplot(data = out, aes(y = Rd, x = ID, group = ID, colour = SIDE)) +
#   facet_wrap(~TREATMENT) +
#   geom_boxplot()

ggplot(data = out_full, aes(y = Rd, x = SIDE, group = SIDE, colour = TREATMENT)) +
  facet_wrap(~ID) +
  geom_boxplot()

out$SIDE_ID <- paste(out$ID, out$SIDE, sep = "_")

xyplot(Photosynthesis ~ New_var | as.factor(ID),
  group = as.factor(SIDE_ID),
  data = out,
  xlab = bquote(~frac(1, 4) ~ "I" ~ Phi[PSII]),
  ylab = "A"
)

out_full %>% filter(TREATMENT == "A2", ID == 2)


lme_lrc2_AD <- out_AD$model
lme_lrc2_AB <- out_AB$model

# LRC2_F_newer <- LRC2_F_new %>%
# left_join(out$coefficients, by = c("ID"))

plot_data <- out_AD$data %>%
  dplyr::mutate(prediction = predict(lme_lrc2_AD))

plot_data_AB <- out_AB$data %>%
  dplyr::mutate(prediction = predict(lme_lrc2_AB))

ggplot2::ggplot(data = plot_data, aes(x = New_var, y = Photosynthesis)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_line(aes(y = prediction), lty = 1, color = "steelblue", size = 1.4) +
  labs(
    title = bquote(~ "Regression for " ~ R[d] ~ "," ~ S ~ "(Side AD)"),
    x = bquote(~frac(1, 4) ~ "I" ~ Phi[PSII]),
    y = "A"
  )

ggplot2::ggplot(data = plot_data, aes(x = New_var, y = Photosynthesis)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = bquote(~ "Regression for " ~ R[d] ~ "," ~ S),
    x = bquote(~frac(1, 4) ~ "I" ~ Phi[PSII]),
    y = "A"
  )

ggplot2::ggplot(data = plot_data_AB, aes(x = New_var, y = Photosynthesis)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_line(aes(y = prediction), lty = 1, color = "steelblue", size = 1.4)+
  labs(
    title = bquote(~ "Regression for " ~ R[d] ~ "," ~ S ~ "(Side AB)"),
    x = bquote(~frac(1, 4) ~ "I" ~ Phi[PSII]),
    y = "A"
  )

ggplot2::ggplot(data = plot_data_AB, aes(x = New_var, y = Photosynthesis)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Regression for Rd, S",
    x = bquote(~frac(1, 4) ~ "I" ~ Phi[PSII]),
    y = "A"
  )

# Collect the data on the estimates of the parameters
S_rd_coeff <- out_full %>%
  select(ID, S, Rd, TREATMENT, SIDE) %>%
  group_by(S, ID, Rd, TREATMENT, SIDE) %>%
  distinct(S)

coef_est <- S_rd_coeff %>%
  group_by(TREATMENT, SIDE) %>%
  summarise(Rd_est = mean(Rd), S_est = mean(S), var_S = sd(S), var_Rd = sd(Rd))

ggplot(coef_est, aes(x = SIDE, colour = TREATMENT)) +
  geom_errorbar(aes(ymax = Rd_est + var_Rd, ymin = Rd_est - var_Rd),
    position = "dodge"
  ) +
  labs(
    title = bquote(~ "Variability of " ~ R[d] ~ " across effects"),
    y = "Rd",
    x = "Side"
  )

ggplot(coef_est, aes(x = SIDE, colour = TREATMENT)) +
  geom_errorbar(aes(ymax = S_est + var_S, ymin = S_est - var_S),
    position = "dodge"
  ) +
  labs(
    title = "Estimate of S",
    y = "S",
    x = "Side"
  )


point_est_AD <- out_AD$model$coefficients$fixed %>% t() %>% as.data.frame()
point_est_AD$SIDE <- "AD"
point_est_AD$TREATMENT <- "A1"

point_est_AB <- out_AB$model$coefficients$fixed %>% t() %>% as.data.frame()
point_est_AB$SIDE <- "AB"
point_est_AB$TREATMENT <- "A1"

point_est_ADA2 <- out_ADA2$model$coefficients$fixed %>% t() %>% as.data.frame()
point_est_ADA2$SIDE <- "AD"
point_est_ADA2$TREATMENT <- "A2"

point_est_ABA2 <- out_ABA2$model$coefficients$fixed %>% t() %>% as.data.frame()
point_est_ABA2$SIDE <- "AB"
point_est_ABA2$TREATMENT <- "A2"

std_errors_AD <- diag(sqrt(out_AD$summary$varFix)) %>% t() %>% as.data.frame()
std_errors_AD$SIDE <- "AD"
std_errors_AD$TREATMENT <- "A1"

std_errors_AB <- diag(sqrt(out_AB$summary$varFix)) %>% t() %>% as.data.frame()
std_errors_AB$SIDE <- "AB"
std_errors_AB$TREATMENT <- "A1"

std_errors_ADA2 <- diag(sqrt(out_ADA2$summary$varFix)) %>% t() %>% as.data.frame()
std_errors_ADA2$SIDE <- "AD"
std_errors_ADA2$TREATMENT <- "A2"

std_errors_ABA2 <- diag(sqrt(out_ABA2$summary$varFix)) %>% t() %>% as.data.frame()
std_errors_ABA2$SIDE <- "AB"
std_errors_ABA2$TREATMENT <- "A2"

est_df <- bind_rows(point_est_AB, point_est_AD, point_est_ABA2, point_est_ADA2)
est_df$Rd <- est_df$`(Intercept)`
est_df$S <- est_df$New_var
est_df <- est_df %>%
  select(SIDE, TREATMENT, Rd, S)

error_df <- bind_rows(std_errors_AD, std_errors_AB, std_errors_ADA2, std_errors_ABA2) %>%
  mutate(S_error = New_var)
error_df$Rd_error <- error_df$`(Intercept)`

error_df <- dplyr::select(error_df, SIDE, TREATMENT, Rd_error, S_error)

est_df <- est_df %>%
  left_join(error_df, by = c("SIDE", "TREATMENT"))


ggplot(est_df, aes(x = SIDE, colour = TREATMENT)) +
  geom_errorbar(aes(ymax = Rd + Rd_error, ymin = Rd - Rd_error),
    position = "dodge"
  ) +
  labs(
    title = bquote(~ "LME model estimate of " ~ R[d]), 
    y = bquote(~ R[d]), 
    x = "Side"
  )

ggplot(est_df, aes(x = SIDE, colour = TREATMENT)) +
  geom_errorbar(aes(ymax = S + S_error, ymin = S - S_error),
    position = "dodge"
  ) +
  labs(
    title = "LME model estimate of S",
    y = "S", x = "Side"
  )

treat_data <- out_AD$data %>%
  group_by(ID, TREATMENT) %>%
  distinct(ID)

rd_s_group_est <- bind_rows(point_est_AD, point_est_AB, point_est_ADA2, point_est_ABA2)


# === Calculate J ==============================================================

LRC2_F_newer <- LRC2_F_new %>%
  filter(SIDE != "ADAB")

LRC2_F_newer$SIDE <- droplevels(LRC2_F_newer$SIDE)

LRC2_F_newest <- LRC2_F_newer %>%
  left_join(S_rd_coeff, by = c("ID", "TREATMENT", "SIDE"))

J_data <- calc_J(LRC2_F_newest, lump_var = S, i_inc_var = PAR, phi_psII_var = PhiPS2)

head(J_data)

out_coefficients <- S_rd_coeff

total_data <- total_data %>%
  left_join(out_coefficients, by = c("ID"))

# CO2_grouped_full <- groupedData(Photosynthesis ~ 1 | ID,
#   outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
#   inner = ~Ci + PAR + SIDE,
#   data = total_data
# )

# === Sco ======================================================================

# ACI1 <- read.table("A-Ci1_2.csv", sep = ",", header = TRUE, na.strings = "")
data_h <- ACI1 %>% dplyr::filter(Ci < 200) %>% filter(SIDE == "AD")
data_h2 <- ACI1 %>% filter(SIDE == "AD")


ggplot2::ggplot(data = data_h, aes(y = Photo, x = Ci)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Ci per individual",
    x = bquote(~C[i]),
    y = "A"
  )

ggplot2::ggplot(data = data_h2, aes(y = Photo, x = Ci)) +
  facet_wrap(~ID) +
  geom_point() +
  geom_smooth() +
  labs(
    title = "Ci per individual",
    x = bquote(~C[i]),
    y = "A"
  )

# # out_coefficients <- bind_rows(out_AB$coefficients, out_AD$coefficients)
# out_coefficients$SIDE <- as.factor(c(
#   rep("AB", nrow(out_AB$coefficients)),
#   rep("AD", nrow(out_AD$coefficients))
# ))
#
# out_coefficients$SIDE <- out_coefficients$SIDE %>%
#   add_level("ADAB")


ID_quo <- quo(ID)
out_coefficients$SIDE <- add_level(out_coefficients$SIDE, "ADAB")

ACI1_RD <- left_join(ACI1_new, out_coefficients,
  by = c("ID", "SIDE", "TREATMENT")
) %>%
  filter(ID != 1)

# # filter(SIDE = "AB")
# ACI1_RD_AD <- ACI1_new %>%
#   left_join(ACI1_new, out_coefficients, by = c(quo_name(ID_quo) = quo_name(ID_quo),
#             quo_name(SIDE_quo) = quo_name(SIDE_quo))) %>%
# filter(ID != 1)

names(ACI1_RD)
unique(ACI1_RD$ID)

unique(ACI2_new$ID)
ACI2_2 <- ACI2_new %>% dplyr::filter(ID != 1) %>% na.omit()
unique(ACI2_2$ID)

ACI1_RD_AD <- ACI1_RD %>%
  filter(SIDE == "AD")

ACI2_2_AD <- ACI2_2 %>%
  filter(SIDE == "AD")

unique(ACI1_RD_AD$ID)
unique(ACI2_2_AD$ID)

ACI1_RD_AD <- ACI1_RD_AD %>% dplyr::filter(ID != 8)

summary(ACI1_RD_AD)
summary(ACI2_2_AD)

data_sco_AD <- calc_Sco(ACI1_RD_AD, ACI2_2_AD,
  Ci_threshold = 200,
  id_var = ID,
  A_var = Photosynthesis,
  other_var = TREATMENT,
  control = lmeControl(
    opt = "nlminb",
    maxIter = 1000,
    msMaxIter = 1000,
    msVerbose = F,
    tolerance = 1e-4,
    pnlsTol = 1e-10,
    pnlsMaxIter = 500,
    niterEM = 1000,
    msMaxEval = 1000,
    sing.tol = 1e-20,
    lower = list(Intercept = 0)
  )
)

head(data_sco_AD$S_co)
head(data_sco_AD)

data_sco_AD <- data_sco_AD %>%
  dplyr::mutate(TREATMENT = TREATMENT.x, SIDE = "AD") %>%
  dplyr::select(-TREATMENT.x, -TREATMENT1, -TREATMENT.y)

ggplot2::ggplot(data = data_sco_AD, aes(x = S_co)) +
  geom_density() +
  facet_wrap(~ID)

ggplot2::ggplot(data = data_sco_AD, aes(y = S_co, x = ID, group = ID, colour = TREATMENT)) +
  # facet_wrap(~ID) +
  geom_boxplot()

# sco_data_AD <- data_sco_AD %>%
#   group_by(ID) %>%
#   filter(!(abs(S_co - median(S_co)) > 2 * sd(S_co)))
# 
# # ggplot2::ggplot(data = sco_data_AD, aes(x = S_co)) +
# #   geom_density() +
# #   facet_wrap(~ID)
# #
# # ggplot2::ggplot(data = sco_data_AD, aes(y = S_co, group = ID, colour = SIDE)) +
# #   geom_boxplot()
# 
# 
# ggplot2::ggplot(data = sco_data_AD, aes(y = S_co, x = ID, group = ID, colour = TREATMENT)) +
#   facet_wrap(~ID) +
#   geom_boxplot()

# sco_data_2 <- sco_data_AD %>%
#   group_by(ID) %>%
#   filter(!(abs(S_co - median(S_co)) > 2 * sd(S_co)))

# ggplot2::ggplot(data = sco_data_2, aes(x = S_co)) +
#   geom_density() +
#   facet_wrap(~ID)

# ggplot2::ggplot(data = sco_data_2, aes(y = S_co, x = ID, group = ID, colour = TREATMENT)) +
#   facet_wrap(~ID) +
#   geom_boxplot()


ACI1_RD_AB <- ACI1_RD %>%
  filter(SIDE == "AB")

ACI2_2_AB <- ACI2_2 %>%
  filter(SIDE == "AB")

unique(ACI1_RD_AB$ID)
unique(ACI2_2_AB$ID)

# ACI1_RD_AB <- ACI1_RD_AB %>% dplyr::filter(ID != 8)

summary(ACI1_RD_AB)
summary(ACI2_2_AB)

data_sco_AB <- calc_Sco(ACI1_RD_AB, ACI2_2_AB,
  Ci_threshold = 200,
  id_var = ID,
  A_var = Photosynthesis,
  other_var = TREATMENT,
  control = lmeControl(
    opt = "nlminb",
    maxIter = 1000,
    msMaxIter = 1000,
    msVerbose = F,
    tolerance = 1e-4,
    pnlsTol = 1e-10,
    pnlsMaxIter = 500,
    niterEM = 1000,
    msMaxEval = 1000,
    sing.tol = 1e-20,
    lower = list(Intercept = 0)
  )
)
data_sco_AB <- data_sco_AB %>%
  dplyr::mutate(TREATMENT = TREATMENT.x, SIDE = "AB") %>%
  dplyr::select(-TREATMENT.x, -TREATMENT1, -TREATMENT.y)

ggplot2::ggplot(data = data_sco_AB, aes(y = S_co, x = ID, group = ID, colour = TREATMENT)) +
  # facet_wrap(~ID) +
  geom_boxplot()

# sco_data_2_AB <- data_sco_AB %>%
#   group_by(ID) %>%
#   filter(!(abs(S_co - median(S_co)) > 2 * sd(S_co)))

ggplot2::ggplot(data = sco_data_2_AB, aes(y = S_co, x = ID, group = ID, colour = TREATMENT)) +
  facet_wrap(~ID) +
  geom_boxplot()

sco_data <- bind_rows(data_sco_AB, data_sco_AD)

S_co_est_df <- sco_data %>%
  group_by(TREATMENT, SIDE) %>%
  summarise(S_co_mean = mean(S_co))

S_co_est <- mean(sco_data$S_co)
names(sco_data)
sco_data_rel <- sco_data %>%
  dplyr::select(ID, SIDE, S_co) %>%
  dplyr::group_by(ID, SIDE) %>%
  dplyr::summarise(mean_sco = mean(S_co)) %>%
  dplyr::mutate(S_co = mean_sco) %>%
  dplyr::select(-mean_sco)

if (group_by_var == "ID") {
  out_coefficients_all <- left_join(out_coefficients, sco_data_rel, by = c("ID", "SIDE"))
} else {
  out_coefficients_all <- left_join(out_coefficients, S_co_est_df, by = c("TREATMENT", "SIDE"))
}
out_coefficients_all <- left_join(out_coefficients, sco_data_rel, by = c("ID", "SIDE"))



ACI1_RD <- dplyr::left_join(ACI1_new, out_coefficients_all, by = c("ID", "SIDE", "TREATMENT"))
ACI1_with_j <- calc_J(ACI1_RD, lump_var = S, i_inc_var = PAR, phi_psII_var = PhiPS2)
summary(ACI1_with_j$J)

LRC1_RD <- dplyr::left_join(LRC1_new, out_coefficients_all, by = c("ID", "SIDE", "TREATMENT"))
# LRC1_with_j <- calc_J(LRC1_RD, lump_var = S, i_inc_var = PAR, phi_psII_var = PhiPS2)

J_max <- ACI1_with_j %>% 
  group_by(ID, TREATMENT, SIDE) %>% 
  summarise(max_J = max(J), max_PAR = max(PAR))

J_max$J <- J_max$max_J

J_max %>% na.omit() %>%  left_join(select(ACI1_with_j, J, PAR), by = "J")

ggplot(J_max, aes(y = max_J, x = SIDE, group = SIDE, colour = TREATMENT)) +
  geom_point() +
  facet_wrap(~ID) +
  labs(
    title = bquote(~J[1300]),
    y = "J", x = "Side"
  )

ACI1_final <- na.omit(ACI1_with_j)
ACI1_final %>%
  dplyr::group_by(ID) %>%
  dplyr::count()
names(ACI1_final)


ACI1_final <- ACI1_final %>% filter(ID != 1)

# === Calculate G_star =========================================================

ACI1_final <- calc_G_star(ACI1_final, S_co_value = S_co_est) %>%
  mutate(G_star1 = G_star)

ACI1_final <- calc_G_star(ACI1_final)

# LRC1_g_star <- calc_G_star(LRC1_RD, S_co_value = S_co_est)
# calc_J(LRC1_g_star, lump_var = S, i_inc_var = PAR, phi_psII_var = PhiPS2)


# === Calculate g_m ============================================================

ACI1_final %>% select(J) %>% glimpse()


# gm estimate using specific values of Ci
ACI1_final_gm <- calc_gm(ACI1_final %>% filter(Ci > 250 & Ci < 300),
  A_par = Photosynthesis,
  Ci_par = Ci,
  Rd_par = Rd,
  G_star_par = G_star,
  J_par = J
)

summary(ACI1_final_gm)

odd_gm <- ACI1_final_gm %>%
  filter(gm < 0.025 | gm > 0.35)
# glimpse()
summary(odd_gm)
head(odd_gm)
nrow(odd_gm)


summary(ACI1_final$gm)
gm_value_ACI1 <- median(ACI1_final_gm$gm) #[(ACI1_final_gm$gm > 0.05) & (ACI1_final_gm$gm < 1.0)])

# this is strangely low, use a vlaue of 0.29

ACI1_final <- ACI1_final %>%
  dplyr::mutate(gm_est = gm_value_ACI1 * 10) %>%
  dplyr::mutate(Cc = Ci - Photosynthesis / gm_est)

# === Data wrangling ===========================================================

head(ACI1_final %>% select(Rd, S, G_star, Cc, Ci, gm_est, Gstar))
head(ACI1_final)

ggplot(
  data = ACI1_final[select_AB, ],
  aes(x = Ci, y = Photosynthesis, colour = TREATMENT, shape = SIDE)
) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ID)


# ACI1_final <- ACI1_final %>% dplyr::mutate(Kc = 1.8 * Kc )

select_AD <- ACI1_final$SIDE == "AD"
head(ACI1_final)
ACI1_final <- ACI1_final %>%
  dplyr::mutate(K_c = Kc * 1.65)

J_est_treat <- ACI1_final %>%
  group_by(TREATMENT, SIDE) %>%
  summarise(J_est = mean(J))

ACI1_final <- ACI1_final %>% left_join(J_est_treat, by = c("TREATMENT", "SIDE"))

ACI1_final$ID <- as.factor(ACI1_final$ID)
summary(ACI1_final$ID)
grouped_ACI <- groupedData(Photosynthesis ~ 1 | ID,
  outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
  inner = ~SIDE,
  data = ACI1_final
)
summary(ACI1_final)

ggplot(data = ACI1_final, aes(y = J_est, x = SIDE, group = SIDE, colour = TREATMENT)) +
  geom_boxplot() +
  facet_wrap(~ID)


ACI1_final$K_o <- ACI1_final$Ko * 0.8

# === Indices for subgroups=====================================================
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

names(ACI1_final)
head(ACI1_final)

grouped_ACI <- groupedData(Photosynthesis ~ 1 | ID,
  outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
  inner = ~SIDE,
  data = ACI1_final
)


grouped_ACI$SIDE <- as.factor(grouped_ACI$SIDE)
summary(grouped_ACI)



# === NLS models ===============================================================
nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star, Kc, K_o, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)



nls_model_A1_AD <- nls(Photosynthesis ~ FvCB(Cc, G_star, Kc, K_o, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AD,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AD <- nls(Photosynthesis ~ FvCB(Cc, G_star, Kc, K_o, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AD,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# === Jmax used =====

jmax_i_a1_ab <- max(ACI1_final[select_A1_AB,]$J)

nls_model_A1_AB_jmax <- nls(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, S, theta),
                       data = ACI1_final,
                       subset = select_A1_AB,
                       start = c(
                         Vcmax = c(60),
                         Jmax = jmax_i_a1_ab,
                         theta = 0.8
                       ),
                       control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

jmax_i_a1_ad <- max(ACI1_final[select_A1_AD,]$J)

nls_model_A1_AD_jmax <- nls(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, S, theta),
                            data = ACI1_final,
                            subset = select_A1_AD,
                            start = c(
                              Vcmax = c(40),
                              Jmax = jmax_i_a1_ad,
                              theta = 0.8
                            ),
                            control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

jmax_i_a2_ab <- max(ACI1_final[select_A2_AB,]$J)

nls_model_A2_AB_jmax <- nls(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, S, theta),
                            data = ACI1_final,
                            subset = select_A2_AB,
                            start = c(
                              Vcmax = c(40),
                              Jmax = jmax_i_a2_ab,
                              theta = 0.8
                            ),
                            control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

jmax_i_a2_ad <- max(ACI1_final[select_A2_AD,]$J)

nls_model_A1_AD_jmax <- nls(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, S, theta),
                            data = ACI1_final,
                            subset = select_A2_AD,
                            start = c(
                              Vcmax = c(80),
                              Jmax = jmax_i_a2_ad,
                              theta = 0.8
                            ),
                            control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# === Try different G_star =======

nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star1, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star1, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# Try different K_c

nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# Both
nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star1, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star1, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# Try Ci
nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)



nls_model_A1_AD <- nls(Photosynthesis ~ FvCB(Ci, G_star, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AD,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AD <- nls(Photosynthesis ~ FvCB(Ci, G_star, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AD,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)


# Try different G_star

nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star1, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star1, Kc, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# Try different K_c

nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# Both
nls_model_A1_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star1, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A1_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Ci, G_star1, K_c, Ko, O, Vcmax, J, Rd),
  data = ACI1_final,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(40)
  ),
  control = list(maxiter = 15000, minFactor = 1e-20, printEval = F, tol = 1e-3)
)

# === NLME model ===============================================================

model_1 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ 1,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(40)
  ),
  subset = select_A1_AB,
  method = "REML",
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = T,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    sing.tol = 1e-45
  )
)


model_2 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ 1,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(40)
  ),
  subset = select_A1_AD,
  method = "REML",
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = T,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0)
    ),
    sing.tol = 1e-45
  )
)

for (i in 10:90) {
  model_3 <- attempt::try_catch(
    expr = nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
      data = grouped_ACI,
      fixed = Vcmax ~ 1,
      random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
      start = c(
        Vcmax = c(i)
      ),
      subset = select_A2_AB,
      method = "REML",
      control = nlmeControl(
        opt = "nlminb",
        # upper = c(
        #   Vcmax = c(500, 500)
        # ),

        maxIter = 500,
        msMaxIter = 500,
        msVerbose = F,
        msTol = 1e-2,
        pnlsTol = 1e-30,
        pnlsMaxIter = 500,
        niterEM = 500,
        msMaxEval = 500,
        minScale = 1e-20,
        returnObject = T,
        lower = c(
          Vcmax = c(0)
        ),
        sing.tol = 1e-50
      )
    ),
    .e = function(a) return(NULL),
    .w = function(a) return(NULL)
  )
  if (!is.null(model_3)) {
    print(paste("works for", i))
    break
  }
}

model_4 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ 1,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(40)
  ),
  subset = select_A2_AD,
  method = "REML",
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = T,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    iter.max = 500,
    sing.tol = 1e-45
  )
)

model_1$coefficients # A1 AB
model_2$coefficients # A1 AD
model_3$coefficients # A2 AB
model_4$coefficients # A2 AD

# ===  NLME model for Side == AD ===============================================
model_5 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ TREATMENT,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(42, -7)
  ),
  subset = select_AD,
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = T,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    iter.max = 500,
    sing.tol = 1e-45
  )
)

model_5_reml <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ TREATMENT,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(42, -4)
  ),
  subset = select_AD,
  method = "REML",
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = F,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    iter.max = 500,
    sing.tol = 1e-45
  )
)

model_5a <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ 1,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(49)
  ),
  subset = select_AD,
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = T,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    iter.max = 500,
    sing.tol = 1e-45
  )
)

model_5b <- nlme(Photosynthesis ~ FvCB(Cc, G_star, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ TREATMENT,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(48, -4)
  ),
  subset = select_AD,
  method = "REML",
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = F,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    iter.max = 500,
    sing.tol = 1e-45
  )
)

summary(model_5)
summary(model_5a)

nlme_5_fixed <- summary(model_5)$tTable
nlme_5_random <- summary(model_5)$coefficients$random

# Compare estimates for Treatments
fixed_comparison <- as.data.frame(nlme_5_fixed)
fixed_comparison$Feature <- rownames(fixed_comparison)
fixed_comparison$Feature <- c(
  rep("Vcmax", 2)
)
# head(grouped_ACI)
fixed_comparison$Treatment <- rep(c("A1", "A2"), nrow(fixed_comparison) / 2)


comparison_nlme_plot <- fixed_comparison %>%
  mutate(Estimate = ifelse(Treatment == "A1", Value, lag(Value) + Value))
# select(Estimate, Std_error, Feature, Treatment)

ggplot(comparison_nlme_plot, aes(x = Treatment, colour = Treatment)) +
  geom_errorbar(aes(ymax = Estimate + Std.Error, ymin = Estimate - Std.Error),
    position = "dodge"
  ) +
  labs(
    title = "NLME model estimate for Side AD",
    y = "Vcmax", x = "Treatment"
  )

# === NLME model for SIde == AB ================================================




model_6 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ TREATMENT,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(49, -4)
  ),
  subset = select_AB,
  method = "REML",
  control = nlmeControl(
    opt = "nlminb",
    # upper = c(
    #   Vcmax = c(500, 500)
    # ),

    maxIter = 500,
    msMaxIter = 500,
    msVerbose = F,
    msTol = 1e-2,
    pnlsTol = 1e-30,
    pnlsMaxIter = 500,
    niterEM = 500,
    msMaxEval = 500,
    minScale = 1e-20,
    returnObject = T,
    lower = c(
      Vcmax = c(0, 0)
    ),
    iter.max = 500,
    sing.tol = 1e-45
  )
)

for (i in 10:90) {
  model_7 <- attempt::try_catch(
    expr = nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
      data = grouped_ACI,
      fixed = Vcmax ~ TREATMENT,
      random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
      start = c(
        Vcmax = c(i, -(i - 42))
      ),
      subset = select_AB,
      method = "REML",
      control = nlmeControl(
        opt = "nlminb",
        # upper = c(
        #   Vcmax = c(500, 500)
        # ),

        maxIter = 500,
        msMaxIter = 500,
        msVerbose = F,
        msTol = 1e-2,
        pnlsTol = 1e-30,
        pnlsMaxIter = 500,
        niterEM = 500,
        msMaxEval = 500,
        minScale = 1e-20,
        returnObject = T,
        lower = c(
          Vcmax = c(0, 0)
        ),
        iter.max = 500,
        sing.tol = 1e-45
      )
    ),
    .e = function(a) return(NULL),
    .w = function(a) return(NULL)
  )
  if (!is.null(model_7)) {
    print(paste("works for", i))
    break
  }
}

summary(ACI1_final[select_A2_AB, ])

ACI1_gm_out <- ACI1_final[ACI1_final$gm > 0, ]

# === Indices for subgroups=====================================================

select_AD_gm <- ACI1_gm_out$SIDE == "AD"
select_AB_gm <- ACI1_gm_out$SIDE == "AB"

select_A1_AD_gm <- as.logical((ACI1_gm_out$TREATMENT == "A1") * (ACI1_gm_out$SIDE == "AD"))
select_A1_AB_gm <- as.logical((ACI1_gm_out$TREATMENT == "A1") * (ACI1_gm_out$SIDE == "AB"))

select_A2_AD_gm <- as.logical((ACI1_gm_out$TREATMENT == "A2") * (ACI1_gm_out$SIDE == "AD"))
select_A2_AB_gm <- as.logical((ACI1_gm_out$TREATMENT == "A2") * (ACI1_gm_out$SIDE == "AB"))

grouped_ACI_gm <- groupedData(Photosynthesis ~ 1 | ID,
  outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
  inner = ~SIDE,
  data = ACI1_gm_out
)

grouped_ACI_gm$K_o <- grouped_ACI_gm$Ko * 0.8

for (i in 10:90) {
  model_8 <- attempt::try_catch(
    expr = nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
      data = grouped_ACI_gm,
      fixed = Vcmax ~ 1,
      random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
      start = c(
        Vcmax = c(i)
      ),
      subset = select_A2_AB_gm,
      control = nlmeControl(
        opt = "nlminb",
        # upper = c(
        #   Vcmax = c(500, 500)
        # ),

        maxIter = 500,
        msMaxIter = 500,
        msVerbose = F,
        msTol = 1e-2,
        pnlsTol = 1e-30,
        pnlsMaxIter = 500,
        niterEM = 500,
        msMaxEval = 500,
        minScale = 1e-20,
        returnObject = T,
        lower = c(
          Vcmax = c(0)
        ),
        sing.tol = 1e-50
      )
    ),
    .e = function(a) return(NULL),
    .w = function(a) return(NULL)
  )
  if (!is.null(model_8)) {
    print(paste("works for", i))
    break
  }
}

model_8 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI_gm,
  fixed = Vcmax ~ 1,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(48)
  ),
  subset = select_A2_AB_gm,
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 500,
    msMaxIter = 500,
    minScale = 1e-02,
    tolerance = 1e-01,
    msVerbose = T,
    niterEM = 500,
    pnlsTol = 1e-03,
    msTol = 1e-01,
    returnObject = T,
    opt = "nlminb",
    msMaxEval = 500,
    lower = c(
      Vcmax = c(0)
    ),
    sing.tol = 1e-20
  )
)

nls_model_A2_AB <- nls(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI_gm,
  subset = select_A2_AB,
  start = c(
    Vcmax = c(50)
  ),
  control = list(maxiter = 25000, minFactor = 1e-20, printEval = F, tol = 1e-02)
)


model_9 <- nlme(Photosynthesis ~ FvCB(Cc, G_star1, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ TREATMENT,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(49, -6)
  ),
  # subset = select_A2_AB_gm,
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 500,
    msMaxIter = 500,
    minScale = 1e-02,
    tolerance = 1e-01,
    msVerbose = T,
    niterEM = 500,
    pnlsTol = 1e-03,
    msTol = 1e-01,
    returnObject = T,
    opt = "nlminb",
    msMaxEval = 500,
    lower = c(
      Vcmax = c(0)
    ),
    sing.tol = 1e-20
  )
)


model_9 <- nlme(Photosynthesis ~ FvCB(Cc, G_star, Kc, K_o, O, Vcmax, J, Rd),
  data = grouped_ACI,
  fixed = Vcmax ~ TREATMENT,
  random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
  start = c(
    Vcmax = c(49, -6)
  ),
  # subset = select_A2_AB_gm,
  method = "REML",
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 500,
    msMaxIter = 500,
    minScale = 1e-02,
    tolerance = 1e-01,
    msVerbose = F,
    niterEM = 500,
    pnlsTol = 1e-03,
    msTol = 1e-01,
    returnObject = T,
    opt = "nlminb",
    msMaxEval = 500,
    lower = c(
      Vcmax = c(0)
    ),
    sing.tol = 1e-20
  )
)

# === Try unique S_c/o per plant ===============================================


ACI1_diff_sco <- ACI1_final %>%
  # left_join(S_co_est_df, by = c("TREATMENT", "SIDE")) %>%
  dplyr::mutate(G_star = 0.5 * O / S_co_mean)

summary(ACI1_diff_sco)

ACI1_diff_sco <- calc_gm(ACI1_diff_sco,
  A_par = Photosynthesis,
  Ci_par = Ci,
  Rd_par = Rd,
  G_star_par = G_star,
  J_par = J
)

odd_gm_sco <- ACI1_diff_sco %>%
  filter(gm < 0.025 | gm > 0.35)
# glimpse()
summary(odd_gm_sco)
head(odd_gm_sco)
nrow(odd_gm_sco)
gm_value_sco <- median(ACI1_diff_sco$gm[(ACI1_diff_sco$gm > 0.005) & (ACI1_diff_sco$gm < 0.40)])

ACI1_diff_sco <- ACI1_diff_sco %>%
  dplyr::mutate(gm_est = gm_value_sco) %>%
  dplyr::mutate(Cc = Ci - Photosynthesis / 0.16)

summary(ACI1_diff_sco)

grouped_g_star <- groupedData(Photosynthesis ~ 1 | ID,
                           outer = ~TREATMENT, # !!random_effects[[n_random_effects]]
                           inner = ~SIDE,
                           data = ACI1_diff_sco
)
grouped_g_star$K_o <- grouped_g_star$Ko * 0.8

head(ACI1_diff_sco)

model_g_star <- nlme(Photosynthesis ~ FvCB(Cc, G_star, Kc, K_o, O, Vcmax, J, Rd),
                     data = grouped_g_star,
                     fixed = Vcmax ~ TREATMENT,
                     random = pdDiag(Vcmax ~ ID), #  %in% TREATMENT
                     start = c(
                       Vcmax = c(49, -6)
                     ),
                     # subset = select_AD,
                     method = "REML",
                     control = nlmeControl(
                       opt = "nlminb",
                       # upper = c(
                       #   Vcmax = c(500, 500)
                       # ),
                       
                       maxIter = 1500,
                       msMaxIter = 500,
                       msVerbose = F,
                       msTol = 1e-2,
                       pnlsTol = 1e-30,
                       pnlsMaxIter = 500,
                       niterEM = 500,
                       msMaxEval = 500,
                       minScale = 1e-20,
                       returnObject = T,
                       lower = c(
                         Vcmax = c(0, 0)
                       ),
                       iter.max = 500,
                       sing.tol = 1e-45
                     )
)


# === Alternative gm ===========================================================

ACI1_gm <- calc_gm(ACI1_final,
                      A_par = Photosynthesis,
                      Ci_par = Ci,
                      Rd_par = Rd,
                      G_star_par = G_star1,
                      J_par = J
)
summary(ACI1_gm)

odd_gm <- ACI1_gm %>%
  filter(gm < 0.025 | gm > 0.35)
# glimpse()
summary(odd_gm)
head(odd_gm)
nrow(odd_gm)


summary(ACI1_gm$gm)
gm_value_ACI1_gm <- median(ACI1_gm$gm[(ACI1_gm$gm > 0.05) & (ACI1_gm$gm < 1.0)])

ACI1_gm <- ACI1_gm %>%
  dplyr::mutate(gm_est = gm_value_ACI1) %>%
  dplyr::mutate(Cc = Ci * gm_est)


plot_pres_data <- ACI1_final %>% filter(ID == 3, SIDE == "AD")

ggplot(plot_pres_data, aes(x = Ci, y = Photosynthesis)) +
  geom_point() +
  geom_smooth(se = F) +
  labs(
    title = "ACI curve",
    x = bquote(~C[i]),
    y = "A"
  )

plot_pres_data_LRC <- LRC1_new %>% filter(ID == 10, SIDE == "AD")

ggplot(plot_pres_data_LRC, aes(x = PAR, y = Photosynthesis)) +
  geom_point() +
  geom_smooth(se = F, method = "loess") +
  labs(
    title = "Light response curve",
    x = "PAR",
    y = "A"
  )  


# === NLME with Jmax ===========================================================



nlme_model_jmax <- nlme(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, S, theta),
                data = grouped_ACI,
                fixed = list(Vcmax ~ TREATMENT, Jmax ~ TREATMENT, theta ~ TREATMENT),
                random = pdDiag(list(Vcmax ~ ID, JMAX ~ 1, theta ~ 1)), #  %in% TREATMENT
                start = c(
                  Vcmax = c(50, 0),
                  Jmax = c(jmax_i_a1_ab, jmax_i_a2_ab - jmax_i_a1_ab),
                  theta = c(2, 0)
                ),
                subset = select_AB,
                method = "REML",
                control = nlmeControl(
                  opt = "nlminb",
                  # upper = c(
                  #   Vcmax = c(500, 500)
                  # ),
                  
                  maxIter = 500,
                  msMaxIter = 500,
                  msVerbose = F,
                  msTol = 1e-2,
                  pnlsTol = 1e-30,
                  pnlsMaxIter = 500,
                  niterEM = 500,
                  msMaxEval = 500,
                  minScale = 1e-20,
                  returnObject = T,
                  lower = c(
                    Vcmax = c(0, 0)
                  ),
                  iter.max = 500,
                  sing.tol = 1e-45
                )
)

s_a1_ab <- median(grouped_ACI$S[select_A1_AB])

nlme_model_jmax_a1_ab <- nlme(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, alpha, theta),
                        data = grouped_ACI,
                        fixed = list(Vcmax ~ 1, Jmax ~ 1, alpha ~ 1, theta ~ 1),
                        random = pdDiag(list(Vcmax ~ ID, JMAX ~ 1, alpha ~ 1, theta ~ 1)), #  %in% TREATMENT
                        start = c(
                          Vcmax = c(40),
                          Jmax = c(jmax_i_a1_ab),
                          theta = c(2),
                          alpha = c(s_a1_ab)
                        ),
                        subset = select_A1_AB,
                        method = "REML",
                        control = nlmeControl(
                          opt = "nlminb",
                          # upper = c(
                          #   Vcmax = c(500, 500)
                          # ),
                          
                          maxIter = 500,
                          msMaxIter = 500,
                          msVerbose = F,
                          msTol = 1e-2,
                          pnlsTol = 1e-30,
                          pnlsMaxIter = 500,
                          niterEM = 500,
                          msMaxEval = 500,
                          minScale = 1e-20,
                          returnObject = T,
                          lower = c(
                            Vcmax = c(0, 0)
                          ),
                          iter.max = 500,
                          sing.tol = 1e-45
                        )
)


nlme_model_jmax_a1_ab <- nlme(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd, PAR, alpha, theta),
                              data = grouped_ACI,
                              fixed = list(Vcmax ~ 1, Jmax ~ 1, alpha ~ 1, theta ~ 1),
                              random = pdDiag(list(Vcmax ~ 1, JMAX ~ 1, alpha ~ 1, theta ~ 1)), #  %in% TREATMENT
                              start = c(
                                Vcmax = c(80),
                                Jmax = c(jmax_i_a1_ab),
                                theta = c(1),
                                alpha = c(s_a1_ab)
                              ),
                              subset = select_A1_AB,
                              method = "REML",
                              control = nlmeControl(
                                opt = "nlminb",
                                # upper = c(
                                #   Vcmax = c(500, 500)
                                # ),
                                
                                maxIter = 500,
                                msMaxIter = 500,
                                msVerbose = F,
                                msTol = 1e-2,
                                pnlsTol = 1e-30,
                                pnlsMaxIter = 500,
                                niterEM = 500,
                                msMaxEval = 500,
                                minScale = 1e-20,
                                returnObject = T,
                                lower = c(
                                  Vcmax = c(0, 0)
                                ),
                                iter.max = 500,
                                sing.tol = 1e-45
                              )
)

rd_a1_ab <- median(grouped_ACI$Rd[select_A1_AB])

nlme_model_jmax_a1_ab <- nlme(Photosynthesis ~ FvCB_jmax(Cc, G_star, Kc, Ko, O, Vcmax, Jmax, Rd_var, PAR, alpha, theta),
                              data = grouped_ACI,
                              fixed = list(Vcmax ~ 1, Jmax ~ 1, Rd_var ~ 1, alpha ~ 1, theta ~ 1),
                              random = pdDiag(list(Vcmax ~ 1, JMAX ~ 1, Rd_var ~ 1, alpha ~ 1, theta ~ 1)), #  %in% TREATMENT
                              start = c(
                                Vcmax = c(40),
                                Jmax = c(jmax_i_a1_ab),
                                Rd_var = 1,
                                theta = c(-1),
                                alpha = c(s_a1_ab)
                              ),
                              subset = select_A1_AB,
                              method = "REML",
                              control = nlmeControl(
                                opt = "nlminb",
                                # upper = c(
                                #   Vcmax = c(500, 500)
                                # ),
                                
                                maxIter = 500,
                                msMaxIter = 500,
                                msVerbose = F,
                                msTol = 1e-2,
                                pnlsTol = 1e-30,
                                pnlsMaxIter = 500,
                                niterEM = 500,
                                msMaxEval = 500,
                                minScale = 1e-20,
                                returnObject = T,
                                lower = c(
                                  Vcmax = c(0, 0)
                                ),
                                iter.max = 500,
                                sing.tol = 1e-45
                              )
)

