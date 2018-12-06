
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


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
  numerator <- C_c - G_star
  denom <- C_c + K_c * (1 + O / K_o)
  return(V_cmax * (numerator / denom) - R_d)
}

# Check Rubisco_limited function
# little_data <- data.frame(V_cmax = 1,
#                           C_c = 1,
#                           G_star = c(3, 4),
#                           K_c = 4,
#                           O = c(5, 8),
#                           K_o = 3,
#                           R_d = 1
#                           )
# Rubisco_limited(little_data$V_cmax,
#                 little_data$C_c,
#                 little_data$G_star,
#                 little_data$K_c,
#                 little_data$O,
#                 little_data$K_o,
#                 little_data$R_d
#                 )

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
  p2 <- ((non_rectangular_hyperbola(PAR, J_max, alpha, theta) / 4) * (C_c - G_star)
         / (C_c  + 2 * G_star)
  )
  ifelse(p1 < p2, p1, p2) - R_d
}

#' @title TPU limited curve
#'
#' @description Returns the expected CO2 assimilation when the rate is limited
#' by Triose Phosphate Use (TPU).
#' @param TPU is the rate of use of triose phosphates but can also be any export
#' of carbon from the Calvin cycle including direct use of photorespiratory
#' glycine or serine.
#' @param R_d (
#' %(\eqn{R_d})%
#' ) is respiratory \eqn{CO_2} release by methods other than by
#' photorespiration and is presumed to be primarily mitochondiral respiraion.
#' @return Nothing at all.
#' @examples
#' TPU_limited(<examples>)
#' More examples and whatnot
TPU_limited <- function(TPU, R_d) {
  return(3 * TPU - R_d)
}

# Pseudo code for fitting aci curve - use in nlme?
#' @title ACI Curve
#' @description Fits the ACI curve using FvCB model with the option to include
#' the TPU limiting curve. Currently not working, merely a placeholder.
#' @param data A dataframe containing all the relevant parameters.
#' @param TPU Boolean indicating the inclusion of the TPU limiting curve
#' (default is FALSE)
aci_curve <- function(data,
                      include_TPU_curve = FALSE,
                      C_c_denom_multiplier = 4,
                      G_star_denom_multiplier = 8,
                      V_cmax_var = V_cmax,
                      C_c_var = C_c,
                      G_star_var = G_star,
                      K_c_var = K_c,
                      O_var = O,
                      K_o_var = K_o,
                      R_d_var = R_d,
                      J_var = J,
                      TPU_var = TPU) {

  # This isn't right.

  # Need to check if the inputs are numbers maube?
  # As in, J = 4 rather than a column
  rubisco_params <- enquos(
    V_cmax_var,
    C_c_var,
    G_star_var,
    K_c_var,
    O_var,
    K_o_var,
    R_d_var
  )

  rubp_params <- enquos(J_var, C_c_var, G_star_var, R_d_var)

  data %<>%
    mutate(
      Rubisco_curve = Rubisco_limited(!!!rubisco_params),
      RuBP_curve = RuBP_limited(!!!rubp_params)
    )

  # rubisco_curve <- Rubisco_limited(!!! rubisco_params)
  # rubp_curve <- RuBP_limited(!!! rubp_params)

  if (include_TPU_curve) {
    tpu_params <- enquos(TPU_var, R_d_var)

    data %<>%
      mutate(TPU_curve = TPU_limited(!!!tpu_params))
  } else {
    data %<>%
      mutate(TPU_curve = rep(Inf, nrow(data)))
  }

  data %<>%
    mutate(fitted_values = pmin(Rubisco_curve, RuBP_curve, TPU_curve))

  # fitted_values <- pmin(rubisco_curve, rubp_curve, tpu_curve)

  # return(fitted_values)
}

# Pseudo code for fitting aci curve in nlme
#' @title ACI Formula
#' @description Fits the ACI curve using FvCB model with the option to include
#' the TPU limiting curve. Currently not working, merely a placeholder.
#' @param TPU Boolean indicating the inclusion of the TPU limiting curve
#' (default is FALSE)
aci_formula <- function(V_cmax, C_c, G_star, K_c, O, K_o, R_d, J,
                        TPU_var = NA,
                        C_c_denom_multiplier = 4,
                        G_star_denom_multiplier = 8
                        ) {

  Rubisco_curve <- Rubisco_limited(V_cmax, C_c, G_star, K_c, O, K_o, R_d)
  RuBP_curve <- RuBP_limited(J, C_c, G_star, R_d,
                             C_c_denom_multiplier = C_c_denom_multiplier,
                             G_star_denom_multiplier = G_star_denom_multiplier
                             )

  if (! is.na(TPU_var)) {
    TPU_curve <- TPU_limited(TPU, R_d)

  } else {
    TPU_curve <- rep(Inf, nrow(data))
  }

  fitted_values <- pmin(Rubisco_curve, RuBP_curve, TPU_curve)

}




# This is Gerrit's model pretty much, placeholder
# Consider scaling option (Bates and Pinheiro reccomend it)
# Add options for covariance matrix, other nlme parameters (maybe as ...)
#' @title Non-Linear Mixed Effects ACI model
#' @description Fits a non-linear mixed effects model using \code{nlme} package.
#' Currently not working, decent placeholder for final version.
#' @param data The dataframe containing photosynthesis data
#' @return An nlme model
nlme_aci <- function(data, fixed, random,
                     V_cmax, C_c, G_star, K_c, O, K_o, R_d, J,
                     TPU_var = NA,
                     C_c_denom_multiplier = 4,
                     G_star_denom_multiplier = 8,
                     starting_values = NA,
                     cov_matrix = pdSymm,
                     C_c_denom_multiplier = 4,
                     G_star_denom_multiplier = 8,
                     ...) {

  if (is.na(starting_values)){
    # Something to find starting values: ols?
    # some sub-function
    y <- T # to avoid throwing an error before I add this
  }


  model <- nlme(y ~ aci_formula(V_cmax, C_c, G_star, K_c, O, K_o, R_d, J,
                                TPU_var = NA,
                                C_c_denom_multiplier = C_c_denom_multiplier,
                                G_star_denom_multiplier = G_star_denom_multiplier
                                ),
    data = data,
    fixed = fixed,
    random = cov_matrix(random),
    start = starting_values,
    ...
  )
  return(model)
}


find_optimal_start <- function(data, fixed, random,
                               V_cmax, C_c, G_star, K_c, O, K_o, R_d, J,
                               TPU_var = NA,
                               C_c_denom_multiplier = 4,
                               G_star_denom_multiplier = 8,
                               starting_values = NA,
                               cov_matrix = pdSymm,
                               ...){
  num_starts <-length(random)

  if (is.na(starting_values)){
    starting_values = list(V_cmax = 100, R_d = 2.0, J = 100)
  }

  # use the nls model as a starting point
  nls_model <- nls(y ~ FvCB(x, PAR, c1, c2, Vcmax, Jmax, alpha, theta, Rd),
      data = CO2light,
      subset = select.high,
      start = starting_values
  )


  coef <- summary(nls_model)$coefficients
  coef_names <- rownames(coef)
  coef <- as.tibble(coef) %>%
    add_column(Coefficient = coef_names) %>%
    dplyr::select(Coefficient, Estimate)

  # use this as the starting values for for NLME?


}
