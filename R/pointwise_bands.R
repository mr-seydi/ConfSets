#' @title Pointwise confidence bands for functional intraclass correlation coefficient, ICC, in Test-Retest Reliability based on F distribution
#'
#' @description This function calculates the pointwise confidence bands based on F distribution for functional ICC in Test-Retest setting. Note that pointwise F bands do not provide simultaneous coverage.
#'
#'@param test A matrix of test measurements with rows representing individuals and columns representing evaluation points.
#'@param retest A matrix of retest measurements with rows representing individuals and columns representing evaluation points.
#'@param individuals A vector of individual identifiers. It must be the combined vector of individual identifiers for test and retest. If NULL, the function assumes that the individuals in the test and retest matrices are the same.
#'@param designation A character string specifying the designation of the ICC. Possible values are "ICC(A,1)" and "ICC(A,k)" Possible designation is "ICC(A,1)", "ICC(A,k)" based on the paper: McGraw KO, Wong SP. Forming inferences about some intraclass correlation coefficients. Psychological methods. 1996 Mar;1(1):30. The default value is "ICC(A,1)".
#'@param alpha A numeric value specifying the significance level for the ICC calculation. The default value is 0.05.
#'
#'@returns A list containing the pointwise ICC values, integrated ICC value, lower and upper pointwise_F bands.
#'
#'@examples
#'# The package contains the datasets, TEST_data, and RETEST_data which are the
#'#test and retest data, respectively.
#'
#' output_balanced_Fbands <- pointwise_F_bands_ICC(test = TEST_data, retest = RETEST_data,
#' individuals = NULL, designation = "ICC(A,1)", alpha = 0.05)
#'
#' # The output_balanced_ICCA1 contains the pointwise and integrated ICC values.
#' output_balanced_Fbands$pointwise_ICC
#' output_balanced_Fbands$integrated_ICC
#' output_balanced_Fbands$lower_F
#' output_balanced_Fbands$upper_F
#'
#' #Consider the scenario with unbalanced case, in total we have 28 individuals,
#' #and for example individual 2 has 1 test and 2 retest trials, individual 3 has
#' #3 test and 2 retest trials, and the other individuals have one trial per session.
#' #Then the individuals vector should be like
#' # c(1:2, rep(3, 3), 4:28, 1, rep(2,2), rep(3, 2), 4:28).
#'
#' output_unbalanced_Fbands <- pointwise_F_bands_ICC(test = TEST_data, retest = RETEST_data,
#' individuals = c(1:2, rep(3, 3), 4:28, 1, rep(2,2), rep(3, 2), 4:28),
#'  designation = "ICC(A,1)", alpha = 0.05)
#'
#' output_unbalanced_Fbands$pointwise_ICC
#' output_unbalanced_Fbands$integrated_ICC
#' output_unbalanced_Fbands$lower_F
#' output_unbalanced_Fbands$upper_F
#'
#'
#'@export
pointwise_F_bands_ICC <- function(test, retest, individuals = NULL, designation = "ICC(A,1)", alpha = 0.05) {


  # Possible values for alpha
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("The significance level must be a number between 0 and 1.")
  }

  # possible designation is "ICC(A,1)", "ICC(A,k)""
  # Check if the designation is valid
  if (!(designation %in% c("ICC(A,1)", "ICC(A,k)"))) {
    stop("Invalid designation. Please specify 'ICC(A,1)' or 'ICC(A,k)'.")
  }

  # Check if the dimensions of the test and retest matrices are the same
  if (!identical(dim(test)[2], dim(retest)[2])) {
    stop("The Number of evaluation points of the test and retest matrices must be the same.")
  }

  if(!is.matrix(test) || !is.matrix(retest)){
    stop("The input data must be a matrix.")
  }

  if(designation == "ICC(A,k)" && dim(test)[1] != dim(retest)[1]){
    stop("The number of individuals in the test and retest matrices must be the same for ICC(A,k).")
  }

  if (is.null(individuals)) {
    individual_test_retest = c(1:dim(test)[1],1:dim(retest)[1])
  } else if ( !(is.null(individuals)) && designation == "ICC(A,k)") {
    stop("The 'individuals' argument is only applicable for ICC(A,1).")
  } else {
    individual_test_retest = individuals
  }


  data = rbind(test, retest)
  t = dim(data)[2] # Number of evaluation points
  N = dim(data)[1] # Number of observed curves in total

  n = length(levels(factor(individual_test_retest))) # Number of individuals

  test_retest_index = factor(c(rep(1, dim(test)[1]), rep(2, dim(retest)[1])))
  # Initialize variables
  nik = table(individual_test_retest, test_retest_index)
  Ni = rowSums(nik)
  Nk = colSums(nik)
  K = length(Nk)
  ks = (N - sum(Ni^2) / N) / (n - 1)
  kt = (N - sum(Nk^2) / N) / (K - 1)


  # Initialize vectors for pointwise ICCs
  pointwise_ICC <- numeric(t)

  # Initialize vectors for lower and upper F bands
  low_F = up_F = numeric(t)

  # Loop for each evaluation point to calculate pointwise ICC
  for (point in 1:t) {
    data.all = data.frame(
      response = c(test[, point], retest[, point]),
      individual = factor(individual_test_retest),
      measurement = test_retest_index
    )

    anova <- stats::aov(response ~ individual + measurement, data = data.all)
    MS <- summary(anova)[[1]]$'Mean Sq'
    MSs <- MS[1]
    MSt <- MS[2]
    MSe <- MS[3]


    # Your formula for pointwise ICC
    if (designation == "ICC(A,1)") {
      # This formula based on paper:
      #Schelin L, Pini A, Markström JL, Häger CK. Test-retest reliability of
      #entire time-series data from hip, knee and ankle kinematics and kinetics
      #during one-leg hops for distance: Analyses using integrated pointwise indices.
      #Journal of biomechanics. 2021 Jul 19;124:110546.
      # Which can be applied to unbalance cases
      pointwise_ICC[point] <- (kt*(MSs-MSe)) / (ks*MSt + kt*MSs + (ks*kt-ks-kt)*MSe)

      a = ks*pointwise_ICC[point]/(kt*(1-pointwise_ICC[point]))
      b = 1 + ks*pointwise_ICC[point]*(kt-1)/(kt*(1-pointwise_ICC[point]))
      nu = (a*MSt+b*MSe)^2 / ((a*MSt)^2/(K-1) + (b*MSe)^2/(N-n-K+1 ))
      Fstar1 = stats::qf(1-alpha/2,n-1,nu)
      Fstar2 = stats::qf(1-alpha/2,nu,n-1)
      low_F[point] = kt*(MSs-Fstar1*MSe)/(kt*MSs + Fstar1*(K*MSt+(kt*ks-kt-ks)*MSe))
      up_F[point] = kt*(Fstar2*MSs-MSe)/((ks*MSt+(kt*ks-kt-ks)*MSe) + Fstar2*kt*MSs)

    } else if (designation == "ICC(A,k)") {

      pointwise_ICC[point] <- (MSs - MSe) / (MSs + (MSt - MSe) / kt)

      a = ks*pointwise_ICC[point]/(kt*(1-pointwise_ICC[point]))
      b = 1 + ks*pointwise_ICC[point]*(kt-1)/(kt*(1-pointwise_ICC[point]))
      nu = (a*MSt+b*MSe)^2 / ((a*MSt)^2/(K-1) + (b*MSe)^2/(N-n-K+1 ))
      Fstar1 = stats::qf(1-alpha/2,n-1,nu)
      Fstar2 = stats::qf(1-alpha/2,nu,n-1)
      low_F[point] = kt*(MSs-Fstar1*MSe)/(kt*MSs + Fstar1*(MSt-MSe))
      up_F[point] = kt*(Fstar2*MSs-MSe)/(MSt- MSe + Fstar2*kt*MSs)

    }


  }


  # Integrated ICC as the mean of pointwise ICCs
  integrated_ICC <- mean(pointwise_ICC)

  # Return the pointwise and integrated ICCs
  return(list(
    pointwise_ICC = pointwise_ICC,
    integrated_ICC = integrated_ICC,
    lower_F = low_F,
    upper_F = up_F
  ))

}
