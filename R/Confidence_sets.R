#' @title Confidence sets for functional intraclass correlation coefficient, ICC, in Test-Retest Reliability
#'
#' @description This function calculates the confidence sets for functional ICC in Test-Retest setting. All Confidence sets explained in Seydi MR, Pini A, Pataky TC, Schelin L. Confidence sets for intraclass correlation coefficients in test–retest curve measurements. Journal of Biomechanics. 2024 Aug 1;173:112232. Note that pointwise F bands do not provide simultaneous coverage. Furthermore, this confidence bands have not investigated for unbalanced cases so far.
#'
#'@param test A matrix of test measurements with rows representing individuals and columns representing evaluation points.
#'@param retest A matrix of retest measurements with rows representing individuals and columns representing evaluation points.
#'@param individuals A vector of individual identifiers. It must be the combined vector of individual identifiers for test and retest. If NULL, the function assumes that the individuals in the test and retest matrices are the same.
#'@param type A character vector specifying the type of confidence sets to be calculated. Possible values are "Rank" for ICC curve confidence bands, and "Normal" for integrated ICC confidence interval. The user can choose one or both. The default value is both "Rank" and "Normal".
#'@param n_boot A numeric value specifying the number of bootstrap samples. The default value is 1000.
#'@param designation A character string specifying the designation of the ICC. Possible values are "ICC(A,1)" and "ICC(A,k)". Possible designation is "ICC(A,1)", "ICC(A,k)" based on the paper: McGraw KO, Wong SP. Forming inferences about some intraclass correlation coefficients. Psychological methods. 1996 Mar;1(1):30. The default value is "ICC(A,1)".
#'@param alpha A numeric value specifying the significance level for the ICC calculation. The default value is 0.05.
#'
#'@returns A list containing the pointwise ICC values, integrated ICC value, and the confidence sets.
#'
#'
#'@examples
#'# The package contains the datasets, TEST_data, and RETEST_data which are the
#'#test and retest data, respectively.
#'
#'
#' set.seed(123)
#' confsets_outputs <- Conf_sets(test = TEST_data, retest = RETEST_data,
#'                              type = c("Rank", "Normal"), n_boot=1000,
#'                              designation = "ICC(A,1)", alpha = 0.05)
#'
#' confsets_outputs$pointwise_ICC
#' confsets_outputs$integrated_ICC
#' confsets_outputs$CSets
#' confsets_outputs$CSets$Rank
#' confsets_outputs$CSets$Normal
#'
#' #Consider the scenario with unbalanced case, in total we have 28 individuals,
#' #and for example the individuals number 2 has attended 1 test and 2 retest trials,
#' #number 3 has attended 3 test and 2 retest trials,
#' #the other individuals have attended only one trials.
#' #Then the individuals vector should be like
#' # c(1:2, rep(3, 3), 4:28, 1, rep(2,2), rep(3, 2), 4:28).
#'
#' set.seed(1234)
#' confsets_outputs <- Conf_sets(test = TEST_data, retest = RETEST_data,
#'                              individuals = c(1:2, rep(3, 3), 4:28, 1, rep(2,2),
#'                              rep(3, 2), 4:28),
#'                              type = c("Rank", "Normal"), n_boot=1000,
#'                              designation = "ICC(A,1)", alpha = 0.05)
#'
#' confsets_outputs$pointwise_ICC
#' confsets_outputs$integrated_ICC
#' confsets_outputs$CSets
#' confsets_outputs$CSets$Rank
#' confsets_outputs$CSets$Normal
#'
#'@export
Conf_sets <- function(test, retest, individuals = NULL, type = c("Rank", "Normal"),
                      n_boot=1000, designation = "ICC(A,1)", alpha = 0.05){

  # Valid types for simultaneous bands
  valid_types <- c("Rank", "Normal")

  # Check if all specified types are valid
  if (!all(type %in% valid_types)) {
    stop(paste("Invalid type specified. Choose from:", paste(valid_types, collapse=", ")))
  }

  # Possible values for n_boot
  if (!is.numeric(n_boot) || n_boot <= 0) {
    stop("The number of bootstrap samples must be a positive integer.")
  }

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


  # # FOR now it must be balanced case
  # if(dim(test)[1] != dim(retest)[1]){
  #   stop("The number of individuals in the test and retest matrices must be the same.")
  # }


  if (is.null(individuals)) {
    individual_test_retest = c(1:dim(test)[1],1:dim(retest)[1])
  } else if ( !(is.null(individuals)) && designation == "ICC(A,k)") {
    stop("The 'individuals' argument is only applicable to ICC(A,1).")
  } else if (length(individuals) != dim(test)[1] + dim(retest)[1]) {
    stop("'individuals' must have length nrow(test) + nrow(retest).")
  } else {
    individual_test_retest = individuals
  }


  E_ICC <- functional_ICC(test=test, retest=retest, individuals = individual_test_retest, designation = designation, alpha = alpha)
  E_integrated_ICC <- E_ICC$integrated_ICC
  E_pointwise_ICC <- E_ICC$pointwise_ICC


  t = dim(test)[2] # Number of evaluation points

  # # Since we are in balance case scenarios 1:dim(test)[1] and 1:dim(retest)[1] are the same
  # n = length(levels(factor(individual_test_retest))) # Number of individuals
  # boot_index = boot_function_index(data=1:n, number_of_resample = n_boot)


  boot_index <- bootstrap_indices(test=test, retest=retest, individuals = individual_test_retest, n_boot = n_boot)


  ICC.Pointwise.boot = matrix(NA,t,n_boot)
  ICC.boot=numeric()

  pb <- progress::progress_bar$new(total = n_boot)
  pb$tick(0)

  for (i in 1:n_boot) {

    pb$tick()

    m1 = test[boot_index[[i]]$test_idx, , drop = FALSE]
    m2 = retest[boot_index[[i]]$retest_idx, , drop = FALSE]
    ICC_out <- functional_ICC(test=m1, retest=m2, individuals = boot_index[[i]]$individuals, designation = designation, alpha = alpha)
    ICC.Pointwise.boot[,i] = ICC_out$pointwise_ICC
    ICC.boot[i] = ICC_out$integrated_ICC
  }


  # Create curve set for simultaneous bands
  curve_set_input <- list(r=1:t, obs=ICC.Pointwise.boot)
  boot_set <- GET::create_curve_set(curve_set_input)

  # Generate the requested bands
  CSets <- list()

  if ("Rank" %in% type) {
    rank_env <- GET::central_region(boot_set, coverage=1-alpha, type="rank", central="mean")
    CSets$Rank <- list(upper=rank_env$hi, lower=rank_env$lo)
  }

  # if ("ERL" %in% type) {
  #   erl_env <- GET::central_region(boot_set, coverage=1-alpha, type="erl", central="mean")
  #   CSets$ERL <- list(upper=erl_env$hi, lower=erl_env$lo)
  # }

  if ("Normal" %in% type) {
    se_integrated_ICC <- stats::sd(ICC.boot)
    z_value <- stats::qnorm(1 - alpha / 2)
    CSets$Normal <- c(E_integrated_ICC - z_value * se_integrated_ICC,
                      min(E_integrated_ICC + z_value * se_integrated_ICC, 1))
  }

  # if ("Pointwise_F" %in% type) {
  #   CSets$Pointwise_F <- list(upper=E_ICC$upper_F, lower=E_ICC$lower_F)
  # }

  # Return results as a list
  return(list(
    pointwise_ICC = E_pointwise_ICC,
    integrated_ICC = E_integrated_ICC,
    CSets = CSets
  ))


}



