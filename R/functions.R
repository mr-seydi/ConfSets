# This function is used to find the indices for the bootstrapping
# boot_function_index <- function(data,number_of_resample){
#   matrix_data=matrix(NA,length(data),number_of_resample)
#   for (i in 1:number_of_resample) {
#     data_sample=base::sample(data,size = length(data),replace = T)
#     matrix_data[,i]=data_sample
#   }
#   return(matrix_data)
# }


bootstrap_indices <- function(test, retest, individuals, n_boot) {

  n_test   <- nrow(test)
  n_retest <- nrow(retest)

  ids_test   <- individuals[seq_len(n_test)]
  ids_retest <- individuals[n_test + seq_len(n_retest)]

  # Row indices of each individual's trials, per session (looked up by id)
  rows_test   <- split(seq_len(n_test),   as.character(ids_test))
  rows_retest <- split(seq_len(n_retest), as.character(ids_retest))

  unique_ids <- as.character(unique(individuals))
  n_ids <- length(unique_ids)

  # Resample the elements of x with replacement (safe when length(x) == 1)
  resample <- function(x) x[sample.int(length(x), length(x), replace = TRUE)]

  lapply(seq_len(n_boot), function(b) {

    # Stage 1: draw individuals with replacement
    boot_ids <- sample(unique_ids, size = n_ids, replace = TRUE)

    # Stage 2: for each draw, resample trials within each session
    t_idx <- lapply(boot_ids, function(g) resample(rows_test[[g]]))
    r_idx <- lapply(boot_ids, function(g) resample(rows_retest[[g]]))

    # New id = draw position i, shared by the test and retest rows of that draw
    new_test_ids   <- rep(seq_along(t_idx), lengths(t_idx))
    new_retest_ids <- rep(seq_along(r_idx), lengths(r_idx))

    list(
      test_idx    = unlist(t_idx),
      retest_idx  = unlist(r_idx),
      individuals = c(new_test_ids, new_retest_ids)

    )
  })
}




##############################################

# # Set seed for reproducibility
# set.seed(123)
#
# # ---------------------------
# # Create the Original Data Sets
# # ---------------------------
# # Test session: individuals and their trial observations
# test_ids <- c(1:2, rep(3, 4), 4:28)
# # Retest session: individuals and their trial observations
# retest_ids <- c(1:2, rep(3, 2), 4:28)
#
# # Create a data frame for the test session.
# # We use dplyr to number trials within each individual.
# library(dplyr)
# test_df <- data.frame(id = test_ids, session = "test") %>%
#   group_by(id) %>%
#   mutate(trial = row_number()) %>%
#   ungroup()
#
# # Create a data frame for the retest session.
# retest_df <- data.frame(id = retest_ids, session = "retest") %>%
#   group_by(id) %>%
#   mutate(trial = row_number()) %>%
#   ungroup()
#
# # For demonstration, let’s simulate a measurement variable
# # (e.g., some outcome measured at each trial)
# # Here we simply use a random normal value.
# test_df <- test_df %>% mutate(value = rnorm(n()))
# retest_df <- retest_df %>% mutate(value = rnorm(n()))
#
# # Check the structure (optional)
# print(test_df)
# print(retest_df)
#
# # ---------------------------
# # Nested (Two-Stage) Bootstrap
# # ---------------------------
# # Step 1: Resample individuals (clusters) with replacement.
# # Get the unique individual IDs from the data.
# unique_ids <- sort(unique(c(test_df$id, retest_df$id)))  # should be 1:28
# n_ids <- length(unique_ids)  # number of unique individuals (28 in our case)
#
# # For one bootstrap replicate, sample individuals with replacement.
# boot_ids <- sample(unique_ids, size = n_ids, replace = TRUE)
# cat("Bootstrapped individual IDs:\n")
# print(boot_ids)
#

#This code will fail when an individual ID is selected more than once

# # Step 2: For each selected individual, resample trials within each session.
# # We will do this separately for test and retest sessions.
#
# # Bootstrap the test session
# boot_test <- do.call(rbind, lapply(boot_ids, function(this_id) {
#   # Extract all test trials for this individual.
#   dt <- test_df %>% filter(id == this_id)
#   n_trials <- nrow(dt)
#   # Resample the trials with replacement. If only one trial, it will be re-selected.
#   dt_boot <- dt[sample(seq_len(n_trials), size = n_trials, replace = TRUE), ]
#   # Optionally, add a column indicating the bootstrapped cluster replicate.
#   dt_boot$boot_id <- this_id
#   return(dt_boot)
# }))
#
# # Bootstrap the retest session
# boot_retest <- do.call(rbind, lapply(boot_ids, function(this_id) {
#   dt <- retest_df %>% filter(id == this_id)
#   n_trials <- nrow(dt)
#   dt_boot <- dt[sample(seq_len(n_trials), size = n_trials, replace = TRUE), ]
#   dt_boot$boot_id <- this_id
#   return(dt_boot)
# }))
#
# # Display the bootstrapped data (first few rows)
# cat("Bootstrapped Test Data (first few rows):\n")
# print(head(boot_test))
# cat("Bootstrapped Retest Data (first few rows):\n")
# print(head(boot_retest))
#
# # ---------------------------
# # Next Steps
# # ---------------------------
# # At this point, you have two bootstrapped data sets (boot_test and boot_retest)
# # that respect both the clustering (individual-level) and the trial-level variability.
# # You would then combine or analyze these bootstrapped data sets with the appropriate
# # ICC estimator for multiple trials per session.
#
# # For example, if you have a function compute_icc(test_data, retest_data) that calculates
# # the ICC for the given sessions, you could call:
# # icc_boot <- compute_icc(boot_test, boot_retest)
# #
# # To perform a full bootstrap, you would wrap the above procedure in a loop (or use replicate())
# # and store the ICC estimates over many replicates (e.g., 1000 times) to form a distribution,
# # from which you can calculate confidence intervals.
#
# # Here’s a skeleton for performing many bootstrap replicates:
# n_boot <- 1000
# icc_estimates <- numeric(n_boot)
# for (i in 1:n_boot) {
#   boot_ids <- sample(unique_ids, size = n_ids, replace = TRUE)
#
#   boot_test <- do.call(rbind, lapply(boot_ids, function(this_id) {
#     dt <- test_df %>% filter(id == this_id)
#     dt[sample(seq_len(nrow(dt)), size = nrow(dt), replace = TRUE), ]
#   }))
#
#   boot_retest <- do.call(rbind, lapply(boot_ids, function(this_id) {
#     dt <- retest_df %>% filter(id == this_id)
#     dt[sample(seq_len(nrow(dt)), size = nrow(dt), replace = TRUE), ]
#   }))
#
#   # Here, you would compute the ICC based on boot_test and boot_retest.
#   # For demonstration, we just compute a dummy statistic (mean difference between sessions).
#   # Replace this with your ICC estimation function.
#   test_mean <- mean(boot_test$value)
#   retest_mean <- mean(boot_retest$value)
#   icc_estimates[i] <- test_mean - retest_mean
# }
#
# # Now, you can use the distribution of icc_estimates to calculate confidence intervals, etc.
# cat("Bootstrap estimates (first 10 replicates):\n")
# print(head(icc_estimates, 10))



# ####Corected########
# n_boot <- 1000
# icc_estimates <- numeric(n_boot)
#
# for (b in 1:n_boot) {
#   # 1. Sample the individual IDs with replacement
#   boot_ids <- sample(unique_ids, size = n_ids, replace = TRUE)
#
#   # 2. Build the bootstrapped sessions simultaneously to keep them aligned
#   boot_test_list <- list()
#   boot_retest_list <- list()
#
#   for (i in seq_along(boot_ids)) {
#     this_id <- boot_ids[i]
#
#     # --- Process Test Session ---
#     dt_t <- test_df %>% filter(id == this_id)
#     n_t <- nrow(dt_t)
#     boot_t <- dt_t[sample(seq_len(n_t), size = n_t, replace = TRUE), ]
#     boot_t$id <- i  # OVERWRITE with unique clone ID (1 to 28)
#     boot_test_list[[i]] <- boot_t
#
#     # --- Process Retest Session ---
#     dt_r <- retest_df %>% filter(id == this_id)
#     n_r <- nrow(dt_r)
#     boot_r <- dt_r[sample(seq_len(n_r), size = n_r, replace = TRUE), ]
#     boot_r$id <- i  # OVERWRITE with the exact same unique clone ID
#     boot_retest_list[[i]] <- boot_r
#   }
#
#   # Combine individual data frames back together
#   boot_test <- do.call(rbind, boot_test_list)
#   boot_retest <- do.call(rbind, boot_retest_list)
#
#   # 3. Compute your ACTUAL ICC here using the corrected 'id' column
#   # icc_estimates[b] <- your_icc_function(boot_test, boot_retest)
#
#   # For demonstration, the dummy statistic still works exactly the same:
#   icc_estimates[b] <- mean(boot_test$value) - mean(boot_retest$value)
# }
