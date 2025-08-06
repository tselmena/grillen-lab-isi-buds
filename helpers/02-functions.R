# This script contains functions to process the data, 
# calculate sample sizes & aucs, and run uni-test logistic regression
# 1) get_data_frames()
# 2) run_weekly_models() 
# 3) get_data_frames_learn()
# 4) get_sample_size()
# 5) compute_auc()
# 6) get_w48_aucs

### 1) get_data_frames() #######################################################
## Description: Creates a df corresponding to a cognitive test;
## Input: 
# - data
# - test_code 
# - target_weeks_w_0 
# - window_weeks, baseline_window_weeks
# - outcome_data, outcome_var
## Output: 
# - BID
# - CDPOS_W240 (response)
# - Baseline score and each followup score at the 5 target weeks
# - Raw delta values
# - Z-score delta values

get_data_frames <- function(data, test_code, target_weeks_w_0 , window_weeks, baseline_window_weeks, outcome_data, outcome_var) {
  # test_data is ADQS_raw; utilizes given week variable, as opposed to manually 
  # calculating it as done with CDR and C3
  test_data <- data |>
    filter(toupper(QSTESTCD) == test_code, TX == "Placebo") |>
    select(BID, WEEK = ADURW, SCORE = QSSTRESN) |>
    # remove any NA results
    filter(!is.na(WEEK), !is.na(SCORE))
  
  if (nrow(test_data) == 0) {
    # return NULL if no data for this test code
    return(NULL)
  }
  
  closest <- test_data |>
    # Create a temporary list-column with the absolute difference to each target week
    mutate(
      tmp = map(WEEK, ~ abs(.x - target_weeks_w_0 )),
      # Find the target week with the minimum difference
      target_week = map_int(tmp, ~ target_weeks_w_0 [which.min(.x)])
    ) |>
    select(-tmp) |>
    
    # Filter to keep only visits within the specified window of the assigned target week
    ### Uses 14 weeks before baseline to capture data for ADL and CFI ### 
    filter(abs(WEEK - target_week) <= if_else(target_week <= 0, baseline_window_weeks, window_weeks)) |>
    
    group_by(BID, target_week) |>
    # For each subject and target week, keep the visit with the smallest time difference
    slice_min(order_by = abs(WEEK - target_week), n = 1, with_ties = FALSE) |>
    ungroup()
  
  baseline_scores <- closest |>
    filter(target_week <= 0) |>
    group_by(BID) |>
    # If multiple baseline candidates, take the one with the latest date
    slice_max(order_by = target_week, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(BID, baseline_score = SCORE)
  
  wide_followup_scores <- closest |>
    filter(target_week > 0) |>
    select(BID, target_week, SCORE) |>
    pivot_wider(
      id_cols = BID,
      names_from = target_week,
      values_from = SCORE,
      names_prefix = paste0(test_code, "_W")
    )
  
  # Check if there are any followup scores to process
  if (ncol(wide_followup_scores) <= 1) {
    return(NULL)
  }
  
  wide_scores_with_deltas <- wide_followup_scores |>
    inner_join(baseline_scores, by = "BID") |>
    mutate(across(
      .cols = starts_with(paste0(test_code, "_W")),
      .fns = ~ .x - baseline_score,
      .names = "delta_{.col}"
    ))
  
  model_data <- wide_scores_with_deltas |>
    inner_join(
      outcome_data |> select(BID, all_of(outcome_var)),
      by = "BID"
    ) |>
    filter(!is.na(.data[[outcome_var]]))
  
  # standardize deltas for comparison (z-scores)
  model_data <- model_data |>
    mutate(across(
      .cols = starts_with("delta_"),
      .fns = ~ as.numeric(scale(.x)),
      .names = "z_{.col}"
    ))
  
  return(model_data)
}

### 2) run_weekly_models() #####################################################

## Description: Runs uni-test logistic model for all 12 tests at a specified week, 
## option to include baseline as covariate

## Input: 
# - data_list,
# - week_num 
# - outcome_var (default = "CDPOS_W240") 
# - include_baseline (default = FALSE)

## Output: 12 x 11 dataframe (12 tests, 11 columns)
# - test 
# - odds_ratio, 
# - SE, test statistic, p.value, CI lower, CI higher 
# - n

run_weekly_models <- function(data_list, week_num, outcome_var = "CDPOS_W240", include_baseline = FALSE) {
  # define the search pattern for the week-specific column
  week_pattern <- paste0("_W", week_num, "$")
  map_dfr(data_list, ~{
    # find the predictor column for the specified week
    z_delta_col <- str_subset(names(.x), paste0("^z_delta_.*", week_pattern))
    # conditionally select data and build the model formula
    if (include_baseline) {
      # with baseline
      model_data <- .x |>
        ungroup() |>
        select(all_of(outcome_var), baseline_score, all_of(z_delta_col)) |>
        na.omit()
      formula <- as.formula(paste(outcome_var, "~", z_delta_col, "+ baseline_score"))
    } else {
      # without baseline
      model_data <- .x |>
        ungroup() |>
        select(all_of(outcome_var), all_of(z_delta_col)) |>
        na.omit()
      formula <- as.formula(paste(outcome_var, "~", z_delta_col))
    }
    
    model <- glm(formula, data = model_data, family = binomial())
    
    # tidy the results 
    tidy_results <- tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term != "(Intercept)")
    
    test_code <- str_extract(z_delta_col, "(?<=z_delta_).*(?=_W)")
    
    # baseline was included, rename its term to be specific
    if (include_baseline && "baseline_score" %in% tidy_results$term) {
      tidy_results <- tidy_results |>
        mutate(term = if_else(
          term == "baseline_score",
          paste0("baseline_score_", test_code),
          term
        )
        )
    }
    # add observation count and rename odds ratio column
    tidy_results |>
      mutate(n_obs = nobs(model)) |>
      rename(odds_ratio = estimate) |> 
      mutate(baseline = ifelse(include_baseline, 1, 0), 
             test_code = test_code)
  }, .id = "test")
}

### 3) get_data_frames_learn() #################################################

## Description: Modified get_data_frames() function to work on LEARN cohort, 
# largely same utility, but no outcome variable

## Input: Returns a df corresponding to a cognitive test
# - data
# - test_code 
# - target_weeks_w_0 
# - window_weeks, baseline_window_weeks

## Output: 
# - BID
# - Baseline score and each followup score at the 5 target weeks
# - Raw delta values
# - Z-score delta values

get_data_frames_learn <- function(data, test_code, target_weeks_w_0 , window_weeks, baseline_window_weeks) {
  test_data <- data |>
    filter(toupper(QSTESTCD) == test_code) |>
    select(BID, WEEK = ADURW, SCORE = QSSTRESN) |>
    filter(!is.na(WEEK), !is.na(SCORE))
  
  if (nrow(test_data) == 0) {
    # Return NULL if no data for this test code
    return(NULL)
  }
  
  closest <- test_data |>
    # Create a temporary list-column with the absolute difference to each target week
    mutate(
      tmp = map(WEEK, ~ abs(.x - target_weeks_w_0 )),
      # Find the target week with the minimum difference
      target_week = map_int(tmp, ~ target_weeks_w_0 [which.min(.x)])
    ) |>
    select(-tmp) |>
    # Filter to keep only visits within the specified window of the assigned target week
    filter(abs(WEEK - target_week) <= if_else(target_week <= 0, baseline_window_weeks, window_weeks)) |>
    group_by(BID, target_week) |>
    # For each subject and target week, keep the visit with the smallest time difference
    slice_min(order_by = abs(WEEK - target_week), n = 1, with_ties = FALSE) |>
    ungroup()
  
  baseline_scores <- closest |>
    filter(target_week <= 0) |>
    group_by(BID) |>
    # If multiple baseline candidates, take the one with the latest date
    slice_max(order_by = target_week, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(BID, baseline_score = SCORE)
  
  wide_followup_scores <- closest |>
    filter(target_week > 0) |>
    select(BID, target_week, SCORE) |>
    pivot_wider(
      id_cols = BID,
      names_from = target_week,
      values_from = SCORE,
      names_prefix = paste0(test_code, "_W")
    )
  
  wide_scores_with_deltas <- wide_followup_scores |>
    inner_join(baseline_scores, by = "BID") |>
    mutate(across(
      .cols = starts_with(paste0(test_code, "_W")),
      .fns = ~ .x - baseline_score,
      .names = "delta_{.col}"
    ))
  
  model_data <- wide_scores_with_deltas |>
    mutate(across(
      .cols = starts_with("delta_"),
      .fns = ~ as.numeric(scale(.x)),
      .names = "z_{.col}"
    ))
  
  return(model_data)
}

### 4) get_sample_size() #######################################################

# Description: Calculates the required sample size needed using specified power, 
# significance level, variance, and mean difference

# Input: A df containing variance and mean, and statistical parameters
# - data: A df that must contain 'var' and 'mean' (delta) columns.
# - alpha
# - beta

# Output: A tibble with the required sample size for each test
# - test_name
# - n: rounded up to the nearest integer
get_sample_size <- function(data, alpha = 0.05, beta = 0.2) {
  sigma_squared <- data$var
  delta <- data$mean
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(1 - beta)
  n <- ( (z_alpha + z_beta)^2 * 2 * sigma_squared ) / delta^2
  sample_size <- 
    tibble(
      test_name = data$test_name, 
      n = ceiling(n)
    ) |> 
    arrange(n)
  sample_size
}

### 5) compute_auc() ###########################################################

## Description: Calculates the AUC given true outcomes and predicted scores (Ian's code)

## Input: vectors of true and predicted values
# - y_true: numeric vector of true binary outcomes (must be 0/1)
# - y_score: numeric vector of predicted scores or probabilities
# - n_cutpoints: number of thresholds to use for approximating the ROC curve

## Output:
# - a single numeric value representing the calculated AUC

compute_auc <- function(y_true, y_score, n_cutpoints = 100) {
  if (!all(y_true %in% c(0, 1))) stop("y_true must be binary (0/1)")
  
  cutpt <- seq(max(y_score), min(y_score), length.out = n_cutpoints)
  sensitivity <- specificity <- numeric(length(cutpt))
  
  for (i in seq_along(cutpt)) {
    prog <- y_score[y_true == 1]
    noprog <- y_score[y_true == 0]
    
    sensitivity[i] <- sum(prog >= cutpt[i]) / length(prog)
    specificity[i] <- sum(noprog < cutpt[i]) / length(noprog)
  }
  
  x <- 1 - specificity
  y <- sensitivity
  auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  return(auc)
}

### 6) get_w48_aucs ############################################################

## Description: Calculates and compares the AUC for two logistic regression models
# using the Week 48 delta score: one with a baseline covariate and one without

## Input: A dataframe for a single cognitive test
# - df: The input df, which must contain the outcome variable,
# baseline_score, and the z-delta column for Week 48
# - outcome_var: A string specifying the name of the outcome column (CDPOS_W240)

## Output: A one-row tibble containing the AUC for both models
# - with_baseline_auc
# - without_baseline_auc

get_w48_aucs <- function(df, outcome_var = "CDPOS_W240") {
  z_delta_col <- str_subset(names(df), "^z_delta_.*_W48$")
  
  model_data_base <- df |>
    ungroup() |>
    select(all_of(outcome_var), baseline_score, all_of(z_delta_col)) |>
    na.omit()
  
  glm_base <- glm(
    paste(outcome_var, "~", z_delta_col, "+ baseline_score"),
    data = model_data_base, family = binomial()
  )
  
  y_true_base <- model_data_base[[outcome_var]]
  y_score_base <- predict(glm_base, type = "response")
  auc_with_baseline <- compute_auc(y_true_base, y_score_base)
  
  model_data_no_base <- df |>
    ungroup() |>
    select(all_of(outcome_var), all_of(z_delta_col)) |>
    na.omit()
  
  glm_no_base <- glm(
    paste(outcome_var, "~", z_delta_col),
    data = model_data_no_base, family = binomial()
  )
  
  y_true_no_base <- model_data_no_base[[outcome_var]]
  y_score_no_base <- predict(glm_no_base, type = "response")
  auc_without_baseline <- compute_auc(y_true_no_base, y_score_no_base)
  
  tibble(
    with_baseline_auc = auc_with_baseline,
    without_baseline_auc = auc_without_baseline
  )
}
