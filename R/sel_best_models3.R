#' Select best models
#'
#' @description
#' This function selects the best models according to user-defined criteria, evaluating statistical significance (partial ROC), predictive ability (omission rates), and model complexity (AIC).
#'
#' @param calibration_results an object of class `calibration_results` returned
#' by the \code{\link{calibration}}() function. Default is NULL.
#' @param cand_models (data.frame) a summary of the evaluation metrics for each
#' candidate model. In the output of the \code{\link{calibration}}(), this
#' data.frame is located in `$calibration_results$Summary`. Default is NULL.
#' @param model_type (character) model type, either "glm" or "glmnet".
#' @param test_concave (logical) whether to remove candidate models presenting
#' concave curves. Default is TRUE.
#' @param omrat_threshold (numeric) the maximum omission rate a candidate model
#' can have to be considered a best model. Default is 10. This value must match
#' one of the values specified in `omrat` in \code{\link{calibration}}().
#' @param allow_tolerance (logical) whether to allow selection of models with
#' minimum values of omission rates even if their omission rate surpasses the
#' `omrat_threshold`. This is only applicable if all candidate models have
#' omission rates higher than the `omrat_threshold`. Default is TRUE.
#' @param tolerance (numeric) The value added to the minimum omission rate if it
#' exceeds the `omrat_threshold`. If `allow_tolerance = TRUE`, selected models
#' will have an omission rate equal to or less than the minimum rate plus this
#' tolerance. Default is 0.01.
#' @param AIC (character) the type of AIC to be calculated: "ws" for AIC
#' proposed by Warren and Seifert (2011), or "nk" for AIC proposed by Ninomiya
#' and Kawano (2016). This is only applicable if model_type = "glmnet".
#' Default is "ws". See References for details.
#' @param significance (numeric) the significance level to select models
#' based on the partial ROC (pROC). Default is 0.05. See Details.
#' @param delta_aic (numeric) the value of delta AIC used as a threshold to
#' select models. Default is 2.
#' @param verbose (logical) whether to display messages during processing.
#' Default is TRUE.
#'
#' @details
#' Partial ROC is calculated following Peterson et al.
#' (2008; http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008).
#'
#' @return
#' If calibration_results is provided, it returns a new calibration_results with
#' the new selected models and summary. If calibration_results is NULL, it
#' returns a list containing the following elements:
#' - selected_models: data frame with the ID and the summary of evaluation
#' metrics for the selected models.
#' - summary: A list containing the delta AIC values for model selection, and
#' the ID values of models that failed to fit, had concave curves,
#' non-significant pROC values, omission rates above the threshold, delta AIC
#' values above the threshold, and the selected models.
#'
#' @export
#' @references
#' Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in generalized
#' linear models." (2016): 2537-2560.
#'
#' Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent:
#' the importance of model complexity and the performance of model selection
#' criteria. Ecological applications, 21(2), 335-342.
#'
#' @examples
#' # Import example of calibration results (output of calibration function)
#' ## GLM
#' data("calib_results_glm", package = "kuenm2")
#'
#' #Select new best models based on another value of omission rate
#' new_best_model <- sel_best_models3(cand_models = calib_results_glm$calibration_results$Summary,
#'                                   model_type = "glm",
#'                                   test_concave = TRUE,
#'                                   omrat_threshold = 5,
#'                                   allow_tolerance = TRUE,
#'                                   tolerance = 0.01,
#'                                   AIC = "ws",
#'                                   significance = 0.05,
#'                                   delta_aic = 10, #Higher value of delta AIC
#'                                   verbose = TRUE)
#' #Compare with best models selected previously, with omission rate of 10 and delta AIC of 2
#' calib_results_glm$summary$Selected #Models 1, 2 and 5 selected
#' new_best_model$summary$Selected #Models 1 and 5 selected
#' #Replace selected models in calib_results
#' calib_results_glm$selected_models <- new_best_model$cand_final
#' calib_results_glm$summary <- new_best_model$summary
#'
sel_best_models3 <- function(data,
                            calibration_results = NULL,
                            cand_models = NULL,
                            model_type = c("glmnet", "glm"),
                            test_concave = TRUE,
                            omrat_threshold = 10,
                            allow_tolerance = TRUE,
                            tolerance = 0.01,
                            AIC = "ws",
                            significance = 0.05,
                            delta_aic = 2,
                            addsamplestobackground = TRUE,
                            weights = NULL,
                            verbose = TRUE) {
  #Check data
  if(is.null(calibration_results) & is.null(cand_models)){
    stop("You must specified calibration_results or cand_models")
  }

  if(!is.null(calibration_results)){
    cand_models <- calibration_results$calibration_results$Summary
  }

  # Adjust AIC column based on model type
  if (model_type == "glmnet") {
    # Remove the unused AIC column in glmnet
    if (AIC == "nk") {
      AIC <- "AIC_nk"
      cand_models$AIC_ws <- NULL
    } else if (AIC == "ws") {
      AIC <- "AIC_ws"
      cand_models$AIC_nk <- NULL
    } else {
      stop("Unsupported AIC type. Please use 'nk' or 'ws'.")
    }
  } else if (model_type == "glm") {
    AIC <- "AIC" # For glm models, we only use a single AIC column
  } else {
    stop("Unsupported model type. Please use 'glmnet' or 'glm'.")
  }

  # Omission rate column name
  om_thr <- paste0("Omission_rate_at_", omrat_threshold, ".mean")

  #proc-pval columns
  #proc_pval <- paste0("pval_pROC_at_", omrat_threshold, ".mean")

  any_bad <- TRUE #To initiate looping
  id_to_remove <- 0

  while(any_bad){
    #Remove bad models
  cand_models <- cand_models[!(cand_models$ID %in% id_to_remove),]


  # Log the number of models being filtered
  if (verbose) {
    message("\nFiltering ", nrow(cand_models), " models")
  }

  # Remove models with errors
  na_models <- cand_models[is.na(cand_models$is_concave), "ID"]
  if (verbose) {
    message("Removing ", length(na_models), " model(s) because they failed to fit")
  }
  cand_models <- cand_models[!is.na(cand_models$is_concave), ]

  # Remove concave curves if test_concave is TRUE
  if (test_concave) {
    concave_models <- cand_models[cand_models$is_concave, "ID"]
    if (verbose) {
      message("Removing ", length(concave_models), " model(s) with concave curves")
    }
    cand_models <- cand_models[!cand_models$is_concave, ]
  } else {
    concave_models <- 0
  }

  # Subset models by omission rate
  high_omr <- cand_models[cand_models[, om_thr] > omrat_threshold / 100, "ID"]
  cand_om <- cand_models[cand_models[, om_thr] <= omrat_threshold / 100, ]
  if (verbose) {
    message(nrow(cand_om), " models were selected with omission rate below ", omrat_threshold, "%")
  }

  # Stop if no models meet the omission rate threshold and allow_tolerance is FALSE
  if (nrow(cand_om) == 0 & !allow_tolerance) {
    stop("There are no models with values of omission rate below ", omrat_threshold, "%. Try with allow_tolerance = TRUE.")
  }

  # Apply tolerance if no models meet the omission rate threshold and allow_tolerance is TRUE
  if (nrow(cand_om) == 0 & allow_tolerance) {
    min_thr <- min(cand_models[, om_thr])
    cand_om <- subset(cand_models, cand_models[, om_thr] <= min_thr + tolerance)
    high_omr <- cand_models[cand_models[, om_thr] > min_thr + tolerance, "ID"]
    if (verbose) {
      message("Minimum value of omission rate (", round(min_thr * 100, 1), "%) is above the selected threshold (", omrat_threshold, "%).\nApplying tolerance and selecting ", nrow(cand_om), " models with omission rate <", round(min_thr * 100 + tolerance, 1), "%")
    }
  }

  # Calculate delta AIC and select models based on delta AIC
  cand_om$dAIC <- cand_om[, AIC] - min(cand_om[, AIC])
  high_aic <- cand_om[cand_om$dAIC > delta_aic, "ID"]
  cand_final <- cand_om[cand_om$dAIC <= delta_aic, ]

  if (verbose) {
    message("Selecting ", nrow(cand_final), " final model(s) with delta AIC <", delta_aic)
  }

  # Validate pROC
  if (verbose) {
    message("Validating pROC...")
  }
  proc_values <- validate_proc(cand_final, data, omrat_threshold,
                               addsamplestobackground, weights,
                               model_type)

  # Create a copy of cand_final to keep the original data unchanged
  cand_final_updated <- cand_final

  # Identify matching rows in cand_final
  match_idx <- match(cand_final_updated$ID, proc_values$ID)

  # Replace NA values in cand_final_updated with corresponding values from proc_values
  for (col in names(proc_values)[-ncol(proc_values)]) {  # Ignoring the ID column
    na_idx <- is.na(cand_final_updated[[col]])  # Identifying NAs in cand_final_updated column
    cand_final_updated[[col]][na_idx] <- proc_values[[col]][match_idx][na_idx]
  }

  # Check if p_value is non-significative
  p_value_omr <- cand_final_updated[,paste0("pval_pROC_at_",
                                            omrat_threshold, ".mean")]
  any_bad <- any(p_value_omr > significance)
  #Get models to remove, if necessary
  id_to_remove <- cand_final_updated$ID[p_value_omr > significance]

  if(any_bad & verbose){
    message("Models with non-significant pROC values were identified. Re-selecting models...")
  }

  }

  # Final results
  sel_res <- list(cand_final = cand_final_updated,
                  summary = list("delta_AIC" = delta_aic,
                                 "omission_rate_thr" = omrat_threshold,
                                 "Errors" = na_models,
                                 "Concave" = concave_models,
                                 "High_omission_rate" = high_omr,
                                 "High_AIC" = high_aic,
                                 "Selected" = cand_final_updated$ID))

  #Update calibration_results if necessary
  if(!is.null(calibration_results)){
    calibration_results$selected_models <- sel_res$cand_final
    calibration_results$summary <- sel_res$summary
    calibration_results$omission_rate <- omrat_threshold
    return(calibration_results)
  } else {
    return(sel_res)}
}


validate_proc <- function(cand_final, data, omission_rate,
                          addsamplestobackground, weights,
                          model_type) {
  #Validate pROC
  res_proc <- lapply(1:nrow(cand_final), function(x){
    grid_x <- cand_final[x,] # Get i candidate model
    #Get formula and reg
    formula_x <- as.formula(grid_x$Formulas)
    reg_x <- grid_x$regm

    # Get background index
    bgind <- which(data$calibration_data$pr_bg == 0)

    # Fit models using k-fold cross-validation
    mods <- try(lapply(1:length(data$kfolds), function(i) {
      notrain <- -data$kfolds[[i]]
      data_i <- data$calibration_data[notrain,]

      # Set weights per k-fold
      if (!is.null(weights)){
        weights_i <- weights[notrain]
      } else {
        weights_i <- NULL
      }

      if (model_type == "glmnet") {
        # Run glmnet model
        mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights_i, calculate_AIC = FALSE)
      } else {
        # Run glm model
        mod_i <- glm_mx(formula = formula_x,
                        family = binomial(link = "cloglog"),
                        data = data_i, weights = weights_i)
      }

      # Predict model
      pred_i <- if (model_type == "glmnet") {
        as.numeric(predict.glmnet_mx(object = mod_i,
                                     newdata = data$calibration_data,
                                     clamp = FALSE, type = "cloglog"))
      } else if (model_type == "glm") {
        enmpa::predict_glm(model = mod_i,
                           newdata = data$calibration_data,
                           type = "response")
      }

      # Extract suitability in train and test points
      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

      #Proc
      proc_i <- lapply(omission_rate, function(omr){
        proc_omr <- enmpa::proc_enm(test_prediction = suit_val_eval,
                                    prediction = pred_i,
                                    threshold = omr)$pROC_summary
        names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                             paste0("pval_pROC_at_", omr))
        return(proc_omr)
      })
      proc_i <- unlist(proc_i)

      # Save metrics in a dataframe
      df_proc <-  if (model_type == "glmnet") {
        data.frame(Replicate = i,
                   t(proc_i),
                   row.names = NULL)
      } else if (model_type == "glm") {
        data.frame(Replicate = i,
                   t(proc_i),
                   row.names = NULL)
      }
      return(df_proc)
    }), silent = TRUE)
    #Get summary
    proc_df <- do.call("rbind", mods)

    #Get means and sd
    means <- sapply(proc_df[, -1], mean)  # Excluindo a coluna Replicate
    sds <- sapply(proc_df[, -1], sd)

    #Create new dataframe
    proc_df <- data.frame(
      t(c(means, sds))
    )
    #Rename
    names(proc_df) <- c(paste0(names(means), ".mean"), paste0(names(sds), ".sd"))

    #Append model ID
    proc_df$ID <- grid_x$ID

    return(proc_df)
  })

  #Bind results
  return(do.call("rbind", res_proc))
}
