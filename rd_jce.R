library(rdrobust)


###################
### Data
###################


# Let rd_data be dataset with the following columns: 
# - record_id: Patient identifier.
# - t: Number of weeks from cutoff (running variable for RD analysis).
# - in_range: Outcome variable (e.g., Time in Range).
# - prev_week_contact: Binary indicator for whether the patient was contacted in the previous week.
# - prev_2w_TIR: Time in Range two weeks prior.
# - Additional demographic information (for covariate balance checks)



###################
### RD Analysis
###################

#-------------------------
# Example for single run
#-------------------------


# Example of single unadjusted analysis with uniform kernel, clustered standard errors
rd = rdrobust(y = rd_data$in_range, 
              x = rd_data$t, 
              kernel = "uniform", 
              cluster = rd_data$record_id)
summary(rd)

# Covariate adjustment example
# Create columns for fixed effects for each patient using dummy variables
dummies = model.matrix(~ factor(record_id), data=rd_data)
dummies = dummies[,-1] #remove intercept
fe_controls = as.data.frame(dummies)

# Combine fixed effects with additional covariates
all_controls = cbind(dummies, rd_data %>% select(prev_week_contact, prev_2w_tir))

# Run RD analysis with covariate adjustments
rd_fe = rdrobust(y = rd_data$in_range,  
                 x = rd_data$t, 
                 kernel = "uniform", 
                 cluster = rd_data$record_id, 
                 covs = all_controls)
summary(rd_fe)


#------------------------
# Multiple Runs
#------------------------


# Function to extract results from an 'rdrobust' object
# Arguments:
# - rd_obj: Output of rdrobust function
# - model_name: Name of the model being evaluated
# - control_type: Type of controls used in the model
# - bw_type: Type of bandwidth used in the model
# Returns: Data Frame with model results

extract_rd_info <- function(rd_obj, model_name, control_type, bw_type) {
    # Extracting main results from rdrobust object
    est <- rd_obj$coef[1]        # Point estimate
    p_val <- rd_obj$pv[3]        # p-value for the treatment effect
    ci_lower <- rd_obj$ci[3,1]   # Lower bound of CI
    ci_upper <- rd_obj$ci[3,2]   # Upper bound of CI
    
    # Bandwidths (left and right)
    bw_left <- rd_obj$bws[1,1]
    bw_right <- rd_obj$bws[1,2]
    
    # Include control_type and bw_type in the returned data
    res <- data.frame(
        model = model_name,
        control_type = control_type,
        bw_type = bw_type,
        est = est,
        p_val = p_val,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        bw_left = bw_left,
        bw_right = bw_right,
        stringsAsFactors = FALSE
    )
    
    return(res)
}


# We want to run the analyses for two bandwidth types:
# MSE optimal and CE optimal bandwidths

#Define bandwidth types
bw_types <- list(
    default_bw = list(bwselect = "mserd", h = NULL),
    cerrd_bw   = list(bwselect = "cerrd", h = NULL)
)

# We want the analysis for 3 models
# 1) Unadjusted model, 
# 2) Covariate adjusted model that drops missing values
# 3) Covariate adjusted model with imputed missing values

#Define control types
control_types <- list(
    no_control = NULL,
    original_controls = all_controls_original,
    imputed_controls = all_controls_imputed
)


# Run analyses for all combinations of control and bandwidth types

results_list <- list()
counter <- 1

for (ctrl_name in names(control_types)) {
    covs_data <- control_types[[ctrl_name]]
    
    for (bw_name in names(bw_types)) {
        bw_choice <- bw_types[[bw_name]]
        
        # Construct a descriptive model name
        model_name <- paste(ctrl_name, bw_name, sep = "_")
        
        # Run rdrobust
        rd_obj <- rdrobust(y = rd_data$in_range, 
                           x = rd_data$t, 
                           kernel = "uniform",
                           cluster = rd_data$record_id,
                           covs = covs_data,
                           bwselect = bw_choice$bwselect,
                           h = bw_choice$h)
        print(summary(rd_obj))
        print(rd_obj$coef)
        print(rd_obj$se)
        print(rd_obj$pv)
        print(rd_obj$ci)
        
        # Extract info
        res <- extract_rd_info(rd_obj, model_name, control_type = ctrl_name, bw_type = bw_name)
        print(res)
        results_list[[counter]] <- res
        counter <- counter + 1
    }
}

# Combine all results into a single data frame
final_results <- do.call(rbind, results_list)


#-----------------------------------------
# Covariate Balance Check
#-----------------------------------------

# Check covariate balance using the CE optimal bandwidth
# - Loop through each covariate of interest
# - Test the discontinuity in each covariate at the cutoff 

# Loop through each covariate and run rdrobust
# let 'covariates' be a list of the names of the columns you want to check

for (covariate in covariates) {
    outcome <- rd_pts[[covariate]]
    result <- rdrobust(y = outcome, x = rd_data$t, 
                       kernel = "uniform", 
                       bwselect = "cerrd", 
                       cluster = rd_data$record_id)
    
    point_estimate <- result$coef[1]
    conf_int <- result$ci[3, ]
    p_value <- result$pv[3]
    
    # Store the results in the results data frame
    results <- rbind(results, data.frame(model = covariate,
                                         est = point_estimate,
                                         ci_lower = conf_int[1],
                                         ci_upper = conf_int[2],
                                         p_val = p_value))
}

# Print the summary table
print(results)
