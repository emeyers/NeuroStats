#' A one-way ANOVA for binned data
#' 
#' An implementation for a one-way ANOVA test on neural data
#' 
#' @param binned_data This argument takes a binned version of a spike train data set. 
#' 
#' @param labels This argument refers to the groups across which the ANOVA is being run. For example, if the ANOVA test is attempting to assess whether the means are different categorical values of "column_3", the argument should be set as labels = column_a. Crucially, it takes as a value a column within binned_data.
#' 
#' @param include_site_info Default is TRUE. This includes the SiteID column of the binned_data, if available. 
#' 
#' @return It returns a table of values that are results from the ANOVA test. 



# the ANOVA test
#' @export
run_ANOVA <- function(binned_data, labels, include_site_info = TRUE) {
  
  # load the binned data (works for both a string file name and a data frame)
  binned_data <- NeuroDecodeR:::check_and_load_binned_data(binned_data)
  
  long_data <- binned_data |>
    rename(labels = all_of(paste0("labels.", labels))) |>
    #select("siteID", paste0("labels.", labels), starts_with("time")) |>
    pivot_longer(starts_with("time"),
                 names_to = "time_period",
                 values_to = "activity") 
  
  anova_results <- long_data |>
    group_by(siteID) |>
    mutate(df1 = length(unique(labels)) - 1) |>
    group_by(siteID, time_period) |>
    mutate(df2 = n() - df1 - 1) |>
    group_by(siteID, labels, time_period) |>
    mutate(group_mean = mean(activity)) |>
    group_by(siteID, time_period) |>
    mutate(overall_mean = mean(activity)) |>
    mutate(group_deviation_squared = (group_mean - overall_mean)^2) |>
    group_by(siteID, time_period) |>
    mutate(residual_squared = (group_mean - activity)^2) |>
    group_by(siteID, time_period) |>
    summarize(df1 = mean(df1),
              df2 = mean(df2),
              SS_group = sum(group_deviation_squared),
              SS_residual = sum(residual_squared),
              F_stat = (1/df1 * SS_group)/(1/df2 * SS_residual),
              p_val = pf(F_stat, df1, df2, lower.tail = FALSE), .groups = "drop")
  
  
  if (include_site_info) {
    
    site_info <- select(binned_data, "siteID", starts_with("site_info")) |>
      distinct()
    
    # For each site, there should only be one unique value for rows for the
    # site_info columns. (i.e., each row of site_info contains redundant
    # information). If this is not the case, send a warning and only use the first 
    # entry for each site. 
    if (max(table(site_info$siteID)) > 1) {
      
      warning(paste("A least one of binned data's site_info columns contains different values",
                    "in different rows for the same siteID. Returning only the site_info value",
                    "for the first row for each site."))
      
      site_info <- site_info |>
        group_by(siteID) |>
        slice(1)
      
    }
    
    
    anova_results <- anova_results |>
      left_join(site_info, by = "siteID") |>
      select("siteID", starts_with("site_info"), everything())
    
  }
  
  
  
  anova_results
  
  
}

devtools::document()

?NeuroStats::run_ANOVA
