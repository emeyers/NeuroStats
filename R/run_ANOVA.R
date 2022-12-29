

library(dplyr)
library(tidyr)
library(ggplot2)
library(NeuroDecodeR)



# Add roxygen documentation, what is listed below could be an example

binned_data_name <- "../inst/extdata/ZD_150bins_50sampled.Rda"

labels <- "stimulus_ID"

anova_results <- run_ANOVA(binned_data_name, labels)


# Add this as a generic plot function...

significance_level <- .0001

anova_results |> 
  mutate(stat_sig = p_val < significance_level) |> 
  group_by(time_period) |>
  summarize(percent_selective = 100 * mean(stat_sig)) |>
  mutate(Time  = NeuroDecodeR:::get_center_bin_time(time_period)) |>
  ggplot(aes(Time, percent_selective)) +
  geom_point() + 
  geom_line() + 
  ylab(paste("Percent of sites selective   (p < ", significance_level, ")")) +
  geom_hline(yintercept = 100 * significance_level, col = "red") +
  theme_bw()


# 
# # test that my anova gives identical results to R's anova function
# 
# iSite = 2
# time_period_to_test <- "time.-100_49"
# 
# load(binned_data_name)
# one_site_data <- binned_data |>
#   filter(siteID == iSite) |>
#   rename(labels = all_of(paste0("labels.", labels))) |>
#   rename(activity = time_period_to_test) 
# 
# r_anova <- anova(lm(activity ~ labels, data = one_site_data))
# 
# my_anova <- filter(anova_results, siteID == iSite,  time_period == time_period_to_test)
# 
# 
# # check if R's anova is giving the same results as my anova 
# # (to within very small rounding error)
# (r_anova$`Pr(>F)`[1] - my_anova$p_val) < 10^-10



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

