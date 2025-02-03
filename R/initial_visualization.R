load("creation_VISp_200bins_25sampled.Rda")

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(NeuroDecodeR)

binned_data_name <- binned_data
labels <- "natural_scene_stimulus_id"

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

anova_results <- run_ANOVA(binned_data_name, labels)

## Based on the above anova_results data frame, choosing the time bin where each 
## neuron is most selective. In cases where there are multiple "most selective"
## time bins for each neuron, this picks one at random.

time_bin_selectivity <- anova_results %>%
  select(siteID, time_period, p_val, F_stat)

best_time_df <- time_bin_selectivity %>%
  group_by(siteID) %>%
  filter(p_val == min(p_val)) %>%
  slice_sample(n = 1)

## Now, I will create the first visualization, which is a scatterplot ranking
## stimuli for a given neuron.
ID_index <- sample(1:length(best_time_df$siteID), 1)
selected_siteID <- best_time_df$siteID[ID_index]
best_time <- best_time_df$time_period[best_time_df$siteID == selected_siteID]

single_neuron_df <- binned_data %>%
  filter(siteID == selected_siteID, !is.na(!!sym(best_time))) %>%
  group_by(labels.natural_scene_stimulus_id) %>%
  summarize(mean_fr = mean(!!sym(best_time), na.rm = TRUE), 
            sd_fr = sd(!!sym(best_time), na.rm = TRUE))

single_neuron_df <- single_neuron_df %>%
  mutate(labels.natural_scene_stimulus_id = reorder(labels.natural_scene_stimulus_id, -mean_fr))

v1.1 = ggplot(data = single_neuron_df, aes(x = labels.natural_scene_stimulus_id, y = mean_fr))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_fr - sd_fr, ymax = mean_fr + sd_fr), width = 0.2, color = "red")+
  theme_minimal() +
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Stimulus", y = "Mean Firing Rate", title = "Firing Rate by Stimulus",
       subtitle = paste(paste("siteID:", selected_siteID, sep = " "), paste("time bin:", best_time, sep = " "), sep = ", "))

## Now, doing this for just the first 25 and last 25 trials (to see if that has an impact).
trial_cutoff <- 25

single_first_section <- binned_data %>%
  filter(trial_number %in% c(1:trial_cutoff)) %>%
  filter(siteID == selected_siteID, !is.na(!!sym(best_time))) %>%
  group_by(labels.natural_scene_stimulus_id) %>%
  summarize(mean_fr = mean(!!sym(best_time), na.rm = TRUE), 
            sd_fr = sd(!!sym(best_time), na.rm = TRUE))

single_first_section <- single_first_section %>%
  mutate(labels.natural_scene_stimulus_id = reorder(labels.natural_scene_stimulus_id, -mean_fr))

v1.2 = ggplot(data = single_first_section, aes(x = labels.natural_scene_stimulus_id, y = mean_fr))+
  geom_point()+
  #geom_errorbar(aes(ymin = mean_fr - sd_fr, ymax = mean_fr + sd_fr), width = 0.2, color = "red")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Stimulus", y = "Mean Firing Rate", 
       title = paste("Firing Rate by Stimulus,", "First", trial_cutoff, "Trials", sep = " "),
       subtitle = paste(paste("siteID:", selected_siteID, sep = " "), 
                        paste("time bin:", best_time, sep = " "), sep = ", "))


single_first_section <- binned_data %>%
  filter(trial_number %in% c(1:trial_cutoff)) %>%
  filter(siteID == selected_siteID, !is.na(!!sym(best_time))) %>%
  group_by(labels.natural_scene_stimulus_id) %>%
  summarize(mean_fr = mean(!!sym(best_time), na.rm = TRUE), 
            sd_fr = sd(!!sym(best_time), na.rm = TRUE))

last_part <- length(unique(binned_data$trial_number)) - trial_cutoff
single_last_section <- binned_data %>%
  filter(trial_number %in% c(last_part:unique(length(binned_data$trial_number)))) %>%
  filter(siteID == selected_siteID, !is.na(!!sym(best_time))) %>%
  group_by(labels.natural_scene_stimulus_id) %>%
  summarize(mean_fr = mean(!!sym(best_time), na.rm = TRUE), 
            sd_fr = sd(!!sym(best_time), na.rm = TRUE))

single_last_section <- single_last_section %>%
  mutate(labels.natural_scene_stimulus_id = reorder(labels.natural_scene_stimulus_id, -mean_fr))

v1.3 = ggplot(data = single_last_section, aes(x = labels.natural_scene_stimulus_id, y = mean_fr))+
  geom_point()+
  #geom_errorbar(aes(ymin = mean_fr - sd_fr, ymax = mean_fr + sd_fr), width = 0.2, color = "red")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Stimulus", y = "Mean Firing Rate", 
       title = paste("Firing Rate by Stimulus,", "Last", trial_cutoff, "Trials", sep = " "),
       subtitle = paste(paste("siteID:", selected_siteID, sep = " "), 
                        paste("time bin:", best_time, sep = " "), sep = ", "))

v1.1
v1.2
v1.3

## Now I make the second set of visualizations, where I focus on looking across
## all neurons to assess the performance on each stimuli.

## It was challenging to determine the best firing rate to use for this analysis,
## since this normally varies on a neuron-by-neuron basis. Thus, in order to 
## do this, I decided to just see which time bin appears the most frequently as
## the best time bin for each neuron. This is not ideal, so I'm looking forward
## to working through alternatives.

time_bin_selectivity <- anova_results %>%
  select(siteID, time_period, p_val, F_stat)

best_time_df <- time_bin_selectivity %>%
  group_by(siteID) %>%
  filter(p_val == min(p_val))

best_time <- names(which.max(table(best_time_df$time_period)))
best_time

## It looks like the most common time bin cited as the most selective was "time.
## 75_275". This passes the initial gut check, as this is relatively close to 
## stimulus onset. 

overall_dist_df <- binned_data %>%
  filter(!is.na(!!sym(best_time))) %>%
  group_by(labels.natural_scene_stimulus_id, siteID) %>%
  summarize(mean_fr = mean(!!sym(best_time), na.rm = TRUE), 
            sd_fr = sd(!!sym(best_time), na.rm = TRUE))

## To remove outliers:
df_no_outliers <- overall_dist_df %>%
  group_by(labels.natural_scene_stimulus_id) %>%
  mutate(
    Q1 = quantile(mean_fr, 0.25),
    Q3 = quantile(mean_fr, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(mean_fr >= lower_bound & mean_fr <= upper_bound) %>%
  select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)

five_no_summary <- df_no_outliers %>%
  group_by(labels.natural_scene_stimulus_id) %>%
  summarize(min = min(mean_fr),
            Q1 = quantile(mean_fr, 0.25),
            median = median(mean_fr),
            Q3 = quantile(mean_fr, 0.75),
            max = max(mean_fr))

five_no_summary <- five_no_summary %>%
  mutate(labels.natural_scene_stimulus_id = reorder(labels.natural_scene_stimulus_id, -median))

ggplot(data = five_no_summary, aes(x = labels.natural_scene_stimulus_id, 
                                   ymin = min, lower = Q1, middle = median,
                                   upper = Q3, ymax = max))+
  geom_boxplot(stat = "identity")+
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Stimulus Labels", y = "Median Firing Rate", title = "Median Firing Rate by Stimulus")
  
