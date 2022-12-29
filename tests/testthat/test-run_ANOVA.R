test_that("run_ANOVA works", {
  
  binned_data_name <- "ZD_150bins_50sampled.Rda"
  
  labels <- "stimulus_ID"
  
  anova_results <- run_ANOVA(binned_data_name, labels)
  
  
  # test that my anova gives identical results to R's anova function
  
  iSite = 2
  time_period_to_test <- "time.-100_49"
  
  load(binned_data_name)
  one_site_data <- binned_data |>
    filter(siteID == iSite) |>
    rename(labels = all_of(paste0("labels.", labels))) |>
    rename(activity = all_of(time_period_to_test))
  
  r_anova <- anova(lm(activity ~ labels, data = one_site_data))
  
  my_anova <- filter(anova_results, siteID == iSite,  time_period == time_period_to_test)
  
  
  # check if R's anova is giving the same results as my anova 
  # (to within very small rounding error)
 
  expect_equal(r_anova$`Pr(>F)`[1], my_anova$p_val)
  
  
})
