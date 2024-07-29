
"Performing Cohen's f^2 for each gridline and point across all the temportal time-points (P30, P60, P120, P150)"

## Load the necessary libraries
library(effectsize)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(effects)


## Load the MRDS data
mrds_data 

## convert data time-points into factors
mrds_data$age <- factor(mrds_data$age, levels = c("30", "60", "120", "150"))

## Function to calculate the f^2

lmer_w_f2 <- function(vertice) {
  metrics <- c("FApar", "FAperp", "MDpar", "MDperp") ## add diffusion metrics 
  
  f2_values <- numeric(length(metrics))
  
  for (i in seq_along(metrics)) {
    m <- metrics[i]
    form <- as.formula(paste(m, "~ age + grp + age*grp + (1 | animalID)")) 
    myMod <- lmer(form, data = vertice)
    t <- contrast(emmeans(myMod, ~grp | age), method = "pairwise")
    ## Both tratio and dferr needs to be modify for each time-point (1:P30, 2:P60, 3:P120, 4:P150)
    tratio <- summary(t)$t.ratio[1] 
    dferr <- summary(t)$df[1]
    f2_value <- t_to_f2(t = tratio, df_error = dferr, paired = TRUE)
    f2_values[i] <- f2_value$Cohens_f2_partial
  }
  
  return(f2_values)
}


## Make each stream/gridline and point unique. This is necessary for the follow iterations
streams <- unique(mrds_data$stream)
points <- unique(mrds_data$point)

## Create an empty list to store the results
results_list <- list()

## Nested loop to iterate between stream/gridlines and points
for (s in seq_along(streams)) {
  for (p in seq_along(points)) {
    stream <- streams[s]
    point <- points[p]
    vertice <- subset(mrds_data, stream == streams[s] & point == points[p])
    f2_test <- lmer_w_f2(vertice)
    
    result_entry <- list(
      stream = stream,
      point = point,
      metrics = c("FApar", "FAperp", "MDpar", "MDperp"),
      f2_values = f2_test
    )
    
    results_list[[paste0("stream_", stream, "_point_", point)]] <- result_entry
  }
}

print(results_list)

## Convert the matrix into dataframe
Cohens_f2_results <- as.data.frame(do.call(rbind, results_list))










