"In this code, we perform a Linear Mixed-effects Model (LMM) analysis and simultaneously calculate Cohen's f² effect size. 
 The workflow consists of the following key components:

    a) Functions: Two custom functions are used to extract Cohen's f² effect sizes and p-values from post hoc tests.

    b) Data iteration: A nested loop iterates through each stream/grid and data point (referred to as a vertex) 
    to apply the LMM and extract the statistics.

    c) Results compilation: The extracted results are then organized and converted into a dataframe for further visualization"



## Load the necessary libraries
library(effectsize)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(effects)


## Load the MRDS data
mrds <- read_csv("MRDS_database")

## convert data time-points into factors
mrds$age <- factor(mrds$age, levels = c("30", "60", "120", "150"))

## This function is to extract the pvalues from the LMER and post-hoc tests
lmer_ph <- function(vertice, index) {
  
  metrics <- c("FApar", "FAperp", "MDpar", "MDperp")
  
  ph_values <- numeric(length(metrics))
  
  for (i in seq_along(metrics)) {
    m <- metrics[i]
    form <- as.formula(paste(m, "~ age + grp + age*grp + (1 | animalID)")) 
    myMod <- lmer(form, data = vertice)
    
    emms <- emmeans(myMod, ~grp | age) %>% pairs(., adjust = "tukey")
    pvalue <- summary(emms)$p.value[index] 
    ph_values[i] <- pvalue
    
  }
  return(ph_values)
}
## Function to calculate the f^2

lmer_w_f2 <- function(vertice, index) {
  metrics <- c("FApar", "FAperp", "MDpar", "MDperp") ## add diffusion metrics 
  
  f2_values <- numeric(length(metrics))
  
  
  for (i in seq_along(metrics)) {
    m <- metrics[i]
    form <- as.formula(paste(m, "~ age + grp + age*grp + (1 | animalID)")) 
    myMod <- lmer(form, data = vertice)
    t <- contrast(emmeans(myMod, ~grp | age), method = "pairwise")
    
    tratio <- summary(t)$t.ratio[index] 
    dferr <- summary(t)$df[index]
    f2_value <- t_to_f2(t = tratio, df_error = dferr, paired = TRUE)
    f2_values[i] <- f2_value$Cohens_f2_partial
    
    
  }
  
  return(f2_values)
}



## Make each stream/gridline and point unique. This is necessary for the follow iterations
streams <- unique(mrds$stream)
points <- unique(mrds$point)

## Create an empty list to store the results
results_list <- list()

index=1:4 #indexing the postnatal days [1=P30, 2=P60, 3=P120, 4=150]

## Nested loop to iterate between stream/gridlines and points
for (s in seq_along(streams)) {
  for (p in seq_along(points)) {
    stream <- streams[s]
    point <- points[p]
    vertice <- subset(mrds, stream == streams[s] & point == points[p])
    
    for (i in index) {
      
      
      f2_test <- lmer_w_f2(vertice, i)
      pvalue <- lmer_ph(vertice, i)
      
      result_entry <- list(
        stream = stream,
        point = point,
        metrics = c("FApar", "FAperp", "MDpar", "MDperp"),
        f2_values = f2_test,
        pval = pvalue,
        index = i
        
      )
      
      results_list[[paste0("stream_", stream, "_point_", point, "index_", i)]] <- result_entry
      
    }
  }
}

print(results_list)


## Convert the matrix into dataframe
lmer_and_f2_df <- as.data.frame(do.call(rbind, results_list))

## Unnest the columns containing nested data
lmer_and_f2_df <- unnest(lmer_and_f2_df, cols = 1:6)




