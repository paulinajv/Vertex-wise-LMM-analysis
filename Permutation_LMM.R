" 
    In this code, we will run a 5,000 Linear Mixed Effect Model (LMM) tests (vertex by vertex) by shuffling the labels each time.

    Currently, the code is designed to analyze one metric at a time. To select the metric you want to measure, 
    you need to modify the variable on line 25. The code can be further improved by implementing a loop to iterate 
    through all metrics automatically, which is a pending task.

    This script will perform the permutation LMM for a single metric across all temporal time points simultaneously. 
    As a result, you will obtain four dataframes (P30, P60, P120, and P150), each containing the p-values from the permutations.

    Please note that this process is computationally intensive and time-consuming. It is recommended to run it in the background.
        
"

## Loading the required packages
library(dplyr)
library(tidyverse)
library(lme4)
library(pbkrtest)
library(effects)
library(emmeans)
library(readr)


set.seed(0002)

## Converting the temporal time points to factors
mrds$age <- factor(mrds$age, levels = c("30", "60", "120", "150"))

## Select the metric to test. 
metric <- mrds$FApar # suing FApar as an example. But it could be either: FApar, FAperp, MDpar, MDperp, FA, MD, RD, AD)

## Select the number of permutations
num_perm <- 5000

## Function to perform LMM and extract p-values
lmerFun <- function(vertices, index) {
  Modl <- lmer(metric ~ age + grp + age * grp + (1 | animalID), data = vertices)
  emms <- emmeans(Modl, ~grp | age) %>% pairs(., adjust = "tukey")
  pvals <- summary(emms)$p.value[index]
  
}


## Create a matrix to store results
streams <- unique(mrds$stream)
points <- unique(mrds$point)

## Assigning the index for p-value extraction
index <- 1:4

## Create a multidimensional matrix to store the results
results <- array(NA, dim = c(length(streams), length(points), num_permutations, length(index)))

## Assign the points and stream/gridlines as a column names in the matrix
colnames(results) <- points
rownames(results) <- streams

## Perform the nested loop to iterate through each set pf stream/gridlines and points
for (s in 1:length(streams)) {
  for (p in 1:length(points)) {
    
    # Subset the data
    vertices <- subset(mrds, stream == streams[s] & point == points[p])
    
    # Loop through permutations
    for (nperm in 1:num_permutations) {
      # Shuffle the response variable
      vertices$FApar <- sample(vertices$FApar) ## Don't forget to add the metric
      # Perform LMM and extract p-values
      lmerTest <- lmerFun(vertices, index)
      # Store the results
      results[s, p, nperm, ] <- lmerTest
      
    }
    
  }
  
}

## Extract the results from the dimensions
dim1 <- results[, , , 1]
dim2 <- results[, , , 2]
dim3 <- results[, , , 3]
dim4 <- results[, , , 4]

## Convert the matrices to dataframes
FApar_5kperm_30 <- as.data.frame(as.table(dim1))
FApar_5kperm_60 <- as.data.frame(as.table(dim2))
FApar_5kperm_120 <- as.data.frame(as.table(dim3))
FApar_5kperm_150 <- as.data.frame(as.table(dim4))



