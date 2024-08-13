" 
    In this code, we will run permutation LMM test for each vertice. In general, the framework 
    is to suffle the labels of each metric an re-run the LMM, and so on for 5k times. Then, the
    pvalue from each posthoc test will be saved in the matrix. The resulting matrices will be
    used to calculate the cluster inference (see the code named Cluster_inference.R)
    
    
   IMPORTANT NOTE:

   1)  Currently, the code is designed to analyze one metric at a time, mostly to keep control of 
       the data, however it could be improved by implementing an extra nested loop. 
       In the meanwhile, to select the metric you want to measure, you need to modify the variable on line 40. 
    
   2)  This script will perform the LMM for a single metric across all temporal time points simultaneously. 
       As a result, you will obtain four dataframes (P30, P60, P120, and P150) with their corresponding p-values.
    
  
   3)  Please keep in mind that this script is computationally intensive and time-consuming. 
       It is recommended to run it in the background if possible.
        
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
metric <- mrds$FApar # using FApar as an example. Here change to the desired metric (FApar, FAperp, FAparpar, FAparperp)

## Use a function to run the LMM for each vertice and extract the p-values of P30, P60, P120 and P150
lmerFun <- function(vertices, index) {
  
  Modl <- lmer(metric ~ age + grp + age*grp + (1 | animalID),
               data = vertices)
  
  emms <- emmeans(Modl, ~grp | age) %>% pairs(., adjust = "tukey")
  
  pvals <- summary(emms)$p.value[index] # the index will iterate the p-values of each temporal time point
  
}

## Select the number of permutations
num_perm <- 5000

## Create a matrix to store results
streams <- unique(mrds$stream)
points <- unique(mrds$point)

## Assigning the index for p-value extraction
index <- 1:4

## Create a multidimensional matrix to store the results
results <- array(NA, dim = c(length(streams), length(points), num_perm, length(index)))

## Assign the points and stream/gridlines as a column names in the matrix
colnames(results) <- points
rownames(results) <- streams

## Perform the nested loop to iterate through each set pf stream/gridlines and points
for (s in 1:length(streams)) {
  for (p in 1:length(points)) {
    
    # Subset the data
    vertices <- subset(mrds, stream == streams[s] & point == points[p])
    
    # Loop through permutations
    for (nperm in 1:num_perm) {
      # Shuffle the response variable
      vertices$metric <- sample(vertices$metric)
      # Perform LMM and extract p-values
      lmerTest <- lmerFun(vertices, index)
      # Store the results
      results[s, p, nperm, ] <- lmerTest
      
    }
    
  }
  
}

## Extract the results from the dimensions
dim1 <- results[, , , 1] # P30
dim2 <- results[, , , 2] #P60
dim3 <- results[, , , 3] #P120
dim4 <- results[, , , 4] #150

## Convert the matrices to dataframes
FApar_5kperm_30 <- as.data.frame(as.table(dim1))
FApar_5kperm_60 <- as.data.frame(as.table(dim2))
FApar_5kperm_120 <- as.data.frame(as.table(dim3))
FApar_5kperm_150 <- as.data.frame(as.table(dim4))



