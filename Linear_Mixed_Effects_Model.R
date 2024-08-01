" 
    In the first part of the code, we will run a Linear Mixed Effect Model (LMM) for each vertice. Then, the second part
    include mostly the same code but including the permutations for each vertice. In general, the framework is to shuffle
    the labels and re-run the LMM, and so on for the entire number of permutations. Then, the resulting matrices will be used 
    to calculate the cluster inference (see the code named Cluster_inference.R)
    
   IMPORTANT NOTE:

    Currently, the code is designed to analyze one metric at a time. To select the metric you want to measure, 
    you need to modify the variable on line 40. The code can be further improved by implementing a loop to iterate 
    through all metrics automatically, which is a pending task.

    This script will perform the LMM for a single metric across all temporal time points simultaneously. 
    As a result, you will obtain four dataframes (P30, P60, P120, and P150) with their corresponding p-values.
    
    
  
    Please keep in mind that the second part of the code (permutations) is computationally intensive and time-consuming. 
    It is recommended to run it in the background.
        
"

## Loading the required packages
library(dplyr)
library(tidyverse)
library(lme4)
library(pbkrtest)
library(effects)
library(emmeans)
library(readr)


" Here we first run the LMM without permutation. 
This is the results we will use as a baseline for the cluster inference"

## Converting the temporal time points to factors
mrds$age <- factor(mrds$age, levels = c("30", "60", "120", "150"))

## Select the metric to test. 
metric <- mrds$FApar # using FApar as an example. But it could be either: FApar, FAperp, FAparpar, FAparperp, FA, MD, RD, AD)

## Use a function to run the LMM for each vertice and extract the p-values of P30, P60, P120 and P150
lmerFun <- function(vertices, index) {
  
  Modl <- lmer(metric ~ age + grp + age*grp + (1 | animalID),
               data = vertices)
  
  emms <- emmeans(Modl, ~grp | age) %>% pairs(., adjust = "tukey")
  
  pvals <- summary(emms)$p.value[index] # the index will iterate the p-values of each temporal time point
  
}


## Make stream/gridlines and points as an unique index
streams <- unique(mrds$stream)
points <- unique(mrds$point)

## Create a multidimensional matrix to store the results, where columns are points, rows are streams/gridlines and 
## dimension z are the [index] which are the pvalues of each temporal time point P30, P60, P120 and P150
results_lmer <- array(NA, dim = c(length(streams), length(points), length(index)))

##  Assign names to columns and rows
colnames(results_lmer) <- points
rownames(results_lmer) <- streams

## Start the iterations through each stream/gridline and point (and age/index)
for (s in 1:length(streams)) {
  for (p in 1:length(points)) {
    
    stream <- streams[i]
    point <- points[j]
    vertices <- mrds[mrds$stream == streams[i] & mrds$point == points[j], ]
    
      for (i in 1:index) {
        
        single_lmer <- lmerFun(vertices, index)
        
        results_lmer[s, p, i, ] <- single_lmer
        
      }
   }
}


## Extract the results from the dimensions 
dim1 <- results[, , , 1]
dim2 <- results[, , , 2]
dim3 <- results[, , , 3]
dim4 <- results[, , , 4]

## Convert the matrices to dataframes
FApar_30 <- as.data.frame(as.table(dim1))
FApar_60 <- as.data.frame(as.table(dim2))
FApar_120 <- as.data.frame(as.table(dim3))
FApar_150 <- as.data.frame(as.table(dim4))


"Here starts the code to run the permutations"


set.seed(0002)


## Select the number of permutations
num_perm <- 5000

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
      vertices$FApar <- sample(vertices$FApar)
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



