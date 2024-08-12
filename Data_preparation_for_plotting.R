" This code is a mid-step to prepare the data for plotting. 

  Here we aimed to:

  a) Get the gridlines coordinates (taken from <tckconvert>)  into our enviroment and extract the 
     stream and point variables based on the filename.

  b) Merge the dataframe <lmer_and_f2_df> with the gridlines coordinates as a single
     dataframe. This merge is based on matching the stream and points observations between 
     the two datasets. "


## Load the libraries
library(tidyverse)
library(dplyr)
library(readr)


lmer_and_f2_df ## This is the final output from the LMM_f2_analysis.R 

## Load my gridlines data
gridlines <-  read_csv("/Users/paulinav/Downloads/Gridlines_coordinates.txt")

## Create a new variable stream and point based from the filename
gridlines$stream <- substring(gridlines$filename, first = 16, last = 17)
gridlines$point <- substring(gridlines$filename, first = 23, last = 24)

## Removing the first zeros of stream column and shift the start count from 0 to 1
gridlines$stream <- ifelse(grepl("^0[1-9]$", gridlines$stream), 
                           sub("^0", "", gridlines$stream), 
                           gridlines$stream) %>% as.numeric() +1


## Merge the analysis results to the coordinates by matching the stream and point columns

# First, lest convert all variables as factors
lmer_and_f2_df[, c("stream", "point")] <- lapply(lmer_and_f2_df[, c("stream", "point")], as.factor)
gridlines[, c("stream", "point")] <- lapply(gridlines[, c("stream", "point")], as.factor)

# Second, merge the two dataframes according to the stream and points references
data_to_plot <- left_join(gridlines, lmer_and_f2_df, by = c("stream", "point"))

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 