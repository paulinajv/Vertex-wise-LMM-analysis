"This script contains two different plotting sections:

      1) Visualization of hemisferes highlithing the p-value and effect size per vertice

      2) Visualization of each vertice the predictor effect from the LMM"


## Load the libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(pbkrtest)
library(effects)
library(emmeans)
library(gridExtra)
library(readr)
library(scales)
library(viridis)
library(paletteer)


"Section 1: vertex-wise visualization" 

## Load the data
mrds_data 
cohens_data 

## Merge the two dataframes into one for data handling
mrds_f2_pval <- cbind(mrds_data, cohens_data)

## Subsetting the data into two categories 
my_datapval <- mrds_f2_pval[, colnames(mrds_f2_pval[9:24])] 
my_dataf2 <- mrds_f2_pval[, colnames(mrds_f2_pval[25:40])]

## Assing the colnames of each category
columns_pval <- colnames(my_datapval)
columns_f2 <- colnames(my_dataf2) 

## Create an empty list to store the results
plot_list <- list()


## Iterate over column indices
for (i in seq_along(columns_f2)) {
  ## Get column names based on index
  cols <- columns_f2[i]
  colval <- columns_pval[i]
  
  p <- ggplot(data = mrds_f2_pval) +  
    geom_point(mapping = aes(x = V1, y = V3,
                             color = ifelse(get(colval) > 0.051, NA, get(colval)),
                             fill = get(colval),
                             size = get(cols)),
               stroke = 1.5, 
               alpha = .8) +
    
    scale_size_area(breaks = c(0,0.03,0.06,0.1,0.12,0.15,0.18)) +
    scale_size(range = c(1,4)) +
    
    scale_color_gradient2(low = "#FDE333FF", 
                          mid = "#F69822FF",  
                          high = "#DD4A60FF", 
                          midpoint = 0.03) +
    theme_minimal() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "white"),
          plot.background = element_rect(fill = "black"),
          legend.title = element_text(color = "white", size=4),
          legend.text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          title = element_text(color = "white")
    ) +
    
    ggtitle(paste0("Vertex-wise analysis", cols)) 
  
## Assign ggtitle and save plot (.svg) into current directory
  plot_results <- p
  plot_list[[cols]] <- plot_results
  svg(paste0("Grays_f2andPval_", cols, ".svg"), width = 8, height = 6)
  print(plot_list[[cols]])
  dev.off()
}


"Section 2: single predictor effect

      note: here we only perform a single predictor effect plot after
      choosing the vertex from the PDF that will represent the data"


## We first generate a PDF with all the predictors (in total 500 plots) 

## Load the data
data

## Function to call the command "predictorEffect"
lmerFun <- function(vertices) {
  ## Change de metric
  Modl <- lmer(FApar ~ age + grp + age*grp + (1 | animalID),
               data = vertices)
  
  p <- plot(predictorEffect("age", Modl),
            lines = list(multiline = TRUE, col = c("#339999", "#333333"), lwd = 3),
            axes = list(transform = "trans"),
            confint = list(style = "auto"), alpha = 1,
            main = paste("Stream", stream, "- Point", point)
  )
}

## Make unique values the stream/gridline and points
streams <- unique(data$stream)
points <- unique(data$point)

pdf(file = "./")
# Nested loop to perform the permutation test
for (i in 1:length(streams)) {
  for (j in 1:length(points)) {
    
    stream <- streams[i]
    point <- points[j]
    vertices <- data[data$stream == streams[i] & data$point == points[j], ]
    plots  <- lmerFun(vertices)
    grid.arrange(plots)
    
  }
}

## PDF file was saved into your current directory



## To plot a single predictor effect as a .svg file 

## Select the stream/gridline and point you want to plot
filter_data <- data %>% filter(point == "4" & stream == "32")

## Convert temporal time points into factors
filter_data$age <- factor(filter_data$age, levels = c("30", "60", "120", "150"))

## Fit the LMM of each metric
myMod <- lmer(FApar ~ age + grp + age*grp + (1 | animalID),
              data = filter_data)
## Calculate the post hocs
emm <- emmeans(myMod, ~grp | age) %>% pairs(., adjust = "tukey")

## Ready to plot the filter data
p <- plot(predictorEffect("age", myMod),
     lines=list(multiline=T, col = c("#FFB90F", "#EEE9E9"), lwd=3),
     axes=list(transform = "trans"),
     confint=list(style="auto"), alpha=1)


## Save the plot as a .svg file in your working directory
svg("predictor", width = 8, height = 6)
print(p)
dev.off()

















