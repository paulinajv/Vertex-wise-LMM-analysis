"This script contains two different plotting sections:

      1) Data handling
      
      2) Visualization of hemisferes highlithing the p-value and effect size per vertice

      3) Visualization of each vertice the predictor effect from the LMM"


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



" Section 1: data handling"

## Here we're going to modify the <data_to_plot> to be organized by columns instead of rows

## Load the final output of Data_preparation_for_plotting.R
data_to_plot <- read_csv("data_to_plot")

## Replacing the index number to the corresponding temporal time point (age)
data_to_plot$index <- replace(data_to_plot$index, data_to_plot$index %in% 1, 30)
data_to_plot$index <- replace(data_to_plot$index, data_to_plot$index %in% 2, 60)
data_to_plot$index <- replace(data_to_plot$index, data_to_plot$index %in% 3, 120)
data_to_plot$index <- replace(data_to_plot$index, data_to_plot$index %in% 4, 150)

## Changing the long format to wider (column-wise)
reshaped_data <- data_to_plot %>%
  pivot_wider(names_from = metrics, values_from = c(f2_values, pval)) %>%
  arrange(filename, V1, V2, V3, stream, point, index)

## Subseting the main dataframe for each index and adding it as a sufix in the column names
df30 <- reshaped_data %>% subset(index == 30) %>% 
  rename_with(~ paste0(., "_30"), .cols =8:15)
df60 <- reshaped_data %>% subset(index == 60) %>% 
  rename_with(~ paste0(., "_60"), .cols =8:15)
df120 <- reshaped_data %>% subset(index == 120) %>% 
  rename_with(~ paste0(., "_120"), .cols =8:15)
df150 <- reshaped_data %>% subset(index == 150) %>% 
  rename_with(~ paste0(., "_150"), .cols =8:15)

## Merge the dataframes by matching the other columns as reference
df <- df30 %>%
  left_join(df60, by = c("filename", "V1", "V2", "V3", "stream", "point")) %>%
  left_join(df120, by = c("filename", "V1", "V2", "V3", "stream", "point")) %>%
  left_join(df150, by = c("filename", "V1", "V2", "V3", "stream", "point"))



"Section 2: vertex-wise visualization" 

## Subset the data
my_datapval <- df %>% select(contains("pval"))
my_dataf2 <- df %>% select(contains("f2"))

## Assign the column names of each category
columns_pval <- colnames(my_datapval)
columns_f2 <- colnames(my_dataf2) 


## Create an empty list to store the results
plot_list <- list()


## Iterate over column indices
for (i in seq_along(columns_f2)) {
  ## Get column names based on index
  cols <- columns_f2[i]
  colval <- columns_pval[i]
  
  p <- ggplot(data = df) +  
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
  svg(paste0("Grays_f2_pval_", cols, ".svg"), width = 8, height = 6)
  print(plot_list[[cols]])
  dev.off()
}


"Section 3: vertex-wise predictor effect"


## We first generate a PDF with all the predictors (in total 500 plots) 

## Load the data
mrds <- read_csv("MRDS_database.txt")

## Function to call the command "predictorEffect"
lmerFun <- function(vertices) {
  ## Change de metric
  Modl <- lmer(FApar ~ age + group + age*group + (1 | animalID),
               data = vertices)
  
  p <- plot(predictorEffect("age", Modl),
            lines = list(multiline = TRUE, col = c("#339999", "#333333"), lwd = 3),
            axes = list(transform = "trans"),
            confint = list(style = "auto"), alpha = 1,
            main = paste("Stream", stream, "- Point", point)
  )
  return(p)
}

## Make unique values the stream/gridline and points
streams <- unique(mrds$stream)
points <- unique(mrds$point)

pdf(file = "plots.pdf")
# Nested loop to perform the permutation test
for (i in 1:length(streams)) {
  for (j in 1:length(points)) {
    
    stream <- streams[i]
    point <- points[j]
    vertices <- mrds[mrds$stream == streams[i] & mrds$point == points[j], ]
    plots  <- lmerFun(vertices)
    grid.arrange(plots)
    
  }
}

dev.off() ## PDF file was saved into your current directory


"note: here we only perform a single predictor effect plot after
      choosing the vertex from the PDF that will represent the data"

## To plot a single predictor effect as a .svg file 

## Select the stream/gridline and point you want to plot
filter_data <- mrds %>% filter(point == "4" & stream == "32")

## Convert temporal time points into factors
filter_data$age <- factor(filter_data$age, levels = c("30", "60", "120", "150"))

## Fit the LMM of each metric
myMod <- lmer(FApar ~ age + group + age*group + (1 | animalID),
              data = filter_data)
## Calculate the post hocs
emm <- emmeans(myMod, ~group | age) %>% pairs(., adjust = "tukey")

## Ready to plot the filter data
p <- plot(predictorEffect("age", myMod),
     lines=list(multiline=T, col = c("#FFB90F", "#EEE9E9"), lwd=3),
     axes=list(transform = "trans"),
     confint=list(style="auto"), alpha=1)


## Save the plot as a .svg file in your working directory
svg("predictor", width = 8, height = 6)
print(p)
dev.off()
getwd()

