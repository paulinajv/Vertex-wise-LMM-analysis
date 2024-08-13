" This code contains the steps we follow to compute the four connectivity and
  cluster inference.

  Here we will process one metric at the time. To change the metric modify the 
  line 34."

## Install the package
install.packages("BiocManager") 
BiocManager::install("EBImage")

## Load the library
library(EBImage)


#list all files from the current directory
list.files()
#crate a list from these files 
list_files <- list.files()
list_files
#create an empy list that will serve as a container to the recieving files
files_all <- list()
files_all
#a loop to read all files
for (i in 1:length(list_files)) {
  full_path <- list_files[i]
  files_all[[i]] <- read.table(full_path)
}

names(files_all) <- list_files
files_all


## Load the permuted LMM data
perm_data <- read.csv("MDpar_1kperm_150") # Change this to the desire metric and age
perm_data

## Change the colnames
names(perm_data) <- c("X", "stream", "point", "df", "pvalue")

## Convert to binary the significant pvalues
threshold = 0.051
perm_data$binary <- ifelse(perm_data$pvalue < threshold, 1, 0)
perm_data

## Create a new matrix for each dataframe
group_wise <- function(group_id) {
  
  
  group_d <- subset(perm_data, df == group_id)
  
  matrices <- matrix(group_d$binary,
                     nrow = length(unique(group_d$stream)), 
                     ncol = length(unique(group_d$point)))
  
  
  i <- Image(data = matrices, dim = c(50,10))
  labl <- bwlabel(i)
  totalconnected <- max(labl)
  tableconn <- table(labl[labl = 0])
  table_df <- as.data.frame(tableconn)
  names(table_df) <- c("Cluster_Num", "Connected")
  
  table_df$Group_ID <- group_id
  
  return(table_df)
  
}


## Create an empty data frame to store the results
clust_results <- data.frame()

## Iterate over unique group IDs in our data frame
  for (group_id in unique(perm_data$df)) {
  
  ## Process each group and append the result to the final dataframe
  clust_results <- rbind(clust_results, group_wise(group_id))
  
      }

## Print the final data frame with results from all groups
print(clust_results)


## Creating a new column with the labels for posterior plot
clust_density <- clust_results %>%
  mutate(label = case_when(
    Connected == 1 ~ "one",
    Connected == 2 ~ "two",
    Connected == 3 ~ "three",
    Connected == 4 ~ "four",
    Connected == 5 ~ "five",
    Connected == 6 ~ "six",
    TRUE ~ "no_connected"  
  ))

## Check the results
head(clust_density)

## Subseting in a new ddataframe without the binary 0. This is just for plotting
clust_analysis <- clust_density[clust_density$Connected < 400, ]
head(clust_analysis)

table(clust_analysis$Connected)

perm_data_final_results <- table(clust_analysis$Connected)


## Plotting a density map of cluster numbers

ggplot(data = clust_analysis, 
       aes(x=Connected, 
           fill=label)) + 
  geom_density(adjust=1.5, alpha=.4)











