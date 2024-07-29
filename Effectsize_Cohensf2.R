## Let's calculate Cohens f2 from our LMM results
library(effectsize)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(effects)


###### COHEN F2 FOR A SINGLE POINT

pvals <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/pvals_from_LMM_dti_r")

mrds_data <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/whole_data_mrdsR")
mrds_data

## getting DTI data ready ###
dti_FA <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/dti_FA_r")
dti_MD <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/dti_MD_r")
dti_RD <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/dti_RD_r")
dti_AD <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/dti_AD_r")

dti_FA$metric <- "FA"
dti_AD$metric <- "AD"
dti_MD$metric <- "MD"
dti_RD$metric <- "RD"

dti_FA <- dti_FA %>% rename(values = FA)
dti_MD <- dti_MD %>% rename(values = MD)
dti_AD <- dti_AD %>% rename(values = AD)
dti_RD <- dti_RD %>% rename(values = RD)

dti_data <- rbind(dti_FA, dti_MD, dti_AD, dti_RD)


##-------------------------###

dt_filter <- dti_data %>% filter(point == "1" & stream == "1" & metric == "FA") %>% pivot_wider(., names_from = metric, values_from = values)
dt_filter$edad <- factor(dt_filter$edad, levels = c("30", "60", "120","150"))

data_filter <- mrds_data %>% filter(point == "1" & stream == "1") 
data_filter$edad <- factor(data_filter$edad, levels = c("30", "60", "120","150"))

#fiteamos el modelo mixto
myMod <- lmer(FA ~ edad + grp + edad*grp + (1 | animalID),
              data = dt_filter)
summary(myMod)

t <- contrast(emmeans(myMod, ~grp | edad), method = "pairwise")
t
tratio <- summary(t)$t.ratio[1]
tratio
dferr <- summary(t)$df[1]
dferr

summary(t)$edad[1]

f2 <- t_to_f2(t = tratio,
        df_error = dferr,
        paired = T)
f2

f2$Cohens_f2_partial

###### Tes the function with one data point

lmer_un_punto <- function(vertices) {
  
  metrics <- c("FApar", "FAperp", "MDpar", "MDperp")
  
  results <- data.frame(variable = character(),
                        t_value = numeric(),
                        p_value = numeric(),
                        stream = numeric(),
                        point = numeric(),
                        edad = numeric(),
                        stringsAsFactors = FALSE)
  
  stream <- unique(vertices$stream)
  point <- unique(vertices$point)
  
  for (m in metrics) {
    
    form <- as.formula(paste(m, "~ edad + grp + edad*grp + (1 | animalID)"))
    
    myMod <- lmer(form, data = vertices)
    #myMod <- lmer(m, paste("~ edad + grp + edad*grp + (1 | animalID)"), data = vertices)
    # Perform contrasts for each model
    t <- contrast(emmeans(myMod, ~grp | edad), method = "pairwise")
    
    # Extract t-ratio and df
    tratio <- summary(t)$t.ratio[1]
    dferr <- summary(t)$df[1]
    
    # Add Edad
    edad <- summary(t)$edad[1]
    
    # Calculate Cohen's f^2
    f2_value <- t_to_f2(t = tratio, df_error = dferr, paired = TRUE)
    
    # Store results in the data frame
    results <- rbind(results, data.frame(
      variable = m,
      edad = edad,
      point = point,
      stream = stream,
      t_value = tratio,
      p_value = summary(t)$p.value[1],
      lmer = f2_value
    ))
  }
  return(results)
}

lmer_un_punto(data_filter)

##-------------------------------------------------------------------------------------------------##
###### PERFORMING COHEN F2 FOR ALL STREAM/POINT PER METRIC AT EACH EDAD

mrds_data$edad <- factor(mrds_data$edad, levels = c("30", "60", "120", "150"))

lmer_w_f2 <- function(vertice) {
  metrics <- c("FApar", "FAperp", "MDpar", "MDperp")
  
  f2_values <- numeric(length(metrics))
  
  for (i in seq_along(metrics)) {
    m <- metrics[i]
    form <- as.formula(paste(m, "~ edad + grp + edad*grp + (1 | animalID)"))
    myMod <- lmer(form, data = vertice)
    t <- contrast(emmeans(myMod, ~grp | edad), method = "pairwise")
    tratio <- summary(t)$t.ratio[2]
    dferr <- summary(t)$df[2]
    f2_value <- t_to_f2(t = tratio, df_error = dferr, paired = TRUE)
    f2_values[i] <- f2_value$Cohens_f2_partial
  }
  
  return(f2_values)
}

streams <- unique(mrds_data$stream)
points <- unique(mrds_data$point)

results_list <- list()

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



data <- as.data.frame(do.call(rbind, results_list))
data

data <- tibble::rownames_to_column(data, "row")

data <- unnest(data, cols = 1:5)


Cohens_f2_results_mrds60 <- data


getwd()
setwd("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/")
write.csv(Cohens_f2_results_mrds60, "Cohens_f2_results_mrds60")


######---------------------------------------------------DTI TEST




####-------------------------------------------------------------------------------####


## in the meanwhile, trying only 30 days for FApar and deal first with ggplot2


dti_FA$edad <- factor(dti_FA$edad, levels = c("30", "60", "120","150"))

fun_single_f2 <- function(vertice) {
  
  myMod <- lmer(FA ~ edad + grp + edad*grp + (1 | animalID), data = vertice)
  
  # Perform contrasts for each model
  t <- contrast(emmeans(myMod, ~grp | edad), method = "pairwise")
  # Extract t-ratio and df
  tratio <- summary(t)$t.ratio[2]
  dferr <- summary(t)$df[2]
  
  # Calculate Cohen's f^2
  f2_value <- t_to_f2(t = tratio, df_error = dferr, paired = TRUE)
  
  f2 <- f2_value$Cohens_f2_partial
  
  return(f2)
  
  return(f2_value)
}
  
#fun_single_f2(data_filter)
#head(data_filter)

streams <- unique(dti_FA$stream)
points <- unique(dti_FA$point) 

results_f02 <- matrix(NA, nrow = length(streams), ncol = length(points))

for (s in 1:length(streams)) {
  
  for (p in 1:length(points)) {
    
    stream <- streams[s]
    point <- points[p]
    
    vertice <- subset(dti_FA, stream == streams[s] & point == points[p])
    
    f2_test <- fun_single_f2(vertice)
    
    results_f02[s,p] <- f2_test

}
}
results_f02


###################### cus we failed to merge the dti data frames with different length. 
#We compute the [...] in one go so it takes less time 


dti_RD$edad <- factor(dti_RD$edad, levels = c("30", "60", "120", "150"))

fun_single_f2 <- function(vertice, indices) {
  
  myMod <- lmer(RD ~ edad + grp + edad * grp + (1 | animalID), data = vertice)
  
  # Perform contrasts for each model
  t <- contrast(emmeans(myMod, ~ grp | edad), method = "pairwise")
  # Extract t-ratio and df
  t_summary <- summary(t)
  tratio <- t_summary$t.ratio[indices]
  dferr <- t_summary$df[indices]
  
  # Calculate Cohen's f^2
  f2_values <- lapply(seq_along(indices), function(i) {
    f2_value <- t_to_f2(t = tratio[i], df_error = dferr[i], paired = TRUE)
    return(f2_value$Cohens_f2_partial)
  })
  
  return(list(tratio = tratio, dferr = dferr, f2 = f2_values))
}

streams <- unique(dti_RD$stream)
points <- unique(dti_RD$point)

indices <- 1:4  # Specify the indices you want to extract

results_f2 <- array(NA, dim = c(length(streams), length(points), length(indices)))

for (s in 1:length(streams)) {
  for (p in 1:length(points)) {
    stream <- streams[s]
    point <- points[p]
    vertice <- subset(dti_RD, stream == streams[s] & point == points[p])
    f2_test <- fun_single_f2(vertice, indices)
    results_f2[s, p, ] <- unlist(f2_test$f2)
  }
}


RD_cohens <- results_f2



dim1 <- RD_cohens[, , 1]
R_30 <- as.data.frame(dim1)

R_30$edad = "30"
R_30$stream = rep(1:50)
R_f2_30 <- pivot_longer(R_30, cols = 1:10, names_to = "point", values_to = "RD_f2_30" )


df_names <- ls(pattern = "_f2_")
df_names

dfs <- mget(df_names)

dit_data <- do.call(cbind, dfs)

# Get the names of the data frames in the global environment
df_names <- ls(pattern = "_f2_")  # Adjust the pattern as needed

# Retrieve the data frames from the global environment
dfs <- mget(df_names)

# Combine the data frames column-wise while preserving column names
combined_df <- dfs[[1]]  # Initialize combined dataframe with the first dataframe
for (i in 2:length(dfs)) {
  combined_df <- cbind(combined_df, dfs[[i]])
  colnames(combined_df)[(ncol(combined_df) - ncol(dfs[[i]]) + 1):ncol(combined_df)] <- colnames(dfs[[i]])
}


dti_data_f2 <- combined_df
getwd()
write.csv(dti_data_f2, "dti_data_f2")


dti_data_f2 <- cbind(new_mix[3:8], combined_df)

dti_data_f2 <- dti_data_f2[, -c(18:20, 22:24, 26:28, 30:32, 34:36)]

## combining the DTI pvals df with the f2 dataframe

pvals_dti <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/pvals_from_LMM_dti_r")

dti_f2_pval <- cbind(pvals_dti, dti_data_f2[7:22])  



#######################
results_f2

pvals <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/pvals_from_LMM_dti_r")
head(pvals)
new_table <- select(pvals, V1, V2, V3, stream, points, FA_EMM_30)
new_table

f2_df <- as.data.frame(results_f2)
f2_df$streams <- rep(1:50)
f2_df_pivot <- pivot_longer(f2_df, cols = 1:10, 
                            names_to = "point", 
                            values_to = "Cohens_f2")
f2_df
f2_df_pivot$point <- sub("V", "", f2_df_pivot$point)
f2_df_pivot

FA2plot <- cbind(new_table, f2_df_pivot)


head(FA2plot)


###########---------------PREPARING FOR PRE-PLOTTING-----------------#############

cohens30 <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/Cohens_f2_results_mrds30")  
cohens60 <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/Cohens_f2_results_mrds60")  
cohens120 <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/Cohens_f2_results_mrds120")  
cohens150 <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/Cohens_f2_results_mrds150")  

##merging dataframes with coordenates 
lmer <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/pvals_from_LMM_mrdsR")
lmer


head(cohens30)
f230 <- pivot_wider(cohens30, id_cols = 2:4, names_from = "metrics", values_from = "f2_values")
f260 <- pivot_wider(cohens60, id_cols = 2:4, names_from = "metrics", values_from = "f2_values")
f2120 <- pivot_wider(cohens120, id_cols = 2:4, names_from = "metrics", values_from = "f2_values")
f2150 <- pivot_wider(cohens150, id_cols = 2:4, names_from = "metrics", values_from = "f2_values")

## adding a suffix
orig_colnames <- colnames(f230)[4:7]
colnames(f230)[4:7] <- paste0(orig_colnames, "_f2_30")
colnames(f260)[4:7] <- paste0(orig_colnames, "_f2_60")
colnames(f2120)[4:7] <- paste0(orig_colnames, "_f2_120")
colnames(f2150)[4:7] <- paste0(orig_colnames, "_f2_150")


whole_f2 <- cbind(f230, f260, f2120, f2150)
whole_f2 <- whole_f2[, -c(8,9,10,15,16,17,22,23,24)]


# Store the last ten rows in a separate dataframe
last_ten_rows <- whole_f2[(nrow(whole_f2) - 9):nrow(whole_f2), ]
# Remove the last ten rows from the original dataframe
whole_f2 <- whole_f2[1:(nrow(whole_f2) - 10), ]
# Concatenate the two dataframes with the last ten rows first
whole_f2 <- rbind(last_ten_rows, whole_f2)

mix <- cbind(lmer, whole_f2)

mix <- mix[, -c(24,25,26)]

getwd()
write.csv(mix, "Plot_data_mix_2paper")

###### ------------------- Plotting 

library(ggplot2)
library(scales)
library(viridis)
library(paletteer)


##### ---- PLOTING ALL ---- ####

getwd()
setwd("/Users/ASUS/OneDrive/Escritorio/")

dtif2 <- read.csv("/Users/ASUS/OneDrive/Escritorio/Paper2/dti_data_f2")
dtipval <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/pvals_from_LMM_dti_r")
dti_f2_pval <- cbind(dtipval, dtif2[8:23])
getwd()
write.csv(dti_f2_pval, "dti_f2_pval")

dti_f2_pval <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/PermLMER_results/dti_f2_pval")

new_mix <- dti_f2_pval

## plotting
my_datapval <- new_mix[, colnames(new_mix[9:24])]
my_dataf2 <- new_mix[, colnames(new_mix[25:40])]


columns_pval <- colnames(my_datapval)
columns_f2 <- colnames(my_dataf2) 
columns_f2
columns_pval   

plot_list <- list()
  
   
   # Iterate over column indices
   for (i in seq_along(columns_f2)) {
     # Get column names based on index
     cols <- columns_f2[i]
     colval <- columns_pval[i]
     
     p <- ggplot(data = new_mix) +  
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
       
       ggtitle(paste0("cohensf2_", cols)) 
     
     # Assign ggtitle and save plot
     plot_results <- p
     plot_list[[cols]] <- plot_results
     svg(paste0("Grays_f2andPval_", cols, ".svg"), width = 8, height = 6)
     print(plot_list[[cols]])
     dev.off()
   }
  
  ############# plotting with the GRAY points and color the pvals only #######
   
   
  #p <-
  ggplot(data = new_mix, 
          mapping= aes(x = V1, y = V3, 
                       fill = RD_EMM_60, 
                       color = ifelse(RD_EMM_60 < 0.051, NA, RD_EMM_60), 
                       size = RD_f2_60)) +
     
     geom_point(stroke = 1.5, alpha = .9) +
     
     scale_color_gradient2(low = "#FDE333FF", 
                           mid = "#F69822FF",  
                           high = "#DD4A60FF", 
                           midpoint = 0.03) +
     
    scale_size_area(breaks = c(0.005, 0.01, 0.03, 0.07)) +
    scale_size(range = c(1, 4)) +
     
     theme_minimal() +
     theme(panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.line = element_line(color = "white"),
           plot.background = element_rect(fill = "black"),
           legend.title = element_text(color = "white", size=4),
           legend.text = element_text(color = "white"),
           axis.title = element_text(color = "white"),
           title = element_text(color = "white")) +
  ggtitle("cohensf2")
   
svg("Grays_f2andPval_RD_f2_60", width = 8, height = 6)
print(p)
dev.off()   


  
###-----------Plott regressions predictors-----------------###
   
mrds <- read.csv("/Users/ASUS/paulinav/MAESTRIA/R_dysplasias/Rdata/results_lmer/whole_data_mrdsR")   
   
d_regrss <- mrds %>% filter(point == "4" & stream == "32")

d_regrss$edad <- factor(d_regrss$edad, levels = c("30", "60", "120", "150"))

d_regrss

myMod <- lmer(FApar ~ edad + grp + edad*grp + (1 | animalID),
      data = d_regrss)
emm <- emmeans(myMod, ~grp | edad) %>% pairs(., adjust = "tukey")
emm

#predict <-
  
  plot(predictorEffect("edad", myMod),
    lines=list(multiline=T, col = c("#FFB90F", "#EEE9E9"), lwd=3),
   # lines=list(multiline=T, col = c("#CC0066", "#660099", "#336666", "#CC6600"), lwd=3),
     axes=list(transform = "trans"),
     confint=list(style="auto"), alpha=1)

  
  
setwd("/Users/ASUS/Desktop/")
svg("Pred_AD", width = 8, height = 6)
print(predict)
dev.off()
  
  
  write.csv(new_mix, "mrds_plot_data")







