# LOAD LIBRARIES ----

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)  
library(combiroc)
library(tibble)
library(gridExtra)
library(pROC)

# Normalise Data ------
#  --------------------- --------------------- ---------------------#

#MinMax Normalisation

minMax <- function(x) {
  res <- (x - min(x)) / (max(x) - min(x))
  return(res)}

data <- add_column(as.data.frame(sapply(data[,c(1:48)], minMax)), data[49])



# Plot codes :-----

### FIGURE 1: BOX AND WHISKER ####
#  --------------------- --------------------- ---------------------# 

# format data levels 
data$Group <- factor(data$Group, levels = c("TB", "ORD"))


# Convert data to long format
long_data <- data %>%
  pivot_longer(cols = 1:9, names_to = "Variable", values_to = "Normailised_MFI")

library(viridis)
scale_fill_viridis(discrete = TRUE)


png("Figure_1_Free_y_Facet_BoxPlot.png", width = 8, height = 8, units = "in", res = 1200)

# Varnames 
long_data

# Rename columns
colnames(long_data)

# Remove underscore in var names 
long_data$Variable
long_data$Variable = gsub("_", "-", long_data$Variable)
long_data$Variable = gsub("IFN-g", "IFN-γ", long_data$Variable)


ggplot(long_data, aes(x = Group, y = Normailised_MFI, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Remove outliers, match fill color
  facet_wrap(~ Variable, scales = "free_y") +  # Facet by variable
  stat_compare_means(
    method = "wilcox.test", 
    p.adjust.method = "fdr", 
    label = "p.format", 
    label.x = 1.5,  
    label.y = 0.75,   
    size = 4.1
  ) +  # Add FDR-adjusted p-values
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    legend.position = "none"  # Remove legend
  ) +
  scale_fill_manual(values = c("TB" = "#E69F00", "ORD" = "#56B4E9")) +  # Custom colors
  labs(title = NULL, x = NULL, y = "Normalised MFI")

dev.off()





### FIGURE 2: FACET ROC CYTOKINES ####
#  --------------------- --------------------- ---------------------# 
# range_TBdx
colnames(data) 
PlotCytokines = c(37, 31, 13, 15, 8, 10, 19, 36, 28, 49) # 49 Group var Group
# VARIABLES: mig ip10 16 6 2ra 1ra gama hgf  MIF

# Subset plot data 
data = range_TBdx[, PlotCytokines]
dim(data)
names(data)

# Remove underscore in var names
colnames(data)[1:9] = gsub("_", "-", colnames(data)[1:9])
colnames(data)[1:9] = gsub("IFN-g", "IFN-γ", colnames(data)[1:9])

# Create a list to store ROC curves
roc_list <- lapply(1:9, function(i) {
  roc_obj <- roc(data$Group, as.numeric(unlist(data[, i])))
  data.frame(
    Specificity = 1 - roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities,
    AUC = round(auc(roc_obj), 2),  # Compute AUC
    Variable = colnames(data[i])
  )
})

# Combine all ROC data into one dataframe
roc_df <- bind_rows(roc_list)

# Extract AUC values for annotation
auc_values <- roc_df %>%
  group_by(Variable) %>%
  summarise(AUC = unique(AUC))

# Gsub var names
names(roc_df)
unique(roc_df$Variable)

# Plot using ggplot2 with facets and AUC annotations

plott.c = ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Variable)) +
  geom_line(size = 1, color = "blue") +  # black ROC curve lines
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), 
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"), 
    strip.background = element_rect(fill = "white", color = "black"),  # white strip background
    strip.text = element_text(color = "black")  # black facet text
  ) +
  labs(title = NULL, x = "1 - Specificity", y = "Sensitivity") +
  geom_text(data = auc_values, aes(x = 0.6, y = 0.2, label = paste("AUC:", AUC)), 
            color = "black", size = 5, inherit.aes = FALSE)

png("Figure_2_Facet_roc_9_cytokines.png", width = 8, height = 8, units = "in", res = 1200)
plott.c
dev.off()



### TABLE 2: ROC details with CI ----
#  --------------------- --------------------- ---------------------# 

# part one: roc plots and roc objects -----

# i. FI Comparison Table

# try p value comparsion 

table(data$Group)

FIdata <- data[, c(2:50)]

names(FIdata)[49] = "Group"


FIdata <- subset(FIdata, FIdata$Group %in% c('0', '1'))
FIdata$Group <- gsub("0", "ORD", FIdata$Group)
FIdata$Group <- gsub("1", "TB", FIdata$Group)

library(dplyr)
library(tidyr)
library(knitr)

FIdata %>%
  gather(variable, value, -Group) %>%
  group_by(variable) %>%
  summarize(p_value = wilcox.test(value ~ Group)$p.value) %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "fdr")) -> result_table

# Print the table using kable
kable(result_table, format = "markdown")

write.csv(result_table, "FI_Table_TBvsORD_FI_comparison.csv")

# Summarise stats 

# Gather the data into long format for easier calculation
FIdata_long <- FIdata %>%
  pivot_longer(cols = -Group, names_to = "cytokine", values_to = "expression")

# Calculate median, lower and upper quartiles (Q1 and Q3) for each cytokine by Group
summary_stats <- FIdata_long %>%
  group_by(Group, cytokine) %>%
  summarise(
    median = median(expression, na.rm = TRUE),
    Q1 = quantile(expression, 0.25, na.rm = TRUE),
    Q3 = quantile(expression, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# View the result
summary_stats

write.csv(summary_stats, "summary_FI_values.csv")

colnames(summary_stats)

# wide 
# Convert the summary_stats to wide format, with 'Group' as separate columns for TB and ORD
summary_stats_wide <- summary_stats %>%
  pivot_wider(
    names_from = Group,       # Group (TB and ORD) will be the new column names
    values_from = c(median, Q1, Q3),  # The values for each cytokine will come from median, Q1, Q3
    names_glue = "{.value}_{Group}"   # This will create columns like 'median_TB', 'Q1_TB', 'Q3_TB' and 'median_ORD', 'Q1_ORD', 'Q3_ORD'
  )

# View the result
summary_stats_wide


# ii. FI ROC table ------
colnames(FIdata)

# Range normalised FIData 

library(pROC)


roclist <- c(1:48)
for (i in roclist){
  rocc = roc(FIdata$Group, as.numeric(unlist(FIdata[,i])))
  plot.roc(rocc, print.auc=TRUE, auc.polygon=TRUE, auc.polygon.col="white", main = colnames(FIdata[i]), print.thres=T)
  #plot.roc(smooth(rocc), add=TRUE, col="blue")
  # rocc_info = data.frame(rocc$thresholds, rocc$sensitivities, rocc$specificities)
  # head(rocc_info)
}
dev.off()


# Create an empty data frame to store AUC values
auc_data <- data.frame(Variable = character(), AUC = numeric(), stringsAsFactors = FALSE)


# Loop through each variable
roclist <- 1:48
for (i in roclist) {
  
  # Calculate ROC and get AUC
  rocc <- roc(FIdata$Group, as.numeric(unlist(FIdata[, i])))
  
  # Get AUC value
  auc_value <- as.numeric(auc(rocc))
  
  # Store AUC value in the data frame
  auc_data <- rbind(auc_data, data.frame(Variable = colnames(FIdata)[i], AUC = auc_value))
}
auc_data
# Close PDF file
dev.off()

# Print the AUC data frame
print(auc_data)

write.csv(auc_data, "FI_AUC_all.csv")

# Merged table ------
auc_data
summary_stats
result_table
summary_stats_wide
names(summary_stats_wide)[1] = "variable"
# 
# Load necessary library
library(dplyr)

# Merging the result_table with summary_stats_wide by 'cytokine'
merged_data1 <- result_table %>%
  left_join(summary_stats_wide, by = "variable")

# View the merged result
names(auc_data)[1] = "variable"
merged_data2 <- merged_data1 %>%
  left_join(auc_data, by = "variable")

# Write CSV
write.csv(merged_data2, "Data_FI_and_AUC_all.csv")


# iii. Add Confidence intervals AUC and SE and SP -----
library(tibble)
library(pROC)

minMax <- function(x) {
  res <- (x - min(x)) / (max(x) - min(x))
  return(res)}

rgFIdata <- add_column(as.data.frame(sapply(FIdata[,c(1:48)], minMax)), FIdata[49])
names(rgFIdata)


# Function to perform ROC analysis and calculate Youden's index, sensitivity, specificity, and AUC CI
roc_analysis <- function(cytokine, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$Group, as.numeric(data[[cytokine]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
  # Compute 95% CI for AUC
  auc_ci <- ci.auc(roc_curve, conf.level = 0.95)  # Bootstrapped CI for AUC
  
  # Compute 95% CI for sensitivity and specificity at Youden's cutoff
  ci_sensitivity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
  ci_specificity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)
  
  # Extract CI values
  sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
  sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
  specificity_lower_ci <- ci_specificity[["specificity"]][1]
  specificity_upper_ci <- ci_specificity[["specificity"]][3]
  
  # Store results in a dataframe
  results_df <- data.frame(
    Cytokine = cytokine,
    AUC = as.numeric(auc(roc_curve)),
    AUC_CI_Lower = auc_ci[1],
    AUC_CI_Upper = auc_ci[3],
    Youden_Cutoff = youden_cutoff,
    Sensitivity = sensitivity,
    Sensitivity_CI_Lower = sensitivity_lower_ci,
    Sensitivity_CI_Upper = sensitivity_upper_ci,
    Specificity = specificity,
    Specificity_CI_Lower = specificity_lower_ci,
    Specificity_CI_Upper = specificity_upper_ci
  )
  
  return(results_df)
}

# Option 1: Raw MFI ----
# List of cytokines (excluding the 'Group' column)
cytokines <- colnames(FIdata)[-which(colnames(FIdata) == "Group")]

# Initialize an empty list
results_listr <- list()

# Loop through each cytokine and perform ROC analysis
for (cytokine in cytokines) {
  results_listr[[cytokine]] <- roc_analysis(cytokine, FIdata)
}

# Combine all results into a single data frame
results_dfrr <- do.call(rbind, results_listr)

# Display the results dataframe
print(results_dfrr)

# Save results as a CSV file
write.csv(results_dfrr, "ROC_Details_RawMFI_withAUC_Cytokines.csv", row.names = FALSE)






# TABLE 3: All Gold combinations ####
#  --------------------- --------------------- ---------------------# 

# results report for specific markers/combinations
allgold <-roc_reports(data, markers_table = tab, case_class = 'A',
                      single_markers =c('MIG'), 
                      selected_combinations = c(34, 23, 119, 29, 114, 134, 309, 43, 153))

# results outputs
allgold$Plot
allgold$Metrics
allgold$Models

write.csv(allgold$Metrics, "Model_Goldcombinations_best_10.csv")

# show markers 
goldcombinations = show_markers(selected_combinations = c(34, 23, 119, 29, 114, 134, 309, 43, 153), markers_table = tab)

write.csv(goldcombinations, "GoldcombinationsModel_markers.csv")


### FIGURE 3A: BEST 3 Gold combinations ####
#  --------------------- --------------------- ---------------------# 
# plot of best 3 actually 
real3res = roc_reports(data, markers_table = tab, case_class = 'A',
                       single_markers =c('MIG'), 
                       selected_combinations = c(153, 309))
real3res$Plot
real3res$Metrics

# manage plot 
Roc_objectn3 = real3res[["Plot"]]
class(Roc_objectn3) # so we can maipulate the roc object 
# saveRDS(Roc_objectn3, "New3SignatureROC_Object.RDS")

# Plot using ggplot2 with facets and AUC annotations
png("Figure_3_BEST3Signatures_roc_cytokines.png", width = 5, height = 4, units = "in", res = 1200)
signature_plot = Roc_objectn3 +
  ggplot2::geom_line(size = 0.5) +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 panel.border = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 0.8, color = "black"),
                 axis.ticks = ggplot2::element_line(size = 0.8))

signature_plot
dev.off()


### FIGURE 3B: Variable importance bar plot #####
#  --------------------- --------------------- ---------------------# 

library(tidyverse)

# List of top 20 models
models <- c(
  "IFN_g-IL_1ra", "IL_16-MIG", "IL_1ra-MIG", "IL_2Ra-MIG", "M_CSF-MIG",
  "IFN_g-IL_16-IL_1ra", "IFN_g-IL_1ra-IL_2Ra", "IL_16-IL_1ra-MIG", "IL_16-IL_2Ra-MIG", "IL_16-M_CSF-MIG",
  "IL_1ra-IL_2Ra-MIG", "IL_1ra-M_CSF-MIG", "IL_2Ra-M_CSF-MIG",
  "IFN_g-IL_16-IL_1ra-IL_2Ra", "IL_16-IL_1ra-IL_2Ra-MIG", "IL_16-IL_1ra-M_CSF-MIG",
  "IL_16-IL_2Ra-M_CSF-MIG", "IL_1ra-IL_2Ra-M_CSF-MIG", "IL_16-IL_1ra-IL_2Ra-M_CSF-MIG"
)

# Convert to a data frame and split into separate variables
var_counts <- str_split(models, "-", simplify = TRUE) %>%
  as.vector() %>%
  table() %>%
  as.data.frame()

var_counts = var_counts[-1,]

# Rename columns
colnames(var_counts) <- c("Variable", "Frequency")

# Remove underscore in var names 
var_counts$Variable
var_counts$Variable = gsub("_", "-", var_counts$Variable)
var_counts$Variable = gsub("IFN-g", "IFN-γ", var_counts$Variable)

# Plot the bar chart
png("Figure 3b. Barplot of most common variables in top20 models.png", width = 8, height = 8, units = "in", res = 1200)

ggplot(var_counts, aes(x = reorder(Variable, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "black") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1),  # Combine size & angle,  
    axis.text.y = element_text(size = 14),  # Increase font size for y-axis values
    axis.title.x = element_text(size = 16), # Increase font size for x-axis label
    axis.title.y = element_text(size = 16),  # Increase font size for y-axis label
    panel.grid = element_blank()
  ) +
  labs(x = "Cytokine", y = "Frequency", 
       title = "Frequency of Variables in Top 20 Combiroc regression Models")
dev.off()








### Figure 4. Treatment response ####
#  --------------------- --------------------- ---------------------# 
# treatment  
range_TxTb.Ord
dtx = range_TxTb.Ord
names(range_TxTb.Ord)

table(range_TxTb.Ord$Timepoint)

dtx <- dtx[, c(15, 28, 37, 13, 33, 49)] # IL2ra,IL_16,MIG, timepoint 
summary(dtx)

# format table to remove month 1 
table(dtx$Timepoint)
dtx = subset(dtx, Timepoint %in% c("Baseline","Month_2","Month_4",  "Month_6","ORD"))
dtx$Timepoint = as.factor(as.character(dtx$Timepoint))
table(dtx$Timepoint)
summary(dtx$Timepoint)
names(dtx)


# Create a data frame in long format (if not already done)
new_dtt_long2 <- dtx %>%
  gather(key = "Cytokine", value = "Value", IL_2Ra, IL_16, MIG, IL_1ra, M_CSF)


# gsub
names(new_dtt_long2)
unique(new_dtt_long2$Cytokine)
new_dtt_long2$Cytokine = gsub("_", "-", new_dtt_long2$Cytokine)
new_dtt_long2$Timepoint = gsub("_", " ", new_dtt_long2$Timepoint)
new_dtt_long2$Timepoint = as.factor(new_dtt_long2$Timepoint)


# Define comparisons
my_comparisons <- list(
  c("Baseline", "Month 2"),
  c("Month 2", "Month 4"),
  c("Month 4", "Month 6"),
  c("Month 6", "ORD")
)
# Create a list to store the individual plots
plots <- list()

# Loop through the cytokines and create a separate plot for each
for (cytokine in unique(new_dtt_long2$Cytokine)) {
  
  # Filter the data for the current cytokine
  cytokine_data <- subset(new_dtt_long2, Cytokine == cytokine)
  
  # Set custom y-axis limits for each cytokine
  y_limits <- switch(cytokine,
                     "IL-2Ra" = c(0, 0.5),
                     "IL-16" = c(0, 0.5),
                     "MIG" = c(0, 0.6),
                     "IL-1ra" = c(0, 0.15),
                     "M-CSF"=c(0, 0.5))
  
  # Generate the plot
  q <- ggplot(cytokine_data, aes(x = Timepoint, y = Value, fill = NA)) +
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.6) +  # No outliers, black outline
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test", 
      p.adjust.method = "fdr",  # Adjust p-values using FDR
      label = "p.signif",  # Show asterisks for significant p-values
      size = 3,  # Font size for asterisks
      label.x = 0.5,  # X position of asterisks
      label.y = 0.4,  # Adjust the y position for p-value labels
      aes(group = Timepoint)  # Ensures comparisons are based on the correct grouping
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(size = 0.8),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black"),
      legend.position = "none",  # No legend
      axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
    ) +
    # scale_fill_manual(values = c("Baseline" = "red", "Month_2" = "blue", "Month_4" = "green", "Month_6" = "purple", "ORD" = "orange")) +  # Custom colors
    labs(title = paste(cytokine), x = NULL, y = "Normalised MFI Expression") +
    coord_cartesian(ylim = y_limits)  # Custom y-axis limits for each cytokine
  
  # Store the plot in the list
  plots[[cytokine]] <- q
}

# Combine the individual plots using grid.arrange with 2 columns
png("Figure4_Trends_Facet_BoxPlot_Cytokines_response to treatment.png", width = 8, height = 12, units = "in", res = 1200)

grid.arrange(plots[["IL-2Ra"]], plots[["IL-16"]], plots[["MIG"]], 
             plots[["IL-1ra"]], plots[["M-CSF"]], ncol = 2)
dev.off()
# Hide NS comparisons ------
# for now manually remove in excel of biorender to save time 

# Add color ----

# Create a list to store the individual plots
plots <- list()

# Loop through the cytokines and create a separate plot for each
for (cytokine in unique(new_dtt_long2$Cytokine)) {
  
  # Filter the data for the current cytokine
  cytokine_data <- subset(new_dtt_long2, Cytokine == cytokine)
  
  # Set custom y-axis limits for each cytokine
  y_limits <- switch(cytokine,
                     "IL-2Ra" = c(0, 0.6),
                     "IL-16" = c(0, 0.5),
                     "MIG" = c(0, 0.6),
                     "IL-1ra" = c(0, 0.6),
                     "M-CSF"=c(0, 0.6))

  
  q <- ggplot(cytokine_data, aes(x = Timepoint, y = Value, fill = Timepoint)) +  # Map fill to Timepoint
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.6) +  # Boxplot with black outlines
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test", 
      p.adjust.method = "fdr",  # Adjust p-values using FDR
      label = "p.format",  # Show asterisks for significant p-values
      size = 3,  # Font size for asterisks
      label.x = 0.5,  # X position of asterisks
      label.y = 0.5,  # Adjust y position for p-value labels
      aes(group = Timepoint)  # Ensures comparisons are based on the correct grouping
    ) + 
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(size = 0.8),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black"),
      legend.position = "none",  # No legend
      axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
    ) +
    scale_fill_manual(values = c(
      "Baseline" = "#E69F00",
      "Month 2" = "#E69F00",
      "Month 4" = "#E69F00",
      "Month 6" = "#E69F00",
      "ORD" = "#56B4E9"
    )) +  # Assign colors to Timepoints
    labs(title = paste(cytokine), x = NULL, y = "Normalised MFI Expression") +
    coord_cartesian(ylim = y_limits)
  
  
  # Store the plot in the list
  plots[[cytokine]] <- q
}

# Combine the individual plots using grid.arrange with 2 columns
png("Figure4_with_color.png", width = 8, height = 12, units = "in", res = 1200)

grid.arrange(plots[["IL-2Ra"]], plots[["IL-16"]], plots[["MIG"]], 
             plots[["IL-1ra"]], plots[["M-CSF"]], ncol = 2)
dev.off()



#  --------------------- --------------------- ---------------------# 
# Three cytokine option: 
# Combine the individual plots using grid.arrange with 2 columns
png("Figure4_three_cytokine_optimal_with_color_pvalue.png", width = 8, height = 8, units = "in", res = 1200)

grid.arrange(plots[["IL-2Ra"]], plots[["M-CSF"]], plots[["MIG"]], ncol = 2)
dev.off()







# TABLE 4: Combiroc tables  #####
#  --------------------- --------------------- ---------------------# 

# A. All signatures unranked -----
# Unranked combinations 
mks
write.csv(mks, "Allcombination_metric_unranked.csv")
#
PathoDat = as.data.frame(FullDx_Data[,c("Sample_ID", "Xpert_Class")])
dim(PathoDat)
write.csv(PathoDat, "SampleIDs_all.csv")

# Table:
library(readxl)
pathdemo <- read_excel("Demobraphy_of_pathogens.xlsx")
table(pathdemo$Pathogens)



# EXTRA Analyses

# 1. Chi square for pathogen ----

# Create a matrix for the observed counts
data_matrix <- matrix(c(33, 14, 
                        26, 9, 
                        86, 62), 
                      nrow = 3, byrow = TRUE)

# Add row and column names
rownames(data_matrix) <- c("Bacteria", "Viruses (+COVID19)", "No pathogen detected")
colnames(data_matrix) <- c("ORD", "TB")

# Perform the chi-square test
chisq_result <- chisq.test(data_matrix)

# Print the results
print(chisq_result)

#
# > print(chisq_result)
# 
# Pearson's Chi-squared test
# 
# data:  data_matrix
# X-squared = 4.4828, df = 2, p-value = 0.1063


# 2. IP 10 comaprisons between MBT and MFI lumniex data -----

# merge teh two datasets by Sample_ID
merged_IP <- merge(MBT_IP, MFI_IP, by = "Sample_ID", all = FALSE)


names(merged_IP)

# simple code 
library(pROC)

# ROC for IP_10
roc_ip10 <- roc(merged_IP$Xpert_Class, merged_IP$IP_10, ci = TRUE)
auc_ip10 <- roc_ip10$auc
ci_ip10 <- ci.auc(roc_ip10)

# ROC for IP.10
roc_ip.10 <- roc(merged_IP$Xpert_Class, merged_IP$IP.10, ci = TRUE)
auc_ip.10 <- roc_ip.10$auc
ci_ip.10 <- ci.auc(roc_ip.10)

# Print AUC with 95% CI
cat("AUC for IP_10:", auc_ip10, " (95% CI:", ci_ip10[1], "-", ci_ip10[3], ")\n")
cat("AUC for IP.10:", auc_ip.10, " (95% CI:", ci_ip.10[1], "-", ci_ip.10[3], ")\n")

# Plot ROC curves
png("Figure6_ROC_IP_10_in_MBT_vs_MFIplot.png", width = 800, height = 600)  # Open a PNG device
plot(roc_ip10, col = "blue", lwd = 2, main = "ROC Curves for IP_10 and IP.10") # MFI is blue 
plot(roc_ip.10, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c("IP_10", "IP.10"), col = c("blue", "red"), lwd = 2)
dev.off() 

#legend("bottomright", legend = c("MFI IP-10", "MBT IP-10"), col = c("blue", "red"), lwd = 2)

# Delong test
roc.test(roc_ip10, roc_ip.10, method = "delong")


# ROC meterics with 95% CI: -----

# Function to perform ROC analysis and calculate Youden's index, sensitivity, specificity, and AUC CI
IP_roc_analysis <- function(cytokine, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$Xpert_Class, as.numeric(data[[cytokine]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
  # Compute 95% CI for AUC
  auc_ci <- ci.auc(roc_curve, conf.level = 0.95)  # Bootstrapped CI for AUC
  
  # Compute 95% CI for sensitivity and specificity at Youden's cutoff
  ci_sensitivity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
  ci_specificity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)
  
  # Extract CI values
  sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
  sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
  specificity_lower_ci <- ci_specificity[["specificity"]][1]
  specificity_upper_ci <- ci_specificity[["specificity"]][3]
  
  # Store results in a dataframe
  results_df <- data.frame(
    Cytokine = cytokine,
    AUC = as.numeric(auc(roc_curve)),
    AUC_CI_Lower = auc_ci[1],
    AUC_CI_Upper = auc_ci[3],
    Youden_Cutoff = youden_cutoff,
    Sensitivity = sensitivity,
    Sensitivity_CI_Lower = sensitivity_lower_ci,
    Sensitivity_CI_Upper = sensitivity_upper_ci,
    Specificity = specificity,
    Specificity_CI_Lower = specificity_lower_ci,
    Specificity_CI_Upper = specificity_upper_ci
  )
  
  return(results_df)
}


# List of cytokines (excluding the 'Group' column)
names(merged_IP)
cytokines <- c("IP.10","IP_10")

# Initialize an empty list
results_ip <- list()

# Loop through each cytokine and perform ROC analysis
for (cytokine in cytokines) {
  results_ip[[cytokine]] <- IP_roc_analysis(cytokine, merged_IP)
}

# Combine all results into a single data frame
results_dfip <- do.call(rbind, results_ip)

# Display the results dataframe
print(results_dfip)

# Save results as a CSV file
write.csv(results_dfip, "ROC_Details_IP10_with_CI_mbt_vs_mfi.csv", row.names = FALSE)





# 3. Paired treatment response pattern -----
#  --------------------- --------------------- ---------------------# 

# bring participant ID 

minMax <- function(x) {
  res <- (x - min(x)) / (max(x) - min(x))
  return(res)}

names(Tx_Od_combine)
treatmentdata <- add_column(as.data.frame(sapply(Tx_Od_combine[,c(2:49)], minMax)), Tx_Od_combine[c(1,50,51)])

rtx = treatmentdata
names(treatmentdata)

table(rtx$Timepoint)

rtx <- rtx[, c(15, 28, 37, 13, 33, 49:51)] # IL2ra,IL_16,MIG, timepoint 
summary(rtx)

# format table to remove month 1 
table(rtx$Timepoint)
rtx = subset(rtx, Timepoint %in% c("Baseline","Month_2","Month_4",  "Month_6","ORD"))
rtx$Timepoint = as.factor(as.character(rtx$Timepoint))
table(rtx$Timepoint)
summary(rtx$Timepoint)
names(rtx)

# gsub timepoint 
rtx$Timepoint = gsub("_", " ", rtx$Timepoint)
rtx$Timepoint = as.factor(rtx$Timepoint)

TreatData = rtx

# Line options :for paired data 

required_timepoints <- c("Baseline", "Month 2", "Month 4", "Month 6")

# Find Participants with data at all required Timepoints


# check sample ID 
TreatData$Sample_ID
TreatData$Sample_ID =gsub("_BL", "", TreatData$Sample_ID)
TreatData$Sample_ID =gsub("_M2", "", TreatData$Sample_ID)
TreatData$Sample_ID =gsub("_M4", "", TreatData$Sample_ID)
TreatData$Sample_ID =gsub("_M6", "", TreatData$Sample_ID)
TreatData$Sample_ID =gsub("BL", "", TreatData$Sample_ID)


participants_with_all_timepoints <- TreatData %>%
  group_by(Sample_ID) %>%
  filter(Timepoint %in% required_timepoints) %>%
  summarize(num_timepoints = n_distinct(Timepoint)) %>%
  filter(num_timepoints == length(required_timepoints)) %>%
  pull(Sample_ID)

# Filter Binary_TB to include only those Participants with data at all Timepoints
filtered_TxData <- TreatData %>%
  filter(Sample_ID %in% participants_with_all_timepoints)

# Print the filtered dataset or use it for further analysis
dim(filtered_TxData)
table(filtered_TxData$Xpert_Class, filtered_TxData$Timepoint)


# N = 51
names(filtered_TxData)
filtered_TxData$Timepoint
table(filtered_TxData$Timepoint)
dim(filtered_TxData)

unique(filtered_TxData$Sample_ID)

filtered_TxData = subset(filtered_TxData, Timepoint%in% 
                           c("Baseline", "Month 2", "Month 4", "Month 6"))

filtered_TxData$Timepoint = as.factor(as.character(filtered_TxData$Timepoint))

filtered_TxData$Timepoint = factor(filtered_TxData$Timepoint,
                                   levels = c("Baseline", "Month 2", "Month 4", "Month 6"))

# vars 
# "IL_2Ra" "IL_16"  "MIG"    "IL_1ra" "M_CSF"
names(filtered_TxData)[1:5] = gsub("_", "-", names(filtered_TxData)[1:5])


# Define the pairs for comparison
my_comparisons <- list(c("Baseline", "Month 2"),
                       c("Month 2", "Month 4"),
                       c("Month 4", "Month 6"))

# Line pplot 
p1 <- ggpaired(filtered_TxData, x = "Timepoint", y = "IL_2Ra",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "IL-2Ra",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco",
               group.by = "Sample_ID") +  # Ensure pairing by Sample_ID
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)

p2 <- ggpaired(filtered_TxData, x = "Timepoint", y = "M_CSF",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "M_CSF",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco",
               group.by = "Sample_ID") +  # Ensure pairing by Sample_ID
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)


p3 <- ggpaired(filtered_TxData, x = "Timepoint", y = "MIG",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "MIG",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco",
               group.by = "Sample_ID") +  # Ensure pairing by Sample_ID
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)


# Arrange the four plots in a 2x2 grid
tiff("Figure_6new_Treatment_paired_n51.tiff", 
     width = 6, height = 12, res=300, unit="in")  # Start PDF device
paired = grid.arrange(p1, p2, p3, ncol = 1, nrow = 3)
paired
dev.off()



# Treatment response ----


##### 1. ROC for baseline vs follow up timepoints with 95% CI  -----


# Function to generate ROC plots for each cytokine with AUC CI
plot_roc_for_cytokine <- function(data, cytokine_name) {
  # Filter for the specific cytokine
  cytokine_data <- data %>% filter(Cytokine == cytokine_name)
  
  # Initialize a list to store ROC objects
  roc_list <- list()
  
  # Define comparisons: Baseline vs. each timepoint
  comparisons <- c("Month 2", "Month 4", "Month 6")
  
  # Generate ROC curves for each comparison
  for (timepoint in comparisons) {
    subset_data <- cytokine_data %>% filter(Timepoint %in% c("Baseline", timepoint))
    
    # Convert Timepoint to a binary factor (Baseline = 0, Other = 1)
    subset_data$BinaryTime <- ifelse(subset_data$Timepoint == "Baseline", 0, 1)
    
    # Compute ROC
    roc_obj <- roc(subset_data$BinaryTime, subset_data$Value)
    
    # Store the ROC object
    roc_list[[timepoint]] <- roc_obj
  }
  
  # Plot ROC curves
  plot(roc_list[[1]], col = "red", lwd = 2, main = paste("ROC for", cytokine_name))
  plot(roc_list[[2]], col = "blue", lwd = 2, add = TRUE)
  plot(roc_list[[3]], col = "green", lwd = 2, add = TRUE)
  
  # Compute CI for AUC
  auc_cis <- sapply(roc_list, function(roc_obj) ci(roc_obj, conf.level = 0.95))
  
  # Add legend with AUC values and CI
  legend("bottomright", 
         legend = paste0(comparisons, ": AUC = ", 
                         round(sapply(roc_list, auc), 2), 
                         " (", round(auc_cis[1, ], 2), "-", round(auc_cis[2, ], 2), ")"), 
         col = c("red", "blue", "green"), lwd = 2)
}

# Apply function to each cytokine
unique_cytokines <- unique(ResponseROC$Cytokine)
for (cyto in unique_cytokines) {
  plot_roc_for_cytokine(ResponseROC, cyto)
}
# all plots manually saved in plot-2025-tretament respnse in pdf



# 2. Treatmment response pattern in slow and fast ----
#  --------------------- --------------------- ---------------------# 

# Slow and Fast 

# Define the list of sample IDs that should be labeled as "slow"
TreatData$Sample_ID
slow_samples <- c("TRI0370", "TRI0569", "TRI0826", "TRI0869", "TRI0874", "TRI0972")

# Add the 'response' column based on Sample_ID
TreatResponse <- TreatData %>%
  mutate(response = ifelse(Sample_ID %in% slow_samples, "slow", "fast"))

# Check the result
table(TreatResponse$response, TreatResponse$Timepoint)

# remove ORD data 

TreatResponseTB = subset(TreatResponse, Timepoint%in% 
                           c("Baseline", "Month 2", "Month 4", "Month 6"))

TreatResponseTB$Timepoint = as.factor(as.character(TreatResponseTB$Timepoint))

TreatResponseTB$Timepoint = factor(TreatResponseTB$Timepoint,
                                   levels = c("Baseline", "Month 2", "Month 4", "Month 6"))

# Reshape the data for plotting
TreatResponseTB_long <- TreatResponseTB %>%
  gather(key = "Cytokine", value = "Expression", IL_2Ra, IL_16, MIG, IL_1ra, M_CSF) %>%
  group_by(Timepoint, response, Cytokine) %>%
  summarise(
    median_expr = median(Expression, na.rm = TRUE),
    IQR_lower = quantile(Expression, 0.25, na.rm = TRUE),
    IQR_upper = quantile(Expression, 0.75, na.rm = TRUE)
  )


# Plot the data
ggplot(TreatResponseTB_long, aes(x = Timepoint, y = median_expr, color = response, group = response)) +
  geom_line(aes(linetype = response), size = 1) +  # Line plot for each group
  geom_errorbar(aes(ymin = IQR_lower, ymax = IQR_upper), width = 0.1) +  # IQR as error bars
  facet_wrap(~ Cytokine, scales = "free_y") +  # Facets for each cytokine
  scale_color_manual(values = c("fast" = "green", "slow" = "red")) +  # Colors for lines
  theme_minimal() + 
  labs(
    title = "Cytokine Expression Trend by Timepoint and Response Group",
    x = "Timepoint",
    y = "Cytokine Expression (Median with IQR)",
    color = "Response Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(size = 0.5, linetype = "dashed")
  )


# Show plot for only IL2RA, M-CSF and MIG 

# Reshape the data for plotting
TreatResponseTB_long <- TreatResponseTB %>%
  gather(key = "Cytokine", value = "Expression", IL_2Ra, IL_16, MIG, IL_1ra, M_CSF) %>%
  group_by(Timepoint, response, Cytokine) %>%
  summarise(
    median_expr = median(Expression, na.rm = TRUE),
    IQR_lower = quantile(Expression, 0.25, na.rm = TRUE),
    IQR_upper = quantile(Expression, 0.75, na.rm = TRUE)
  )


# subset 
unique(TreatResponseTB_long$Cytokine)
TreatResponseTB_long2 = subset(TreatResponseTB_long,Cytokine %in% c("IL_2Ra","MIG","M_CSF"))

# Plot the data
ggplot(TreatResponseTB_long2, aes(x = Timepoint, y = median_expr, color = response, group = response)) +
  geom_line(aes(linetype = response), size = 1) +  # Line plot for each group
  geom_errorbar(aes(ymin = IQR_lower, ymax = IQR_upper), width = 0.1) +  # IQR as error bars
  facet_wrap(~ Cytokine, scales = "free_y") +  # Facets for each cytokine
  scale_color_manual(values = c("fast" = "green", "slow" = "red")) +  # Colors for lines
  theme_minimal() + 
  labs(
    title = "Cytokine Expression Trend by Timepoint and Response Group",
    x = "Timepoint",
    y = "Cytokine Expression (Median with IQR)",
    color = "Response Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(size = 0.5, linetype = "dashed")
  )



# print p value for (wilcox_results for comaprison of baseline vs follow up timepoints) ----

# Reshape the data to long format
TreatResponseTB_long <- TreatResponseTB %>%
  pivot_longer(cols = c("IL_2Ra", "IL_16", "MIG", "IL_1ra", "M_CSF"),
               names_to = "Cytokine",
               values_to = "Expression")

# Initialize an empty dataframe to store results
wilcox_results <- data.frame(Cytokine = character(),
                             Timepoint = character(),
                             p_value = numeric(),
                             test_statistic = numeric(),
                             stringsAsFactors = FALSE)

# Loop through each cytokine and timepoint to perform the Wilcoxon test
for (cytokine in unique(TreatResponseTB_long$Cytokine)) {
  
  for (timepoint in unique(TreatResponseTB_long$Timepoint)) {
    
    # Subset the data for the current cytokine and timepoint
    subset_data <- TreatResponseTB_long %>%
      filter(Cytokine == cytokine, Timepoint == timepoint) %>%
      drop_na()
    
    # Perform the Wilcoxon test comparing 'slow' and 'fast' responses
    if (nrow(subset_data) > 1) {
      test_result <- wilcox.test(Expression ~ response, data = subset_data)
      
      # Store results
      wilcox_results <- wilcox_results %>%
        add_row(
          Cytokine = cytokine,
          Timepoint = timepoint,
          p_value = test_result$p.value,
          test_statistic = test_result$statistic
        )
    }
  }
}

# View the results
print(wilcox_results)

# Write to csv
write.csv(wilcox_results, "TreatmentResponse_wilcox_results.csv")






# EXTRA Analyses : Reporting 95% CI for best combination models -----
#  --------------------- --------------------- ---------------------# 
real3res
View(real3res)
str(real3res)

Model1 = real3res[["Models"]][["Combination 153"]]
Model2 = real3res[["Models"]][["Combination 309"]]
Model3 = real3res[["Models"]][["MIG"]]


#-----------------------------------------------------------------------------------------# 
# MODEL 1 ------
#-----------------------------------------------------------------------------------------# 

# Model 1 ROC
library(pROC)

# Extract the model
Model1 = real3res[["Models"]][["Combination 153"]]

# Get predicted probabilities
pred_probs <- predict(Model1, type = "response")

# Get the actual outcome variable from the dataset
actual_outcome <- real3res[["Models"]][["Combination 153"]][["y"]]

# # Compute ROC curve
# Model1_roc <- roc(actual_outcome, pred_probs)
# Compute AUC with confidence interval
# Model1_auc_ci <- ci.auc(Model1_roc)  
# print(Model1_auc_ci)  # Print AUC CI

# # Plot the ROC curve
# plot(Model1_roc, col = "blue", main = "ROC Curve for Combination 153")
# auc(Model1_roc)  # Get AUC value

# Compute ROC curve
Model1_roc <- roc(actual_outcome, pred_probs, ci = TRUE)


# Model 1 details 

# Compute best threshold for sensitivity and specificity
best_coords <- coords(Model1_roc, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden", transpose = FALSE)

# Extract values for Youden's cutoff, sensitivity, and specificity
youden_cutoff <- as.numeric(best_coords["threshold"])
sensitivity <- as.numeric(best_coords["sensitivity"])
specificity <- as.numeric(best_coords["specificity"])

# Compute 95% CI for AUC
auc_ci <- ci.auc(Model1_roc, conf.level = 0.95)  # Bootstrapped CI for AUC

# Compute 95% CI for sensitivity and specificity at Youden's cutoff
ci_sensitivity <- ci.coords(Model1_roc, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
ci_specificity <- ci.coords(Model1_roc, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)

# Extract CI values
sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
specificity_lower_ci <- ci_specificity[["specificity"]][1]
specificity_upper_ci <- ci_specificity[["specificity"]][3]

# Store results in a dataframe
Model1_df <- data.frame(
  AUC = as.numeric(auc(Model1_roc)),
  AUC_CI_Lower = auc_ci[1],
  AUC_CI_Upper = auc_ci[3],
  Youden_Cutoff = youden_cutoff,
  Sensitivity = sensitivity,
  Sensitivity_CI_Lower = sensitivity_lower_ci,
  Sensitivity_CI_Upper = sensitivity_upper_ci,
  Specificity = specificity,
  Specificity_CI_Lower = specificity_lower_ci,
  Specificity_CI_Upper = specificity_upper_ci
)

print(Model1_df)




#-----------------------------------------------------------------------------------------# 
# MODEL 2 ------
#-----------------------------------------------------------------------------------------# 

# Model 2 ROC
library(pROC)
Model2 = real3res[["Models"]][["Combination 309"]]

# Get predicted probabilities
pred_probs <- predict(Model2, type = "response")

# Get the actual outcome variable from the dataset
actual_outcome <- real3res[["Models"]][["Combination 309"]][["y"]]

# Compute ROC curve
Model1_roc <- roc(actual_outcome, pred_probs,ci = TRUE)



# Model 2 details 

# Compute ROC curve
Model2_roc <- roc(actual_outcome, pred_probs, ci = TRUE)

# Compute best threshold for sensitivity and specificity
best_coords <- coords(Model2_roc, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden", transpose = FALSE)

# Extract values for Youden's cutoff, sensitivity, and specificity
youden_cutoff <- as.numeric(best_coords["threshold"])
sensitivity <- as.numeric(best_coords["sensitivity"])
specificity <- as.numeric(best_coords["specificity"])

# Compute 95% CI for AUC
auc_ci <- ci.auc(Model2_roc, conf.level = 0.95)  # Bootstrapped CI for AUC

# Compute 95% CI for sensitivity and specificity at Youden's cutoff
ci_sensitivity <- ci.coords(Model2_roc, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
ci_specificity <- ci.coords(Model2_roc, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)

# Extract CI values
sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
specificity_lower_ci <- ci_specificity[["specificity"]][1]
specificity_upper_ci <- ci_specificity[["specificity"]][3]

# Store results in a dataframe
Model2_df <- data.frame(
  AUC = as.numeric(auc(Model2_roc)),
  AUC_CI_Lower = auc_ci[1],
  AUC_CI_Upper = auc_ci[3],
  Youden_Cutoff = youden_cutoff,
  Sensitivity = sensitivity,
  Sensitivity_CI_Lower = sensitivity_lower_ci,
  Sensitivity_CI_Upper = sensitivity_upper_ci,
  Specificity = specificity,
  Specificity_CI_Lower = specificity_lower_ci,
  Specificity_CI_Upper = specificity_upper_ci
)

print(Model2_df)






#-----------------------------------------------------------------------------------------# 
# MODEL 3 ------
#-----------------------------------------------------------------------------------------# 


## Model 3 ROC
library(pROC)
Model3 = real3res[["Models"]][["MIG"]]

# Get predicted probabilities
pred_probs <- predict(Model3, type = "response")

# Get the actual outcome variable from the dataset
actual_outcome <- real3res[["Models"]][["MIG"]][["y"]]

# Compute ROC curve
Model3_roc <- roc(actual_outcome, pred_probs,ci = TRUE)


# Model 3 details

# Compute ROC curve
Model3_roc <- roc(actual_outcome, pred_probs, ci = TRUE)

# Compute best threshold for sensitivity and specificity
best_coords <- coords(Model3_roc, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden", transpose = FALSE)

# Extract values for Youden's cutoff, sensitivity, and specificity
youden_cutoff <- as.numeric(best_coords["threshold"])
sensitivity <- as.numeric(best_coords["sensitivity"])
specificity <- as.numeric(best_coords["specificity"])

# Compute 95% CI for AUC
auc_ci <- ci.auc(Model3_roc, conf.level = 0.95)  # Bootstrapped CI for AUC

# Compute 95% CI for sensitivity and specificity at Youden's cutoff
ci_sensitivity <- ci.coords(Model3_roc, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
ci_specificity <- ci.coords(Model3_roc, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)

# Extract CI values
sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
specificity_lower_ci <- ci_specificity[["specificity"]][1]
specificity_upper_ci <- ci_specificity[["specificity"]][3]

# Store results in a dataframe
Model3_df <- data.frame(
  AUC = as.numeric(auc(Model3_roc)),
  AUC_CI_Lower = auc_ci[1],
  AUC_CI_Upper = auc_ci[3],
  Youden_Cutoff = youden_cutoff,
  Sensitivity = sensitivity,
  Sensitivity_CI_Lower = sensitivity_lower_ci,
  Sensitivity_CI_Upper = sensitivity_upper_ci,
  Specificity = specificity,
  Specificity_CI_Lower = specificity_lower_ci,
  Specificity_CI_Upper = specificity_upper_ci
)

print(Model3_df)


# COMBINE IN ONE DF

# Format the output to two decimal places
format_ci <- function(val, ci_lower, ci_upper) {
  paste0(format(round(val, 2), nsmall = 2), " (", format(round(ci_lower, 2), nsmall = 2), "-", format(round(ci_upper, 2), nsmall = 2), ")")
}

# Combine the results
combined_df <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3"),
  AUC = c(format_ci(Model1_df$AUC, Model1_df$AUC_CI_Lower, Model1_df$AUC_CI_Upper),
          format_ci(Model2_df$AUC, Model2_df$AUC_CI_Lower, Model2_df$AUC_CI_Upper),
          format_ci(Model3_df$AUC, Model3_df$AUC_CI_Lower, Model3_df$AUC_CI_Upper)),
  Sensitivity = c(format_ci(Model1_df$Sensitivity, Model1_df$Sensitivity_CI_Lower, Model1_df$Sensitivity_CI_Upper),
                  format_ci(Model2_df$Sensitivity, Model2_df$Sensitivity_CI_Lower, Model2_df$Sensitivity_CI_Upper),
                  format_ci(Model3_df$Sensitivity, Model3_df$Sensitivity_CI_Lower, Model3_df$Sensitivity_CI_Upper)),
  Specificity = c(format_ci(Model1_df$Specificity, Model1_df$Specificity_CI_Lower, Model1_df$Specificity_CI_Upper),
                  format_ci(Model2_df$Specificity, Model2_df$Specificity_CI_Lower, Model2_df$Specificity_CI_Upper),
                  format_ci(Model3_df$Specificity, Model3_df$Specificity_CI_Lower, Model3_df$Specificity_CI_Upper))
)

combined_df

# Save to CSV
write.csv(combined_df, "Signatures_model_results.csv", row.names = FALSE)

# save ----
save.image("Saved_on_090325.RData")
