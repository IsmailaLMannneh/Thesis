# Codes

# Introduction
# This document describes the workflow for POC analyses 

# Load libraries
library(readxl)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(ggplot2)
library(gridExtra)
```

#------------------------------------------------------------------------------#
# PART 2A: SUMMARY STATISTICS 
#------------------------------------------------------------------------------#

```{r datasets}
# Datasets 
# 1. AllData # Category = Group
# 2. AllNeg
# 3. WithinTB
# 4. WithinORD
# 5. XpertBasedData # Cat = GeneXpert_Ultra
# 6. XpertHIV_Neg
```

# All TB vs ORD using Composite ref
```{r summary gene expression}
# Dataset 1: Summarise TB vs ORD in ALl dxn Data --------

# Calculate mean, 25th percentile, 75th percentile, and IQR value
AllData_Summary <- data %>%
  group_by(Group) %>%
  summarise(
    KLF2_mean = mean(KLF2, na.rm = TRUE),
    KLF2_Q1 = quantile(KLF2, 0.25, na.rm = TRUE),
    KLF2_Q3 = quantile(KLF2, 0.75, na.rm = TRUE),
    KLF2_IQR = KLF2_Q3 - KLF2_Q1,
    DUSP3_mean = mean(DUSP3, na.rm = TRUE),
    DUSP3_Q1 = quantile(DUSP3, 0.25, na.rm = TRUE),
    DUSP3_Q3 = quantile(DUSP3, 0.75, na.rm = TRUE),
    DUSP3_IQR = DUSP3_Q3 - DUSP3_Q1,
    GBP5_mean = mean(GBP5, na.rm = TRUE),
    GBP5_Q1 = quantile(GBP5, 0.25, na.rm = TRUE),
    GBP5_Q3 = quantile(GBP5, 0.75, na.rm = TRUE),
    GBP5_IQR = GBP5_Q3 - GBP5_Q1,
    TB_score_mean = mean(TB_score, na.rm = TRUE),
    TB_score_Q1 = quantile(TB_score, 0.25, na.rm = TRUE),
    TB_score_Q3 = quantile(TB_score, 0.75, na.rm = TRUE),
    TB_score_IQR = TB_score_Q3 - TB_score_Q1
  )

# Print the summary statistics
write.csv(AllData_Summary, "Summary1_AllData_Summary.csv")
```


# Summary: TB vs ORD in HIV seronegative patients 
```{r summary hiv seronegatives }
AllNeg_summary <- AllNeg %>%
  group_by(Group) %>%
  summarise(
    KLF2_mean = mean(KLF2, na.rm = TRUE),
    KLF2_Q1 = quantile(KLF2, 0.25, na.rm = TRUE),
    KLF2_Q3 = quantile(KLF2, 0.75, na.rm = TRUE),
    KLF2_IQR = KLF2_Q3 - KLF2_Q1,
    DUSP3_mean = mean(DUSP3, na.rm = TRUE),
    DUSP3_Q1 = quantile(DUSP3, 0.25, na.rm = TRUE),
    DUSP3_Q3 = quantile(DUSP3, 0.75, na.rm = TRUE),
    DUSP3_IQR = DUSP3_Q3 - DUSP3_Q1,
    GBP5_mean = mean(GBP5, na.rm = TRUE),
    GBP5_Q1 = quantile(GBP5, 0.25, na.rm = TRUE),
    GBP5_Q3 = quantile(GBP5, 0.75, na.rm = TRUE),
    GBP5_IQR = GBP5_Q3 - GBP5_Q1,
    TB_score_mean = mean(TB_score, na.rm = TRUE),
    TB_score_Q1 = quantile(TB_score, 0.25, na.rm = TRUE),
    TB_score_Q3 = quantile(TB_score, 0.75, na.rm = TRUE),
    TB_score_IQR = TB_score_Q3 - TB_score_Q1
  )
write.csv(AllNeg_summary, "Summary2_AllNeg_summary.csv")
```


#  Summarise TB vs ORD  based on XpertBasedData
```{r summary based on Xpert}

XpertData_summary <- XpertBasedData %>%
  group_by(GeneXpert_Ultra) %>%
  summarise(
    KLF2_mean = mean(KLF2, na.rm = TRUE),
    KLF2_Q1 = quantile(KLF2, 0.25, na.rm = TRUE),
    KLF2_Q3 = quantile(KLF2, 0.75, na.rm = TRUE),
    KLF2_IQR = KLF2_Q3 - KLF2_Q1,
    DUSP3_mean = mean(DUSP3, na.rm = TRUE),
    DUSP3_Q1 = quantile(DUSP3, 0.25, na.rm = TRUE),
    DUSP3_Q3 = quantile(DUSP3, 0.75, na.rm = TRUE),
    DUSP3_IQR = DUSP3_Q3 - DUSP3_Q1,
    GBP5_mean = mean(GBP5, na.rm = TRUE),
    GBP5_Q1 = quantile(GBP5, 0.25, na.rm = TRUE),
    GBP5_Q3 = quantile(GBP5, 0.75, na.rm = TRUE),
    GBP5_IQR = GBP5_Q3 - GBP5_Q1,
    TB_score_mean = mean(TB_score, na.rm = TRUE),
    TB_score_Q1 = quantile(TB_score, 0.25, na.rm = TRUE),
    TB_score_Q3 = quantile(TB_score, 0.75, na.rm = TRUE),
    TB_score_IQR = TB_score_Q3 - TB_score_Q1
  )
write.csv(XpertData_summary, "Summary3_XpertData_summary.csv") 
```


# Summarise TB vs ORD gene expression in seronegatives in XpertBasedData
```{r}
XpertHIV_Neg_summary <- XpertHIV_Neg %>%
  group_by(GeneXpert_Ultra) %>%
  summarise(
    KLF2_mean = mean(KLF2, na.rm = TRUE),
    KLF2_Q1 = quantile(KLF2, 0.25, na.rm = TRUE),
    KLF2_Q3 = quantile(KLF2, 0.75, na.rm = TRUE),
    KLF2_IQR = KLF2_Q3 - KLF2_Q1,
    DUSP3_mean = mean(DUSP3, na.rm = TRUE),
    DUSP3_Q1 = quantile(DUSP3, 0.25, na.rm = TRUE),
    DUSP3_Q3 = quantile(DUSP3, 0.75, na.rm = TRUE),
    DUSP3_IQR = DUSP3_Q3 - DUSP3_Q1,
    GBP5_mean = mean(GBP5, na.rm = TRUE),
    GBP5_Q1 = quantile(GBP5, 0.25, na.rm = TRUE),
    GBP5_Q3 = quantile(GBP5, 0.75, na.rm = TRUE),
    GBP5_IQR = GBP5_Q3 - GBP5_Q1,
    TB_score_mean = mean(TB_score, na.rm = TRUE),
    TB_score_Q1 = quantile(TB_score, 0.25, na.rm = TRUE),
    TB_score_Q3 = quantile(TB_score, 0.75, na.rm = TRUE),
    TB_score_IQR = TB_score_Q3 - TB_score_Q1
  )

# Print the summary statistics
write.csv(XpertHIV_Neg_summary, "Summary4_XpertHIV_Neg_summary.csv") 
```



#------------------------------------------------------------------------------#
# PART 2B: DEMOGRAPHIC CHARACTERISTICS  
#------------------------------------------------------------------------------#
```{r demo data import}
# Demographic details 
# Set Dir 
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/Cepheid3HR/Task3_3HR_data_cleaning")

# data -----
library(readxl)
CepheidDataset <- read_excel("CepheidAnalysis_DATA1.xlsx")

# select column -----
names(CepheidDataset)
CepheidData = CepheidDataset[, c(1, 53, 3, 6, 7, 9, 10, 11, 12, 16:21, 27, 28, 
                                 30,43, 46, 47,49, 52, 54,58)]

# "How many previous episodes?"  is Have TB history. 
#  "TB_history" = recentTB history (last 12 months)
# reorganise names 

#----------------------------------------------------------- ###########
# split into data and demo 

# A. Data --------------
names(CepheidData)
Ceph_data = CepheidData[, c(1, 2, 3, 5, 21, 16, 18, 10:15, 19, 20, 21, 23:25 )]
names(Ceph_data)

# Clean data 
names(Ceph_data)[1] = "ParticipantID"
names(Ceph_data)[13] = "LOT"
Ceph_data$LOT = as.character(Ceph_data$LOT)
Ceph_data = subset(Ceph_data, enrolled %in% c("Yes"))
Ceph_data <- Ceph_data[Ceph_data$ParticipantID != 'TRI0500', ]
summary(Ceph_data)

# Convert character columns to factors while keeping numeric columns unchanged
Ceph_data <- data.frame(lapply(Ceph_data, function(x) {
  if (is.character(x)) {
    return(as.factor(x))  # Convert character columns to factors
  } else {
    return(x)  # Leave numeric columns unchanged
  }
}))

# Create a composite group 
Ceph_data$Composite_Ref = Ceph_data$CompositeDx

# Replace values 
Ceph_data$Composite_Ref = gsub("CAT A", "Definite TB", Ceph_data$Composite_Ref)
Ceph_data$Composite_Ref = gsub("CAT B", "Definite TB", Ceph_data$Composite_Ref)
Ceph_data$Composite_Ref = gsub("CAT C", "Definite TB", Ceph_data$Composite_Ref)
Ceph_data$Composite_Ref = gsub("CAT E", "Probable TB", Ceph_data$Composite_Ref)
Ceph_data$Composite_Ref = gsub("CAT H", "No TB", Ceph_data$Composite_Ref)
Ceph_data$Composite_Ref = gsub("CAT I", "No TB", Ceph_data$Composite_Ref)
#Ceph_data$Composite_Ref = gsub("CAT U (Enrolled)", "Unknown TB", Ceph_data$Composite_Ref)
Ceph_data$Composite_Ref <- gsub("CAT U \\(Enrolled\\)", "Unknown TB", Ceph_data$Composite_Ref)


# Composite Binary 
Ceph_data$BinaryClass = Ceph_data$Composite_Ref
Ceph_data$BinaryClass = gsub("Definite TB", "TB", Ceph_data$BinaryClass)
Ceph_data$BinaryClass = gsub("Probable TB", "TB", Ceph_data$BinaryClass)
Ceph_data$BinaryClass = gsub("Unknown TB", "TB", Ceph_data$BinaryClass) # FOR ALL COMPARISONS
Ceph_data$BinaryClass = gsub("No TB", "ORD", Ceph_data$BinaryClass)

# # Geneexpert 
# Ceph_data$GeneXpert_Ultra = gsub("MTB detected", "TB", Ceph_data$GeneXpert_Ultra)
# Ceph_data$GeneXpert_Ultra = gsub("MTB not detected", "ORD", Ceph_data$GeneXpert_Ultra)


###### ----------------------------------------------------------- ###########
# B. demo ------
Ceph_Demo = CepheidData[, -c(10:14)]
names(Ceph_Demo)

# Clean data 
names(Ceph_Demo)[1] = "ParticipantID"
Ceph_Demo = subset(Ceph_Demo, enrolled %in% c("Yes"))
Ceph_Demo <- Ceph_Demo[Ceph_Demo$ParticipantID != 'TRI0500', ]
summary(Ceph_Demo)


# Convert character columns to factors while keeping numeric columns unchanged
Ceph_Demo <- data.frame(lapply(Ceph_Demo, function(x) {
  if (is.character(x)) {
    return(as.factor(x))  # Convert character columns to factors
  } else {
    return(x)  # Leave numeric columns unchanged
  }
}))

names(Ceph_Demo)
summary(Ceph_Demo)


# Create a composite group 
Ceph_Demo$Composite_Ref = Ceph_Demo$CompositeDx

# Replace values 
Ceph_Demo$Composite_Ref = gsub("CAT A", "Definite TB", Ceph_Demo$Composite_Ref)
Ceph_Demo$Composite_Ref = gsub("CAT B", "Definite TB", Ceph_Demo$Composite_Ref)
Ceph_Demo$Composite_Ref = gsub("CAT C", "Definite TB", Ceph_Demo$Composite_Ref)
Ceph_Demo$Composite_Ref = gsub("CAT E", "Probable TB", Ceph_Demo$Composite_Ref)
Ceph_Demo$Composite_Ref = gsub("CAT H", "No TB", Ceph_Demo$Composite_Ref)
Ceph_Demo$Composite_Ref = gsub("CAT I", "No TB", Ceph_Demo$Composite_Ref)
#Ceph_Demo$Composite_Ref = gsub("CAT U (Enrolled)", "Unknown TB", Ceph_Demo$Composite_Ref)
Ceph_Demo$Composite_Ref <- gsub("CAT U \\(Enrolled\\)", "Unknown TB", Ceph_Demo$Composite_Ref)

# Composite Binary 
Ceph_Demo$BinaryClass = Ceph_Demo$Composite_Ref
Ceph_Demo$BinaryClass = gsub("Definite TB", "TB", Ceph_Demo$BinaryClass)
Ceph_Demo$BinaryClass = gsub("Probable TB", "TB", Ceph_Demo$BinaryClass)
Ceph_Demo$BinaryClass = gsub("Unknown TB", "TB", Ceph_Demo$BinaryClass)
Ceph_Demo$BinaryClass = gsub("No TB", "ORD", Ceph_Demo$BinaryClass)

# Geneexpert 
# Ceph_Demo$GeneXpert_Ultra = gsub("MTB detected", "TB", Ceph_Demo$GeneXpert_Ultra)
# Ceph_Demo$GeneXpert_Ultra = gsub("MTB not detected", "ORD", Ceph_Demo$GeneXpert_Ultra)

names(Ceph_Demo)

# calculate Age 
Ceph_Demo$Date.of.Birth <- as.Date(Ceph_Demo$Date.of.Birth)

# Calculate age and create a new column 'age'
Ceph_Demo$Age <- floor(as.numeric(difftime(Sys.Date(), Ceph_Demo$Date.of.Birth, units = "weeks")) / 52.25)
summary(Ceph_Demo$Age)
summary(Ceph_Demo)

#save.image("CephDataceaning.RData")
```

 # Demographic 
```{r datasets}
# Define datasets 
names(Ceph_Demo)
dim(Ceph_Demo)
# [1] 300  22

# demography based on Based on Tb categories 
table(Ceph_Demo$Composite_Ref)

# Next demo based on Binary grouping Subset 
# here we exclude unknown TB and possible TB. Our dataset has no possible TB group. 

# Note: Binaryclass is 'group' in POC analyses
AllTBSets = Ceph_Demo
names(AllTBSets)
table(AllTBSets$Composite_Ref)

# For binary classification take all tb except unknown TB 
DefiniteBinarySet = subset(Ceph_Demo, Composite_Ref %in% c("Definite TB","Probable TB", "No TB"))
dim(DefiniteBinarySet)
names(DefiniteBinarySet)
table(DefiniteBinarySet$Composite_Ref)
table(DefiniteBinarySet$BinaryClass) # matched composite 
table(DefiniteBinarySet$GeneXpert_Ultra)

table(DefiniteBinarySet$Smear_grade)
table(DefiniteBinarySet$Smear_grade, DefiniteBinarySet$BinaryClass)
table(DefiniteBinarySet$Xpert_grade, DefiniteBinarySet$BinaryClass)

# table(AllData$Group)
# ORD  TB 
# 185 105 

  #    ORD TB
  # High       0 59
  # Low        0  9
  # Medium     0 26
  # Trace      1  0
  # Very low   0  4

  # ORD  TB
  # 0  185  15
  # 1+   0  26
  # 2+   0  23
  # 3+   0  41

```

# Method 1: Comparisons based on Binary variable (composite reference)

```{r compariosn 1 gtsummary}
library(gtsummary)
# Convert all variables except Age and BMI to factors
AllTBSets = Ceph_Demo

# Remove columns using dplyr::select
AllTBSets <- AllTBSets %>%
  select(-Date.of.Birth, -ParticipantID, -enrolled, -Phase, 
         -CompositeDx, -How.many.previous.episodes., -TTP1, 
         -Composite_Ref, -Composite, -TB_history, -Lot...)
names(AllTBSets)

# Convert all variables except Age and BMI to factors, if they aren't already
AllTBSets <- AllTBSets %>%
  mutate(across(-c(Age, BMI), as.factor))

# Create the summary table, comparing by BinaryClass (TB vs ORD)
gtsummary_table <- AllTBSets %>%
  tbl_summary(
    by = BinaryClass,  # Group by the BinaryClass variable
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",  # Show mean and SD for Age and BMI
      all_categorical() ~ "{n} / {N} ({p}%)"  # Show count and percentage for categorical variables
    ),
    missing = "ifany"  # Handle missing data
  ) %>%
  add_p()  # Add p-values for group comparisons (TB vs ORD)

# Display the summary table
gtsummary_table

# The analyses for demography must show patients and ofcourse based on composite but 
# subsequent analyses would only included definete TB vs ORD
```


# Method 2: Using Two separate tables, one for all comparisons and the other for within Tb comparisons 

```{r comparions}
# Load necessary libraries
library(dplyr)
library(gtsummary)

# Subset dataset for the full comparison between TB and ORD
full_comparison <- Ceph_Demo %>%
  select(-Date.of.Birth, -ParticipantID, -enrolled, -Phase, -CompositeDx, 
        -How.many.previous.episodes., -TTP1, -Composite_Ref, 
        -Xpert_grade, -Liquid_culture1, -Smear_grade, -TB_history)

# Create a dataset for variables that only apply to the TB group
tb_only <- AllTBSets %>%
  filter(BinaryClass == "TB") %>%
  select(Xpert_grade, Liquid_culture1, Smear_grade) 


# Comparison 1. Create the summary table comparing TB vs ORD for common variables
#-------------------------------------------------------------------------------------#
# All comparisons, Note: dx based on composstie collapse as binary 

gtsummary_full <- full_comparison %>%
  tbl_summary(
    by = BinaryClass,  # Group by the BinaryClass variable
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",  # Show mean and SD for continuous variables (Age, BMI)
      all_categorical() ~ "{n} / {N} ({p}%)"  # Show count and percentage for categorical variables
    ),
    missing = "ifany"  # Handle missing data
  ) %>%
  add_p()  # Add p-values for comparison

write.csv(as.data.frame(gtsummary_full), "Table1. Comparison demographic characteristics by composite.csv")


# 2. Create a separate summary table for the TB group only
#----------------------------------------------------------------------#
# Summary table for TB group only
gtsummary_tb_only <- tb_only %>%
  tbl_summary(
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"  # Show count and percentage for categorical variables
    ),
    missing = "ifany")  # Handle missing data

# Print the summary tables
gtsummary_full
gtsummary_tb_only

# save.image("Patients_Demography_done.RData")
# load("Patients_Demography_done.RData")
```



#----------------------------------------------------------------#
# PART 3: GENE EXPRESSION COMPARISON OF 
#----------------------------------------------------------------# 

```{r datasets}
# All comparisons box plots TIFF ------ 
#----------------------------------------------------------------#
# Datasets 
# 1. AllData # Category = Group
# 2. AllNeg
# 3. WithinTB
# 4. WithinORD
# 5. XpertBasedData # Cat = GeneXpert_Ultra
# 6. XpertHIV_Neg
```

# 1. Mean expressions of genes and TB score in All TB vs ORD (CRS ref)
```{r gene expression}
# Panel of four plots 
p1 <- ggplot(AllData, aes(Group, GBP5)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "GBP5")

p2 <- ggplot(AllData, aes(Group, DUSP3)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  #scale_y_continuous(limits = c(16, 23)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 23) +
  labs(title = "DUSP3")

p3 <- ggplot(AllData, aes(Group, KLF2)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "KLF2")

p4 <- ggplot(AllData, aes(Group, TB_score)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  #scale_y_continuous(limits = c(20, 27)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 27) +
  labs(title = "TB Score")

# Arrange the four plots in a 2x2 grid
tiff("Boxplot1.AllData_Mean_expressions_TB_vs_ORD.tiff", 
     width = 10, height = 8, units = "in", res = 300)
AllData_comparison <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off() 
# This is figure 1. 
```


# 2. Mean expressions of genes and TB score in HIV neg TB vs ORD (AllNeg
```{r genes}
# Create individual plots for each variable
ai <- ggplot(AllNeg, aes(Group, GBP5)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = Group), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "GBP5")

aii <- ggplot(AllNeg, aes(Group, DUSP3)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = Group), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 23)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 23) +
  labs(title = "DUSP3")

aiii <- ggplot(AllNeg, aes(Group, KLF2)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = Group), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "KLF2")

aiv <- ggplot(AllNeg, aes(Group, TB_score)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = Group), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(20, 27)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, #label.y = 27) +
  labs(title = "TB Score")

# Arrange the four plots in a 2x2 grid
tiff("Boxplot2.AllNeg_Mean_expressions_TB_vs_ORD.tiff", 
     width = 10, height = 8, units = "in", res = 300)
allneg.comparisons = grid.arrange(ai, aii, aiii, aiv, ncol = 2, nrow = 2)
dev.off()
```


# 3. Mean expressions comparison based on Xpert (XpertBasedData
```{r gens}
names(XpertBasedData) # GeneXpert_Ultra

y1 <- ggplot(XpertBasedData, aes(GeneXpert_Ultra, GBP5)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "GBP5")

y2 <- ggplot(XpertBasedData, aes(GeneXpert_Ultra, DUSP3)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 23)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 23) +
  labs(title = "DUSP3")

y3 <- ggplot(XpertBasedData, aes(GeneXpert_Ultra, KLF2)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "KLF2")

y4 <- ggplot(XpertBasedData, aes(GeneXpert_Ultra, TB_score)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(20, 27)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, #label.y = 27) +
  labs(title = "TB Score")

# Arrange the four plots in a 2x2 grid
tiff("Boxplot3.XpertBasedData_Mean_expressions_TB_vs_ORD.tiff", 
     width = 10, height = 8, units = "in", res = 300)
XpertBasedData.comparisons = grid.arrange(y1, y2, y3, y4, ncol = 2, nrow = 2)
dev.off()
```


# 4. Mean expressions comparison in HIV neg based on Xpert (XpertHIV_Neg) -----
```{r xpert seronegatives}
names(XpertHIV_Neg) # GeneXpert_Ultra

y1 <- ggplot(XpertHIV_Neg, aes(GeneXpert_Ultra, GBP5)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "GBP5")

y2 <- ggplot(XpertHIV_Neg, aes(GeneXpert_Ultra, DUSP3)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 23)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 23) +
  labs(title = "DUSP3")

y3 <- ggplot(XpertHIV_Neg, aes(GeneXpert_Ultra, KLF2)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "KLF2")

y4 <- ggplot(XpertHIV_Neg, aes(GeneXpert_Ultra, TB_score)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = GeneXpert_Ultra), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(20, 27)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, #label.y = 27) +
  labs(title = "TB Score")

# Arrange the four plots in a 2x2 grid
tiff("Boxplot4.XpertHIV_Neg_Mean_expressions_TB_vs_ORD.tiff", 
     width = 10, height = 8, units = "in", res = 300)
XpertHIV_Neg.comparisons = grid.arrange(y1, y2, y3, y4, ncol = 2, nrow = 2)
dev.off()
```




# Other comparisons (based on Composite reference)


# a. Within TB patients: by HIV status 
```{r hiv tb }
# names(WithinTB) # There is no NA 
summary(WithinTB$HIV)

# Create individual plots for each variable
withintb_hivcomparison1 <- ggplot(WithinTB, aes(HIV, GBP5)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25)+ 
  #, label.y = 22) +
  labs(title = "GBP5")

withintb_hivcomparison2 <- ggplot(WithinTB, aes(HIV, DUSP3)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  #scale_y_continuous(limits = c(16, 23)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) +
  #, label.y = 23) +
  labs(title = "DUSP3")

withintb_hivcomparison3 <- ggplot(WithinTB, aes(HIV, KLF2)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) +
  #, label.y = 22) +
  labs(title = "KLF2")


withintb_hivcomparison4 <- ggplot(WithinTB, aes(HIV, TB_score)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_y_continuous(limits = c(-5, 0)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25, label.y = 0) +
  labs(title = "TB Score") 

# Arrange the four plots in a 2x2 grid
tiff("Boxplot5_comparisons withintb_by_HIV_status.tiff", 
     width = 8, height = 8, res=300, unit="in")
withintb_by_HIV_status = grid.arrange(withintb_hivcomparison1, withintb_hivcomparison2,
                                      withintb_hivcomparison3, withintb_hivcomparison4,
                                      ncol = 2, nrow = 2)
dev.off()

#table(WithinTB$HIV, WithinTB$Group) # 103 vs 2
```


# Within TB by HIV plot export (Figure)
```{r plot}
# Figure export
# Arrange the four plots in a 2x2 grid
tiff("Figure_WithinTB_by_HIV_status.tiff", 
     width = 6, height = 6, res=300, unit="in")
ggplot(WithinTB, aes(HIV, TB_score)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_y_continuous(limits = c(-5, 0)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) +
  labs(title = "TB Score") 
dev.off()
```


# b. Within TB seronegatives by HIV 
```{r plots}
levels(WithinORD$HIV) # there a missing hiv datapoint
WithinORD <- WithinORD %>%
  filter(!is.na(HIV))
# Create individual plots for each variable
w1 <- ggplot(WithinORD, aes(HIV, GBP5)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + 
  labs(title = "GBP5")


w2 <- ggplot(WithinORD, aes(HIV, DUSP3)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + 
  labs(title = "DUSP3")

w3 <- ggplot(WithinORD, aes(HIV, KLF2)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + 
  labs(title = "KLF2")


w4 <- ggplot(WithinORD, aes(HIV, TB_score)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + 
  labs(title = "TB_score")
# Arrange the four plots in a 2x2 grid
tiff("Boxplot6_comparisons-WithinORD_by_HIV_status.tiff", 
     width = 10, height = 10, res=300, unit="in")
WithinORD_by_HIV_status = grid.arrange(w1, w2, w3, w4,ncol = 2, nrow = 2)
dev.off()
```


# c. TB score by HIV within ORD patients
```{r ord}
# Plots to display: TB score within ORD patients 
tiff("Figure_WithinORD_by_HIV_status.tiff", 
     width = 6, height = 6, res=300, unit="in")
  ggplot(WithinORD, aes(HIV, TB_score)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = HIV), width = 0.15, size = 4, alpha = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + 
  labs(title = "TB_score")
dev.off()
```


d. Within TB by Xpert grades
```{r xpert grades}
# Define comparisons for Xpert_grade
my_comparisons <- list(
  #c("Trace", "Very low"),
                       c("Low", "Trace/Very low"),
                       c("Medium", "Low"),
                       c("High", "Medium"))

# Note only one trace to replace with verylow then combine categories
WithinTBx2 = WithinTBx
WithinTBx2$Xpert_grade = gsub("Trace", "Very low", WithinTBx2$Xpert_grade)
WithinTBx2$Xpert_grade = gsub("Very low", "Trace/Very low", WithinTBx2$Xpert_grade)
WithinTBx2$Xpert_grade <- factor(WithinTBx2$Xpert_grade, 
                                 levels = c("Trace/Very low", "Low", "Medium", "High"))

                       

# Arrange the four plots in a 2x2 grid
tiff("Boxplot2_comparisons_WithinTB_by_Xpert_grade.tiff", 
     width = 6, height = 6, res=1200, unit="in")

tb_score_plot <- ggplot(WithinTBx2, aes(x = Xpert_grade, y = TB_score)) +
  geom_boxplot(aes(color = Xpert_grade), width = 0.6, size = 1.2) +  # Boxplot
  geom_jitter(aes(color = Xpert_grade), width = 0.15, size = 3, alpha = 0.3) +  # Jitter for individual points
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD")) +  # Custom colors for Xpert_grade
  scale_y_continuous(limits = c(-5, 0)) +  # Set y-axis limits to range from -5 to 0
  theme_minimal() +  # Minimal theme for a clean look
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  labs(title = "TB Score Across Xpert Grades",
       x = "Xpert Grade",
       y = "TB Score") +
  stat_compare_means(comparisons = my_comparisons,  # Add p-values for comparisons
                     method = "wilcox.test",
                     label = "p.signif",  # Show significant p-values
                     label.x = 1.5,  # Adjust label position on x-axis if needed
                     label.y = -0.5,
                     p.adjust.method = "BH")  # Adjust label position on y-axis to avoid overlap

# Display the plot
print(tb_score_plot)
dev.off()
```


# d. Within TB by  smear grade
```{r}
Smear <- read_excel("HR_smear_grades.xlsx")[1:2]
names(Smear)[1] = "ParticipantID"
names(Smear)[2] = "Grade_smear"

# Within TB data 
names(WithinTB)
dim(WithinTB)

WithinTB$ParticipantID
TBsmearData = merge(WithinTB, Smear, by="ParticipantID")

dim(TBsmearData)

# grades 
TBsmearData$Grade_smear
class(TBsmearData$Grade_smear)
TBsmearData$Grade_smear = as.factor(TBsmearData$Grade_smear)
table(TBsmearData$Grade_smear)
# 0 1+ 2+ 3+ 
# 15 26 23 41 

# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Define comparisons 
sm_comparisons <- list(c("0", "1+"),
                       c("1+", "2+"),
                       c("2+", "3+"))

# Arrange the four plots in a 2x2 grid
tiff("Boxplot2_comparisons_WithinTB_by_Smear_grade.tiff", 
     width = 6, height = 6, res=300, unit="in")

tb_smear_plot <- ggplot(TBsmearData, aes(x = Grade_smear, y = TB_score)) +
  geom_boxplot(aes(color = Grade_smear), width = 0.6, size = 1.2) +  # Boxplot
  geom_jitter(aes(color = Grade_smear), width = 0.15, size = 3, alpha = 0.3) +  # Jitter for individual points
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD")) +  # Custom colors for Xpert_grade
  #scale_y_continuous(limits = c(-5, 0)) +  # Set y-axis limits to range from -5 to 0
  theme_minimal() +  # Minimal theme for a clean look
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  labs(title = "TB Score Across varying smear Grades",
       x = "Smear Grade",
       y = "TB Score") +
  stat_compare_means(comparisons = sm_comparisons,  # Add p-values for comparisons
                     method = "wilcox.test",
                     label = "p.signif",  # Show significant p-values
                     label.x = 1.5,  # Adjust label position on x-axis if needed
                     label.y = -0.5,
                     p.adjust.method = "BH")  

# Display the plot
print(tb_smear_plot)
dev.off()

# MORE (added 06022025)


#---------------------------------------------------------------#
# Between Smear pos TB and smear neg tb  
#---------------------------------------------------------------#
names(AllData)
# create smear dataset 

# make a new copy 

BacillaryData = AllData
dim(BacillaryData)

# Smear is not available in Alldata 
# Add Smear
dim(Smear)
names(Smear)
Baci_Data = merge(BacillaryData, Smear, by="ParticipantID")
dim(Baci_Data)

Baci_Data$Grade_smear

# Create a new column 'SmearBinary' based on 'Grade_smear'
Baci_Data <- Baci_Data %>%
  mutate(SmearBinary = ifelse(Grade_smear == "0", "Negative", "Positive"))
Baci_Data = merge(BacillaryData, Smear, by="ParticipantID")

# Ensure SmearBinary column is created correctly
Baci_Data <- Baci_Data %>%
  mutate(SmearBinary = ifelse(Grade_smear == "0", "Negative", "Positive"))

Baci_Data$SmearBinary

# write.csv(Baci_Data, "Smear_Baci_Data.csv")




# BOX plot  
#_____________________________________________________#

# Create a box plot comparing positive and negative smear groups

# Arrange the four plots in a 2x2 grid
tiff("Boxplot_TBscore_in_smaer_pos and negatives.tiff", 
     width = 6, height = 6, res=1200, unit="in")

smearbinary <- ggplot(Baci_Data, aes(x = SmearBinary, y = TB_score)) +
  geom_boxplot(aes(color = SmearBinary), width = 0.6, size = 1.2) + 
  geom_jitter(aes(color = SmearBinary), width = 0.15, size = 3, alpha = 0.3) +  
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +  
  scale_y_continuous(limits = c(-6, 1)) +
  theme_minimal() +  # Minimal theme for a clean look
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  labs(title = "TB Score in  postive and negative smear Grades",
       x = "Smear Grade",
       y = "TB Score") +
  stat_compare_means(
    #comparisons = sm_comparisons,  # Add p-values for comparisons
                     method = "wilcox.test",
                     label = "p.signif",  # Show significant p-values
                     label.x = 1.5,  # Adjust label position on x-axis if needed
                     #label.y = 1.5,
                     p.adjust.method = "BH")

# Display the plot
print(smearbinary)
dev.off()




## ROC plot  
#_____________________________________________________#

# Function to generate ROC plot for a gene
smearBinaryRoc <- function(gene_name, data) {
  rocc <- roc(data$SmearBinary, as.numeric(data[[gene_name]]))  # Generate ROC curve
  
  # Calculate AUC and 95% CI using bootstrap with 5000 iterations
  ci_roc <- ci.auc(rocc, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve
  ggroc(rocc, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = gene_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocc), 2), 
                           "\n95% CI:", 
                           round(ci_roc[1], 2), "-", 
                           round(ci_roc[3], 2)), size = 5)
}



# Step 1: Generate ROC plots for each gene
SmearNg_vs_Pos_roc <- smearBinaryRoc("TB_score", Baci_Data)
# Arrange the four plots in a 2x2 grid
tiff("ROC_between_Smear_Pos_TB_and_ORD.tiff", 
     width = 6, height = 6, res=1200, unit="in")
SmearNg_vs_Pos_roc$
dev.off()

# ROCC details 
# Define the target gene  
target_gene <- "TB_score"  # Replace with the desired gene  

# Perform ROC analysis for the selected gene  
resultsm <- smearbinary(TB_score, Baci_Data)  

# Display the result  
print(result)



# SP at 70% ROC detail






# 
# # plot 
# library(ggplot2)
# library(ggpubr)
# library(dplyr)
# 
# # Ensure SmearBinary column is created correctly
# Baci_Data <- Baci_Data %>%
#   mutate(SmearBinary = ifelse(Grade_smear == "0", "Negative", "Positive"))
# 
# 


#---------------------------------------------------------------#
# Within Smear pos TB only 
#---------------------------------------------------------------#

# Manually edited to Smear_neg_TB.csv contain last col sm-negTB

# load data 
SmearNegTb = read_csv("Smear_neg_TB.csv")
SmearNegTb = SmearNegTb[,-1]
SmearNegTb$Sm_negTB
SmearNegTb = subset(SmearNegTb, Sm_negTB%in% c("TB","ORD"))
SmearNegTb$Sm_negTB = as.factor(as.character(SmearNegTb$Sm_negTB))

# Arrange the four plots in a 2x2 grid
tiff("Boxplot_between_Smear_Neg_TB_and_ORD.tiff", 
     width = 6, height = 6, res=1200, unit="in")

# Box plot with jitter and exact p-value
ggplot(SmearNegTb, aes(x = Sm_negTB, y = TB_score)) +
  geom_boxplot(aes(color = Sm_negTB), width = 0.6, size = 1.2) + 
  geom_jitter(aes(color = Sm_negTB), width = 0.15, size = 3, alpha = 0.3) +  
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +  
  scale_y_continuous(limits = c(-6, 1)) +
  theme_minimal() +  
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),
        axis.ticks = element_line(size = 1)) +
  labs(title = "TB Score in Smear neg TB and ORD patients",
       x = "Smear Grade",
       y = "TB Score") +
  stat_compare_means(
    method = "wilcox.test",  # Use Wilcoxon test
    label = "p.format",  # Display exact p-value instead of asterisks
    p.adjust.method = "BH",
    label.x = 1.5  # Adjust label position if needed
  )

dev.off()


# ROC 
# Function to generate ROC plot for a gene
smearRoc <- function(gene_name, data) {
  rocc <- roc(data$Sm_negTB, as.numeric(data[[gene_name]]))  # Generate ROC curve
  
  # Calculate AUC and 95% CI using bootstrap with 5000 iterations
  ci_roc <- ci.auc(rocc, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve
  ggroc(rocc, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = gene_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocc), 2), 
                           "\n95% CI:", 
                           round(ci_roc[1], 2), "-", 
                           round(ci_roc[3], 2)), size = 5)
}



# Step 1: Generate ROC plots for each gene
names(SmearNegTb)
SmearNgTB_roc <- smearRoc("TB_score", SmearNegTb)
# Arrange the four plots in a 2x2 grid
tiff("ROC_between_Smear_Neg_TB_and_ORD.tiff", 
     width = 6, height = 6, res=1200, unit="in")
SmearNgTB_roc
dev.off()



# Roc details and rocc at SE70

# smaer roc details 

# data Baci_Data gene TB_score


# Function to perform ROC analysis and calculate Youden's index, sensitivity, and specificity
roc_detail_smear <- function(gene_name, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$SmearBinary, as.numeric(data[[gene_name]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
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
    Gene = gene_name,
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

# List of genes to analyze
genes <- c("TB_score")

# Initialize an empty list 
results_list <- list()

# Loop through each gene and perform ROC analysis
for (gene in genes) {
  results_list[[gene]] <- roc_detail_smear(gene, Baci_Data)
}
# Combine all results into a single data frame
results_dfs <- do.call(rbind, results_list)

# Display the results dataframe
print(results_dfs)
write.csv(results_dfs, "1_Roc_details_smear_pos_vs_neg.csv")




# At TPP

# List of genes to analyze
genes <- c("TB_score")

# Initialize a list to store results
results_list <- list()

# Loop through each gene and compute ROC analysis
for (gene_name in genes) {
  
  # Compute ROC curve
  roc_curve <- roc(Baci_Data$SmearBinary, as.numeric(Baci_Data[[gene_name]]))
  
  # Get thresholds, sensitivities, and specificities
  all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
  
  # Find the closest specificity to 70%
  fixed_spec <- 0.70
  closest_idx <- which.min(abs(all_coords$specificity - fixed_spec))
  
  # Get the corresponding threshold and sensitivity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # Compute 95% CI for sensitivity at the fixed specificity using ci.se()
  ci_sens <- ci.se(roc_curve, specificities = fixed_spec, boot.n = 2000)
  sensitivity_ci_lower <- ci_sens[1]
  sensitivity_ci_upper <- ci_sens[3]
  
  # Store results in a data frame
  results_list[[gene_name]] <- data.frame(
    Gene = gene_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Sensitivity_CI_Lower = sensitivity_ci_lower,
    Sensitivity_CI_Upper = sensitivity_ci_upper,
    Specificity = closest_specificity
  )
}

# Combine results into a single data frame
results_sm70 <- do.call(rbind, results_list)

# Display the results
print(results_sm70)

# Save the results to a CSV file if needed
#write.csv(results_df70, "HR_sens_70Sp_with_CI.csv", row.names = FALSE)
write.csv(results_sm70, "Roc_details_for_smear_pos_vs_neg_TPP.csv")


```



# ------------------------------------------------------------------#
# kruskal_test for all comparison within TB and ORD categories 


```{r}
# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(TB_score ~ Grade_smear, data = TBsmearData)
print(kruskal_test_result)

kruskal_xpgrade <- kruskal.test(TB_score ~ Xpert_grade, data = WithinTB)
print(kruskal_xpgrade)
# # 4. With and withoout TB history
# # 5. Look at others demo info.
```



# ------------------------------------------------------------------#
# PART 4 : ROCC analyses
# ------------------------------------------------------------------#
```{r data}
# Datasets 
# 1. AllData # Category = Group
# 2. AllNeg
# 3. WithinTB
# 4. WithinORD
# 5. XpertBasedData # Cat = GeneXpert_Ultra
# 6. XpertHIV_Neg
```



# ROC 1. AllDATA  TB v ORD (compsoite ref)

# define function for rocc
```{r fucntion }
# Function to generate ROC plot for a gene
generate_roc_plot <- function(gene_name, data) {
  rocc <- roc(data$Group, as.numeric(data[[gene_name]]))  # Generate ROC curve
  
  # Calculate AUC and 95% CI using bootstrap with 5000 iterations
  ci_roc <- ci.auc(rocc, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve
  ggroc(rocc, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = gene_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocc), 2), 
                           "\n95% CI:", 
                           round(ci_roc[1], 2), "-", 
                           round(ci_roc[3], 2)), size = 5)
}

```


# Plot rocc  
```{r}
# Step 1: Generate ROC plots for each gene
names(AllData)
roc_gbp5 <- generate_roc_plot("GBP5", AllData)
roc_dusp3 <- generate_roc_plot("DUSP3", AllData)
roc_klf2 <- generate_roc_plot("KLF2", AllData)
roc_tb_score <- generate_roc_plot("TB_score", AllData)

# Step 2: Arrange the ROC plots in a 2x2 grid and save to pdf 
tiff("Roc1_Accuracy_of_genes_and_tbscore_in_alldata_TB_and_ORD.tiff", 
     width = 10, height = 10, units = "in", res=300)  # Start PDF device
Alldata_rocc = grid.arrange(roc_gbp5, roc_dusp3, roc_klf2, roc_tb_score, ncol = 2, nrow = 2)
dev.off() # 
```

# Get roc info 
```{r roc results}
# Function to perform ROC analysis and calculate Youden's index, sensitivity, and specificity
roc_analysis <- function(gene_name, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$Group, as.numeric(data[[gene_name]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
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
    Gene = gene_name,
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

# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize an empty list 
results_list <- list()

# Loop through each gene and perform ROC analysis
for (gene in genes) {
  results_list[[gene]] <- roc_analysis(gene, AllData)
}
# Combine all results into a single data frame
results_df <- do.call(rbind, results_list)

# Display the results dataframe
print(results_df)
write.csv(results_df, "1_Roc_details_AllData_new.csv")


```

# When optimized for TPP : one at a time, check loop
```{r rocc details}
# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize a list to store results
results_list <- list()

# Loop through each gene and compute ROC analysis
for (gene_name in genes) {
  
  # Compute ROC curve
  roc_curve <- roc(AllData$Group, as.numeric(AllData[[gene_name]]))
  
  # Get thresholds, sensitivities, and specificities
  all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
  
  # Find the closest specificity to 70%
  fixed_spec <- 0.70
  closest_idx <- which.min(abs(all_coords$specificity - fixed_spec))
  
  # Get the corresponding threshold and sensitivity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # Compute 95% CI for sensitivity at the fixed specificity using ci.se()
  ci_sens <- ci.se(roc_curve, specificities = fixed_spec, boot.n = 2000)
  sensitivity_ci_lower <- ci_sens[1]
  sensitivity_ci_upper <- ci_sens[3]
  
  # Alternatively, compute 95% CI for sensitivity using ci.coords()
  # Uncomment the following block if you want to use ci.coords instead
  # ci_sensitivity <- ci.coords(roc_curve, x = fixed_spec, input = "specificity", ret = "sensitivity", boot.n = 2000)
  # sensitivity_ci_lower <- ci_sensitivity[1]
  # sensitivity_ci_upper <- ci_sensitivity[3]
  
  # Store results in a data frame
  results_list[[gene_name]] <- data.frame(
    Gene = gene_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Sensitivity_CI_Lower = sensitivity_ci_lower,
    Sensitivity_CI_Upper = sensitivity_ci_upper,
    Specificity = closest_specificity
  )
}

# Combine results into a single data frame
results_df70 <- do.call(rbind, results_list)

# Display the results
print(results_df70)

# Save the results to a CSV file if needed
#write.csv(results_df70, "HR_sens_70Sp_with_CI.csv", row.names = FALSE)
write.csv(results_df70, "Roc_details_AllData_TPP.csv") #redone correctly
```


#### corrected to SE90
```{r set se90}
# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize a list to store results
results_list <- list()

# Loop through each gene and compute ROC analysis
for (gene_name in genes) {
  
  # Compute ROC curve
  roc_curve <- roc(AllData$Group, as.numeric(AllData[[gene_name]]))
  
  # Get thresholds, sensitivities, and specificities
  all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
  
  # Find the closest sensitivity to 90%
  fixed_sens <- 0.90
  closest_idx <- which.min(abs(all_coords$sensitivity - fixed_sens))
  
  # Get the corresponding threshold and specificity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # Compute 95% CI for specificity at the fixed sensitivity using ci.sp()
  ci_spec <- ci.sp(roc_curve, sensitivities = fixed_sens, boot.n = 2000)
  specificity_ci_lower <- ci_spec[1]
  specificity_ci_upper <- ci_spec[3]
  
  # Store results in a data frame
  results_list[[gene_name]] <- data.frame(
    Gene = gene_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Specificity = closest_specificity,
    Specificity_CI_Lower = specificity_ci_lower,
    Specificity_CI_Upper = specificity_ci_upper
  )
}

# Combine results into a single data frame
results_df90 <- do.call(rbind, results_list)

# Display the results
print(results_df90)

# Save the results to a CSV file if needed
write.csv(results_df90, "Roc_details_AllData_TPP_90Se.csv")
```


# ROC 2. All HIV seronegative TB vs ord  ----------

# plot rocc
```{r}
# data 
AllNeg
# Step 1: Generate ROC plots for each gene
AllNegroc_gbp5 <- generate_roc_plot("GBP5", AllNeg)
AllNegroc_dusp3 <- generate_roc_plot("DUSP3", AllNeg)
AllNegroc_klf2 <- generate_roc_plot("KLF2", AllNeg)
AllNegroc_tb_score <- generate_roc_plot("TB_score", AllNeg)

# Step 2: Arrange the ROC plots in a 2x2 grid
tiff("Roc2_Accuracy_of_genes_and_tbscore_in_allneg_TB_and_ORD.tiff", 
    width = 10, height = 10, units = "in", res=300) 
AllNeg_rocc = grid.arrange(AllNegroc_gbp5, AllNegroc_dusp3, AllNegroc_klf2, 
                                  AllNegroc_tb_score, ncol = 2, nrow = 2)
dev.off()
```


# ROCC details 
```{r roc detail}
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize an empty list 
Seroneg_list <- list()

# Loop through each gene and perform ROC analysis
for (gene in genes) {
  Seroneg_list[[gene]] <- roc_analysis(gene, AllNeg)
}
# Combine all results into a single data frame
Neg.results_df <- do.call(rbind, Seroneg_list)

# Display the results dataframe
print(Neg.results_df)
write.csv(Neg.results_df, "1_Roc_details_AllNeg_new.csv")
```

# optimizinf for TPP
```{r}
# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize a list to store results
results_list <- list()

# Loop through each gene and compute ROC analysis
for (gene_name in genes) {
  
  # Compute ROC curve
  roc_curve <- roc(AllNeg$Group, as.numeric(AllNeg[[gene_name]]))
  
  # Get thresholds, sensitivities, and specificities
  all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
  
  # Find the closest sensitivity to 90%
  fixed_sens <- 0.90
  closest_idx <- which.min(abs(all_coords$sensitivity - fixed_sens))
  
  # Get the corresponding threshold and specificity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # Compute 95% CI for specificity at the fixed sensitivity using ci.sp()
  ci_spec <- ci.sp(roc_curve, sensitivities = fixed_sens, boot.n = 2000)
  specificity_ci_lower <- ci_spec[1]
  specificity_ci_upper <- ci_spec[3]
  
  # Store results in a data frame
  results_list[[gene_name]] <- data.frame(
    Gene = gene_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Specificity = closest_specificity,
    Specificity_CI_Lower = specificity_ci_lower,
    Specificity_CI_Upper = specificity_ci_upper
  )
}

# Combine results into a single data frame
results_allneg <- do.call(rbind, results_list)

# Display the results
print(results_allneg)

# Display the results
print(results_allneg)

# Save the results to a CSV file if needed
write.csv(results_allneg, "Roc_details_AllNeg_TPP.csv")
```







#------------------------------------------------------------------#
# Overlay plot
#------------------------------------------------------------------#
```{r}

# Function to generate overlay ROC plot for a gene
generate_overlay_roc_plot <- function(gene_name, data1, data2) {
  
  # Generate ROC curves for both datasets
  roc1 <- roc(data1$Group, as.numeric(data1[[gene_name]]))  # ROC for AllData
  roc2 <- roc(data2$Group, as.numeric(data2[[gene_name]]))  # ROC for AllNeg
  
  # Calculate AUC and 95% CI for both datasets
  ci_roc1 <- ci.auc(roc1, method = "bootstrap", boot.n = 5000)
  ci_roc2 <- ci.auc(roc2, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve of AllData
  plot <- ggroc(roc1, legacy.axes = TRUE, col = "black", size = 1) +
    geom_line(data = ggroc(roc2)$data, aes(x = 1 - specificity, y = sensitivity), 
              color = "#FF5733", size = 1) +  # Overlay ROC for AllNeg
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = gene_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI for both datasets
    annotate("text", x = 0.55, y = 0.20,
             label = paste0("All Definite TB vs ORD AUC: ", round(auc(roc1), 2), 
                            " (CI ", round(ci_roc1[1], 2), "-", round(ci_roc1[3], 2), ")"), 
             size = 4.5, color = "black") +
    
    annotate("text", x = 0.55, y = 0.15,
             label = paste0("HIV Neg TB vs ORD AUC: ", round(auc(roc2), 2), 
                            " (CI ", round(ci_roc2[1], 2), "-", round(ci_roc2[3], 2), ")"), 
             size = 4.5, color = "#FF5733")
  
  return(plot)
}

# Step 1: Generate overlay ROC plots for each gene
overlay_roc_gbp5 <- generate_overlay_roc_plot("GBP5", AllData, AllNeg)
overlay_roc_dusp3 <- generate_overlay_roc_plot("DUSP3", AllData, AllNeg)
overlay_roc_klf2 <- generate_overlay_roc_plot("KLF2", AllData, AllNeg)
overlay_roc_tb_score <- generate_overlay_roc_plot("TB_score", AllData, AllNeg)

# Step 2: Arrange the four ROC plots in a 2x2 grid
tiff("Overlay_Roc_Alldata.and.AllNeg_TB_and_ORD.tiff", 
     width = 10, height = 10, units = "in", res=300)
grid.arrange(overlay_roc_gbp5, overlay_roc_dusp3, overlay_roc_klf2, overlay_roc_tb_score, ncol = 2)
dev.off()

```



# ------------------------------------------------------------------#
# Part 4B. GeneXpert based ROCC analyses
# ------------------------------------------------------------------#
```{r}
#data 
XpertBasedData 
XpertHIV_Neg 


# 1. Roc TB vs ORD in XpertBasedData as group -----

# Change function to add xpert ultra
generate_Xpertroc <- function(gene_name, data) {
  rocc <- roc(data$GeneXpert_Ultra, as.numeric(data[[gene_name]]))  # Generate ROC curve
  
  # Calculate AUC and 95% CI using bootstrap with 5000 iterations
  ci_roc <- ci.auc(rocc, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve
  ggroc(rocc, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = gene_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocc), 2), 
                           "\n95% CI:", 
                           round(ci_roc[1], 2), "-", 
                           round(ci_roc[3], 2)), size = 5)
}

# Step 1: Generate ROC plots for each gene
roc_gbp5 <- generate_Xpertroc("GBP5", XpertBasedData) # 
roc_dusp3 <- generate_Xpertroc("DUSP3", XpertBasedData) #
roc_klf2 <- generate_Xpertroc("KLF2", XpertBasedData) # 
roc_tb_score <- generate_Xpertroc("TB_score", XpertBasedData) # 

# Step 2: Arrange the ROC plots in a 2x2 grid and save to pdf 
tiff("x1.Roc_of_genes_and_tbscore_in_Aaa_xpert_TB_and_ORD.tiff", 
     width = 10, height = 10, units = "in", res=300)  # Start PDF device
plotrocs = grid.arrange(roc_gbp5, roc_dusp3, roc_klf2, roc_tb_score, ncol = 2, nrow = 2)
dev.off() # 


# Function to perform ROC analysis and calculate Youden's index, sensitivity, and specificity
roc_analysis <- function(gene_name, data) {
  # Calculate ROC curve
  Xpertroc_curve <- roc(data$GeneXpert_Ultra, as.numeric(data[[gene_name]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(Xpertroc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
  # Compute 95% CI for sensitivity and specificity at Youden's cutoff
  ci_sensitivity <- ci.coords(Xpertroc_curve, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
  ci_specificity <- ci.coords(Xpertroc_curve, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)
  
  # Extract CI values
  sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
  sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
  specificity_lower_ci <- ci_specificity[["specificity"]][1]
  specificity_upper_ci <- ci_specificity[["specificity"]][3]
  
  # Store results in a dataframe
  results_df <- data.frame(
    Gene = gene_name,
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

# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize an empty list 
Xpertresults_list <- list()

# Loop through each gene and perform ROC analysis
for (gene in genes) {
  Xpertresults_list[[gene]] <- roc_analysis(gene, XpertBasedData)
}

# Combine all results into a single data frame
results_dfxx <- do.call(rbind, Xpertresults_list)

# Display the results dataframe
print(results_dfxx)
write.csv(results_dfxx, "Roc_details_XpertData_new.csv")


# When optimized for TPP : one at a time, check loop ------
# #-------------------------------------------------------------------------#
# 
# # Create a list of variable names
# vars <- c("GBP5", "DUSP3", "KLF2", "TB_score")
# 
# # Define a function to perform the ROC analysis and extract statistics
# perform_xproc <- function(var_name, data, fixed_sens = 0.90) {
#   # Generate ROC curve
#   roc_obj <- roc(data$GeneXpert_Ultra, data[[var_name]])
#   
#   # Get thresholds, sensitivities, and specificities
#   all_coords <- coords(roc_obj, ret = c("threshold", "sensitivity", "specificity"))
#   
#   # Find the closest sensitivity to the fixed_sens value
#   closest_idx <- which.min(abs(all_coords$sensitivity - fixed_sens))
#   
#   # Get the corresponding threshold, sensitivity, and specificity
#   closest_threshold <- all_coords$threshold[closest_idx]
#   closest_sensitivity <- all_coords$sensitivity[closest_idx]
#   closest_specificity <- all_coords$specificity[closest_idx]
#   
#   # CI for the specificity at the closest threshold
#   ci_spec <- ci.sp(roc_obj, thresholds = closest_threshold)
#   
#   # Check if row with name '0.9' exists in ci_spec
#   if ("0.9" %in% rownames(ci_spec)) {
#     row_0_9 <- ci_spec["0.9", ]  # Extract the row where rowname = 0.9
#     
#     # Get column 1 and 3 values at row '0.9'
#     ci_0_9_lower <- row_0_9[1]  # Column 1
#     ci_0_9_upper <- row_0_9[3]  # Column 3
#   } else {
#     ci_0_9_lower <- NA
#     ci_0_9_upper <- NA
#   }
#   
#   # Return the results as a data frame
#   return(data.frame(
#     Variable = var_name,
#     Threshold = closest_threshold,
#     Sensitivity = closest_sensitivity,
#     Specificity = closest_specificity,
#     CI_Lower_0.9 = ci_0_9_lower,
#     CI_Upper_0.9 = ci_0_9_upper
#   ))
# }
# 
# # Initialize an empty data frame to store the results
# results_df <- data.frame()
# 
# # Loop through each variable and perform ROC analysis, collecting results in the data frame
# for (var in vars) {
#   result <- perform_xproc(var, XpertBasedData)
#   results_df <- rbind(results_df, result)
# }
# 
# # Print the final data frame
# print(results_df)

# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize a list to store results
results_list <- list()

# Loop through each gene and compute ROC analysis
for (gene_name in genes) {
  
  # Compute ROC curve
  roc_curve <- roc(XpertBasedData$Group, as.numeric(XpertBasedData[[gene_name]]))
  
  # Get thresholds, sensitivities, and specificities
  all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
  
 # Find the closest sensitivity to 90%
  fixed_sens <- 0.90
  closest_idx <- which.min(abs(all_coords$sensitivity - fixed_sens))
  
  # Get the corresponding threshold and specificity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # Compute 95% CI for specificity at the fixed sensitivity using ci.sp()
  ci_spec <- ci.sp(roc_curve, sensitivities = fixed_sens, boot.n = 2000)
  specificity_ci_lower <- ci_spec[1]
  specificity_ci_upper <- ci_spec[3]
  
  # Store results in a data frame
  results_list[[gene_name]] <- data.frame(
    Gene = gene_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Specificity = closest_specificity,
    Specificity_CI_Lower = specificity_ci_lower,
    Specificity_CI_Upper = specificity_ci_upper
  )
}


# Combine results into a single data frame
results_xnegSE90 <- do.call(rbind, results_list)

# Display the results
print(results_xnegSE90)
write.csv(results_xnegSE90, "Roc_details_XpertData_TPP.csv")




# 2. Roc TB vs ORD in XpertNeg as group ------


# Step 1: Generate ROC plots for each gene
roc_gbp5 <- generate_Xpertroc("GBP5", XpertHIV_Neg) # 
roc_dusp3 <- generate_Xpertroc("DUSP3", XpertHIV_Neg) #
roc_klf2 <- generate_Xpertroc("KLF2", XpertHIV_Neg) # 
roc_tb_score <- generate_Xpertroc("TB_score", XpertHIV_Neg) # 

# Step 2: Arrange the ROC plots in a 2x2 grid and save to pdf 
tiff("x1.Roc_of_genes_and_tbscore_in_XpertNeg_TB_and_ORD.tiff", 
     width = 10, height = 10, units = "in", res=300)  # Start PDF device
plotrocs = grid.arrange(roc_gbp5, roc_dusp3, roc_klf2, roc_tb_score, ncol = 2, nrow = 2)
dev.off() # 

#-------------------------------------#
# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize an empty list 
Xpertresults_list <- list()

# Loop through each gene and perform ROC analysis
for (gene in genes) {
  Xpertresults_list[[gene]] <- roc_analysis(gene, XpertHIV_Neg)
}

# Combine all results into a single data frame
results_dfxx <- do.call(rbind, Xpertresults_list)

# Display the results dataframe
print(results_dfxx)
write.csv(results_dfxx, "Roc_details_XpertNeg_new.csv")



#-------------------------------------------------------------------------#
# When optimized for TPP : one at a time, check loop -----

# List of genes to analyze
genes <- c("GBP5", "DUSP3", "KLF2", "TB_score")

# Initialize a list to store results
results_list <- list()

# Loop through each gene and compute ROC analysis
for (gene_name in genes) {
  
  # Compute ROC curve
  roc_curve <- roc(XpertHIV_Neg$Group, as.numeric(XpertHIV_Neg[[gene_name]]))
  
  # Get thresholds, sensitivities, and specificities
  all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
  
  # Find the closest sensitivity to 90%
  fixed_sens <- 0.90
  closest_idx <- which.min(abs(all_coords$sensitivity - fixed_sens))
  
  # Get the corresponding threshold and specificity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # Compute 95% CI for specificity at the fixed sensitivity using ci.sp()
  ci_spec <- ci.sp(roc_curve, sensitivities = fixed_sens, boot.n = 2000)
  specificity_ci_lower <- ci_spec[1]
  specificity_ci_upper <- ci_spec[3]
  
  # Store results in a data frame
  results_list[[gene_name]] <- data.frame(
    Gene = gene_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Specificity = closest_specificity,
    Specificity_CI_Lower = specificity_ci_lower,
    Specificity_CI_Upper = specificity_ci_upper
  )
}


# Combine results into a single data frame
results_xneg <- do.call(rbind, results_list)

# Display the results
print(results_xneg)

# Save the results to a CSV file if needed
write.csv(results_xneg, "Roc_details_XpertNeg_TPP.csv", row.names = FALSE)

#-------------------------------------------------------------------------#


#-------------------------------------------------------------------------#
# Overlay plot -----
#-------------------------------------------------------------------------#

# Function to generate overlay ROC plot for a gene
generate_overlay_roc_plotX <- function(gene_name, data1, data2) {
  
  # Generate ROC curves for both datasets
  roc1 <- roc(data1$GeneXpert_Ultra, as.numeric(data1[[gene_name]]))  # ROC for XpertBasedData
  roc2 <- roc(data2$GeneXpert_Ultra, as.numeric(data2[[gene_name]]))  # ROC for XpertHIV_Neg
  
  # Calculate AUC and 95% CI for both datasets
  ci_roc1 <- ci.auc(roc1, method = "bootstrap", boot.n = 5000)
  ci_roc2 <- ci.auc(roc2, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve of AllData
  plot <- ggroc(roc1, legacy.axes = TRUE, col = "black", size = 1) +
    geom_line(data = ggroc(roc2)$data, aes(x = 1 - specificity, y = sensitivity), color = "#FF5733", size = 1) +  # Overlay ROC for AllNeg
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = gene_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI for both datasets
    annotate("text", x = 0.55, y = 0.20,
             label = paste0("GeneXpert TB vs ORD AUC: ", round(auc(roc1), 2), 
                            " (CI ", round(ci_roc1[1], 2), "-", round(ci_roc1[3], 2), ")"), 
             size = 4.5, color = "black") +
    
    annotate("text", x = 0.55, y = 0.15,
             label = paste0("HIV- GeneXpert TB vs ORD AUC: ", round(auc(roc2), 2), 
                            " (CI ", round(ci_roc2[1], 2), "-", round(ci_roc2[3], 2), ")"), 
             size = 4.5, color = "#FF5733")
  
  return(plot)
}

# Step 1: Generate overlay ROC plots for each gene
overlay_rocx_gbp5 <- generate_overlay_roc_plotX("GBP5", XpertBasedData, XpertHIV_Neg)
overlay_rocx_dusp3 <- generate_overlay_roc_plotX("DUSP3", XpertBasedData, XpertHIV_Neg)
overlay_rocx_klf2 <- generate_overlay_roc_plotX("KLF2", XpertBasedData, XpertHIV_Neg)
overlay_rocx_tb_score <- generate_overlay_roc_plotX("TB_score", XpertBasedData, XpertHIV_Neg)

# Step 2: Arrange the four ROC plots in a 2x2 grid
# Step 2: Arrange the four ROC plots in a 2x2 grid
tiff("Overlay_Roc_GeneXpert.and.XpertHIVNeg_TB_and_ORD.tiff", 
     width = 12, height = 12, units = "in", res=300)
grid.arrange(overlay_rocx_gbp5, overlay_rocx_dusp3, overlay_rocx_klf2, overlay_rocx_tb_score, ncol = 2)
dev.off()

```






#---------------------------------------------------------------------------#
# # Part 4C. ROC AUC Comparisons using (Delong Test
#---------------------------------------------------------------------------#

# Perform DeLong's test AllData ------
```{r delong alldata}
delong_gbp5 <- roc.test(roc(AllData$Group, AllData$GBP5), roc(AllNeg$Group, AllNeg$GBP5))
delong_dusp3 <- roc.test(roc(AllData$Group, AllData$DUSP3), roc(AllNeg$Group, AllNeg$DUSP3))
delong_klf2 <- roc.test(roc(AllData$Group, AllData$KLF2), roc(AllNeg$Group, AllNeg$KLF2))
delong_tb_score <- roc.test(roc(AllData$Group, AllData$TB_score), roc(AllNeg$Group, AllNeg$TB_score)) 
# Print results
print(delong_gbp5)
print(delong_dusp3)
print(delong_klf2)
print(delong_tb_score)
```


# Perform DeLong's test  GeneXpertData
```{r}
# Perform DeLong's test to compare the two ROC curves for GBP5
delong_gbp5x <- roc.test(roc(XpertBasedData$GeneXpert_Ultra, XpertBasedData$GBP5), 
                        roc(XpertHIV_Neg$GeneXpert_Ultra, XpertHIV_Neg$GBP5))
delong_dusp3x <- roc.test(roc(XpertBasedData$GeneXpert_Ultra, XpertBasedData$DUSP3), 
                         roc(XpertHIV_Neg$GeneXpert_Ultra, XpertHIV_Neg$DUSP3))
delong_klf2x <- roc.test(roc(XpertBasedData$GeneXpert_Ultra, XpertBasedData$KLF2), 
                        roc(XpertHIV_Neg$GeneXpert_Ultra, XpertHIV_Neg$KLF2))
delong_tb_scorex <- roc.test(roc(XpertBasedData$GeneXpert_Ultra, XpertBasedData$TB_score), 
                            roc(XpertHIV_Neg$GeneXpert_Ultra, XpertHIV_Neg$TB_score))

# Print results
print(delong_gbp5x)
print(delong_dusp3x)
print(delong_klf2x)
print(delong_tb_scorex)

```


# Figure 5. POC 
# Two overlay plots for TB score against Xpert and Composite ------
```{r plot}
# Get the two plots 
overlay_rocx_tb_score
tiff("Overlay_Rocs_against_composite_and_GeneXpert.tiff", 
     width = 8, height = 16, units = "in", res=300)
grid.arrange(overlay_roc_tb_score, overlay_rocx_tb_score, ncol = 1, nrow=2)
dev.off()
```

```{r plot in cols}
# Get the two plots 
overlay_rocx_tb_score
tiff("Overlay_Rocs_composite_and_GeneXpert_in_columns.tiff", 
     width = 11, height = 6, units = "in", res=300)
grid.arrange(overlay_roc_tb_score, overlay_rocx_tb_score, ncol = 2, nrow=1)
dev.off()

overlay_rocx_tb_score
pdf("Overlay_Rocs_composite_and_GeneXpert_in_columns.pdf", 
     width = 11, height = 6)
grid.arrange(overlay_roc_tb_score, overlay_rocx_tb_score, ncol = 2, nrow=1)
dev.off()

```




# ------------------------------------------------------------------#
# PART 5: Effects of treatment on POC disease score ----
# ------------------------------------------------------------------#
```{r treatment}
# Corrected

# Set Dir
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/Cepheid3HR/Task6_Data_analyses2_Treatment/")
# rm(list=ls())
# graphics.off()
# first exported without LDA total
TxDataset <- read_excel("CepheidTR_Aug2024.xlsx")

# Data with LDA total
TxDataset2 <- read_excel("CepheidTR2_Aug2024.xlsx")
colnames(TxDataset2)
data2id = TxDataset2[, c(1,12)]
names(data2id)[1] = "ParticipantID"


# rem unwanted colms
TxData = TxDataset[, -c(2,4,5,13, 17:ncol(TxDataset))]
names(TxData)
summary(TxData)

# renames cols 
names(TxData)[1] = "ParticipantID"
names(TxData)[2] = "Timepoint"
names(TxData)[3] = "GeneXpert_status"
names(TxData)[9] = "LOT"
names(TxData)[11] = "GeneXpertUltra"
names(TxData)[12] = "GeneXpert_grade"


# do they match perfectly 
data2id$ParticipantID
TxData$ParticipantID

# Check if all ParticipantIDs match
all_match <- all(data2id$ParticipantID %in% TxData$ParticipantID) && 
  all(TxData$ParticipantID %in% data2id$ParticipantID)

if (all_match) {
  print("All ParticipantIDs match perfectly between data2id and TxData.")
} else {
  print("ParticipantIDs do not match perfectly between data2id and TxData.")
}

# Merge datasets by index not common col
rownames(data2id)
rownames(TxData)

TxData <- merge(data2id[2], TxData, by = "row.names") %>%
  select(-Row.names)  # Remove the Row.names column if not needed
names(TxData)
# rearrage col 1 and 2
TxData = TxData[, c(2,1, 3:ncol(TxData))]

# rename cols 
TxData <- TxData %>%
  rename(TB_score = "New LDA Total:",
         GBP5 = "Ct: CF6 target", 
         DUSP3 = "Ct: FAM target",
         KLF2 = "Ct: CF3 target")


# replace cat for groups column
TxData$Group = gsub("CAT A", "Definite TB", TxData$Group)
TxData$Group = gsub("CAT B", "Definite TB", TxData$Group)
TxData$Group = gsub("CAT C", "Definite TB", TxData$Group)
TxData$Group = gsub("CAT E", "Probable TB", TxData$Group)
TxData$Group = gsub("CAT H", "No TB", TxData$Group)
TxData$Group = gsub("CAT I", "No TB", TxData$Group)
#TxData$Group = gsub("CAT U (Enrolled)", "Unknown TB", TxData$Group)
TxData$Group <- gsub("CAT U \\(Enrolled\\)", "Unknown TB", TxData$Group)

TxData$TB_score


#---------------------------------------------------------------------#
# # Creat Binary group for Tb ORD only 


# 1. Binary: is auto classification grouped def and prob TB as one 
TxData$Binary = TxData$Group
TxData$Binary = gsub("Definite TB", "TB", TxData$Binary)
TxData$Binary = gsub("Probable TB", "TB", TxData$Binary)
#TxData$Binary = gsub("Unknown TB", "TB", TxData$Binary)
TxData$Binary = gsub("No TB", "ORD", TxData$Binary)

# replace genexpert values 
TxData$GeneXpertUltra
TxData$GeneXpertUltra = gsub("MTB not detected", "ORD", TxData$GeneXpertUltra)
TxData$GeneXpertUltra = gsub("MTB detected", "TB", TxData$GeneXpertUltra)
table(TxData$Timepoint)

# Fill Group value (based on Binary composite) ------
TxData$Binary
TxData <- TxData %>%
  group_by(ParticipantID) %>%
  mutate(Binary = Binary[Timepoint == "Baseline"][1]) %>%
  ungroup()

# Fill missing 'Xpert for the other timepoints 
TxData <- TxData %>%
  group_by(ParticipantID) %>%
  mutate(GeneXpertUltra = GeneXpertUltra[Timepoint == "Baseline"][1]) %>%
  ungroup()


# Check the filled dataset
table(TxData$Binary, TxData$Timepoint)
table(TxData$GeneXpertUltra, TxData$Timepoint) # perfect


# Now remove screening failure and unknown TB 
# Basline data TxData_comp


#-----------------------------------------------------------------#
# Dataset 1. Binary dataset Alldata -----

# based on Composite : Binary 
Binary_TB = subset(TxData, Binary %in% c("TB")) # composite 
Binary_TB$Binary = as.factor(as.character(Binary_TB$Binary))
names(Binary_TB)
Binary_TB <- Binary_TB %>%
  mutate(across(c(3,10:14), as.factor),   
         across(c(2,5:9), as.numeric))
table(Binary_TB$Binary, Binary_TB$Timepoint)
summary(Binary_TB$TB_score)
summary(Binary_TB)
Binary_TB = Binary_TB[!is.na(Binary_TB$TB_score),]
summary(Binary_TB)
table(Binary_TB$Binary, Binary_TB$Timepoint)


#  2. HIV neg Binary data 

# AllnegTx = subset(Binary_TB, HIV %in% c("Negative"))
# AllnegTx$Negative = as.factor(as.character(AllnegTx$Negative))

# Dataset 2. Based on Xpert : ------
Xpert_TB = subset(TxData, GeneXpertUltra %in% c("TB"))
Xpert_TB$GeneXpertUltra = as.factor(as.character(Xpert_TB$GeneXpertUltra))
Xpert_TB <- Xpert_TB %>%
  mutate(across(c(3,10:14), as.factor),   
         across(c(2, 5:9), as.numeric))

table(Xpert_TB$GeneXpertUltra, Xpert_TB$Timepoint)
```


# Comparisons
```{r comparisons}
# Datasets 
# 1. Binary_TB
# 2. Xpert_TB


# ----------------------------------------------------------------------------#
# A. Based on BINARY GROUPING) -------

# Define a function to calculate p-values for a specific gene
calculate_p_values <- function(data, gene, comparisons) {
  p_values <- sapply(comparisons, function(pair) {
    formula <- as.formula(paste(gene, "~ Timepoint"))
    test <- wilcox.test(formula, data = filter(data, Timepoint %in% pair))
    return(test$p.value)
  })
  return(p_values)
}

# Apply BH correction and filter significant comparisons
get_significant_comparisons <- function(p_values, comparisons, alpha = 0.05) {
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  significant_comps <- comparisons[which(adjusted_p_values < alpha)]
  return(significant_comps)
}

# Function to plot each gene with its specific significant comparisons
plot_gene <- function(data, gene, comparisons, label_x, label_y) 
  {
  p_values <- calculate_p_values(data, gene, comparisons)
  significant_comparisons <- get_significant_comparisons(p_values, comparisons)
  
  ggplot(data, aes(Timepoint, .data[[gene]])) +
    geom_boxplot(aes(color = Timepoint), width = 0.6, size = 1.2) +
    #scale_y_continuous(limits = c(16, 28)) +
    scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")) +  
    # Updated color palette
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.8, color = "black"),
          axis.ticks = element_line(size = 0.8)) +
    stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons, 
                       label = "p.signif", 
                       label.x = label_x[1:length(significant_comparisons)],  # Custom x positions
                       label.y = label_y[1:length(significant_comparisons)]) +  # Custom y positions
    labs(title = gene, x = "Timepoint", y = paste(gene, "Ct"))
}

# Define comparison ----
my_comparisons <- list(c("Baseline", "Month 2"), 
                       c("Month 2", "Month 4"), 
                       c("Month 4", "Month 6"))

# Define custom x and y positions for labels based on significant comparisons
label_x_positions <- c(1.25, 2.25, 3.25)
label_y_positions <- c(27, 27, 27)      

# Plot genes
t1 <- plot_gene(Binary_TB, "GBP5", my_comparisons, label_x_positions, label_y_positions)
t2 <- plot_gene(Binary_TB, "DUSP3", my_comparisons, label_x_positions, label_y_positions)
t3 <- plot_gene(Binary_TB, "KLF2", my_comparisons, label_x_positions, label_y_positions)

# adjust for plot 4
label_x_positions <- c(1.25, 2.25, 3.25)
label_y_positions <- c(0,0,0) 

t4 <- plot_gene(Binary_TB, "TB_score", my_comparisons, label_x_positions, label_y_positions)

# Arrange the four plots in a 2x2 grid
tiff("A1_Treatment_CompositeBinary_over_time.tiff", width = 6, height = 6, res=300, unit="in")  # Start PDF device
comp_Binary_TR <- grid.arrange(t1, t2, t3, t4, ncol = 2, nrow = 2)
dev.off()




################################## Unpaired line ##################################
ggpaired(Binary_TB, x = "Timepoint", y = "TB_score",
         color = "Timepoint",
         xlab = "Timepoint",
         ylab = "TB_score Ct",
         line.color = "gray", 
         width = 0.6,
         line.size = 0.4,
         palette = "jco") +
  stat_compare_means(paired = FALSE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)


####################################################################

# --------------------------------------------------------------------#
# check --------
kruskal_test <- kruskal.test(TB_score ~ Timepoint, data = Binary_TB)

# Print the result of the Kruskal-Wallis test
kruskal_test 


# --------------------------------------------------------------------#
# Plot only TB score ----
#------
# Arrange the four plots in a 2x2 grid
tiff("TBscore_CompositeBinary_over_time.tiff", width = 6, height = 6, res=300, unit="in") 
t4
dev.off()
```


# Line options :for paired data (n=61)
```{r List of required Timepoints}

required_timepoints <- c("Baseline", "Month 2", "Month 4", "Month 6")

# Find Participants with data at all required Timepoints
participants_with_all_timepoints <- Binary_TB %>%
  group_by(ParticipantID) %>%
  filter(Timepoint %in% required_timepoints) %>%
  summarize(num_timepoints = n_distinct(Timepoint)) %>%
  filter(num_timepoints == length(required_timepoints)) %>%
  pull(ParticipantID)

# Filter Binary_TB to include only those Participants with data at all Timepoints
filtered_Binary_TB <- Binary_TB %>%
  filter(ParticipantID %in% participants_with_all_timepoints)

# Print the filtered dataset or use it for further analysis
dim(filtered_Binary_TB)
table(filtered_Binary_TB$Binary, filtered_Binary_TB$Timepoint)

filtered_Binary_TB
# remove zero value
#filtered_Binary_TB <- filtered_Binary_TB[filtered_Binary_TB$GBP5 != 0, ]
filtered_Binary_TB <- filtered_Binary_TB[filtered_Binary_TB$ParticipantID != "TRI1502", ]
dim(filtered_Binary_TB)
table(filtered_Binary_TB$Binary, filtered_Binary_TB$Timepoint)

# Line chat 
# Define the pairs for comparison
my_comparisons <- list(c("Baseline", "Month 2"),
                       c("Month 2", "Month 4"),
                       c("Month 4", "Month 6"))

c.1 = ggpaired(filtered_Binary_TB, x = "Timepoint", y = "GBP5",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "GBP5 Ct",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco") +
  #geom_jitter(aes(color = Binary), width = 0.15, size = 4, alpha = 0.3) +
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)

c.2 = ggpaired(filtered_Binary_TB, x = "Timepoint", y = "DUSP3",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "DUSP3 Ct",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco") +
  #geom_jitter(aes(color = Binary), width = 0.15, size = 4, alpha = 0.3) +
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)

c.3 = ggpaired(filtered_Binary_TB, x = "Timepoint", y = "KLF2",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "KLF2 Ct",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco") +
  #geom_jitter(aes(color = Binary), width = 0.15, size = 4, alpha = 0.3) +
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)


c.4 = ggpaired(filtered_Binary_TB, x = "Timepoint", y = "TB_score",
               color = "Timepoint",
               xlab = "Timepoint",
               ylab = "TB_score Ct",
               line.color = "gray", 
               width = 0.6,
               line.size = 0.4,
               palette = "jco") +
  #geom_jitter(aes(color = Binary), width = 0.15, size = 4, alpha = 0.3) +
  stat_compare_means(paired = TRUE,
                     method = "wilcox.test", 
                     comparisons = my_comparisons)

# Arrange the four plots in a 2x2 grid
tiff("A2_Treatment_Binary_over_time_all_paired_n61.tiff", 
     width = 10, height = 10, res=300, unit="in")  # Start PDF device
Binary_over_time_all_paired = grid.arrange(c.1, c.2, c.3, c.4, ncol = 2, nrow = 2)
dev.off()
```



# More on EOT: Question does EOT Tb score differ from ORD levels 
```{r}
TxData$Binary
names(TxData)
# 
TxORDscore = TxData[TxData$Binary=="ORD", c("TB_score", "Binary")] # ord =185
TxM6Score = Binary_TB[Binary_TB$Timepoint == "Month 6", c("TB_score", "Binary")] #m6=75
dim(TxM6Score)

# row bind 
Eotx = rbind(TxORDscore, TxM6Score)



#---------------------------------------------------------------------#
# Compare and plot 

# Create the box plot with p-value
tiff("EoTx_TB_vs_ORD_BinaryTB.tiff", 
     width = 6, height = 6, res=300, unit="in") 

EOTplot <- ggplot(Eotx, aes(Binary, TB_score)) +
  geom_boxplot(aes(color = Binary), width = 0.6, size = 1.2) +
  geom_jitter(aes(color = Binary), width = 0.15, size = 4, alpha = 0.3) + 
  #scale_y_continuous(limits = c(16, 22)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.25) + #, label.y = 22) +
  labs(title = "TB_score at End of tretament vs ORD levels")

print(EOTplot)
dev.off()

```



# EOT roc curve
```{r}
# data Eotx contains TB at M6 and ORD at baseline 
# fucntion: generate_roc_plot("gene", data)
# grouping var = "Group"
# rename accordingly
names(Eotx)[2] = "Group"
# save to tiff 
tiff("Roc_tbscore_at_EOT.tiff", width = 10, height = 10, units = "in", res=300)  
# plot rocc 
Month6roc <- generate_roc_plot("TB_score", Eotx)
Month6roc
dev.off() 
```



```{r save}
# save 
#save.image("C3HR_dataanalyses_Complete.RData")
#load("C3HR_dataanalyses_Complete.RData")
```





#----------------------------------------------------------------------------#
# PART 2: MBT 
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
# SECTION A. PROTEIN COMPARISONS
#----------------------------------------------------------------------------#
```{r clear Cephied enviroment}
# rm(list = ls())
```

```{r set directory}
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task5_Data_analyses/")
```

```{r libraries}
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pROC)
library(tidyr)
library(gridExtra)
```

# Load datasets
```{r load data}
# Load or rerun script 
# load("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task3_MBT_data_cleaning/MBTDataceaning_step1.RData")
```

```{r fefine Data sets }
# All timepoint and basline dataset 
All_timepoint # All dataset useful for tretament monitoring (all data)
Baseline_all_Data # baseline TB and ORD

# Based on Comp binary -------
TB_All_timepoint_comp # TB, All timepoints but within TB binary dxn
TB_baseline_comp # TB, baseline, based on composite binary dx 
# ORD data 
ORD_Alltime_comp # ORD, All timepoints; ord defined by Binary (definite only)
ORD_baseline_comp # ORD, baseline dataset

# Based on Xpert ------
Baseline_All_data_xpert  # TB, All timepints, defined Xpert, no NA
TB_baseline_Xpert # TB, Baselien def by Xpert 

ORD_Alltime_Xpert # All timepints, ORD defined Xpert, no NA
ORD_baseline_Xpert # baselien TB data, def by Xpert
```

# ------------------------------------------------------------------#
#  Section B. Proteins comparisons
# ------------------------------------------------------------------#
# 1. Proteins between TB and ORD (binary comparison) using ref standard
```{r comparions composite ref}
names(Baseline_all_Data)
# Remove rows with NA in the SAA1 column, since they are matching samples 
Baseline_all_Data <- Baseline_all_Data[!is.na(Baseline_all_Data$SAA1), ]

# does Binary have NA
Baseline_all_Data <- Baseline_all_Data[!is.na(Baseline_all_Data$Binary), ] #n=300
# this where it changes: dim changes here 
Baseline_all_Data = subset(Baseline_all_Data, Binary %in% c("TB", "ORD"))
Baseline_all_Data$Binary = as.factor(as.character(Baseline_all_Data$Binary))

p1 <- ggplot(Baseline_all_Data, aes(Binary, SAA1)) +
  geom_boxplot(aes(color = Binary), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 12)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25) + #, label.y = 9) +
  labs(title = "SAA1")

p2 <- ggplot(Baseline_all_Data, aes(Binary, CRP)) +
  geom_boxplot(aes(color = Binary), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25) + #, label.y = 3.5) +
  labs(title = "CRP")

p3 <- ggplot(Baseline_all_Data, aes(Binary, IP.10)) +
  geom_boxplot(aes(color = Binary), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25) + #label.y = 3.5) +
  labs(title = "IP-10")

# # Box plot 
# # Arrange the four plots in a 2x2 grid
# tiff("1.TBvORD_preotein_comparison_by_composite_ref.tiff", 
#      width = 10, height = 10, units = "in", res = 300) 
# TBvORD_preotein_comparison_composite = grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)
# dev.off() # 
```


# plot option 2: facet plot : cancelled out
```{r facet}
# # Load necessary libraries
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # Reshape data to long format
# Baseline_long <- Baseline_all_Data %>%
#   pivot_longer(
#     cols = c(SAA1, CRP, IP.10), # Columns to pivot
#     names_to = "Protein",        # Name of new key column
#     values_to = "Value"          # Name of new value column
#   )
# 
# # Create the combined plot with facets and custom y-axis limits
# p_combined <- ggplot(Baseline_long, aes(x = Binary, y = Value)) +
#   geom_boxplot(aes(color = Binary), width = 0.6, size = 1.2, outlier.shape = NA) +  # Exclude outliers
#   scale_color_manual(values = c("#00AFBB", "#E7B800")) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(size = 0.8, color = "black"),
#         axis.ticks = element_line(size = 0.8),
#         strip.text = element_text(size = 12)) +
#   stat_compare_means(method = "wilcox.test", 
#                      label.x = 1.25, label.y.npc = "top") + # Adjust `label.y.npc` for dynamic positioning
#   labs(title = "Distribution of Protein in TB and ORD patients") +
#   facet_wrap(~Protein, scales = "free_y") +
#   theme(strip.text = element_text(size = 12)) #+  # Customize facet label text size
#   #coord_cartesian(ylim = c(0, 3)) # Default y-axis limit for all plots
# 
# # Adjust y-axis limits manually after plotting
# p_combined <- p_combined +
#   facet_wrap(~Protein, scales = "free_y") +
#   theme(
#     strip.text.x = element_text(size = 12),
#     panel.spacing = unit(1, "lines")
#   ) +
#   scale_y_continuous(
#     limits = c(0, ifelse(Baseline_long$Protein == "SAA1", 3,
#                          ifelse(Baseline_long$Protein == "CRP", 3, 20)))
#   )
# 
# # Print the plot
# # Corrected tiff function call
# tiff("Distribution_of_Protein_in_TB_and_ORD_patients_based_on_composite.tiff", 
#      width = 12, height = 6, units = "in", res = 300)
# # Your plotting code goes here
# print(p_combined) # Strange CRP higher in ORD 
# dev.off()
```


# 2. Protein expressions based on Xpert grouping 
```{r xpert comparison}
names(Baseline_All_data_xpert) # original dataset
# Remove rows with NA in the SAA1 column, since they are matching samples 
Baseline_All <- Baseline_All_data_xpert[!is.na(Baseline_All_data_xpert$SAA1), ]

# does GeneXpert_Ultra have NA
Baseline_All <- Baseline_All[!is.na(Baseline_All$GeneXpert_Ultra), ] #n=298
Baseline_All$GeneXpert_Ultra = as.factor(as.character(Baseline_All$GeneXpert_Ultra))

q1 <- ggplot(Baseline_All, aes(GeneXpert_Ultra, SAA1)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25, label.y = 9) +
  labs(title = "SAA1")

q2 <- ggplot(Baseline_All, aes(GeneXpert_Ultra, CRP)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 4)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25, label.y = 3.5) +
  labs(title = "CRP")

q3 <- ggplot(Baseline_All, aes(GeneXpert_Ultra, IP.10)) +
  geom_boxplot(aes(color = GeneXpert_Ultra), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 4)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25, label.y = 3.5) +
  labs(title = "IP-10")

# # Arrange the four plots in a 2x2 grid
# pdf("2.TBvORD_preotein_comparison_by_Xpert_ref.pdf", width = 10, height = 10) 
# TBvORD_preotein_comparison_by_Xpert_ref = grid.arrange(q1, q2, q3, ncol = 2, nrow = 2)
# dev.off() # CRP higher in ORD observed
```


# 3. Comparison across different xpert grades 
```{r xpert grades}
# Should be within TB comparison TB_baseline_Xpert
TB_baseline_Xpert
#clean it 
# Remove rows with NA in the SAA1 column, since they are matching samples 
TB_baseline <- TB_baseline_Xpert[!is.na(TB_baseline_Xpert$SAA1), ]
TB_baseline <- TB_baseline[!is.na(TB_baseline$GeneXpert_Ultra), ] #n=114

# define levels 
TB_baseline$Xpert_grade <- factor(TB_baseline$Xpert_grade,
                                  levels = c("Trace", "Very low","Low", "Medium", "High"))

# Define the comparisons for Xpert_grade
xpert_comparisons <- list(c("Trace", "Very low"),
                          c("Very low", "Low"),
                          c("Low", "Medium"),
                          c("Medium", "High"))

# Function to calculate p-values and apply BH correction for each gene
calculate_xpert_pvalues <- function(data, protein, comparisons) {
  p_values <- sapply(comparisons, function(pair) {
    formula <- as.formula(paste(protein, "~ Xpert_grade"))
    test <- wilcox.test(formula, data = filter(data, Xpert_grade %in% pair))
    return(test$p.value)
  })
  
  # Apply BH correction
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  
  # Filter significant comparisons
  significant_comps <- comparisons[which(adjusted_p_values < 0.05)]
  return(list(p_values = adjusted_p_values, significant_comps = significant_comps))
}

# Function to plot each protein across Xpert_grade
plot_gene_xpert <- function(data, protein, comparisons, label_x, label_y) {
  pval_result <- calculate_xpert_pvalues(data, protein, comparisons)
  
  ggplot(data, aes(Xpert_grade, .data[[protein]])) +  # Xpert_grade on x-axis
    geom_boxplot(aes(color = Xpert_grade), width = 0.6, size = 1.2) +
    geom_jitter(aes(color = Xpert_grade), width = 0.2, size = 3, alpha = 0.5) +
    scale_y_continuous(limits = c(0, 5)) +
    scale_color_manual(values = c("#D62728", "#9467BD", 
                                  "#FF7F0E", "#2CA02C","#1F77B4")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.8, color = "black"),
          axis.ticks = element_line(size = 0.8)) +
    stat_compare_means(method = "wilcox.test", comparisons = pval_result$significant_comps,
                       label = "p.signif", 
                       label.x = label_x[1:length(pval_result$significant_comps)],
                       label.y = label_y[1:length(pval_result$significant_comps)]) +
    labs(title = protein, x = "Xpert Grade", y = paste(protein, "Expression"))
}

# Custom x and y positions for p-value labels
label_x_positions <- c(1.25, 2.25, 3.25)
label_y_positions <- c(4.8, 3, 3)

# Plot each gene
g1 <- plot_gene_xpert(TB_baseline, "SAA1", xpert_comparisons, label_x_positions, label_y_positions)
g2 <- plot_gene_xpert(TB_baseline, "CRP", xpert_comparisons, label_x_positions, label_y_positions)
g3 <- plot_gene_xpert(TB_baseline, "IP.10", xpert_comparisons, label_x_positions, label_y_positions)

# Arrange the four plots in a 2x2 grid
pdf("3_Protein_comparisons_across_xpert_grade.pdf", width = 10, height = 10)
Protein_comparisons_across_xpert_grade <- grid.arrange(g1, g2, g3, ncol = 2, nrow = 2)
dev.off() # No differnce 
```


# 4. Comparison across different smear grades 
```{r smear grades}
# clean NA 
TB_smear <- TB_baseline_Xpert[!is.na(TB_baseline_Xpert$SAA1), ] # data
TB_smear <- TB_smear[!is.na(TB_smear$GeneXpert_Ultra), ] #no missing xpert
TB_smear <- TB_smear[!is.na(TB_smear$Smear_grade), ] #no missing smear grade n100

# gsub
class(TB_smear$Smear_grade)
# Replace strings in TB_smear$Smear_grade
TB_smear$Smear_grade <- gsub("Negative", "0", TB_smear$Smear_grade)
TB_smear$Smear_grade <- gsub("1\\+ = less than 1 AFB per field", "1+", TB_smear$Smear_grade)
TB_smear$Smear_grade <- gsub("2\\+ = 1 to 10 AFB per field", "2+", TB_smear$Smear_grade)
TB_smear$Smear_grade <- gsub("3\\+ = more than 10 AFB per field", "3+", TB_smear$Smear_grade)

table(TB_smear$Smear_grade)

TB_smear$Smear_grade <- factor(TB_smear$Smear_grade,
                                  levels = c("0", "1+", "2+", "3+"))


# create new comparions but function above remain valid 
# Define the comparisons for Smear_grade
smear_comparisons <- list(c("0", "1+"), 
                          c("1+", "2+"), 
                          c("2+", "3+"))


# Define new functions 
# Function to calculate p-values and apply BH correction for each gene
calculate_xpert_pvalues <- function(data, protein, comparisons) {
  p_values <- sapply(comparisons, function(pair) {
    formula <- as.formula(paste(protein, "~ Smear_grade"))
    test <- wilcox.test(formula, data = filter(data, Smear_grade %in% pair))
    return(test$p.value)
  })
  
  # Apply BH correction
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  
  # Filter significant comparisons
  significant_comps <- comparisons[which(adjusted_p_values < 0.05)]
  return(list(p_values = adjusted_p_values, significant_comps = significant_comps))
}

# Function to plot each protein across Smear_grade
plot_gene_xpert <- function(data, protein, comparisons, label_x, label_y) {
  pval_result <- calculate_xpert_pvalues(data, protein, comparisons)
  
  ggplot(data, aes(Smear_grade, .data[[protein]])) +  # Smear_grade on x-axis
    geom_boxplot(aes(color = Smear_grade), width = 0.6, size = 1.2) +
    geom_jitter(aes(color = Smear_grade), width = 0.2, size = 3, alpha = 0.5) +
    scale_y_continuous(limits = c(0, 5)) +
    scale_color_manual(values = c("#D62728", "#9467BD", 
                                  "#FF7F0E", "#2CA02C","#1F77B4")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.8, color = "black"),
          axis.ticks = element_line(size = 0.8)) +
    stat_compare_means(method = "wilcox.test", comparisons = pval_result$significant_comps,
                       label = "p.signif", 
                       label.x = label_x[1:length(pval_result$significant_comps)],
                       label.y = label_y[1:length(pval_result$significant_comps)]) +
    labs(title = protein, x = "Xpert Grade", y = paste(protein, "Expression"))
}

# Custom x and y positions for p-value labels
label_x_positions <- c(1.25, 2.25, 3.25)
label_y_positions <- c(4.8, 3, 3)

#plot data 
# Plot each gene
s1 <- plot_gene_xpert(TB_smear, "SAA1", smear_comparisons, label_x_positions, label_y_positions)
s2 <- plot_gene_xpert(TB_smear, "CRP", smear_comparisons, label_x_positions, label_y_positions)
s3 <- plot_gene_xpert(TB_smear, "IP.10", smear_comparisons, label_x_positions, label_y_positions)

# Arrange the four plots in a 2x2 grid
pdf("4_Protein_comparisons_across_smear_grade.pdf", width = 10, height = 10)
Protein_comparisons_across_smear_grade <- grid.arrange(s1, s2, s3, ncol = 2, nrow = 2)
dev.off() # CRP 0 significantly higher than other grades 
```


# 5. Comparison across HIV status 
```{r baseline HIV impoteted}
Baseline_All
HIVcols <- readRDS("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task5_Data_analyses/HIV_columns_from_C3HR.RData")
names(HIVcols)
MBT_HIVData = merge(Baseline_all_Data, HIVcols, by="ParticipantID")
names(MBT_HIVData)
# ensure no NA 
MBT_HIVData = MBT_HIVData[!is.na(MBT_HIVData$HIV),]
MBT_HIVData$HIV = gsub("Positive, confirmed with second Rapid", "Positive", MBT_HIVData$HIV, ignore.case = TRUE)
table(MBT_HIVData$HIV)
MBT_HIVData$HIV = as.factor(MBT_HIVData$HIV)

# Within TB
WithinTB_Mbt = subset(MBT_HIVData, Binary %in% c("TB"))

# plots 
h1 <- ggplot(WithinTB_Mbt, aes(HIV, SAA1)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 12)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25) + #, label.y = 9) +
  labs(title = "SAA1")

h2 <- ggplot(WithinTB_Mbt, aes(HIV, CRP)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25) + #, label.y = 3.5) +
  labs(title = "CRP")

h3 <- ggplot(WithinTB_Mbt, aes(HIV, IP.10)) +
  geom_boxplot(aes(color = HIV), width = 0.6, size = 1.2) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.25) + #label.y = 3.5) +
  labs(title = "IP-10")

# Box plot 
# Arrange the four plots in a 2x2 grid
tiff("WithinTB_preotein_comparison_accross_HIV_status.tiff", 
     width = 10, height = 10, units = "in", res = 300) 
proteinHIV = grid.arrange(h1, h2, h3, ncol = 2, nrow = 2)
dev.off() # 
```




# ------------------------------------------------------------------#
#  Section B. SUMMARY STAISTICS
# ------------------------------------------------------------------#



#  A. TB vs ORD based on Binary composite 
```{r summary based on composite ref}
# Table of comparison
library(gtsummary)
dim(Baseline_all_Data) # 290  16

# Create summary table with means, Q1, Q3, and p-values
summary_table_gtsummary <- Baseline_all_Data %>%
  select(SAA1, CRP, IP.10, Binary) %>%
  tbl_summary(
    by = Binary, 
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})" # Shows Mean (Q1, Q3)
    ),
    digits = all_continuous() ~ 2 # Rounds numbers to 2 decimal places
  ) %>%
  add_p(test = all_continuous() ~ "wilcox.test") # Adds p-values from Wilcoxon test

# Print the gtsummary table
summary_table_gtsummary
write.csv(summary_table_gtsummary, "Protein distribution.csv")
```


```{r Alternatively Summary}
# Baseline_all_Data
# Calculate mean, 25th percentile, 75th percentile, and IQR value
names(Baseline_all_Data)
Bl.summary <- Baseline_all_Data %>%
  group_by(Binary) %>%
  summarise(
    SAA1_mean = median(SAA1, na.rm = TRUE),
    SAA1_Q1 = quantile(SAA1, 0.25, na.rm = TRUE),
    SAA1_Q3 = quantile(SAA1, 0.75, na.rm = TRUE),
    CRP_mean = median(CRP, na.rm = TRUE),
    CRP_Q1 = quantile(CRP, 0.25, na.rm = TRUE),
    CRP_Q3 = quantile(CRP, 0.75, na.rm = TRUE),
    IP.10_mean = median(IP.10, na.rm = TRUE),
    IP.10_Q1 = quantile(IP.10, 0.25, na.rm = TRUE),
    IP.10_Q3 = quantile(IP.10, 0.75, na.rm = TRUE),)

# Print the summary statistics
write.csv(Bl.summary.xpert, "Bl.expression.csv") # perfect 
```



# ------------------------------------------------------------------#
#  Section C. ROC Analyses ----
# ------------------------------------------------------------------#

# ROC 1. preins based on composite refernce 
```{r function}
# data is Baseline_all_Data (clean already)

# Function to generate ROC plot for a gene 
generate_roc_plot <- function(protein_name, data) {
  rocc <- roc(data$Binary, as.numeric(data[[protein_name]]))  # Generate ROC curve gene_name
  
  # Calculate AUC and 95% CI using bootstrap with 5000 iterations
  ci_roc <- ci.auc(rocc, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve
  ggroc(rocc, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "blue") +  # Add diagonal line
    labs(title = protein_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(size = 1.2, color = "black"),
          axis.ticks = element_line(size = 1)) +
    
    # Annotate with AUC and CI
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocc), 2), 
                           "\n95% CI:", 
                           round(ci_roc[1], 2), "-", 
                           round(ci_roc[3], 2)), size = 5)
}
```

```{r plot and save}
# Step 1: Generate ROC plots for each gene
roc_a <- generate_roc_plot("SAA1", Baseline_all_Data)
roc_b <- generate_roc_plot("CRP", Baseline_all_Data)
roc_c <- generate_roc_plot("IP.10", Baseline_all_Data)

# Step 2: Arrange the ROC plots in a 2x2 grid and save to pdf 
pdf("6. ROCC of protein markers based on composite binary.pdf", width = 10, height = 10)  # Start PDF device
ROCC_proteins_based_on_binary = grid.arrange(roc_a, roc_b, roc_c, ncol = 2, nrow = 2)
dev.off() # Close PDF device
```



# Roc details 
```{r rocc details}
# Function to perform ROC analysis and calculate Youden's index, sensitivity, and specificity
roc_details <- function(gene_name, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$Binary, as.numeric(data[[gene_name]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
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
    Gene = gene_name,
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

# List of genes to analyze
genes <- c("SAA1", "CRP", "IP.10")

# Initialize an empty list 
protein_list1 <- list()

# Loop through each gene and perform ROC analysis
for (i in genes) {
  protein_list1[[i]] <- roc_details(i, Baseline_all_Data)
}
# Combine all results into a single data frame
protein_df1 <- do.call(rbind, protein_list1)

# Display the results dataframe
print(protein_df1)
write.csv(protein_df1, "1_Roc_protein_df_BL.all.csv")
```


# ROC TPP benchmarking
```{r tpp}
# Create a list of variable names
vars <- c("SAA1", "CRP", "IP.10")

# Define a function to perform the ROC analysis and extract statistics
perform_proc <- function(var_name, data, fixed_spec = 0.70) {
  # Generate ROC curve
  roc_obj <- roc(data$Binary, data[[var_name]])
  all_coords <- coords(roc_obj, ret = c("threshold", "sensitivity", "specificity"))
  closest_idx <- which.min(abs(all_coords$specificity - fixed_spec))

  # Get the corresponding threshold, sensitivity, and specificity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]

  # CI for the specificity at the closest threshold
  ci_sens <- ci.se(roc_obj, thresholds = closest_threshold)

  # Check if row with name '0.9' exists in ci_spec
  if ("0.7" %in% rownames(ci_sens)) {
    row_0_7 <- ci_sens["0.7", ]

    # Get column 1 and 3 values at row '0.9'
    ci_0_7_lower <- row_0_7[1]  # Column 1
    ci_0_7_upper <- row_0_7[3]  # Column 3
  } else {
    ci_0_7_lower <- NA
    ci_0_7_upper <- NA
  }

  # Return the results as a data frame
  return(data.frame(
    Variable = var_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Specificity = closest_specificity,
    CI_Lower_0.7 = ci_0_7_lower,
    CI_Upper_0.7 = ci_0_7_upper
  ))
}

# Initialize 
results_dfa <- data.frame()
# Loop 
for (var in vars) {
  result <- perform_proc(var, Baseline_all_Data)
  results_dfa <- rbind(results_dfa, result)
}
# Print the final data frame
print(results_dfa)
write.csv(results_dfa, "Roc_baselineBinary_TPP.csv")
```



# excluding HIV positives
```{r seronegatives roc}
# ROCCs
HIVnegIDs <- readRDS("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task5_Data_analyses/AllSeronegativesIDs.RDS")
dim(Baseline_all_Data)
dim(HIVnegIDs)
SeroNegMBTData = merge(Baseline_all_Data, HIVnegIDs, by="ParticipantID")

# 
dim(SeroNegMBTData)

# Step 1: Generate ROC plots for each gene
nroc_a <- generate_roc_plot("SAA1", SeroNegMBTData) # 0.86 (0.81-0.90)
nroc_b <- generate_roc_plot("CRP", SeroNegMBTData) # 0.85 (0.81-0.89)
nroc_c <- generate_roc_plot("IP.10", SeroNegMBTData) # 0.88 (0.82-0.91)
```

# Roc details 
```{r roc detail}
# Function to perform ROC analysis and calculate Youden's index, sensitivity, and specificity
roc_details <- function(gene_name, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$Binary, as.numeric(data[[gene_name]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
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
    Gene = gene_name,
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

# List of genes to analyze
genes <- c("SAA1", "CRP", "IP.10")

# Initialize an empty list 
protein_list <- list()

# Loop through each gene and perform ROC analysis
for (i in genes) {
  protein_list[[i]] <- roc_details(i, SeroNegMBTData)
}
# Combine all results into a single data frame
protein_df <- do.call(rbind, protein_list)

# Display the results dataframe
print(protein_df)
write.csv(protein_df, "1_Roc_protein_df_Seroneg.csv")
```


# TPP ---- # Not done because HIV numbers are very small to influence performance

# Overlay Composite Binary rocc 
```{r overlap plots}
# Binary 
# Function to generate overlay ROC plot for a gene with DeLong test
Binary_overlay <- function(protein_name, data1, data2) {
  
  # Generate ROC curves for both datasets
  roc1 <- roc(data1$Binary, as.numeric(data1[[protein_name]]))  # ROC for AllData
  roc2 <- roc(data2$Binary, as.numeric(data2[[protein_name]]))  # ROC for AllNeg
  
  # Calculate AUC and 95% CI for both datasets
  ci_roc1 <- ci.auc(roc1, method = "bootstrap", boot.n = 5000)
  ci_roc2 <- ci.auc(roc2, method = "bootstrap", boot.n = 5000)
  
  # Perform DeLong test to compare the ROC curves
  delong_test <- roc.test(roc1, roc2, method = "delong")
  delong_p_value <- round(delong_test$p.value, 3)  # Round the p-value to 3 decimal places
  
  # Create ggplot object for ROC curve of AllData
  plot <- ggroc(roc1, legacy.axes = TRUE, col = "black", size = 1) +
    geom_line(data = ggroc(roc2)$data, aes(x = 1 - specificity, y = sensitivity), 
              color = "#FF5733", size = 1) +  # Overlay ROC for AllNeg
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = protein_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI for both datasets
    annotate("text", x = 0.55, y = 0.20,
             label = paste0("AUC: ", round(auc(roc1), 2), 
                            " (CI ", round(ci_roc1[1], 2), "-", round(ci_roc1[3], 2), ")"), 
             size = 4.5, color = "black") +
    
    annotate("text", x = 0.55, y = 0.15,
             label = paste0("AUC: ", round(auc(roc2), 2), 
                            " (CI ", round(ci_roc2[1], 2), "-", round(ci_roc2[3], 2), ")"), 
             size = 4.5, color = "#FF5733") +
    
    # Annotate with DeLong p-value
    annotate("text", x = 0.55, y = 0.10,
             label = paste0("DeLong p = ", delong_p_value), 
             size = 4.5, color = "blue")  # You can change color if needed
  
  return(plot)
}

# Generate overlay ROC plots for each protein
overlay_1 <- Binary_overlay("SAA1", Baseline_all_Data, SeroNegMBTData)
overlay_2 <- Binary_overlay("CRP", Baseline_all_Data, SeroNegMBTData)
overlay_3 <- Binary_overlay("IP.10", Baseline_all_Data, SeroNegMBTData)

# Arrange the plots in a grid and save to TIFF file
tiff("Overlay_Roc_MBT_allvsSeroNeg_TB_and_ORD.tiff", 
     width = 10, height = 10, units = "in", res=300)
grid.arrange(overlay_1, overlay_2, overlay_3, ncol = 2)
dev.off()
```

#  ROC 2. ROCC comparison based on Xpert ----
```{r roc xpert diagnosis}
# 1. Summary 
#------------------------------------------------------------#
# Xpert expression summary 

# Create summary table with means, Q1, Q3, and p-values
summary_XpertData <- Baseline_All %>%
  select(SAA1, CRP, IP.10, GeneXpert_Ultra) %>%
  tbl_summary(
    by = GeneXpert_Ultra, 
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})" # Shows Mean (Q1, Q3)
    ),
    digits = all_continuous() ~ 2 # Rounds numbers to 2 decimal places
  ) %>%
  add_p(test = all_continuous() ~ "wilcox.test") # Adds p-values from Wilcoxon test

# Print the gtsummary table
summary_XpertData
write.csv(summary_XpertData, "summary_XpertData.csv")



# 2. ROCC plots 
#------------------------------------------------------------#

# Function to generate ROC plot for a gene
generate_roc_plot <- function(protein_name, data) {
  rocc <- roc(data$GeneXpert_Ultra, as.numeric(data[[protein_name]]))  # Generate ROC curve gene_name
  
  # Calculate AUC and 95% CI using bootstrap with 5000 iterations
  ci_roc <- ci.auc(rocc, method = "bootstrap", boot.n = 5000)
  
  # Create ggplot object for ROC curve
  ggroc(rocc, legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", color = "blue") +  # Add diagonal line
    labs(title = protein_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(size = 1.2, color = "black"),
          axis.ticks = element_line(size = 1)) +
    
    # Annotate with AUC and CI
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocc), 2), 
                           "\n95% CI:", 
                           round(ci_roc[1], 2), "-", 
                           round(ci_roc[3], 2)), size = 5)
}

# Step 1: Generate ROC plots for each gene
roc_1 <- generate_roc_plot("SAA1", Baseline_All)
roc_2 <- generate_roc_plot("CRP", Baseline_All)
roc_3 <- generate_roc_plot("IP.10", Baseline_All)

# Step 2: Arrange the ROC plots in a 2x2 grid and save to pdf 
pdf("5. ROCC of protein markers based on Xpert.pdf", width = 10, height = 10)  # Start PDF device
ROCC_proteins_based_on_Xpert = grid.arrange(roc_1, roc_2, roc_3, ncol = 2, nrow = 2)
dev.off() # Close PDF device



# 3. Roc details 
#------------------------------------------------------------------------- #

# Function to perform ROC analysis and calculate Youden's index, sensitivity, and specificity
roc_Xpert_details <- function(gene_name, data) {
  # Calculate ROC curve
  roc_curve <- roc(data$GeneXpert_Ultra, as.numeric(data[[gene_name]]))
  
  # Calculate Youden's index and optimal cutoff
  youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  
  # Extract values for Youden's cutoff, sensitivity, and specificity
  youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
  sensitivity <- as.numeric(youden_coords["sensitivity"])
  specificity <- as.numeric(youden_coords["specificity"])
  
  # Compute 95% CI for sensitivity and specificity at Youden's cutoff
  ci_sensitivity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
  ci_specificity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)
  
  # Extract CI values
  sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
  sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
  specificity_lower_ci <- ci_specificity[["specificity"]][1]
  specificity_upper_ci <- ci_specificity[["specificity"]][3]
  
  # Store results in a dataframe
  results_dfX <- data.frame(
    Gene = gene_name,
    Youden_Cutoff = youden_cutoff,
    Sensitivity = sensitivity,
    Sensitivity_CI_Lower = sensitivity_lower_ci,
    Sensitivity_CI_Upper = sensitivity_upper_ci,
    Specificity = specificity,
    Specificity_CI_Lower = specificity_lower_ci,
    Specificity_CI_Upper = specificity_upper_ci
  )
  
  return(results_dfX)
}

# List of genes to analyze
genes <- c("SAA1", "CRP", "IP.10")

# Initialize an empty list 
protein_listx <- list()

# Loop through each gene and perform ROC analysis
for (i in genes) {
  protein_listx[[i]] <- roc_Xpert_details(i, Baseline_All)
}
# Combine all results into a single data frame
protein_dfx <- do.call(rbind, protein_listx)

# Display the results dataframe
print(protein_dfx)
write.csv(protein_dfx, "Roc_protein_dfxpert.csv")


# 4. TPP 
#------------------------------------------------------------------------- #
# Create a list of variable names
proteins <- c("SAA1", "CRP", "IP.10")

# Define a function to perform the ROC analysis and extract statistics
XpertTPP <- function(protein_name, data, fixed_spec = 0.70) {
  # Generate ROC curve
  roc_obj <- roc(data$GeneXpert_Ultra, data[[protein_name]])
  all_coords <- coords(roc_obj, ret = c("threshold", "sensitivity", "specificity"))
  closest_idx <- which.min(abs(all_coords$specificity - fixed_spec))
  
  # Get the corresponding threshold, sensitivity, and specificity
  closest_threshold <- all_coords$threshold[closest_idx]
  closest_sensitivity <- all_coords$sensitivity[closest_idx]
  closest_specificity <- all_coords$specificity[closest_idx]
  
  # CI for the specificity at the closest threshold
  ci_sens <- ci.se(roc_obj, thresholds = closest_threshold)
  
  # Check if row with name '0.9' exists in ci_spec
  if ("0.7" %in% rownames(ci_sens)) {
    row_0_7 <- ci_sens["0.7", ]
    
    # Get column 1 and 3 values at row '0.9'
    ci_0_7_lower <- row_0_7[1]  # Column 1
    ci_0_7_upper <- row_0_7[3]  # Column 3
  } else {
    ci_0_7_lower <- NA
    ci_0_7_upper <- NA
  }
  
  # Return the results as a data frame
  return(data.frame(
    Variable = protein_name,
    Threshold = closest_threshold,
    Sensitivity = closest_sensitivity,
    Specificity = closest_specificity,
    CI_Lower_0.7 = ci_0_7_lower,
    CI_Upper_0.7 = ci_0_7_upper
  ))
}

# Initialize an empty list 
resultx <- list()

# Loop through each gene and perform ROC analysis
for (i in proteins) {
  resultx[[i]] <- XpertTPP(i, Baseline_All)
}
# Combine all results into a single data frame
resultxdf <- do.call(rbind, resultx)

# Display the results dataframe
print(resultxdf)
write.csv(resultxdf, "Roc_protein_TPPdfxpert.csv")



# Overlay GeneXpert ROCC 
# ------------------------------------------------------------------#
# GeneXpert_Ultra 

# Function to generate overlay ROC plot for a gene with DeLong test
GeneXpert_overlay <- function(protein_name, data1, data2) {
  
  # Generate ROC curves for both datasets
  roc1 <- roc(data1$GeneXpert_Ultra, as.numeric(data1[[protein_name]]))  # ROC for AllData
  roc2 <- roc(data2$GeneXpert_Ultra, as.numeric(data2[[protein_name]]))  # ROC for AllNeg
  
  # Calculate AUC and 95% CI for both datasets
  ci_roc1 <- ci.auc(roc1, method = "bootstrap", boot.n = 5000)
  ci_roc2 <- ci.auc(roc2, method = "bootstrap", boot.n = 5000)
  
  # Perform DeLong test to compare the ROC curves
  delong_test <- roc.test(roc1, roc2, method = "delong")
  delong_p_value <- round(delong_test$p.value, 3)  # Round the p-value to 3 decimal places
  
  # Create ggplot object for ROC curve of AllData
  plot <- ggroc(roc1, legacy.axes = TRUE, col = "black", size = 1) +
    geom_line(data = ggroc(roc2)$data, aes(x = 1 - specificity, y = sensitivity), 
              color = "#FF5733", size = 1) +  # Overlay ROC for AllNeg
    geom_abline(linetype = "dashed", color = "gray") +  # Add diagonal line
    labs(title = protein_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)) +
    # Annotate with AUC and CI for both datasets
    annotate("text", x = 0.55, y = 0.20,
             label = paste0("AUC: ", round(auc(roc1), 2), 
                            " (CI ", round(ci_roc1[1], 2), "-", round(ci_roc1[3], 2), ")"), 
             size = 4.5, color = "black") +
    
    annotate("text", x = 0.55, y = 0.15,
             label = paste0("AUC: ", round(auc(roc2), 2), 
                            " (CI ", round(ci_roc2[1], 2), "-", round(ci_roc2[3], 2), ")"), 
             size = 4.5, color = "#FF5733") +
    
    # Annotate with DeLong p-value
    annotate("text", x = 0.55, y = 0.10,
             label = paste0("DeLong p = ", delong_p_value), 
             size = 4.5, color = "blue")  # You can change color if needed
  
  return(plot)
}

# Generate overlay ROC plots for each protein
overlay_a <- GeneXpert_overlay("SAA1", Baseline_All, SeroNegMBTData)
overlay_b <- GeneXpert_overlay("CRP", Baseline_All, SeroNegMBTData)
overlay_c <- GeneXpert_overlay("IP.10", Baseline_All, SeroNegMBTData)

# Arrange the plots in a grid and save to TIFF file
tiff("Overlay_Roc_MBT_Xpert.allvsSeroNeg_TB_and_ORD.tiff", 
     width = 10, height = 10, units = "in", res=300)
grid.arrange(overlay_a, overlay_b, overlay_c, ncol = 2)
dev.off()
```



#------------------------------------------------------------------------- #
# Part 3: Treatment response
#------------------------------------------------------------------------- #
```{r part tretament to review}
names(Baseline_all_Data)
class(Baseline_all_Data$Timepoint)
unique(Baseline_all_Data$Timepoint)

# stats 
# Timepoint 
Time1 <- kruskal.test(SAA1 ~ Timepoint, data = All_timepoint) #  p-value = 6.389e-15
Time2 <- kruskal.test(CRP ~ Timepoint, data = All_timepoint) # p-value = 0.02165
Time3 <- kruskal.test(IP.10 ~ Timepoint, data = All_timepoint) # p-value = 0.0003319

# Smear grade ????
smear1 <- kruskal.test(SAA1 ~ Smear_grade, data = TB_smear) 
smear2 <- kruskal.test(CRP ~ Smear_grade, data = TB_smear) 
smear3 <- kruskal.test(IP.10 ~ Smear_grade, data = TB_smear) # ns

# Xpert grade 
xp1 <- kruskal.test(SAA1 ~ Xpert_grade, data = TB_baseline)
xp2 <- kruskal.test(CRP ~ Xpert_grade, data = TB_baseline)
xp3 <- kruskal.test(IP.10 ~ Xpert_grade, data = TB_baseline) # ns


# 
# Plots 
#------------------------------------------------------------------------- #

#save.image("StopMBT180924.RData")

# tx data start at All_timepoint
table(All_timepoint$Timepoint)
names(All_timepoint)
dim(All_timepoint)
names(All_timepoint)
table(All_timepoint$Timepoint)

# Propagate TB status to the data 
RxData = All_timepoint
RxData <- RxData %>%
  group_by(ParticipantID) %>%
  mutate(Binary = dplyr::first(Binary[Timepoint == "Baseline"], default = NA)) %>%
  ungroup() 
dim(RxData)
table(RxData$Timepoint)

# Remove NA in protein Data 
RxData = RxData[!is.na(RxData$SAA1), ]
RxData <- RxData[!is.na(RxData$Binary), ]

# subset only TB 
RxDataTb = subset(RxData, Binary %in% c("TB"))
names(RxDataTb)
table(RxDataTb$Timepoint) # Perfect 



#-----------------------------------------------------------------------#
# Create summary table with means, Q1, Q3, and p-values
summary_table_time <- RxDataTb %>%
  select(SAA1, CRP, IP.10, Timepoint) %>%
  tbl_summary(
    by = Timepoint, 
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})" # Shows Mean (Q1, Q3)
    ),
    digits = all_continuous() ~ 2 # Rounds numbers to 2 decimal places
  ) %>%
  add_p(test = all_continuous() ~ "wilcox.test") # Adds p-values from Wilcoxon test

# Print the gtsummary table
summary_table_time
#write.csv(summary_table_time, "Protein by timepoint.csv")


#-----------------------------------------------------------------------#
# plot 

# Function to plot boxplots and perform pairwise comparisons with separate y-axis limits
plot_protein_over_time <- function(data, protein, y_limits) {
  
  # Define pairwise comparisons (adjust these if you have different timepoints)
  my_comparisons <- list(c("Baseline", "Month 4"),
                         c("Baseline", "Month 6"),
                         c("Month 4", "Month 6"))
  
  # Custom y positions for p-value labels based on the IQR of the protein
  label_y_positions <- c(y_limits[2] - 0.2, 
                         y_limits[2] - 0.1, 
                         y_limits[2] - 0.05) 
  
  # Boxplot with pairwise comparisons and custom p-value label positions
  p <- ggplot(data, aes(x = Timepoint, y = .data[[protein]])) +
    geom_boxplot(aes(color = Timepoint), width = 0.6, size = 1.2) +
    scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C")) +  
    scale_y_continuous(limits = y_limits) + 
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.8, color = "black"),
          axis.ticks = element_line(size = 0.8)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif",
                       label.y = label_y_positions) + 
    labs(title = paste(protein, "Expression Over Time"), x = "Timepoint", y = paste(protein, "Expression"))
  return(p)
}

# Plot 
SAA_filtered = RxDataTb[!RxDataTb$SAA1==47.630,] # remove extreme value 
dim(SAA_filtered)

r1 = plot_protein_over_time(SAA_filtered, "SAA1", c(0, 2))
r2 = plot_protein_over_time(RxDataTb, "CRP", c(0, 2))
r3 = plot_protein_over_time(RxDataTb, "IP.10", c(0, 3))

# Arrange the plots in a grid
tiff("MBT_Binary_over_time.tiff", width = 10, height = 10, res = 300, unit = "in")
Binary_TR <- grid.arrange(r1, r2, r3, ncol = 2, nrow = 2)
dev.off()


save.image("Stop_21092024.RData")

# Load or rerun script 
#load("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task3_MBT_data_cleaning/Stop_21092024.RData") 
```




# SECTION 2: MBT MODEL 


#---------------------------------------------------------------------------
# SET DIRECTORY
#---------------------------------------------------------------------------

```{r}
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task5_Data_analyses/")
```

```{r load enviroment}
#LOAD DATA ENVIRONMENT FROM RDA
#---------------------------------------------------------------------------
load("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task3_MBT_data_cleaning/MBT_analyses_part2.RData") 
```


```{r LOAD LIBRARIES}
#----------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pROC)
library(tidyr)
library(gridExtra)
library(caret)
library(DMwR)        # For SMOTE
library(randomForest) # For Random Forest
library(xgboost)     # For XGBoost
library(e1071)       # For SVM
library(rpart)       # For Decision Tree
library(ROSE)

# load("MBT_analyses_part1.RData") 
```


# DEFINE DATASETS 
```{r datasets}

# data based on Xpert 
Baseline_All

# data based on comp binary BL
Baseline_all_Data # no NA 

# Xpert based
Baseline_All_data_xpert

# Within TB only BL
TB_baseline_Xpert
# Within ORD ony BL
ORD_baseline_Xpert
```



#------------------------------------------------------------------------------------------------#
#                                           MODELLING   
#------------------------------------------------------------------------------------------------#
```{r IDENTIFY DATA}
Baseline_all_Data

Dt = Baseline_all_Data
dim(Dt)
names(Dt)
names(Dt)[16] = "Group"  # Binary renamed as group

# SUBSET REQUIRED COLUMNS  
table(Dt$Group)
names(Dt)

#subset only the quant data 
Dt = Dt[, c("Group", "SAA1","CRP", "IP.10")]


# FORMAT DATA 
# Create model with caret 
table(Dt$Group)
Dt$Group = as.factor(as.character(Dt$Group))
str(Dt)
summary(Dt)
# Reorder levels of Group
Dt$Group <- factor(Dt$Group, levels = c("TB", "ORD"))

# Check for missing values
sum(is.na(Dt))
sum(is.na(Dt$Group))
```



```{r SPLIT DATA and set train control} 
# ------------------------------------------------------------------#
# 1. Split data ----
set.seed(123)  # Ensure reproducibility
train_index <- createDataPartition(Dt$Group, p = 0.7, list = FALSE)
train_data <- Dt[train_index, ]
test_data <- Dt[-train_index, ]
dim(train_data)
dim(test_data)


# SET TRAINING CONTROL  
# ------------------------------------------------------------------#
# Set up 10-fold cross-validation
train_control <- trainControl(method = "cv", number = 10, 
                              classProbs = TRUE, 
                              summaryFunction = twoClassSummary)

dim(train_data) # 211
dim(test_data) # 89
```


# Train models 
```{r train and select method}
# i. Random Forest 
set.seed(123)
rf_model <- train(Group ~ SAA1 + CRP + IP.10, 
                  data = train_data, 
                  method = "rf", 
                  trControl = train_control, 
                  metric = "ROC")
rf_model$results

# ii. XGB
set.seed(123)
xgb_model <- train(Group ~ SAA1 + CRP + IP.10, 
                   data = train_data, 
                   method = "xgbTree", 
                   trControl = train_control, 
                   metric = "ROC")

# iii. SVM
set.seed(123)
svm_model <- train(Group ~ SAA1 + CRP + IP.10, 
                   data = train_data, 
                   method = "svmRadial", 
                   trControl = train_control, 
                   metric = "ROC")

# iv. Decision trees
set.seed(123)
dt_model <- train(Group ~ SAA1 + CRP + IP.10, 
                  data = train_data, 
                  method = "rpart", 
                  trControl = train_control, 
                  metric = "ROC")


# 17 dec we decided to add a few linear models 

# Regular logistic regression model
glm_model <- train(Group ~ SAA1 + CRP + IP.10, 
                   data = train_data, 
                   method = "glm", 
                   family = "binomial", 
                   trControl = train_control, 
                   metric = "ROC")

glmnet_model <- train(Group ~ SAA1 + CRP + IP.10, 
                   data = train_data, 
                   method = "glmnet", 
                   family = "binomial", 
                   trControl = train_control, 
                   metric = "ROC")

# Train a linear XGBoost model 
xgb_linear <- train(
  Group ~ SAA1 + CRP + IP.10,  
  data = train_data,
  method = "xgbLinear",                 # Use XGBoost linear booster
  trControl = train_control,            # Cross-validation setup
  metric = "ROC",                       # Optimize for ROC AUC
  verbose = FALSE)

# Ridge Regression model 
set.seed(123)
ridge_model <- train(
  Group ~ SAA1 + CRP + IP.10,
  data = train_data,
  method = "glmnet",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(alpha = 0,  # Ridge regression
                         lambda = seq(0.0001, 1, length = 100)) # Range of penalties
)
print(ridge_model)

library(caret)
set.seed(123)

# Train Elastic Net model
elastic_net_model <- train(
  Group ~ SAA1 + CRP + IP.10,  # Formula
  data = train_data,
  method = "glmnet",                    # Use glmnet for Elastic Net
  trControl = train_control,            # Cross-validation setup
  metric = "ROC",                       # Optimize for ROC AUC
  tuneGrid = expand.grid(
    alpha = seq(0, 1, by = 0.1),        # Range of alpha (0 = Ridge, 1 = Lasso)
    lambda = seq(0.0001, 1, length = 50) # Range of lambda (penalty strength)
  )
)

# Print the model summary
print(elastic_net_model)

# Plot to visualize optimal alpha and lambda
plot(elastic_net_model)



# 4. Compare models in training set 
results <- resamples(list(RandomForest = rf_model, 
                          Logistic = glm_model,
                          Regression = glmnet_model,
                          XGBoost = xgb_model, 
                          XG_Linear = xgb_linear,
                          Ridge_regression = ridge_model,
                          EN = elastic_net_model,
                          SVM = svm_model, 
                          DecisionTree = dt_model))

# Summary of model performances
summary(results)
Rankings_df <- as.data.frame(summary(results)$statistics)
write.csv(Rankings_df, "modelrankings.csv")

# Boxplot to visualize comparison
bwplot(results)
```


# Method 2: fit a glm model 

```{r logistic regression}
# Regular logistic regression model
set.seed(123)
glm_model <- train(Group ~ SAA1 + CRP + IP.10, 
                   data = train_data, 
                   method = "glm", 
                   family = "binomial", 
                   trControl = train_control, 
                   metric = "ROC")
```



```{r eval training}
# turns out we overide the model to the model$finalmodel value which still works

# Note set up of the model 
# 10-fold cross-validation
# train_control <- trainControl(method = "cv", number = 10, 
#                               classProbs = TRUE, 
#                               summaryFunction = twoClassSummary)
# dim(train_data) # 211
# dim(test_data) # 89
# # Regular logistic regression model
# set.seed(123)
# model.glm <- train(Group ~ SAA1 + CRP + IP.10, 
#                    data = train_data, 
#                    method = "glm", 
#                    family = "binomial", 
#                    trControl = train_control, 
#                    metric = "ROC")

model.glm$results

# Class prediction
pred_trdata <- predict(model.glm, newdata = train_data) 
pred_trdata <- factor(pred_trdata, levels = levels(train_data$Group))

# probability
pred_train_prob <- predict(model.glm, newdata = train_data, type = "prob")


# Compute the confusion matrix: predict class vs true class
cmtrain <- confusionMatrix(pred_trdata, train_data$Group)
print(cmtrain)


# ROC 
# positive class 
# Extract probability of the positive class 
tb_predictions <- pred_train_prob[, "TB"]  # Extracts only the "TB" probability

# Compute ROC curve
trainroc_glm <- roc(train_data$Group, tb_predictions)

# Plot ROC curve
plot(trainroc_glm)

# Print AUC
auc(trainroc_glm)

# ROC plot 
trainplot = plot(trainroc_glm, col = "red", lwd = 2, main = "ROC Curves for Logistic regression model", 
                    legacy.axes = TRUE) +
  
  # Add AUC values
  legend("bottomright", legend = c(
    paste("Regression Model: AUC =", round(auc(trainroc_glm), 2))))
```


# PREDICTIONS; Test set 
```{r predict}
# Class prediction
glm_predictions <- predict(glm_model, newdata = test_data) # 89
class(glm_predictions)

# probability
#glm_prob <- predict(glm_model, newdata = test_data, type = "prob")[, "TB"] #89. didnt work second time

# resolve 
glm_prob <- predict(glm_model, newdata = test_data, type = "prob")

glm_prob = glm_prob[, "TB"]

# Confusion matrix to assess performance
cm = confusionMatrix(glm_predictions, test_data$Group) #

# ROC 
glm_roc <- roc(test_data$Group, glm_prob)

# ROC plot 
# glm_rocc = plot(glm_roc, col = "red", lwd = 2, main = "ROC Curves for Logistic regression model", legacy.axes = TRUE) +
# 
# # Add AUC values
# legend("bottomright", legend = c(
#   paste("Regression Model: AUC =", round(auc(glm_roc), 2))))

# amended
# Plot the ROC curve
plot(glm_roc, col = "red", lwd = 2, 
     main = "ROC Curve for Logistic Regression Model", 
     legacy.axes = TRUE)

# Add AUC to legend
legend("bottomright", 
       legend = paste("Regression Model: AUC =", round(auc(glm_roc), 2)), 
       col = "red", lwd = 2)

# save model
#saveRDS(glm_roc, "Glm_Model_Object.RDS")
```


# Prediction in the training set 
```{r predict}
# Class prediction
glm_train <- predict(glm_model, newdata = train_data) # prob
# Predict probabilities for the "TB" class
glm_train_prob <- predict(glm_model, newdata = train_data, type = "response")
cm_train = confusionMatrix(glm_train, train_data$Group) #

# ROC 
glm_roc_train <- roc(train_data$Group, glm_train_prob)

# ROC plot 
glm_train_roc = plot(glm_roc_train, col = "red", lwd = 2, main = "ROC Curves for Logistic regression model in training set", legacy.axes = TRUE) +

# Add AUC values
legend("bottomright", legend = c(
  paste("Regression Model: AUC =", round(auc(glm_roc), 2))))
```


# ROC CURVE 

# linear roc curve not smoothed 
```{r}

set.seed(1234)
library(pROC)

# Define the output TIFF file with specified resolution and size
tiff("MBT_rocplot.tiff", width = 6, height = 6, units = "in", res = 300)

# Create ROC object
roc_curve <- roc(test_data$Group, glm_prob)

# Calculate AUC and its confidence interval
auc_value <- auc(roc_curve)
auc_int <- ci.auc(roc_curve) 
# Calculate confidence interval for AUC
lower_ci <- auc_int[1]
upper_ci <- auc_int[3]

# Plot ROC curve
plot(roc_curve, main = NULL, col = "black", lwd = 2)

# Add shaded confidence interval
plot(ci_roc, type = "shape", col = rgb(0.1, 0.4, 0.8, 0.2), add = TRUE)

# Add legend with AUC and confidence interval
legend("bottomright", legend = c(
  paste("Regression Model: AUC =", round(auc_value, 2), 
        "[", round(lower_ci, 2), "-", round(upper_ci, 2), "]")
))

# Close the graphics device to save the file
dev.off()
```


# Get plot without shading 
```{r}
# Plot ROC curve
mbt_plot_true = plot(roc_curve, main = NULL, col = "black", lwd = 2)
# Add shaded confidence interval : un-comment to add CI shading
# plot(ci_roc, type = "shape", col = rgb(0.1, 0.4, 0.8, 0.2), add = TRUE)

# Add legend with AUC and confidence interval
legend("bottomright", legend = c(
  paste("Regression Model: AUC =", round(auc_value, 2), 
        "[", round(lower_ci, 2), "-", round(upper_ci, 2), "]")
))
mbt_plot_true
```


# GET PERFORMANCE METRICS OF THE MODEL  
```{r metrics}
# Get metrics for all thresholds
roc_metrics <- coords(roc_curve, x = "all", ret = c("threshold", "sensitivity", "specificity"))

# View metrics
print(roc_metrics)

# Example: Get metrics for a specific threshold (e.g., Youden's J index optimal threshold)
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
print(optimal_coords)
```

# Get CI values 
```{r ci roc}
# Create ROC object
roc_curve <- roc(test_data$Group, glm_prob)

# Calculate Youden's index and optimal cutoff
youden_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Extract values for Youden's cutoff, sensitivity, and specificity
youden_cutoff <- as.numeric(youden_coords["threshold"])  # Ensure this is numeric
sensitivity <- as.numeric(youden_coords["sensitivity"])
specificity <- as.numeric(youden_coords["specificity"])

# Compute 95% CI for sensitivity and specificity at Youden's cutoff
ci_sensitivity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "sensitivity", boot.n = 5000)
ci_specificity <- ci.coords(roc_curve, x = youden_cutoff, input = "threshold", ret = "specificity", boot.n = 5000)

# Extract CI values
sensitivity_lower_ci <- ci_sensitivity[["sensitivity"]][1]
sensitivity_upper_ci <- ci_sensitivity[["sensitivity"]][3]
specificity_lower_ci <- ci_specificity[["specificity"]][1]
specificity_upper_ci <- ci_specificity[["specificity"]][3]

# Store results in a dataframe
results_dfr <- data.frame(
  Youden_Cutoff = youden_cutoff,
  Sensitivity = sensitivity,
  Sensitivity_CI_Lower = sensitivity_lower_ci,
  Sensitivity_CI_Upper = sensitivity_upper_ci,
  Specificity = specificity,
  Specificity_CI_Lower = specificity_lower_ci,
  Specificity_CI_Upper = specificity_upper_ci
)

print(results_dfr)

# Save the data frame as a CSV file
write.csv(results_dfr, "MBT_metrics.csv", row.names = FALSE)

# Print a message confirming the file was saved
cat("Metrics saved as 'MBT_metrics.csv'.\n")
```



# At TPP 
```{r tpp}
# # Create ROC object
# roc_curve <- roc(test_data$Group, glm_prob)
# 
# # Get thresholds, sensitivities, and specificities
# all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))
# 
# # Find the closest specificity to 70%
# fixed_spec <- 0.70
# closest_idx <- which.min(abs(all_coords$specificity - fixed_spec))
# 
# # Get the corresponding threshold and sensitivity
# closest_threshold <- all_coords$threshold[closest_idx]
# closest_sensitivity <- all_coords$sensitivity[closest_idx]
# closest_specificity <- all_coords$specificity[closest_idx]
# 
# # Compute 95% CI for sensitivity at the fixed specificity using ci.se()
# ci_sens <- ci.se(roc_curve, specificities = fixed_spec, boot.n = 2000)
# sensitivity_ci_lower <- ci_sens[1]
# sensitivity_ci_upper <- ci_sens[3]
# 
# # Store results in a data frame
# mbt_tpp <- data.frame(
#   Threshold = closest_threshold,
#   Sensitivity = closest_sensitivity,
#   Sensitivity_CI_Lower = sensitivity_ci_lower,
#   Sensitivity_CI_Upper = sensitivity_ci_upper,
#   Specificity = closest_specificity)
# 
# # print results 
# print(mbt_tpp)
# 
# # Save the results to a CSV file if needed
# write.csv(mbt_tpp, "Roc_details_MBT_TPP.csv") #redone correctly
```

```{r tpp at 90SE}
# Create ROC object
roc_curve <- roc(test_data$Group, glm_prob)

# Get thresholds, sensitivities, and specificities
all_coords <- coords(roc_curve, ret = c("threshold", "sensitivity", "specificity"))

# Find the closest sensitivity to 90%
fixed_sens <- 0.90
closest_idx <- which.min(abs(all_coords$sensitivity - fixed_sens))

# Get the corresponding threshold and specificity
closest_threshold <- all_coords$threshold[closest_idx]
closest_sensitivity <- all_coords$sensitivity[closest_idx]
closest_specificity <- all_coords$specificity[closest_idx]

# Compute 95% CI for sensitivity at the fixed sensitivity using ci.se()
ci_spec <- ci.sp(roc_curve, sensitivities = fixed_sens, boot.n = 2000)
spec_ci_lower <- ci_spec[1]
spec_ci_upper <- ci_spec[3]

# Store results in a data frame
mbt_tpp <- data.frame(
  Threshold = closest_threshold,
  Specificity = closest_specificity,
  spec_ci_lower = spec_ci_lower,
  spec_ci_upper = spec_ci_upper,
  Sensitivity = closest_sensitivity)

# Print results
print(mbt_tpp)

# Save the results to a CSV file if needed
write.csv(mbt_tpp, "Roc_details_MBT_TPP.csv")
```


# save 
```{r save objects}
# Load r environemnt from part 2 

#save.image("MBT_analyses_FINAl_stop_on_120125.RData")

# load("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task5_Data_analyses/MBT_analyses_FINAl_stop_on_120125.RData")

```



# Before Overlay. Can we check performnce in HIV negatives only 
Baseline_all_Data

names(Baseline_all_Data)

```{r}
# HIV seroneg Mbt 
# excluding hiv pos 
HIVnegIDs <- readRDS("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task5_Data_analyses/AllSeronegativesIDs.RDS")
dim(Baseline_all_Data)
dim(HIVnegIDs)
SeroNegMBTData = merge(Baseline_all_Data, HIVnegIDs, by="ParticipantID")

# 
dim(SeroNegMBTData)

Seronegdata = SeroNegMBTData
names(Seronegdata)
names(Seronegdata)[16] = "Group"  # Binary renamed as group

#subset only the quant data 
Seronegdata = Seronegdata[, c("Group", "SAA1","CRP", "IP.10")]
table(Seronegdata$Group)
Seronegdata$Group = as.factor(as.character(Seronegdata$Group))
# Reorder levels of Group
Seronegdata$Group <- factor(Seronegdata$Group, levels = c("TB", "ORD"))

# Check for missing values
sum(is.na(Seronegdata))
sum(is.na(Seronegdata$Group))
```


# SPLIT DATA 
```{r}
# 1. Split data ----
set.seed(123)  # Ensure reproducibility
new.train_index <- createDataPartition(Seronegdata$Group, p = 0.7, list = FALSE)
new.train_data <- Seronegdata[new.train_index, ]
new.test_data <- Seronegdata[-new.train_index, ]

dim(new.train_data)
dim(test_data)
```


# Test model 
```{r}
#  prediction
model.glm
pred_SN <- predict(model.glm, newdata = new.test_data)
pred_SN <- factor(pred_SN, levels = levels(new.test_data$Group))
# probability
pred_SNpro <- predict(model.glm, newdata = new.test_data, type = "prob")

# Compute the confusion matrix
cm.sn <- confusionMatrix(pred_SN, new.test_data$Group)
print(cm.sn)
```

# Roc plot 
```{r}
# positive class 
# Extract probability of the positive class 
pred_sn_tb <- pred_SNpro[, "TB"]  
# Compute ROC curve
roc_seroneg <- roc(new.test_data$Group, pred_sn_tb)

# Plot ROC curve
plot(roc_seroneg)
auc(roc_seroneg)
plot(roc_seroneg)


# ROC plot 
rocc_seronegs = plot(roc_seroneg, col = "red", lwd = 2, main = "ROC Curves for Logistic regression model",legacy.axes = TRUE) +
  
  # Add AUC values
  legend("bottomright", legend = c(
    paste("Regression Model: AUC =", round(auc(roc_seroneg), 2))))
```



# Combined mbt plots 
```{r mbt plots}
# function to genrate overlay roc then plot the overlap roc 

# Function to overlay two ROC curves and perform DeLong test
Overlay_ROC <- function(roc1, roc2, label1 = "ROC 1", label2 = "ROC 2") {
  
  # Calculate AUC and 95% CI for both ROC curves
  ci_roc1 <- ci.auc(roc1, method = "bootstrap", boot.n = 5000)
  ci_roc2 <- ci.auc(roc2, method = "bootstrap", boot.n = 5000)

  # Perform DeLong test
  delong_test <- roc.test(roc1, roc2, method = "delong")
  delong_p_value <- round(delong_test$p.value, 3)  # Round to 3 decimals

  # Create ggplot object for ROC curves
  plot <- ggroc(roc1, legacy.axes = TRUE, col = "black", size = 1) +
    geom_line(data = ggroc(roc2)$data, aes(x = 1 - specificity, y = sensitivity), 
              color = "#FF5733", size = 1) +  # Overlay second ROC curve
    geom_abline(linetype = "dashed", color = "gray") +  # Diagonal reference line
    labs(title = "Overlayed ROC Curves", x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1)
    ) +
    
    # Annotate with AUC and CI for both ROC curves
    annotate("text", x = 0.55, y = 0.20, 
             label = paste0(label1, " AUC: ", round(auc(roc1), 2), 
                            " (CI ", round(ci_roc1[1], 2), "-", round(ci_roc1[3], 2), ")"), 
             size = 4.5, color = "black") +
    
    annotate("text", x = 0.55, y = 0.15, 
             label = paste0(label2, " AUC: ", round(auc(roc2), 2), 
                            " (CI ", round(ci_roc2[1], 2), "-", round(ci_roc2[3], 2), ")"), 
             size = 4.5, color = "#FF5733") +
    
    # Annotate with DeLong p-value
    annotate("text", x = 0.55, y = 0.10, 
             label = paste0("DeLong p = ", delong_p_value), 
             size = 4.5, color = "blue")

  return(plot)
}

# Generate ROC curves
roc_seroneg <- roc(new.test_data$Group, pred_sn_tb)  # First ROC curve
# Second ROC curve
glm_roc <- roc(test_data$Group, glm_prob)

# Call function to plot the overlayed ROC curves
# Arrange the plots in a grid and save to TIFF file
tiff("Overlay_MBTserostatus_rocs.tiff", 
     width = 10, height = 10, units = "in", res=300)
Overlay_MBTserostatus_roc = Overlay_ROC(roc_seroneg, roc_curve, 
                                        label1 = "Seroneg ROC", label2 = "GLM ROC")
Overlay_MBTserostatus_roc
dev.off()
```

#-----------------------------------------------------------------------------------# 
# Overlay all BMT plots 
```{r overlay Mbt plots}
mbtroc1= overlay_1
mbtroc2= overlay_2
mbtroc3= overlay_3
mbtroc4 = Overlay_MBTserostatus_roc

# Arrange the plots in a grid and save to TIFF file
tiff("Overlay_all_MBT_plots.tiff", 
     width = 10, height = 10, units = "in", res=300)
grid.arrange(mbtroc1, mbtroc2, mbtroc3, mbtroc4, ncol = 2)
dev.off()

# for report only 
# Arrange the plots in a grid and save to TIFF file
png("Overlay_all_MBT_plots.png", 
     width = 10, height = 10, units = "in", res=300)
grid.arrange(mbtroc1, mbtroc2, mbtroc3, mbtroc4, ncol = 2)
dev.off()
```




#------------------------------------------------------------------------------------------------------------#
# SECTION 4: PATHOGEN IMPACT ON POC PERFORMANCE 
#------------------------------------------------------------------------------------------------------------#


# I. CEPHEID HOST RESPONSE
#---------------------------------------------#
```{r load libraries and set wd}
library(readxl)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(janitor)
library(ggpubr)
library(gridExtra)

# Set dir 
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/Pathogen_Data")
```

# Define datsets
```{r define datasets}                        
# 1. Pathogen Data 
PathData <- read_excel("Pathogen_Data_09092024.xlsx")
PathDat = PathData[, c(1, 31:33)] # colnames: Mixed" "VIRUS_only" "BACTERIA_only"  
# rename column
names(PathDat)[1] = "ParticipantID"

# 2. Merge with Cepheid Data for 3 gens comparisons
  # Merge with Data binary for diagnostic groups : colname Group
Composite_data <- readRDS("Composite_data_table.RDS")
names(Composite_data)

MergedpathDatasets = merge(PathDat, Composite_data, by="ParticipantID")
dim(MergedpathDatasets)
names(MergedpathDatasets)
# Select columns
MergedpathData = MergedpathDatasets[, c(1:4,23,21)] # data with group var
names(MergedpathData)
dim(MergedpathData)

# 3. MBT: we will apply this to the test set in MBT analyses 

# 4: COVID data 
names(Ceph_Demo)
Covid_data = Ceph_Demo[, c("ParticipantID", "Combined_COVID")]
names(Covid_data)[2] = "COVID"

# Merger MergedpathData and COVID data
Cephpathcovid = merge(MergedpathData, Covid_data, by="ParticipantID")
names(Cephpathcovid)
dim(Cephpathcovid)
#[1] 282   7

# Safe as csv
# write.csv(Cephpathcovid, "pathogenHRcovid.csv")

# load modified data 
Cephpath_complete = read.csv("pathogenHRcovid_amended.csv")
names(Cephpath_complete)
Pathogen_subset = Cephpath_complete[, c(1,6, 12, 5)]
dim(Pathogen_subset) # [1] 282   4

# Merge with DefiniteBinarySet to remove unknown TBs
names(DefiniteBinarySet)
table(DefiniteBinarySet$Composite_Ref)
table(DefiniteBinarySet$BinaryClass)
UnknowTB.ids = DefiniteBinarySet[, c("ParticipantID", "BinaryClass")]
# Merge step
Pathogen_final = merge(Pathogen_subset, UnknowTB.ids, by="ParticipantID")
dim(Pathogen_final) # [1] 272   5
# this mean the 10, unknown TB are excluded from the analyses 

# format data 
Pathogen_final$Pathogens = as.factor(Pathogen_final$Pathogens)
# table(Pathogen_final$Group) # ORD  TB 176  96

# Reorder factor levels
# gsub the pathogen groups 
Pathogen_final$Pathogens = gsub("_", "+", Pathogen_final$Pathogens)

Pathogen_final$Pathogens <- factor(Pathogen_final$Pathogens, 
                                   levels = c("TB", "ORD", "ORD+Bacteria", "ORD+Virus", "ORD+Mixed"))
```


# compare TB score across pathogen groups 
```{r compare}
# Comparisonsof expressions 
library(ggplot2)
library(ggpubr)

# Define comparisons
comparisons <- list(
  c("ORD", "ORD+Bacteria"),
  c("ORD", "ORD+Mixed"),
  c("ORD", "ORD+Virus"),
  c("ORD", "TB"),
  c("ORD+Bacteria", "ORD+Mixed"),
  c("ORD+Bacteria", "ORD+Virus"),
  c("ORD+Bacteria", "TB"),
  c("ORD+Mixed", "ORD+Virus"),
  c("ORD+Mixed", "TB"),
  c("ORD+Virus", "TB"))

# Plot with reordered levels
patho_plot <- ggplot(Pathogen_final, aes(x = Pathogens, y = TB_score)) + 
  geom_boxplot(aes(color = Pathogens), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2", "#FC4E07")) +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)
  ) +  
  stat_compare_means(comparisons = comparisons, 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE    # Hide non-significant results
                     #label = "p.signif"
                     ) +  
  labs(title = "TB Score Expression Across Pathogen Groups", x = "Pathogen Group", y = "TB Score")

print(patho_plot)

# To hide non-significant comparisons: since we are unable to complete this in ggplot, we 
# do this retrospectively nby removing the comparisons we already know are not significant 

# Define comparisons
comparisons <- list(
  # c("ORD", "ORD+Bacteria"),
  # c("ORD", "ORD+Mixed"),
  # c("ORD", "ORD+Virus"),
  c("ORD", "TB"),
  # c("ORD+Bacteria", "ORD+Mixed"),
  # c("ORD+Bacteria", "ORD+Virus"),
  c("ORD+Bacteria", "TB"),
 # c("ORD+Mixed", "ORD+Virus"),
  c("ORD+Mixed", "TB"),
  c("ORD+Virus", "TB"))

# Plot with reordered levels
patho_plot_hide.ns <- ggplot(Pathogen_final, aes(x = Pathogens, y = TB_score)) + 
  geom_boxplot(aes(color = Pathogens), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2", "#FC4E07")) +  
  scale_y_continuous(limits = c(-5, 5)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(comparisons = comparisons, 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE    # Hide non-significant results
                     #label = "p.signif"
                     ) +  
  labs(
    #title = "TB Score Expression Across Pathogen Groups", 
    x = "Pathogen Group", y = "TB Score")
print(patho_plot_hide.ns)
png("TB_Score_Across_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
patho_plot_hide.ns
dev.off()
```


# Idea 2. Plots within TB and Within ORD 

# Cepheid
```{r within ORD}

# subset only ord 
Pathogen_ORd = subset(Pathogen_final, Pathogens %in% 
                        c("ORD", "ORD+Bacteria", "ORD+Virus", "ORD+Mixed"))
# levels 
Pathogen_ORd$Pathogens <- as.factor(as.character(Pathogen_ORd$Pathogens, 
                                   levels = c("ORD", "ORD+Bacteria", "ORD+Virus", "ORD+Mixed")))
# Plot 
# Define comparisons
comps <- list(
  c("ORD+Bacteria", "ORD"),
  c("ORD+Mixed", "ORD"),
  c("ORD+Virus", "ORD"))

# Plot with reordered levels
HR_patho_within_ORD <- ggplot(Pathogen_ORd, aes(x = Pathogens, y = TB_score)) + 
  geom_boxplot(aes(color = Pathogens), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2")) +  
  scale_y_continuous(limits = c(-5, 2.5)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(
    #comparisons = comps, 
                     method = "kruskal.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +  
  labs(
    x = NULL, y = "TB Score")
print(HR_patho_within_ORD)

png("TB_Score_withinORD_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
HR_patho_within_ORD
dev.off()



# Within TB -----

# load amended pathogenHRcovid_amended2 which contains a new column called Pathogen_group for within TB groups 

# load modified data 
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/Pathogen_Data")
HR_patho_within_tb = read.csv("pathogenHRcovid_amended2.csv")
names(HR_patho_within_tb)

# select cols
HR_patho_within_tb = HR_patho_within_tb[, c(1,6, 13, 5)]
dim(Pathogen_subset) # [1] 282   4

# Merge with DefiniteBinarySet to remove unknown TBs
names(DefiniteBinarySet)
#UnknowTB.ids = DefiniteBinarySet[, c("ParticipantID", "BinaryClass")]

# Merge step
HR_patho_within_tb = merge(HR_patho_within_tb, UnknowTB.ids, by="ParticipantID")
dim(HR_patho_within_tb) # [1] 272   5

# format data 
HR_patho_within_tb$Pathogen_group = gsub("_", "+", HR_patho_within_tb$Pathogen_group)

# subset only ord 
HR_patho_within_tb = subset(HR_patho_within_tb, Pathogen_group %in% 
                        c("TB+only", "TB+Bacteria", "TB+Virus", "TB+Mixed"))

HR_patho_within_tb$Pathogen_group <- factor(HR_patho_within_tb$Pathogen_group, 
                              levels = c("TB+only", "TB+Bacteria", "TB+Virus", "TB+Mixed"))

# Plot 
# Define comparisons
compp <- list(
  c("TB+Bacteria", "TB+only"),
  c("TB+Mixed", "TB+only"),
  c("TB+Virus", "TB+only"))

# Plot with reordered levels
HR_patho_within_tb_plot <- ggplot(HR_patho_within_tb, aes(x = Pathogen_group, y = TB_score)) + 
  geom_boxplot(aes(color = Pathogen_group), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2")) +  
  scale_y_continuous(limits = c(-5, 0)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(
    #comparisons = compp, 
                     method = "kruskal.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +  
  labs(x = NULL,
       y = "TB Score")
print(HR_patho_within_tb_plot)

png("TB_Score_withinTB_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
HR_patho_within_tb_plot
dev.off()
```

# MBT 
```{r within mbt}

# Within ord : we already had TB vs all ords, subset to remove TB group and plot.
MBTpathogensData

MBTpathogensData$Pathogens
# subset
MBTpathogensData_ord = subset(MBTpathogensData, Pathogens %in% 
                              c("ORD" ,"ORD+Bacteria", "ORD+Virus","ORD+Mixed"))

# control levels
MBTpathogensData_ord$Pathogens = as.factor(as.character(MBTpathogensData_ord$Pathogens))
MBTpathogensData_ord$Pathogens <- factor(MBTpathogensData_ord$Pathogens, 
                              levels = c("ORD" ,"ORD+Bacteria", "ORD+Virus","ORD+Mixed"))


names(MBTpathogensData_ord)

# Plot with reordered levels
mbt_Pathogen_within_ord_plot <- ggplot(MBTpathogensData_ord, aes(x = Pathogens, y = MBT_Score)) + 
  geom_boxplot(aes(color = Pathogens), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2", "#FC4E07")) +  
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(
    #comparisons = comparisons, 
                     method = "kruskal.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE    # Hide non-significant results
                     #label = "p.signif"
                     ) +  
  labs(
    #title = "TB Score Expression Across Pathogen Groups", 
    x = "Pathogen Group", y = "MBT Score")
print(mbt_Pathogen_within_ord_plot)
png("MBT_Score_withinORD_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
mbt_Pathogen_within_ord_plot
dev.off()



# within TB 

# merge with pathogen data by sample id 
MBTpathogensData.tb = merge(MBTpathogensData, HR_patho_within_tb, by="ParticipantID")
names(MBTpathogensData.tb)
dim(MBTpathogensData.tb)

MBTpathogensData.tb = MBTpathogensData.tb[, c(1,3,13, 14)]

# plot
path_mbt_tb_plot <- ggplot(MBTpathogensData.tb, aes(x = Pathogen_group, y = MBT_Score)) + 
  geom_boxplot(aes(color = Pathogen_group), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2")) +  
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(
    #comparisons = compp, 
                     method = "kruskal.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +  
  labs(x = NULL,
       y = "MBT Score")
print(path_mbt_tb_plot)

png("MBT_Score_withinTB_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
path_mbt_tb_plot
dev.off()


# Within TB append mixed to have mixed infection plot : manually

# create dataset with 4 mixed infection for TB
Dt
subset_dt <- Dt[Dt$ParticipantID %in% c("TRI0973", "TRI1477", "TRI0100", "TRI0578"), ]
names(subset_dt)
subset_dt_pred <- predict(model.glm, newdata = subset_dt, type = "prob")
subset_dt_pred = subset_dt_pred[, c("TB")]
# append the two below to the previously plotted TB data to include mixed infection
# TRI0100	0.4194828	TB+Mixed	TB
# TRI0578	0.7563323	TB+Mixed	TB

# load data 
MBTpathogensData.tb_plusmix <- read_csv("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/Cepheid3HR/Task3_3HR_data_cleaning/Redone with LDA total/MBTpathogensData_with_manuallyappended_2mixedTB.csv")

# plot data 
# order
MBTpathogensData.tb_plusmix$Pathogen_group = factor(MBTpathogensData.tb_plusmix$Pathogen_group, 
                                                    levels = c("TB+only","TB+Bacteria", "TB+Virus", "TB+Mixed"))
MBTpathogensData.tb_plusmix <- ggplot(MBTpathogensData.tb_plusmix, aes(x = Pathogen_group, y = MBT_Score)) + 
  geom_boxplot(aes(color = Pathogen_group), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2")) +  
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(
    #comparisons = compp, 
                     method = "kruskal.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +  
  labs(x = NULL,
       y = "MBT Score")
print(path_mbt_tb_plot)

png("MBT_Score_withinTB_including_mixed_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
MBTpathogensData.tb_plusmix
dev.off()

```


# compare with ROC curves 

```{r rocc}
# subset data for each binary comprison with TB 
names(Pathogen_final)
Pathogen_final$Pathogens #TB ORD ORD+Bacteria ORD+Virus ORD+Mixed

# Subsets 
ORD_rocc = subset(Pathogen_final, Pathogens %in% c("TB", "ORD"))
ORD_rocc$Pathogens = as.factor(as.character(ORD_rocc$Pathogens))

ORDbac_rocc = subset(Pathogen_final, Pathogens %in% c("TB", "ORD+Bacteria"))
ORDbac_rocc$Pathogens = as.factor(as.character(ORDbac_rocc$Pathogens))

ORDvir_rocc = subset(Pathogen_final, Pathogens %in% c("TB", "ORD+Virus"))
ORDvir_rocc$Pathogens = as.factor(as.character(ORDvir_rocc$Pathogens))

ORDmix_rocc = subset(Pathogen_final, Pathogens %in% c("TB", "ORD+Mixed"))
ORDmix_rocc$Pathogens = as.factor(as.character(ORDmix_rocc$Pathogens))

# Open a TIFF device
tiff("Pathogen_roc_plots_with_auc.tiff", width = 6.5, height = 6, units = "in", res = 300)

# Plot the first ROC curve
ORD_roc=plot.roc(ORD_rocc$Pathogens, ORD_rocc$TB_score, col = "#00AFBB")
# Calculate AUC and CI for the first ROC curve
roc_auc_ORD <- auc(ORD_roc)
ci_auc_ORD <- ci.auc(ORD_roc)

# Add AUC and CI to the plot (adjust x and y to avoid overlap)
text(0.4, 0.2, label = paste("TB vs ORD  AUC =", round(roc_auc_ORD, 2), "95% CI: ", round(ci_auc_ORD[1], 2), "-", round(ci_auc_ORD[3], 2)), 
     col = "#00AFBB", cex = 0.9)

# Add subsequent ROC curves and their AUCs with CI
ORDbac_roc=lines.roc(ORDbac_rocc$Pathogens, ORDbac_rocc$TB_score, col = "#E7B800")
roc_auc_ORD_bac <- auc(ORDbac_roc)
ci_auc_ORD_bac <- ci.auc(ORDbac_roc)
text(0.4, 0.15, label = paste("TB vs ORD+BAC AUC =", round(roc_auc_ORD_bac, 2), "95% CI: ", round(ci_auc_ORD_bac[1], 2), "-", round(ci_auc_ORD_bac[3], 2)), 
     col = "#E7B800", cex = 0.9)

ORDvir_roc=lines.roc(ORDvir_rocc$Pathogens, ORDvir_rocc$TB_score, col = "#0073C2")
roc_auc_ORD_vir <- auc(ORDvir_roc)
ci_auc_ORD_vir <- ci.auc(ORDvir_roc)
text(0.4, 0.1, label = paste("TB vs ORD+VIR AUC =", round(roc_auc_ORD_vir, 2), "95% CI: ", round(ci_auc_ORD_vir[1], 2), "-", round(ci_auc_ORD_vir[3], 2)), 
     col = "#0073C2", cex = 0.9)

ORDmix_roc=lines.roc(ORDmix_rocc$Pathogens, ORDmix_rocc$TB_score, col = "#FC4E07")
roc_auc_ORD_mix <- auc(ORDmix_roc)
ci_auc_ORD_mix <- ci.auc(ORDmix_roc)
text(0.4, 0.05, label = paste("TB vs ORD+MIX AUC =", round(roc_auc_ORD_mix, 2), "95% CI: ", round(ci_auc_ORD_mix[1], 2), "-", round(ci_auc_ORD_mix[3], 2)), 
     col = "#FC4E07", cex = 0.9)

# Close the TIFF device to save the plot
dev.off()
```


# compare AUCS
```{r delong}
# Perform DeLong test
delong_ORD_vs_ORDbacteria <- roc.test(ORD_roc, ORDbac_roc, method = "delong")
delong_ORD_vs_ORDvirus <- roc.test(ORD_roc, ORDvir_roc, method = "delong")
delong_ORD_vs_ORDmix <- roc.test(ORD_roc, ORDmix_roc, method = "delong")
delong_ORDbacteria_vs_ORDvirus <- roc.test(ORDbac_roc, ORDvir_roc, method = "delong")
delong_ORDbacteria_vs_ORDmix <- roc.test(ORDbac_roc, ORDmix_roc, method = "delong")
delong_ORDvirus_vs_ORDmix <- roc.test(ORDvir_roc, ORDmix_roc, method = "delong")

# Extract p-values
p_value_ORD_vs_ORDbacteria <- delong_ORD_vs_ORDbacteria$p.value
p_value_ORD_vs_ORDbacteria
p_value_ORD_vs_ORDvirus <- delong_ORD_vs_ORDvirus$p.value
p_value_ORD_vs_ORDvirus
p_value_ORD_vs_ORDmix <- delong_ORD_vs_ORDmix$p.value
p_value_ORD_vs_ORDmix
p_value_ORDbacteria_vs_ORDvirus <- delong_ORDbacteria_vs_ORDvirus$p.value
p_value_ORDbacteria_vs_ORDvirus
p_value_ORDbacteria_vs_ORDmix <- delong_ORDbacteria_vs_ORDmix$p.value
p_value_ORDbacteria_vs_ORDmix
p_value_ORDvirus_vs_ORDmix <- delong_ORDvirus_vs_ORDmix$p.value
p_value_ORDvirus_vs_ORDmix

# > p_value_ORD_vs_ORDbacteria
# [1] 0.6613769
# > p_value_ORD_vs_ORDvirus
# [1] 0.8977884
# > p_value_ORD_vs_ORDmix
# [1] 0.486438
# > p_value_ORDbacteria_vs_ORDvirus
# [1] 0.5994183
# > p_value_ORDbacteria_vs_ORDmix
# [1] 0.3271625
# > p_value_ORDvirus_vs_ORDmix
# [1] 0.5939587
```



II. MBT TEST

```{r MBT probability}

# Can we do the same for MBT ??/
# Option 1 is to work with the tests set, which would tally with the auc in the mbt model performce or ORD all. 
# Option 2 is to calculate the probaility for TB for all samples and do the same as done for TB score 

# Option 1: 
# resplit data after editting data to include participant IDs

Dt = Baseline_all_Data


table(Baseline_all_Data$Composite)

# remove unknown TB
# MBT_defn = subset(Baseline_all_Data, Composite %in% c("Definite TB", "No TB", "Probable TB"))
              # lets forget about this for now: get back to it later after finish with Dt

dim(Dt)
names(Dt)
names(Dt)[16] = "Group"  # Binary renamed as group

# SUBSET REQUIRED COLUMNS  
table(Dt$Group)
names(Dt)

#subset only the quant data 
Dt = Dt[, c("ParticipantID", "Group", "SAA1","CRP", "IP.10")]

# FORMAT DATA 
# Create model with caret 
table(Dt$Group)

Dt$Group = as.factor(as.character(Dt$Group))
str(Dt)
summary(Dt)
# Reorder levels of Group
Dt$Group <- factor(Dt$Group, levels = c("TB", "ORD"))

# Check for missing values
sum(is.na(Dt))
sum(is.na(Dt$Group))


# ------------------------------------------------------------------#
# 1. Split data ----
set.seed(123)  # Ensure reproducibility
train_index <- createDataPartition(Dt$Group, p = 0.7, list = FALSE)
train_data <- Dt[train_index, ]
test_data <- Dt[-train_index, ]
dim(train_data)
dim(test_data)

# > dim(train_data)
# [1] 211   5
# > dim(test_data)
# [1] 89  5

# Now that we have athe participant IDs lets merger them



# resolve 
# glm_prob <- predict(glm_model, newdata = test_data, type = "prob")
glm_probabily = predict(glm_model, newdata = test_data, type = "prob")
Mbt_test_df = as.data.frame(glm_probabily)
names(Mbt_test_df)
rownames(Mbt_test_df)
dim(Mbt_test_df) # [1] 89  5
dim(test_data) # [1] 89  2
names(test_data)
rownames(test_data)

# merge data by rownames
merged_mbt <- merge(Mbt_test_df, test_data,by = "row.names", all = TRUE)
rownames(merged_mbt) <- merged_mbt$Row.names
merged_mbt <- merged_mbt[, -1]
names(merged_mbt)
names(merged_mbt)[2] = "MBT_Score"

# Merge with pathogen data: Pathogen_final: note this is only test set and includes other TB
names(Pathogen_final)
dim(Pathogen_final)

# Merge step: when merged it will exclude other TB so n could go down
MBTpathogensData = merge(merged_mbt, Pathogen_final, by="ParticipantID")
dim(MBTpathogensData)
names(MBTpathogensData)

```

```{r comparisons}
# Define comparisons
comparisons <- list(
  # c("ORD", "ORD+Bacteria"),
  # c("ORD", "ORD+Mixed"),
  # c("ORD", "ORD+Virus"),
  c("ORD", "TB"),
  # c("ORD+Bacteria", "ORD+Mixed"),
  # c("ORD+Bacteria", "ORD+Virus"),
  c("ORD+Bacteria", "TB"),
 # c("ORD+Mixed", "ORD+Virus"),
  c("ORD+Mixed", "TB"),
  c("ORD+Virus", "TB"))

# Plot with reordered levels
MBT_patho_plot <- ggplot(MBTpathogensData, aes(x = Pathogens, y = MBT_Score)) + 
  geom_boxplot(aes(color = Pathogens), width = 0.6, size = 1.2) +  
  scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2", "#FC4E07")) +  
  scale_y_continuous(limits = c(0, 1.6)) +
  #use codes to show all comaprisons
    #scale_color_manual(values = c("#6A3D9A", "#00AFBB", "#E7B800", "#0073C2", "#FC4E07")) +  
    # scale_y_continuous(limits = c(0, 2.5)) +
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8)) +  
  stat_compare_means(comparisons = comparisons, 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = TRUE    # Hide non-significant results
                     #label = "p.signif"
                     ) +  
  labs(
    #title = "MBT Score Expression Across Pathogen Groups", 
    x = "Pathogen Group", y = "MBT Score")
print(MBT_patho_plot)

png("MBT_Score_Across_Pathogen_Groups.png",
    width = 8, height = 6, units = "in", res=300)
MBT_patho_plot
dev.off()
```


# compare with ROC curves 
```{r rocc}
# subset data for each binary comprison with TB 
names(MBTpathogensData)
Pathogen_final$Pathogens #TB ORD ORD+Bacteria ORD+Virus ORD+Mixed

# Subsets 
MBT_rocc = subset(MBTpathogensData, Pathogens %in% c("TB", "ORD"))
MBT_rocc$Pathogens = as.factor(as.character(MBT_rocc$Pathogens))

MBTbac_rocc = subset(MBTpathogensData, Pathogens %in% c("TB", "ORD+Bacteria"))
MBTbac_rocc$Pathogens = as.factor(as.character(MBTbac_rocc$Pathogens))

MBTvir_rocc = subset(MBTpathogensData, Pathogens %in% c("TB", "ORD+Virus"))
MBTvir_rocc$Pathogens = as.factor(as.character(MBTvir_rocc$Pathogens))

MBTmix_rocc = subset(MBTpathogensData, Pathogens %in% c("TB", "ORD+Mixed"))
MBTmix_rocc$Pathogens = as.factor(as.character(MBTmix_rocc$Pathogens))

# Open a TIFF device
tiff("MBT_Pathogen_roc_plots_with_auc.tiff", width = 6.5, height = 6, units = "in", res = 300)

# Plot the first ROC curve
MBT_roc=plot.roc(MBT_rocc$Pathogens, MBT_rocc$MBT_Score, col = "#00AFBB")
# Calculate AUC and CI for the first ROC curve
roc_auc_MBT <- auc(MBT_roc)
ci_auc_MBT <- ci.auc(MBT_roc)

# Add AUC and CI to the plot (adjust x and y to avoid overlap)
text(0.4, 0.2, label = paste("TB vs ORD  AUC =", round(roc_auc_ORD, 2), "95% CI: ", round(ci_auc_ORD[1], 2), "-", round(ci_auc_ORD[3], 2)), 
     col = "#00AFBB", cex = 0.9)

# Add subsequent ROC curves and their AUCs with CI 
MBTbac_roc=lines.roc(MBTbac_rocc$Pathogens, MBTbac_rocc$MBT_Score, col = "#E7B800")
roc_auc_MBT_bac <- auc(MBTbac_roc)
ci_auc_MBT_bac <- ci.auc(MBTbac_roc)
text(0.4, 0.15, label = paste("TB vs ORD+BAC AUC =", round(roc_auc_ORD_bac, 2), "95% CI: ", round(ci_auc_ORD_bac[1], 2), "-", round(ci_auc_ORD_bac[3], 2)), 
     col = "#E7B800", cex = 0.9)

MBTvir_roc=lines.roc(MBTvir_rocc$Pathogens, MBTvir_rocc$MBT_Score, col = "#0073C2")
roc_auc_MBT_vir <- auc(MBTvir_roc)
ci_auc_MBT_vir <- ci.auc(MBTvir_roc)
text(0.4, 0.1, label = paste("TB vs ORD+VIR AUC =", round(roc_auc_ORD_vir, 2), "95% CI: ", round(ci_auc_ORD_vir[1], 2), "-", round(ci_auc_ORD_vir[3], 2)), 
     col = "#0073C2", cex = 0.9)

MBTmix_roc=lines.roc(MBTmix_rocc$Pathogens, MBTmix_rocc$MBT_Score, col = "#FC4E07")
roc_auc_MBT_mix <- auc(MBTmix_roc)
ci_auc_MBT_mix <- ci.auc(MBTmix_roc)
text(0.4, 0.05, label = paste("TB vs ORD+MIX AUC =", round(roc_auc_ORD_mix, 2), "95% CI: ", round(ci_auc_ORD_mix[1], 2), "-", round(ci_auc_ORD_mix[3], 2)), 
     col = "#FC4E07", cex = 0.9)

# Close the TIFF device to save the plot
dev.off()
```

# Compare AUCs

# compare AUCS
```{r delong}
# Perform DeLong test
# MBT_roc
# MBTbac_roc
# MBTvir_roc
# MBTmix_roc

m.delong_ORD_vs_ORDbacteria <- roc.test(MBT_roc, MBTbac_roc, method = "delong")
m.delong_ORD_vs_ORDvirus <- roc.test(MBT_roc, MBTvir_roc, method = "delong")
m.delong_ORD_vs_ORDmix <- roc.test(MBT_roc, MBTmix_roc, method = "delong")
m.delong_ORDbacteria_vs_ORDvirus <- roc.test(MBTbac_roc, MBTvir_roc, method = "delong")
m.delong_ORDbacteria_vs_ORDmix <- roc.test(MBTbac_roc, MBTmix_roc, method = "delong")
m.delong_ORDvirus_vs_ORDmix <- roc.test(MBTvir_roc, MBTmix_roc, method = "delong")

# Extract p-values
mp_value_ORD_vs_ORDbacteria <- m.delong_ORD_vs_ORDbacteria$p.value
mp_value_ORD_vs_ORDbacteria
mp_value_ORD_vs_ORDvirus <- m.delong_ORD_vs_ORDvirus$p.value
mp_value_ORD_vs_ORDvirus
mp_value_ORD_vs_ORDmix <- m.delong_ORD_vs_ORDmix$p.value
mp_value_ORD_vs_ORDmix
mp_value_ORDbacteria_vs_ORDvirus <- m.delong_ORDbacteria_vs_ORDvirus$p.value
mp_value_ORDbacteria_vs_ORDvirus
mp_value_ORDbacteria_vs_ORDmix <- m.delong_ORDbacteria_vs_ORDmix$p.value
mp_value_ORDbacteria_vs_ORDmix
mp_value_ORDvirus_vs_ORDmix <- m.delong_ORDvirus_vs_ORDmix$p.value
mp_value_ORDvirus_vs_ORDmix

# resulting p values 
# > # Extract p-values
# > mp_value_ORD_vs_ORDbacteria
# [1] 0.8525191
# > mp_value_ORD_vs_ORDvirus
# [1] 0.7946423
# > mp_value_ORD_vs_ORDmix
# [1] 0.7501673
# > mp_value_ORDbacteria_vs_ORDvirus
# [1] 0.698164
# > mp_value_ORDbacteria_vs_ORDmix
# [1] 0.8828729
# > mp_value_ORDvirus_vs_ORDmix
# [1] 0.6335158
```



# II. MBT and tretament 
#-------------------------------#
```{r data setup}
# Propagate TB status to the data 
RxData = All_timepoint
# RxData <- RxData %>%
#   group_by(ParticipantID) %>%
#   mutate(Binary = dplyr::first(Binary[Timepoint == "Baseline"], default = NA)) %>%
#   ungroup() 
# dim(RxData)
# table(RxData$Timepoint)

# Remove NA in protein Data 
RxData = RxData[!is.na(RxData$SAA1), ]
RxData <- RxData[!is.na(RxData$Binary), ]

# # subset only TB 
# RxDataTb = subset(RxData, Binary %in% c("TB"))
# names(RxDataTb)
# table(RxDataTb$Timepoint) # Perfect 


# # Baseline  Month 4  Month 6 
# 115       59       54 

# Write RxData to csv amend propagate TB staus and load data 
write.csv(RxData, "RxData.csv")
RxData_amended= read_csv("RxData_amended.csv")
# Method 1: Changes in probability 

# Method 2: model diff baseline ord and m6 TB se and sp
# remove unknown TB: We manually added M4 and M4 to labels in composite
MBTrxcombine = subset(RxData_amended, Composite %in% c("Definite TB","No TB", "Probable TB", "M4", "M6"))
dim(MBTrxcombine)
table(MBTrxcombine$Composite)
# Definite TB       No TB Probable TB 
# 103         185           2 
names(MBTrxcombine)
# make timepoint subsets 
table(MBTrxcombine$Timepoint)

# Baseline  ----
MbtBL = subset(MBTrxcombine, Timepoint %in% c("Baseline"))
# subset baseline ord to append to other timepoint for comparison with ord 
MbtBL_ORD = subset(MbtBL, Binary %in% c("ORD"))

# Month 4 -----
MbtM4 = subset(MBTrxcombine, Timepoint %in% c("Month 4"))
# append ord 
MbtM4_all = rbind(MbtM4, MbtBL_ORD)
dim(MbtM4_all) # [1] 244  16

# Month 6 ----
MbtM6 = subset(MBTrxcombine, Timepoint %in% c("Month 6"))
MbtM6_all = rbind(MbtM6, MbtBL_ORD)
dim(MbtM6_all) # [1] 239  16
# 
```


 #-----------------------------------------------------#
```{r MBTR se and sp for bl vs m4 or m6}
# MBT Model predcition of M6 vs BL data 

names(MbtM6_all)
DataMsix = MbtM6_all[, c(1, 3:5, 16)]
names(DataMsix)[5] = "Group"
DataMsix$Group = as.factor(as.character(DataMsix$Group))

# Predicitions 
M6Pred <- predict(glm_model, newdata = DataMsix) # prob
# Predict probabilities for the "TB" class
M6Prob <- predict(glm_model, newdata = DataMsix, type = "prob")
names(M6Prob)
M6Probdata = M6Prob[, c("TB")]
M6_cm = confusionMatrix(M6Pred, DataMsix$Group) #

# ROC 
M6_roc <- roc(DataMsix$Group, M6Probdata)
# Call:
#   roc.default(response = DataMsix$Group, predictor = M6Probdata)
# 
# Data: M6Probdata in 185 controls (DataMsix$Group ORD) < 54 cases (DataMsix$Group TB).
# Area under the curve: 0.5064


#-----------------------------------------------------#
# The same for M4 

# MBT Model predcition of M6 vs BL data 
names(MbtM4_all)
DataMf = MbtM4_all[, c(1, 3:5, 16)]
names(DataMf)[5] = "Group"
DataMf$Group = as.factor(as.character(DataMf$Group))

# Predicitions 
M4Pred <- predict(glm_model, newdata = DataMf) # prob
# Predict probabilities for the "TB" class
M4Prob <- predict(glm_model, newdata = DataMf, type = "prob")
names(M4Prob)
M4Probdata = M4Prob[, c("TB")]
M4_cm = confusionMatrix(M4Pred, DataMf$Group) #

# ROC 
M4_roc <- roc(DataMf$Group, M4Probdata)
M4_roc
# Call:
#   roc.default(response = DataMf$Group, predictor = M4Probdata)
# 
# Data: M4Probdata in 185 controls (DataMf$Group ORD) < 59 cases (DataMf$Group TB).
# Area under the curve: 0.5986


#-----------------------------------------------------#
# Lets do this for BL too
MbtBL
# MBT Model predcition of M6 vs BL data 
names(MbtBL)
DataBL = MbtBL[, c(1, 3:5, 16)]
names(DataBL)[5] = "Group"
DataBL$Group = as.factor(as.character(DataBL$Group))

# Predicitions 
BLPred <- predict(glm_model, newdata = DataBL) # prob
# Predict probabilities for the "TB" class
BLProb <- predict(glm_model, newdata = DataBL, type = "prob")
names(BLProb)
BLProbdata = BLProb[, c("TB")]
BL_cm = confusionMatrix(BLPred, DataBL$Group) #

# ROC 
BL_roc <- roc(DataBL$Group, BLProbdata)
BL_roc
```


# compare Probability scores for TB.
```{r compare prob}
# Baseline-----
# add prob colun to the data, filter out TB 
BLProbdata = as.data.frame(BLProbdata)
names(BLProbdata)[1] = "MBT_Prob.score"
# DataBL
BL_datasetm = merge(DataBL,BLProbdata, by = "row.names") %>%
  select(-Row.names) 
names(BL_datasetm)
dim(BL_datasetm)

BL_data.df = subset(BL_datasetm, Group %in% c("TB"))
BL_data.df$Group = as.factor(as.character(BL_data.df$Group))
dim(BL_data.df)
# Subset Group and prob then replace TB with the timepoint 
BL_data.df = BL_data.df[, c(5:6)]
BL_data.df$Group = gsub("TB", "Baseline", BL_data.df$Group)



# Month 4-----
# M4Probdata DataMf
# add prob colun to the data, filter out TB 
M4Probdata = as.data.frame(M4Probdata)
names(M4Probdata)[1] = "MBT_Prob.score"
# DataBL
DataM4dataset = merge(DataMf,M4Probdata, by = "row.names") %>%
  select(-Row.names) 
names(DataM4dataset)
dim(DataM4dataset)

M4_data.df = subset(DataM4dataset, Group %in% c("TB"))
M4_data.df$Group = as.factor(as.character(M4_data.df$Group))
dim(M4_data.df)
# Subset Group and prob then replace TB with the timepoint 
M4_data.df = M4_data.df[, c(5:6)]
M4_data.df$Group = gsub("TB", "Month_4", M4_data.df$Group)




# Month 6-----
# DataMsix M6Probdata
# add prob colun to the data, filter out TB 
M6Probdata = as.data.frame(M6Probdata)
names(M6Probdata)[1] = "MBT_Prob.score"
# DataBL
DataM6dataset = merge(DataMsix,M6Probdata, by = "row.names") %>%
  select(-Row.names) 
names(DataM6dataset)
dim(DataM6dataset)

M6_data.df = subset(DataM6dataset, Group %in% c("TB"))
M6_data.df$Group = as.factor(as.character(M6_data.df$Group))
dim(M6_data.df)
# Subset Group and prob then replace TB with the timepoint 
M6_data.df = M6_data.df[, c(5:6)]
M6_data.df$Group = gsub("TB", "Month_6", M6_data.df$Group)



# Comapre mean probailities 
# rbind all three 
# Combine the datasets
combined_dff <- rbind(BL_data.df, M4_data.df, M6_data.df)
names(combined_dff)
combined_dff$Group = as.factor(combined_dff$Group)

# compare Groups (timepoint)
# Perform Kruskal-Wallis test
kruskal_mbt_tx <- kruskal.test(MBT_Prob.score ~ Group, data = combined_dff)
print(kruskal_mbt_tx)

summary_mbt <- combined_dff %>%
  group_by(Group) %>%
  summarise(
    Mean_MBT_Prob_Score = mean(MBT_Prob.score, na.rm = TRUE),
    Median_MBT_Prob_Score = median(MBT_Prob.score, na.rm = TRUE),
    Count = n())
print(summary_mbt)


# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Define comparisons
comparizon <- list(c("Baseline", "Month_4"), 
                    c("Baseline", "Month_6"), 
                    c("Month_4", "Month_6"))

# Create the box plot with p-value
#tiff("EoTx_MBT_Score_TB_vs_ORD_BinaryTB.tiff", 
  #   width = 6, height = 6, res=300, unit="in") 


# Plot
pm <- ggplot(combined_dff, aes(x = Group, y = MBT_Prob.score, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, size = 1.2) +  # Boxplot
  geom_jitter(width = 0.15, size = 4, alpha = 0.3) +           # Add jittered points
  stat_compare_means(comparison = comparizon, method = "wilcox.test") +  # Pairwise p-values
  stat_compare_means(label.y = 0.01) +                            # Global p-value on top
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +   # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  labs(title = "MBT Probability Score Across Timepoints",
       y = "MBT_ProbScore", x = "Group")

# Display the plot
print(pm)
#dev.off()
```


# Compare to baseline ORD 
```{r compare eot and ord}
# compare EOT and M6 -------
ORDmbt_baseline = BLProb[, c("TB")]
ORDmbt_baseline = as.data.frame(ORDmbt_baseline)
names(ORDmbt_baseline)

names(ORDmbt_baseline)[1] = "MBT_Prob.score"
# merge with BL data 
OD_datasetm = merge(DataBL,ORDmbt_baseline, by = "row.names") %>%
  select(-Row.names) 
names(OD_datasetm)
dim(OD_datasetm)

OD_data.df = subset(OD_datasetm, Group %in% c("ORD"))
OD_data.df$Group = as.factor(as.character(OD_data.df$Group))
dim(OD_data.df)
names(OD_data.df)
# Subset Group and prob then replace TB with the timepoint 
OD_data.df = OD_data.df[, c(5:6)]

# Mergee with M6 M6_data.df using rbind
EOTm_dff <- rbind(OD_data.df, M6_data.df)
names(EOTm_dff)
EOTm_dff$Group = as.factor(EOTm_dff$Group)

# Plot 
# Create the box plot with p-value
tiff("EoTx_vs_BL_ORD.tiff", 
     width = 6, height = 6, res=300, unit="in") 


# Plot
pt <- ggplot(EOTm_dff, aes(x = Group, y = MBT_Prob.score, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, size = 1.2) + 
  geom_jitter(width = 0.15, size = 4, alpha = 0.3) + 
  stat_compare_means(method = "wilcox.test") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +   # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 0.8)) +
  labs(title = "MBT Probability at EOT vs Baseline ORD",
       y = "MBT Score", x = "Group")

# Display the plot
print(pt)

dev.off()
```

```{r m6 roc}
# Load the pROC package
library(pROC)
M6_roc <- roc(DataMsix$Group, M6Probdata)
M6_ci <- ci.auc(M6_roc)

# Plot the ROC curve
# Create the box plot with p-value
tiff("EoTx_vs_BL_ROCC.tiff", 
     width = 6, height = 6, res=300, unit="in") 

plot(M6_roc, col = "#2C7BB6", lwd = 2, main = "ROC Curve with AUC and CI")

# Add the AUC with CI as text on the plot
auc_text <- paste0("AUC = ", round(auc(M6_roc), 3),
                   " (95% CI: ", round(M6_ci[1], 3), "-",
                   round(M6_ci[3], 3), ")")

legend("bottomright", legend = auc_text, col = "#2C7BB6", lwd = 2, bty = "n")

dev.off()
```



```{r save}
# Last saved: 05-Feb-2025
# save.image("POC_Analyses_ALL_Objects_01.02.25.RData")
# At this point all analyses were complete, pending review of codes #
```


###### Luminex
```{r luminex ip10 comparsions }

# Load or rerun script 
load("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Triage POC paper/MBT/Task3_MBT_data_cleaning/MBTDataceaning_step1.RData")
names(Baseline_all_Data)
dim(Baseline_all_Data)

# Subset IP 10 data 
MBT_IP = Baseline_all_Data[, c(1,5,16)]
dim(MBT_IP)

# Save data 
write.csv(MBT_IP, "MBT_IP.csv") # for luminex analyses 

```
