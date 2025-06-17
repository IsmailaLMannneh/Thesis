# Chapter 4: 
# PART 1 (ROSALIND)----
# data Normalizarion and DEG analyses done in ROSALIND cloud interface 
#---------------------------------------------------#


# PART 2: DEMOGRAPHIC CHARACTERISTICS ------
#---------------------------------------------------#
### DEMOGRAPHIC Chracteristics comparisons 
library(readr)
library(readxl)

#Read in the files 
setwd("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Nanostring paper/Demography")
SampList <- read_excel("Baseline_samplesList.xlsx")
names(SampList)
DemoTab <- read_excel("Table_D.xlsx") #
names(DemoTab)
DemoTab <- DemoTab[, c(1, 4, 5, 8:11)]
names(DemoTab)[1] = "PaticipantID"

# Merge with with sample IDs 
# find common column 
commonID <- intersect(names(DemoTab), names(SampList))
Merged_Tab <- merge(DemoTab, SampList, by=commonID) 
names(Merged_Tab)
dim(Merged_Tab)

#format data 
Merged_Tab$Age <- round(Merged_Tab$Age, digits = 0)
Merged_Tab$Sex <- as.factor(Merged_Tab$Sex)
Merged_Tab$Response <- as.factor(Merged_Tab$Response)
Merged_Tab$Outcome <- as.factor(Merged_Tab$Outcome)

summary(Merged_Tab$Age)
summary(Merged_Tab$Sex)
summary(Merged_Tab$Response) 
summary(Merged_Tab$Outcome)


# STATISTICS 

#Chi Square and Fisher exact test  

# sex
table(Merged_Tab$Sex, Merged_Tab$Outcome)
table(Merged_Tab$Sex, Merged_Tab$Response)
prop.table(table(Merged_Tab$Sex, Merged_Tab$Outcome), 1) # as percentage 
chisq.test(Merged_Tab$Sex, Merged_Tab$Response, correct=FALSE)

# Make dataframe 
your_data <- data.frame(
  group = c("Cured", "Failure", "Relapse"),
  proportion_male = c(0.728, 0.714, 0.66))

# Chi-squared test for independence
chi_square_result <- chisq.test(matrix(c(your_data$proportion_male, 1 - your_data$proportion_male), ncol = 2))

# Age 
library(dplyr)
summary_stats <- Merged_Tab %>%
  group_by(Outcome) %>%
  summarise(
    median_age = median(Age),
    lower_quartile = quantile(Age, 0.25),
    upper_quartile = quantile(Age, 0.75))

print(summary_stats)

# Kruskal-Wallis test 
kruskal.test(Age ~ Outcome, data = Merged_Tab)
boxplot(Age ~ Outcome, data = Merged_Tab, main = "Age distribution by Group", xlab = "Group", ylab = "Age")


# HIV
Hiv <- subset(Merged_Tab, Merged_Tab$HIV %in% c('Pos'))
table(Hiv$Outcome)

# TB history
history <- subset(Merged_Tab, Merged_Tab$TB_History %in% c('Yes'))
table(history$Outcome)

#### Part 2: Signatures and feature selections





# PART 3: GENE SIGNATURES -------

# A. TREATMENT RESPONSE -----

# # read in data -----
dataa <- readRDS("M2_data.RDS") # M2_Response_data_allfeats.RDS
a <- dataa[, -ncol(dataa)]
a <- a %>% mutate_all(as.numeric)

# feature selection by filteration
a_pro <- preProcess(a,method = c("corr", "nzv"),cutoff = 0.8)
a_pro1 <- predict(a_pro,a)
dim(a_pro1)
y <- as.factor(dataa$Response) # still index matched


datacomb <- cbind(y, a_pro1)
dim(datacomb)

# Split data
set.seed(123) # Set seed for reproducibility
trainIn <- createDataPartition(datacomb$y, p = .8,
                               list = FALSE,
                               times = 1)
training <- datacomb[trainIn,]
testing  <- datacomb[-trainIn,]

table(training$y)
table(testing$y)

# Set up your training control with upsampling the smaller category
set.seed(123)
train_control <- trainControl(method = "cv",  # Cross-validation method
                              number = 5)

model <- train(y~.,  # Specify your formula
               method = "rf",
               data=training, # Specify your model
               trControl = train_control,ntrees=500) # Pass the training control
summary(model)
model$results
model$finalModel


featz = varImp(model)

# check trianing
predt <- predict(model, newdata = training)
conf_predt <- confusionMatrix(reference=training$y, data = predt,positive = "Slow")
print(conf_predt)
accuracy <- conf_predt$overall["Accuracy"]
print(accuracy)

# predict hold out
predict_test <- predict(model, newdata = testing)
conf_test <- confusionMatrix(reference = testing$y, data = predict_test,positive = "Slow")
print(conf_test)
accuracy <- conf_test$overall["Accuracy"] #
print(accuracy) # 0.5455



# load in Data
dataa <- readRDS("M2_Data.RDS") # M2_Response_data_allfeats.RDS
dim(dataa)
names(dataa)

# remove response variable (last col)
a <- dataa[, -ncol(dataa)]

# make all predictors numeric
a <- a %>% mutate_all(as.numeric)

# feature selection by filteration
a_pro <- preProcess(a,method = c("corr", "nzv"),cutoff = 0.8)

# filtered data
a_pro1 <- predict(a_pro,a)
dim(a_pro1)
# [1]  59 291

y <- as.factor(dataa$Response) # still index matched
datacomb <- cbind(y, a_pro1) # honestly no need to combine them into a df
names(datacomb)[-1] # name of predictors

# setseed
set.seed(123)

# Initialize variables outside the loop
auc_values <- numeric(500)
roc_list <- list()
importance <- matrix(0, nrow = 500, ncol = ncol(a_pro1))  # Assuming ncol(training) is the number of features
dim(importance)
accuracy <- rep(NA, 500)

# Loop for 500 iterations
for (i in 1:500) {
  # Split data into training and testing sets
  trainIn <- createDataPartition(datacomb$y, p = .8, list = FALSE, times = 1)
  training <- datacomb[trainIn,]
  testing  <- datacomb[-trainIn,]
  
  #dim(training)
  #importance <- matrix(0, nrow = 500, ncol = ncol(training))
  
  # Train the model
  train_control <- trainControl(method = "cv", number = 5)
  model <- train(y ~ ., data = training, method = "rf", trControl = train_control, ntrees = 500)
  
  # Predict probabilities on testing data
  pred_prob <- predict(model, newdata = testing, type = "prob")
  pred_y <- predict(model, newdata = testing)
  
  # Collect feature importance if 'features' is defined
  feature_import <- varImp(model)
  features <- rownames(feature_import$importance)[which(feature_import$importance>=0.50)]
  importance[i,] <- (names(a_pro1)%in%features)*1
  
  # Calculate accuracy
  cm <- confusionMatrix(pred_y, testing$y, positive = "Slow")
  accuracy[i] <- cm$overall["Accuracy"]
  
  # Calculate ROC curve
  roc_data <- roc(testing$y, pred_prob[, "Slow"])
  
  # Store AUC value
  auc_values[i] <- auc(roc_data)
  
  # Store ROC curve data
  roc_list[[i]] <- roc_data
}


# Find the index of models with accuracy > 0.70
Model_index <- which(accuracy >= 0.00)

# Get the features in the top n 29 models given by # which(accuracy >= 0.70)

# Total number of top models
total_top_models <- length(Model_index)

# Count the occurrence of each feature in the top models
feature_counts <- colSums(importance[Model_index, ])

# Calculate the percentage of occurrence for each feature
feature_percentages <- (feature_counts / total_top_models) * 100

# Get the features with 100% occurrence.
selected_features <- which(feature_percentages == 100.00000) # lets start with 70
length(selected_features) # n =157

# Plot
# Plot the occurrences of each feature
barplot(feature_percentages,
        main = "Occurrences of Features in Models with Accuracy > 0.7",
        xlab = "Features",
        ylab = "Occurrences",
        col = "skyblue",
        las = 2)  # Rotate x-axis labels vertically if necessary


# Subset the Data for these features using the feature matrix (a_pro1)
subset_data_feats <- a_pro1[, selected_features] # n=157, not y is still y and index remain unchanged
names(subset_data_feats)
dim(subset_data_feats)



# Built refined model

library(caret)
subset_data_set = cbind(y, subset_data_feats)

# Step 1: Split data into training and testing sets
set.seed(123)  # Set seed for reproducibility
trainIndexes <- createDataPartition(subset_data_set$y, p = 0.8, list = FALSE, times = 1)
trainingdata <- subset_data_set[trainIndexes, ]
testingdata <- subset_data_set[-trainIndexes, ]

# Step 2: Train a machine learning model
model2 <- train(y ~ ., data = trainingdata, method = "rf", ntree=500)

# Step 3: Evaluate the model's performance
prediction1 <- predict(model2, newdata = testingdata, type = "prob")
prediction1b <- predict(model2, newdata = testingdata)

# Calculate accuracy
conf_matrix = confusionMatrix(prediction1b, testingdata$y, positive = "Slow")
print(conf_matrix)
Sub_accuracy <- conf_matrix$overall["Accuracy"]
print(Sub_accuracy)


# Calculate ROC curve
Sub_roc_data <- roc(testingdata$y, prediction1[, "Slow"])
# Roc AUC value
Subs_AUC <- auc(Sub_roc_data)
print(Subs_AUC)

# Print results for test set ------
# > print(Sub_accuracy)
# Accuracy
# 0.8181818
# > print(Subs_AUC)
# Area under the curve: 0.9333


# Save R Data -------
# #save.image("Month2_modelling.RData")



# Predict baseline and week 2 using this model ----

# 1. Predict Baseline ------

# load data
dd = readRDS("Month2_modelling.RData")

# remove response variable (last col)
d <- dd[, -ncol(dd)] # all predictors

# make all predictors numeric
d <- d %>% mutate_all(as.numeric) # make numeric 
names(d)

#response 
y <- as.factor(dd$Response) # make and format y

# filtered data 
d1 <- predict(a_pro,d)
dim(d1)
names(d1)
# [1]  59 291


feats <- d1[, selected_features] # n=157, note y index  unchanged
names(feats)

#now combine 
combine_M2 = cbind(y, feats)

summary(combine_M2)

# rename 
validationdata = combine_M2
names(validationdata)

# Step 3: Evaluate the model's performance
Predictprob <- predict(model2, newdata = validationdata, type = "prob")
predictclass <- predict(model2, newdata = validationdata) 

# Calculate accuracy
cfm = confusionMatrix(predictclass, validationdata$y, positive = "Slow")
print(cfm) # AUC
M2_accuracy <- cfm$overall["Accuracy"]
print(M2_accuracy)

# Calculate ROC curve
rocval <- roc(validationdata$y, Predictprob[, "Slow"])
# Roc AUC value
M2_AUC <- auc(rocval) 
print(M2_AUC) #Area under the curve: 0.4851

# saveRDS(rocval, "BLData.RDS")
# saveRDS(M2_AUC, "BL.AUC.RDS")



# Week 2 Predictions -----
week2data = readRDS("W2_Response_data_allfeats.RDS")
dd = week2data

# remove response variable (last col)
d <- dd[, -ncol(dd)] # all predictors

# make all predictors numeric
d <- d %>% mutate_all(as.numeric) # make numeric 
names(d)

#response 
y <- as.factor(dd$Response) # make and format y

# filtered data 
d1 <- predict(a_pro,d)
dim(d1)
names(d1)
# [1]  59 291


feats <- d1[, selected_features] # n=157, not y is still y and index remain unchanged
names(feats)

# combine 
combine_M2 = cbind(y, feats)

summary(combine_M2)

# rename 
validationdata = combine_M2
names(validationdata)

# Step 3: Evaluate the model's performance
Predictprob <- predict(model2, newdata = validationdata, type = "prob")
predictclass <- predict(model2, newdata = validationdata) 

# Calculate accuracy
cfm = confusionMatrix(predictclass, validationdata$y, positive = "Slow")
print(cfm) # AUC
M2_accuracy <- cfm$overall["Accuracy"]
print(M2_accuracy)

# Calculate ROC curve
rocval <- roc(validationdata$y, Predictprob[, "Slow"])
# Roc AUC value
M2_AUC <- auc(rocval) 
print(M2_AUC) #Area under the curve: 0.5708

# saveRDS(rocval, "W2Data.RDS")
# saveRDS(M2_AUC, "W2AUC.RDS")



#---------------------------------------------------------------------#
# Part 2: ROC plots 
#---------------------------------------------------------------------#

# Plot ROC curves -----

# Load rocs  
BaselineROCdata = readRDS("BLData.RDS") 
Week2ROCdata = readRDS("W2Data.RDS") 
Month2ROCdata = readRDS("M2Data.RDS") 

# Load roc aucs 
BaseAUC = readRDS("BL.AUC.RDS")
Wk2AUC = readRDS("W2AUC.RDS")
Mt2AUC = readRDS("M2AUC.RDS")

# Plot 
library(pROC)

BaselineROCdata_smooth <- smooth(BaselineROCdata, method = "density")
Week2ROCdata_smooth <- smooth(Week2ROCdata, method = "density")
Month2ROCdata_smooth <- smooth(Month2ROCdata, method = "density")

# Plot the smoothed ROC curves
plot(BaselineROCdata, type = "l", col = "blue", lty = 1, lwd = 2, main = "ROC Curve")
lines(Week2ROCdata, col = "green", lty = 1, lwd = 2)
lines(Month2ROCdata, col = "black", lty = 1, lwd = 2)


# Load roc aucs 
BaseAUC #= readRDS("BL.AUC.RDS")
Wk2AUC #= readRDS("W2AUC.RDS")
Mt2AUC #= readRDS("M2AUC.RDS")

legend("bottomright", legend = c(paste("Baseline AUC =", round(BaseAUC, 3)),
                                 paste("Week 2 AUC =", round(Wk2AUC, 3)),
                                 paste("Month 2 AUC =", round(Mt2AUC, 3))),
       col = c("blue", "green", "black"), lty = 1, cex = 1.3)


# Adjust plot 
#---------------------------------------------------

# Color pallet 
#-----------------
# Baseline = Black
# week 2 = brown
# Month2 = Blue 
# M6 = Green

# Plot the smoothed ROC curves

# initial tiff object 
tiff("TxResponseSlowFast_roc_plot.tiff", 
     width = 6, height = 6, units = "in", res = 300)

plot(BaselineROCdata, type = "l", col = "black", lty = 1, lwd = 2, #main = "ROC Curve",
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(Week2ROCdata, col = "brown", lty = 1, lwd = 2)
lines(Month2ROCdata, col = "blue", lty = 1, lwd = 2)
legend("bottomright", legend = c(paste("Baseline AUC =", round(BaseAUC, 3)),
                                 paste("Week 2 AUC =", round(Wk2AUC, 3)),
                                 paste("Month 2 AUC =", round(Mt2AUC, 3))),
       col = c("black", "brown", "blue"), lty = 1, cex = 1.0)

dev.off()

# save.image("Final_TxResponseSlowFast.RData")

# B. TREATMENT OUTCOME -----

# load libraries
library(dplyr)
library(readxl)
library(caret)
library(readr)
library(tibble)
library(pROC)

# Load data
M6Fail <- read_excel("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Nanostring paper/Data/Outcome_venn_May2025/M6FailureGenes.xlsx")
M6Fail$Gene

M6Recur <- read_excel("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Nanostring paper/Data/Outcome_venn_May2025/M6RecurrentTbGenes.xlsx")
M6Recur$Gene

# look for intersect DE genes 
common_genes = intersect(M6Fail$Gene, M6Recur$Gene)

# Find  genes in Failre only 
Fail_diff_genes <- setdiff(M6Fail$Gene, intersect(M6Fail$Gene, M6Recur$Gene))
Recur_diff_genes <- setdiff(M6Recur$Gene, intersect(M6Fail$Gene, M6Recur$Gene))



# FIND INTERSECTING GENES 

# Define the sets
Fail_only <- length(Fail_diff_genes)
Recur_only <- length(Recur_diff_genes)
both <- length(common_genes)

# Create the Venn diagram with equal-sized circles
venn.plot <- venn.diagram(
  x = list(Failure = M6Fail$Gene, Recurrent_TB = M6Recur$Gene),
  category.names = c("Failure", "Recurrent_TB"),
  filename = NULL)

# Draw the Venn diagram
grid.draw(venn.plot)



# MODEL TX OUTCOME USING M6 DATA
OutcomeDb <- read_excel("~/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/Manuscripts/Nanostring paper/Data/Outcome_venn_May2025/AllData_Outcome.xlsx")
View(OutcomeDb)

# transpose data
Outcome_tran <- t(OutcomeDb)
View(Outcome_tran)
colnames(Outcome_tran) <- Outcome_tran[1, ]
Outcome_tran <- Outcome_tran[-1, ]
colnames(Outcome_tran)

Outcome_tran <- as.data.frame(Outcome_tran)
View(Outcome_tran)
Outcome_tran <- rownames_to_column(Outcome_tran, "Sample_Replacement")
colnames(Outcome_tran)
dim(Outcome_tran)



#Join with attribute file 
Outcome_Attributes <- read_csv("~/OneDrive - London School of Hygiene and Tropical Medicine/Manuscripts/Nanostring paper/Demography/MetaData/Sample_Attribute_all2.csv", show_col_types = FALSE)
colnames(Outcome_Attributes)
colnames(Outcome_Attributes)[2] = "Sample_Replacement"

Attributes <- Outcome_Attributes[, c(2, 5, 7)]
colnames(Attributes)
dim(Attributes)

# merge with metadata 
common_colnm <- intersect(names(Attributes), names(Outcome_tran))
#merge 
Outcome_dd <- merge(Attributes, Outcome_tran, by=common_colnm)
names(Outcome_dd)
dim(Outcome_dd)



#Subset data by timepoint: 
# Month 6 dataset
Outcome_Month_6 <- Outcome_dd[Outcome_dd$Timepoint =="Month_6", ]
Outcome_Month_2 <- Outcome_dd[Outcome_dd$Timepoint =="Month_2", ]
Outcome_Month_4 <- Outcome_dd[Outcome_dd$Timepoint =="Month_4", ]

# inspect month 6 data 
dim(Outcome_Month_6)
names(Outcome_Month_6)

# M6 data for modeling 
data6 = Outcome_Month_6[, -c(1:2)]
names(data6)

# Subset varsiables and recode to good vs poor nomenclature 
common_genes_intersect <- intersect(common_genes, colnames(data6))
finData <- data6[, c(common_genes_intersect, "Outcome"), drop = FALSE]
unique(finData$Outcome)
finData$Outcome = gsub("Cure", "Good", finData$Outcome)
finData$Outcome = gsub("Relapse", "Poor", finData$Outcome)
finData$Outcome = gsub("Failure", "Poor", finData$Outcome)



# Feature elimination and selection 
names(finData)
a <- finData[, -ncol(finData)]
a <- a %>% mutate_all(as.numeric)

# feature selection by filteration 
a_pro <- preProcess(a,method = c("corr", "nzv"),cutoff = 0.8)
a_pro1 <- predict(a_pro,a)
dim(a_pro1)
y <- as.factor(finData$Outcome) # still index matched 
datacomb <- cbind(y, a_pro1)


# # Split data 
library(DMwR)
data1 = datacomb
# Perform SMOTE
oversampled_data <- SMOTE(y ~ ., data = data1, perc.over = 350, k = 5)

set.seed(1234)       

# Initialize variables
auc_values <- numeric(500)
roc_list <- list()
accuracy <- rep(NA, 500)  
importance <- matrix(0, nrow = 500, ncol = ncol(a_pro1))

# Loop for 500 iterations
for (i in 1:500) {
  # Train the model (same as in your previous code)
  # Split data into training and testing sets
  trainIn <- createDataPartition(oversampled_data$y, p = .8, list = FALSE, times = 1)
  training <- oversampled_data[trainIn,]
  testing  <- oversampled_data[-trainIn,]
  
  
  # Train the model
  train_control <- trainControl(method = "cv", number = 5)
  model <- train(y ~ ., data = training, method = "rf", trControl = train_control, ntrees = 500)
  
  # Predict probabilities on testing data
  pred_prob <- predict(model, newdata = testing, type = "prob")
  pred_y <- predict(model, newdata = testing)
  
  # Calculate ROC curve
  roc_data <- roc(testing$y, pred_prob[, "Poor"])
  
  #Plot ROC curve
  #Collect feature importance if 'features' is defined
  feature_import <- varImp(model)
  features <- rownames(feature_import$importance)[which(feature_import$importance>=0.50)]
  importance[i,] <- (names(a_pro1)%in%features)*1
  
  # Calculate accuracy
  cm <- confusionMatrix(pred_y, testing$y, positive = "Poor")
  accuracy[i] <- cm$overall["Accuracy"]
  
  # Calculate ROC curve
  roc_data <- roc(testing$y, pred_prob[, "Poor"])
  
  # Store AUC value
  auc_values[i] <- auc(roc_data)
  
  # Store ROC curve data
  roc_list[[i]] <- roc_data
}


# Summarize AUC values
mean_auc <- mean(auc_values)
sd_auc <- sd(auc_values)
range <- range(auc_values)
IQR <- IQR(auc_values)
lowerQ <- quantile(auc_values)
cat("Mean AUC:", mean_auc, "\n")
cat("Standard Deviation of AUC:", sd_auc, "\n")

# Plot summary ROC curve
plot(roc_list[[1]], type = "n", main = "Summary ROC Curve")
for (i in 1:500) {
  lines(roc_list[[i]], col = i)
}



# Select all models 
Model_index <- which(accuracy >= 0.00)


# Total number of top models
total_top_models <- length(Model_index)

# Count the occurrence of each feature in the top models
feature_counts <- colSums(importance[Model_index, ])

# Calculate the percentage of occurrence for each feature
feature_percentages <- (feature_counts / total_top_models) * 100

# Get the features with 100% occurrence.
selected_features <- which(feature_percentages == 100.00000) # lets start with 70
length(selected_features) # n =157

# Plot
# Plot the occurrences of each feature
barplot(feature_percentages,
        main = "Occurrences of Features in Models with Accuracy > 0.7",
        xlab = "Features",
        ylab = "Occurrences",
        col = "skyblue",
        las = 2)  # Rotate x-axis labels vertically if necessary


# Subset the Data for these features using the feature matrix (a_pro1)
subset_data_feats <- a_pro1[, selected_features] 
featurenames = names(subset_data_feats)



# Construct a refined model: -----

library(caret)
subset_data_set = cbind(y, subset_data_feats)

# Step 1: Split data into training and testing sets
set.seed(123)  # Set seed for reproducibility
trainIndexes <- createDataPartition(subset_data_set$y, p = 0.8, list = FALSE, times = 1)
trainingdata <- subset_data_set[trainIndexes, ]
testingdata <- subset_data_set[-trainIndexes, ]

# Step 2: Train a machine learning model
model2 <- train(y ~ ., data = trainingdata, method = "rf", ntree=500)

# Step 3: Evaluate the model's performance
prediction1 <- predict(model2, newdata = testingdata, type = "prob")
prediction1b <- predict(model2, newdata = testingdata)

# Calculate accuracy
conf_matrix = confusionMatrix(prediction1b, testingdata$y, positive = "Poor")
print(conf_matrix)
Sub_accuracy <- conf_matrix$overall["Accuracy"]
print(Sub_accuracy)


# Calculate ROC curve
Sub_roc_data <- roc(testingdata$y, prediction1[, "Poor"])
# Roc AUC value
Subs_AUC <- auc(Sub_roc_data)
print(Subs_AUC)


# Predict outcome of earlier tretament timepoints 


# Predict Baseline ------
Outcome_Baseline 

# format data 
# Outcome_Month_6 data 
names(Outcome_Baseline)
Outcome_Baseline$Outcome

# M6 data for modeling 
Outcome_Baseline = Outcome_Baseline[, -c(1, 2)]
names(Outcome_Baseline)

Outcome_Baseline$Outcome = gsub("Cure", "Good", Outcome_Baseline$Outcome)
Outcome_Baseline$Outcome = gsub("Relapse", "Poor", Outcome_Baseline$Outcome)
Outcome_Baseline$Outcome = gsub("Failure", "Poor", Outcome_Baseline$Outcome)
dim(Outcome_Baseline)
Outcome_Baseline = Outcome_Baseline[, c(2:774, 1)]

# Outcome_Month_2 

# load data
dd = Outcome_Baseline


names(Outcome_Baseline)
# remove response variable (last col)
d <- dd[, -ncol(dd)] # all predictors

# make all predictors numeric
d <- d %>% mutate_all(as.numeric) # make numeric 
names(d)

#response 
y <- as.factor(dd$Outcome) # make and format y

#now combine 
combine_M2 = cbind(y, d)

# Subset 
combine_M2 = combine_M2[, c("y", featurenames)]

# rename 
validationdata = combine_M2
names(validationdata)

# Step 3: Evaluate the model's performance
Predictprob <- predict(model2, newdata = validationdata, type = "prob")
predictclass <- predict(model2, newdata = validationdata) 

# Calculate accuracy
cfm = confusionMatrix(predictclass, validationdata$y, positive = "Poor")
print(cfm) # AUC
M2_accuracy <- cfm$overall["Accuracy"]
print(M2_accuracy)

# Calculate ROC curve
rocval <- roc(validationdata$y, Predictprob[, "Poor"])
# Roc AUC value
M2_AUC <- auc(rocval) 
print(M2_AUC)


# Month 2 Predictions -----
week2data = Outcome_Month_2 # note its month 2 dont want to bother changing text

# format data 
Outcome_Month_2 = Outcome_Month_2[, -c(1, 2)]

Outcome_Month_2$Outcome = gsub("Cure", "Good", Outcome_Month_2$Outcome)
Outcome_Month_2$Outcome = gsub("Relapse", "Poor", Outcome_Month_2$Outcome)
Outcome_Month_2$Outcome = gsub("Failure", "Poor", Outcome_Month_2$Outcome)
Outcome_Month_2 = Outcome_Month_2[, c(2:774, 1)]

dd = Outcome_Month_2

# remove response variable (last col)
d <- dd[, -ncol(dd)] # all predictors

# make all predictors numeric
d <- d %>% mutate_all(as.numeric) # make numeric 
names(d)

#response 
y <- as.factor(dd$Outcome) # make and format y


#now combine 
combine_M2 = cbind(y, d)

summary(combine_M2)

# rename 
validationdata = combine_M2
names(validationdata)

# Step 3: Evaluate the model's performance
Predictprob <- predict(model2, newdata = validationdata, type = "prob")
predictclass <- predict(model2, newdata = validationdata) 

# Calculate accuracy
cfm = confusionMatrix(predictclass, validationdata$y, positive = "Poor")
print(cfm) # AUC
M2_accuracy <- cfm$overall["Accuracy"]
print(M2_accuracy)

# Calculate ROC curve
rocval <- roc(validationdata$y, Predictprob[, "Poor"])
# Roc AUC value
M2_AUC <- auc(rocval) 
print(M2_AUC)



# Plot ROC curves -----
# Plot 
library(pROC)


# initial tiff object 
png("TxOutcomes_roc_plot.png", 
    width = 6, height = 6, units = "in", res = 300)

# Color pallet #
# Baseline = Black
# week 2 = brown
# Month2 = Blue 
# M6 = Green

# Plot the smoothed ROC curves
plot(BaselineROCdata, type = "l", col = "black", lty = 1, lwd = 2, #main = "ROC Curve", 
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(M6ROCdata, col = "green", lty = 1, lwd = 2)
lines(Month2ROCdata, col = "brown", lty = 1, lwd = 2)


legend("bottomright", legend = c(paste("Month 6 AUC =", sprintf("%.3f", M6AUC)),
                                 paste("Month 2 AUC =", round(Mt2AUC, 3)),
                                 paste("Baseline AUC =", round(BaseAUC, 3))),
       col = c("green", "brown", "black"), lty = 1, cex = 1)

dev.off()

#save.image("TreatmentOutcomemodelling.RData") 






# PART 4: GENE SET ANALYSES WITH DAVID ------
#---------------------------------------------------#

# DAVID 

# Cluster 2 genes
cluster2_genes <- c("ACE", "ACVR1B", "ADAMTS1", "AHR", "AKT1", "ALB", "ALOX5", "AR", "APOE", 
                    "ATP1B1", "BAX", "BCL2", "BRAF", "BRCA1", "CASP3", "CCL2", "CDK2", 
                    "CFTR", "CTNNB1", "CXCL8", "EGFR", "ERBB2", "FAS", "GAPDH", "GSK3B", 
                    "HLA-A", "IGF1", "IKBKG", "IL1B", "IL6", "ITGAV", "JAK2", "MAPK1", 
                    "NFKB1", "PIK3CA", "PTEN", "SRC")

# Cluster 3 genes
cluster3_genes <- c("AKT3", "ALOX12", "ALOX5", "APOBEC3G", "ARRB2", "ATG7", "C5AR1", "CALM1", 
                    "CASP3", "CCL28", "CCR4", "CCR5", "CCRL2", "CD247", "CD3E", "CD3G", 
                    "CD4", "CD40", "CD40LG", "CD8B", "CDK4", "CEACAM3", "CPA3", "CR1", 
                    "CSF2RA", "CSF2RB", "CTSA", "CTSS", "CX3CR1", "CXCL1", "CXCL10", "CXCL2", 
                    "CXCR2", "CXCR6", "DNAJA2", "F5", "FCAR", "FCRL2", "FOS", "FURIN", 
                    "GBP4", "GLA", "GSK3B", "HAVCR2", "HLA-A", "HLA-C", "HLA-DPB1", "HLA-E", 
                    "HMOX1", "HPGD", "IFI16", "IFI6", "IFITM3", "IFNAR2", "IKBKG", "IL12RB1", 
                    "IL16", "IL18", "IL1R2", "IL21R", "IL2RG", "IRF9", "ITGAX", "JAK3", 
                    "KDM6B", "KIR2DL3", "KPNB1", "LAMP2", "LCN2", "LCP1", "LDHB", "LILRB2", 
                    "LYN", "MAP1LC3A", "MAP2K3", "MAPK14", "MT2A", "MVP", "NAE1", "NAMPT", 
                    "NCF2", "NCF4", "NDUFS8", "NFATC3", "NFE2L2", "NFKB2", "NKG7", "NT5E", 
                    "OAS3", "OASL", "PANX1", "PECAM1", "PELI1", "PIK3CB", "PLCG2", "PLEK", 
                    "PRKCQ", "PRKCSH", "PSAP", "PSEN1", "PSMB8", "PSMB9", "PSTPIP1", "PTPN6", 
                    "RAF1", "RBCK1", "RNF31", "S100A12", "SELL", "SERPINA1", "SMAD3", "SP1", 
                    "SPI1", "SPIB", "STAT3", "STAT4", "TBX21", "TCF7", "TCIRG1", "TIFA", 
                    "TLR1", "TLR5", "TLR7", "TNF", "TNFRSF10B", "TNFRSF1A", "TOLLIP", "TRIM21", 
                    "TRIM22", "TRIM25", "GUCY1B1", "JAML", "CD45R0", "FCGR1A/B", "CD45RB")

# Cluster 4 genes
cluster4_genes <- c("AGT", "AKT1", "ALB", "AMPK", "BAX", "BCL2", "CDH1", "CDKN2A", "E2F1", 
                    "EGR1", "FAS", "GATA3", "HMOX1", "HSP90AB1", "IGF1", "IL6", "JAK2", 
                    "KLF4", "MAPK14", "NFKB1", "NRF2", "P53", "PIK3CA", "PTEN", "RELA", 
                    "STAT3", "TGFBR2", "VEGF", "ZEB1")

# Cluster 6 genes (updated)
cluster6_genes <- c("AKT3", "ALOX12", "ALOX5", "APOBEC3G", "ARRB2", "ATG7", "C5AR1", "CALM1", 
                    "CASP3", "CCL28", "CCR4", "CCR5", "CCRL2", "CD247", "CD3E", "CD3G", 
                    "CD4", "CD40", "CD40LG", "CD8B", "CDK4", "CEACAM3", "CPA3", "CR1", 
                    "CSF2RA", "CSF2RB", "CTSA", "CTSS", "CX3CR1", "CXCL1", "CXCL10", "CXCL2", 
                    "CXCR2", "CXCR6", "DNAJA2", "F5", "FCAR", "FCRL2", "FOS", "FURIN", 
                    "GBP4", "GLA", "GSK3B", "HAVCR2", "HLA-A", "HLA-C", "HLA-DPB1", "HLA-E", 
                    "HMOX1", "HPGD", "IFI16", "IFI6", "IFITM3", "IFNAR2", "IKBKG", "IL12RB1", 
                    "IL16", "IL18", "IL1R2", "IL21R", "IL2RG", "IRF9", "ITGAX", "JAK3", 
                    "KDM6B", "KIR2DL3", "KPNB1", "LAMP2", "LCN2", "LCP1", "LDHB", "LILRB2", 
                    "LYN", "MAP1LC3A", "MAP2K3", "MAPK14", "MT2A", "MVP", "NAE1", "NAMPT", 
                    "NCF2", "NCF4", "NDUFS8", "NFATC3", "NFE2L2", "NFKB2", "NKG7", "NT5E", 
                    "OAS3", "OASL", "PANX1", "PECAM1", "PELI1", "PIK3CB", "PLCG2", "PLEK", 
                    "PRKCQ", "PRKCSH", "PSAP", "PSEN1", "PSMB8", "PSMB9", "PSTPIP1", "PTPN6", 
                    "RAF1", "RBCK1", "RNF31", "S100A12", "SELL", "SERPINA1", "SMAD3", "SP1", 
                    "SPI1", "SPIB", "STAT3", "STAT4", "TBX21", "TCF7", "TCIRG1", "TIFA", 
                    "TLR1", "TLR5", "TLR7", "TNF", "TNFRSF10B", "TNFRSF1A", "TOLLIP", "TRIM21", 
                    "TRIM22", "TRIM25", "GUCY1B1", "JAML", "CD45R0", "FCGR1A/B", "CD45RB")

# Combine all clusters
all_genes <- c(cluster2_genes, cluster3_genes, cluster4_genes, cluster6_genes)
all_genes_Cl2.4.6 <- c(cluster2_genes, cluster4_genes, cluster6_genes)

# Save the list to a text file
write.table(all_genes, file = "DAVID_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(all_genes, file = "~/Desktop/DAVID_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_genes, file = "~/Desktop/DAVID_gene_list246.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)




# Load necessary libraries
setwd("~/OneDrive - London School of Hygiene and Tropical Medicine/Manuscripts/Nanostring paper/DAVID")
library(ggplot2)

# Read the DAVID output .txt file
# Adjust `sep` if your file uses different delimiters (e.g., "\t" for tab-separated files)
david_results <- read.table("Top50bp.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preview the first few rows of the data
head(david_results)

# Ensure the columns for biological processes and enrichment scores are correctly identified
# Replace "Term" and "EnrichmentScore" with the actual column names in your file
# Check column names with: colnames(david_results)
david_results$EnrichmentScore <- as.numeric(david_results$Fold.Enrichment)  # Ensure numeric
david_results <- david_results[!is.na(david_results$Fold.Enrichment), ]     # Remove NA values

# Calculate -log10(FDR)
david_results$LogFDR <- -log10(david_results$FDR)

# Remove "GO: 7-digit number~" prefix
david_results$Term <- sub("GO:\\d+~", "", david_results$Term)

# Select the top 10 biological processes
top_bp <- david_results[order(david_results$FDR), ][1:10, ]

# Define the significance threshold (-log10(0.05))
significance_threshold <- -log10(0.05)


# initiate 
png("Top10_DwnregBP.png", width = 10, height = 6, units = "in", res = 300)

# Create the bar plot
ggplot(top_bp, aes(x = reorder(Term, LogFDR), y = LogFDR)) +
  geom_bar(stat = "identity", fill = "green") +
  geom_hline(yintercept = significance_threshold, 
             linetype = "dashed", color = "white", size = 1) +
  coord_flip() +
  labs(
    title = NULL,
    x = NULL,
    y = "-log10 adj p value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1) )

#
dev.off()

# Unregulated cluster 
cluster5_genes <- c("ACSL3", "ACVR1", "ADORA2A", "AHR", "AKT1", "ALOX15", "ALPK1", "ANPEP", 
                    "AP1M1", "APEX1", "APP", "ATG10", "ATG12", "BATF", "BCL2", "BCL3",
                    "BCR", "BECN1", "BLK", "BNIP3", "BPI", "C3AR1", "CARD11", "CARD16",
                    "CASP5", "CCL17", "CCNC", "CCR3", "CCR7", "CD163", "CD2", "CD244", 
                    "CD38", "CD6", "CD69", "CD79B", "CD8A", "CDH1", "CSF1", "CSF1R", 
                    "CTLA4", "CTSL", "CTSW", "CUL1", "CXCL16", "CXCL5", "CXCL8", "CXCR5", 
                    "DDAH2", "DDIT3", "DEFA4", "DHX58", "EGLN1", "EIF2AK2", "EIF2AK3", "EOMES",
                    "EPHX2", "FASLG", "FBXO6", "FOXO1", "FOXP3", "GAB2", "GATA3", "GLB1", 
                    "GSTM4", "GZMB", "HDC", "HLA-DMB", "HLA-DOB", "HMGB1", "ICOS", "ICOSLG", 
                    "IDO1", "IFI27", "IFI35", "IGFBP7", "IKBKB", "IKBKE", "IL10RB", "IL12RB2", 
                    "IL15", "IL18BP", "IL18R1", "IL1R1", "IL1RL1", "IL23A", "IL23R", "IL27RA", 
                    "IL2RA", "IL2RB", "IL3RA", "IL5RA", "IRAK3", "IRF4", "IRF7", "ISG15", 
                    "ITGAM", "ITGB7", "ITK", "ITLN1", "ITPR3", "JUN", "KIR2DL1", "KLRC1", 
                    "KRAS", "LAG3", "LANCL1", "LILRA3", "LTBR", "LTF", "MAFB", "MAP2K7", 
                    "MAPK13", "MAPK8", "MAPK9", "MEFV", "MLKL", "MS4A2", "MS4A4A", "MSRA", 
                    "MTOR", "MX1", "MYD88", "NCR1", "NCR3", "NEO1", "NFAT5", "NFATC1", 
                    "NFKB1", "NGLY1", "NLRP3", "NOTCH1", "NRAS", "NTNG2", "OS9", "OSM", 
                    "P2RX7", "PARP1", "PDHB", "PELI2", "PFKFB3", "PIK3C3", "PIK3R4", "PLAUR", 
                    "PLCG1", "PNOC", "PRDM1", "PRF1", "PRKCA", "PTGER2", "PTGER4", "RASGRP4", 
                    "RELA", "RELB", "RIPK1", "RIPK2", "RNF135", "RPS6KB1", "RUNX3", "SCARB2", 
                    "SH2D1A", "SIGLEC5", "SLC11A1", "SMAD4", "SMAD5", "SOCS1", "SOCS3", "SOD1", 
                    "STAT2", "SUGT1", "TAB1", "TBK1", "TCL1A", "TCN2", "THBS1", "THOP1", 
                    "TIGIT", "TLR9", "TNFRSF17", "TNFRSF25", "TNFRSF4", "TNFRSF9", "TNFSF4", "TRAF2", 
                    "TRAF3", "TRAF6", "TRIM5", "TXK", "TYK2", "ULK1", "ULK2", "VEGFA", 
                    "VRK3", "WIPI1", "XAF1", "FAM30A", "SEM1", "XCL1/2", "SELENOS", "KIR3DL1/2")

write.table(cluster5_genes, file = "Cluster5_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Run DAVID and download output as UpTop50bp.txt

# load data 
david_upreg <- read.table("UpTop50bp.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preview 
head(david_upreg)

# Ensure the columns for biological processes and enrichment scores are correctly identified
david_upreg$EnrichmentScore <- as.numeric(david_upreg$Fold.Enrichment)  # Ensure numeric
david_upreg <- david_upreg[!is.na(david_upreg$Fold.Enrichment), ]     # Remove NA values

# Calculate -log10(FDR)
david_upreg$LogFDR <- -log10(david_upreg$FDR)

# Remove "GO: 7-digit number~" prefix
david_upreg$Term <- sub("KW-\\d{4}~", "", david_upreg$Term)
david_upreg$Term

# Select the top 10 biological processes
top_10bp <- david_upreg[order(david_upreg$FDR), ][1:10, ]

# Define the significance threshold (-log10(0.05))
sign_threshold <- -log10(0.05)


# Inititae tiff 
#tiff("Top10_UpregBP.tiff", width = 7, height = 10, units = "in", res = 300)
png("Top10_UpregBP.png", width = 6, height = 6, units = "in", res = 300)

# Create the bar plot
ggplot(top_10bp, aes(x = reorder(Term, LogFDR), y = LogFDR)) +
  geom_bar(stat = "identity", fill = "purple") +
  geom_hline(yintercept = sign_threshold, 
             linetype = "dashed", color = "white", size = 1) +
  coord_flip() +
  labs(
    title = NULL,
    x = NULL,
    y = "-log10 adj p value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1) )

#
dev.off()




# Run DAVID KEGG for Outcome differences 

# Run DAVID and download output as UpTop50bp.txt

# load data 
KeggAll <- read.table("Goodpoor 14112024/Txt/KEGG_all339.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preview 
head(KeggAll)

# Ensure the columns for biological processes and enrichment scores are correctly identified
KeggAll$EnrichmentScore <- as.numeric(KeggAll$Fold.Enrichment)  # Ensure numeric
KeggAll <- KeggAll[!is.na(KeggAll$Fold.Enrichment), ]     # Remove NA values

# Calculate -log10(FDR)
KeggAll$LogFDR <- -log10(KeggAll$FDR)

# Remove "GO: 7-digit number~" prefix
# Removing everything before and including the colon
KeggAll$Term <- gsub("^[^:]*: ?", "", KeggAll$Term)
# View the result
print(KeggAll$Term)


# Select the top 10 biological processes
top_10Kegg <- KeggAll[order(KeggAll$FDR), ][1:10, ]

# Define the significance threshold (-log10(0.05))
signal_threshold <- -log10(0.05)


# Inititae tiff 
#tiff("Top10_UpregBP.tiff", width = 7, height = 10, units = "in", res = 300)
png("Top10_Kegg_pathwyas.png", width = 6, height = 6, units = "in", res = 300)

# Create the bar plot
ggplot(top_10Kegg, aes(x = reorder(Term, LogFDR), y = LogFDR)) +
  geom_bar(stat = "identity", fill = "#404040") +
  geom_hline(yintercept = signal_threshold, 
             linetype = "dashed", color = "white", size = 1) +
  coord_flip() +
  labs(
    title = NULL,
    x = NULL,
    y = "-log10 adj p value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1) )

#
dev.off()


# 2. Reactome pathways 

# load data 
ReactomeAll <- read.table("Goodpoor 14112024/Txt/Reactome339.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preview 
head(ReactomeAll)

# Ensure the columns for biological processes and enrichment scores are correctly identified
ReactomeAll$EnrichmentScore <- as.numeric(ReactomeAll$Fold.Enrichment)  # Ensure numeric
ReactomeAll <- ReactomeAll[!is.na(ReactomeAll$Fold.Enrichment), ]     # Remove NA values

# Calculate -log10(FDR)
ReactomeAll$LogFDR <- -log10(ReactomeAll$FDR)

# Remove "GO: 7-digit number~" prefix
# Removing everything up to and including the tilde (~)
ReactomeAll$Term <- gsub("^[^~]*~", "", ReactomeAll$Term)

# View the result
print(ReactomeAll$Term)



# Select the top 10 biological processes
top_10Reactome <- ReactomeAll[order(ReactomeAll$FDR), ][1:10, ]

# Define the significance threshold (-log10(0.05))
signal_react <- -log10(0.05)


# Inititae tiff 
#tiff("Top10_UpregBP.tiff", width = 7, height = 10, units = "in", res = 300)
png("Top10_Reactome_pathwyas.png", width = 6, height = 6, units = "in", res = 300)

# Create the bar plot
ggplot(top_10Reactome, aes(x = reorder(Term, LogFDR), y = LogFDR)) +
  geom_bar(stat = "identity", fill = "#404040") +
  geom_hline(yintercept = signal_react, 
             linetype = "dashed", color = "white", size = 1) +
  coord_flip() +
  labs(
    title = NULL,
    x = NULL,
    y = "-log10 adj p value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1) )

#
dev.off()

# 3. Go BP 

# load data 
BP.All <- read.table("Goodpoor 14112024/Txt/BP_339.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preview 
head(BP.All)

# Ensure the columns for biological processes and enrichment scores are correctly identified
BP.All$EnrichmentScore <- as.numeric(BP.All$Fold.Enrichment)  # Ensure numeric
BP.All <- BP.All[!is.na(BP.All$Fold.Enrichment), ]     # Remove NA values

# Calculate -log10(FDR)
BP.All$LogFDR <- -log10(BP.All$FDR)

# Remove "GO: 7-digit number~" prefix
BP.All$Term <- sub("GO:\\d+~", "", BP.All$Term)

# View the result
print(BP.All$Term)



# Select the top 10 biological processes
top_10go <- BP.All[order(BP.All$FDR), ][1:10, ]

# Define the significance threshold (-log10(0.05))
signal_GO <- -log10(0.05)


# Inititae tiff 
#tiff("Top10_UpregBP.tiff", width = 7, height = 10, units = "in", res = 300)
png("Top10_GO_BP_pathwyas.png", width = 6, height = 6, units = "in", res = 300)

# Create the bar plot
ggplot(top_10go, aes(x = reorder(Term, LogFDR), y = LogFDR)) +
  geom_bar(stat = "identity", fill = "#404040") +
  geom_hline(yintercept = signal_GO, 
             linetype = "dashed", color = "white", size = 1) +
  coord_flip() +
  labs(
    title = NULL,
    x = NULL,
    y = "-log10 adj p value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1) )

#
dev.off()




# 4. Molecular fiunctions 

# Top 10 MF


# load data 
MF.All <- read.table("Goodpoor 14112024/Txt/BP_339.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Preview 
head(MF.All)

# Ensure the columns for biological processes and enrichment scores are correctly identified
MF.All$EnrichmentScore <- as.numeric(MF.All$Fold.Enrichment)  # Ensure numeric
MF.All <- MF.All[!is.na(MF.All$Fold.Enrichment), ]     # Remove NA values

# Calculate -log10(FDR)
MF.All$LogFDR <- -log10(MF.All$FDR)

# Remove "GO: 7-digit number~" prefix
MF.All$Term <- sub("GO:\\d+~", "", MF.All$Term)

# View the result
print(MF.All$Term)


# Select the top 10 biological processes
top_10MF <- MF.All[order(MF.All$FDR), ][1:10, ]

# Define the significance threshold (-log10(0.05))
signal_MF <- -log10(0.05)


# Inititae tiff 
#tiff("Top10_UpregBP.tiff", width = 7, height = 10, units = "in", res = 300)
png("Top10_GO_MF_pathwyas.png", width = 6, height = 6, units = "in", res = 300)

# Create the bar plot
ggplot(top_10MF, aes(x = reorder(Term, LogFDR), y = LogFDR)) +
  geom_bar(stat = "identity", fill = "#404040") +
  geom_hline(yintercept = signal_MF, 
             linetype = "dashed", color = "white", size = 1) +
  coord_flip() +
  labs(
    title = NULL,
    x = NULL,
    y = "-log10 adj p value") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1) )

#
dev.off()


# PART 5: THESES PLOTS -----
#---------------------------------------------------------------------#

#-----------------------------------------------------# 
# Figure 1: Volcano Plot
#-----------------------------------------------------# 

# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Set working directory
setwd("~/Library/..")

# Read data
data <- read.delim("Month_2 vs Baseline Comparisons.txt", header = TRUE, stringsAsFactors = FALSE)

names(data)
dim(data) # n127 genes
data = data[-c(2:5)]
names(data)[1] = "Gene"


# Create column for significance
data <- data %>%
  mutate(
    neg_log10_p = -log10(padj),
    significance = case_when(
      padj < 0.05 & log2FoldChange > 0.585 ~ "Up",
      padj < 0.05 & log2FoldChange < -0.585 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Add labels for top genes
top_down <- data %>%
  filter(padj < 0.05, log2FoldChange < 0.585) %>%
  arrange(log2FoldChange) %>%
  slice(1:10)

top_up <- data %>%
  filter(padj < 0.05, log2FoldChange > 0.585) %>%
  arrange(desc(log2FoldChange)) %>%
  slice(1:5)

# Combine the two sets
top_genes_to_label <- bind_rows(top_down, top_up)
top_genes_to_label$Gene


# Define custom colors
volcano_colors <- c("Up" = "green", "Down" = "purple", "NS" = "grey")


volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = volcano_colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    title = "Differential Expression Volcano Plot",
    color = "Regulation"
  ) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  geom_text_repel(
    data = top_genes_to_label,
    aes(label = Gene),
    size = 3,
    max.overlaps = 100
  )


# Save or print plot
print(volcano_plot)
ggsave("Figure1_updated.pdf", volcano_plot, width = 6, height = 6)



#----------------------------------------------------- # 
# Figure 3A Slow vs fast: Heatmap at baseline -------
#----------------------------------------------------- # 

# Load libraries
library(ComplexHeatmap)
library(grid)

# Set working directory
setwd("~/Library/...")

# Read data
data = read.csv("Baseline_Data_Table.csv")
samples = colnames(data)[c(-1, -61)]

# Split colname and convert to data frame to get sample conditions
sampleCond = as.data.frame(t(do.call(data.frame, strsplit(samples, "\\.\\."))))
colnames(sampleCond) = c("sample_id", "condition")

# Fix row cluster names
data$Cluster[data$Cluster == 1] = "Up regulated"
data$Cluster[data$Cluster == 2] = "Down regulated"

# Define custom colors for annotations
row_colors <- c("Up regulated" = "green", "Down regulated" = "purple")
column_colors <- c("Slow responder" = "#FDAE61", "Fast responder" = "blue")

# Extract gene expression matrix
dat = data[, c(-1, -61)]

# Set rownames to gene names 
rownames(dat) <- data[[1]] 

# Convert Cluster to factor with correct levels
data$Cluster <- factor(data$Cluster, levels = c("Up regulated", "Down regulated"))

# Double-check order is preserved
stopifnot(all(rownames(dat) == data[[1]]))  # should return TRUE

# Row annotation
rha <- rowAnnotation(
  Cluster = data$Cluster,
  col = list(Cluster = c("Up regulated" = "green", "Down regulated" = "purple")),
  show_annotation_name = FALSE)


# Create top annotation (column-wise)
ha <- HeatmapAnnotation(
  condition = sampleCond$condition,
  col = list(condition = column_colors),
  gp = gpar(col = "black", lwd = 0.5),
  show_legend = TRUE
)

# Plot heatmap
fig3a = Heatmap(
  dat,
  left_annotation = rha,
  top_annotation = ha,
  row_split = data$Cluster,
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_title = NULL,
  name = "Expression")

print(fig3a)

# save heatmap
pdf("Fig_3A_heatmap.pdf", width = 8, height = 8)
print(fig3a)
dev.off()




#----------------------------------------------------- # 
# Figure 3B Slow vs fast: Heatmap at Month 2 -------
#----------------------------------------------------- # 

# Load libraries
library(ComplexHeatmap)
library(grid)

# Set working directory
setwd("~/Library/...")

# Read data
ddata = read.csv("Month2_Data_Tablecsv.csv")
samples = colnames(ddata)[c(-1, -61)]

# Split colname and convert to data frame to get sample conditions
sampleCond = as.data.frame(t(do.call(data.frame, strsplit(samples, "\\.\\."))))
colnames(sampleCond) = c("sample_id", "condition")

# Rename condition labels
sampleCond$condition = gsub("Slow.Month_2.", "Slow responder", sampleCond$condition)
sampleCond$condition = gsub("Fast.Month_2.", "Fast responder", sampleCond$condition)


# Fix row cluster names
ddata$Cluster[ddata$Cluster == 1] = "Up regulated"

# Define custom colors for annotations
row_colors <- c("Up regulated" = "green", "Down regulated" = "purple")
column_colors <- c("Slow responder" = "#FDAE61", "Fast responder" = "blue")

# Extract gene expression matrix
dat = ddata[, c(-1, -61)]

# Set rownames to gene names 
rownames(dat) <- ddata[[1]] 
rownames(dat)
# [1] "GBP5"  "CASP5" "GBP1" 


# Convert Cluster to factor with correct levels
ddata$Cluster <- factor(ddata$Cluster, levels = c("Up regulated", "Down regulated"))


# Row annotation
rha <- rowAnnotation(
  Cluster = ddata$Cluster,
  col = list(Cluster = c("Up regulated" = "green", "Down regulated" = "purple")),
  show_annotation_name = FALSE)


# Create top annotation (column-wise)
ha <- HeatmapAnnotation(
  condition = sampleCond$condition,
  col = list(condition = column_colors),
  gp = gpar(col = "black", lwd = 0.5),
  show_legend = TRUE
)

# Plot heatmap
fig3b = Heatmap(
  dat,
  left_annotation = rha,
  top_annotation = ha,
  row_split = ddata$Cluster,
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_title = NULL,
  name = "Expression")

print(fig3b)


# save heatmap
pdf("Fig_3B_heatmap.pdf", width = 8, height = 8)
print(fig3b)
dev.off()

# save rds
save.image("Figure_3B_updated090425.Rdata")



#-------------------------------------------------------------------#
# Figure 5 
#-------------------------------------------------------------------#

# load libraries
library(readxl)
library(ggplot2)
library(readr)
library(ggplot2)
library(reshape2)

## Set workign directory
setwd("~/Library/..")

#load cell data 

celldata <- read_excel("CellScoreData.xlsx")
dim(celldata)
names(celldata)
names(celldata)[1] = "Sample_Name"

# load metadata 
metadata <- read_csv("MetaData.csv")
names(metadata)
metadata$Outcome

# subset outcome and sample id only
metta = metadata[, c(1, 4, 6)]
metta$Outcome = as.factor(as.character(metta$Outcome))
names(metta)[1] = "Sample_Name"

# Merge with Celldata
Merged_data = merge(celldata, metta, by= "Sample_Name")
dim(Merged_data)
names(Merged_data)

#  *** Subset baseline data ***
Baseline_data = subset(Merged_data, Timepoint =="Baseline")
names(Baseline_data)
Baseline_data = Baseline_data[c(-1,-8)] # drop sample name column
names(Baseline_data)


# Reshape the data for ggplot
cellplotdata_long <- melt(Baseline_data, id.vars = "Outcome", 
                          variable.name = "Cell_Type", 
                          value.name = "Score")


p5 = ggplot(cellplotdata_long, aes(x = Outcome, y = Score)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.15, size = 1, alpha = 0.5, aes(color = Outcome)) + 
  scale_color_manual(values = c("Cure" = "orange", 
                                "Failure" = "#654321",  
                                "Recurrent" = "blue")) +  
  labs(y = "Cell-type score") +
  facet_wrap(~ Cell_Type, scales = "free_y") + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),  
        panel.grid = element_blank(), 
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = "lightgray", color = "gray"))

print(p5)

# save pdf 
pdf("Fig_5_updated.pdf", width = 8, height = 8)
print(p5)
dev.off()

