#*******************************************************************************
# Project: MUTATED - roboBayes
# Script Purpose: Train Random Forest
# Date: 2023-04-06
# Author: Jenna Abrahamson
#*******************************************************************************
library(terra)
library(data.table)
library(caret)
library(randomForest)
library(doParallel)

# Calculate number of cores
no_cores <- detectCores() - 1

# Create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)


# Load training data set
train_data_path <- "./all_regions.Rdata"
load(train_data_path)

# ~~~~~~~~~~~~~
# Train RF
# ~~~~~~~~~~~~~
global_train <- final_features[,3:115]
train <- na.omit(global_train)

# 10 fold cross validation
control <- trainControl(method='repeatedcv', 
                        number=5, 
                        allowParallel=TRUE,
                        search='grid')

# Create tune grid
mtry <- sqrt(ncol(train))
tunegrid <- expand.grid(.mtry=mtry)

print("Made it to trianing!")

# Train Random Forest
global_RF <- train(as.factor(Class) ~ ., 
                   data = train, 
                   method = "rf", 
                   metric = "Accuracy", 
                   tuneGrid = tunegrid, 
                   trControl = control)

# Get summary stats
print(global_RF)


# Write out model to load later
save(global_RF, file="./models/RF_all.Rdata")

# Done
stopCluster(cl)
