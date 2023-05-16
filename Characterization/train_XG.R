#*******************************************************************************
# Project: MUTATED - roboBayes
# Script Purpose: Train XGBoost
# Date: 2023-03-04
# Author: Jenna Abrahamson
#*******************************************************************************
# Load libraries
library(terra)
library(data.table)
library(caret)
library(xgboost)
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
# Train XGBoost
# ~~~~~~~~~~~~~
global_train <- final_features[,3:115]
train <- na.omit(global_train)

# Specify cross validation method
xgb_trcontrol = trainControl(
  method = "cv",
  number = 5,  
  allowParallel = TRUE,
  verboseIter = FALSE,
  returnData = FALSE
)

# Create grid space to search for best hyperparameters
xgbGrid <- expand.grid(nrounds = c(100,200),  
                       max_depth = c(10, 15, 20, 25),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1
                      )

print("made it to train function!")

# Train XG Boost
global_XG <- train(as.factor(Class) ~ ., 
                   data = train, 
                   trControl = xgb_trcontrol,
                   tuneGrid = xgbGrid,
                   method = "xgbTree")

# Get summary stats
print(global_XG)

# Write out model to load later
save(global_XG, file="./models/globalXG_all.Rdata")

# Done
stopCluster(cl)
