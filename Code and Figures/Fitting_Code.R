##############################
#### Code for simulation #####
##############################


seedValue<-1000
noOfObservations <-200  
noOfIrrelevantFeatures<- 10  
noOfRelevantFeatures<-9  
simulatedDataLabelName <- "class_label"
probability_of_label_1 <-70
noiseAmplitudes<- "0 , .1,10"
correlationDirectionality<- "1,-1,2,-56"  




create_train_test <- function(data, size = 0.8, train = TRUE) {
  n_row = nrow(data)
  total_row = size * n_row
  train_sample <- 1: total_row
  if (train == TRUE) {
    return (data[train_sample, ])
  } else {
    return (data[-train_sample, ])
  }
}







simulatedData<- getSimulatedData(noOfObservations = noOfObservations, 
                                 noOfIrrelevantFeatures = noOfIrrelevantFeatures,
                                 noOfRelevantFeatures = noOfRelevantFeatures,  
                                 simulatedDataLabelName = simulatedDataLabelName,
                                 probability_of_label_1 = probability_of_label_1,
                                 noiseAmplitudes = noiseAmplitudes, 
                                 correlationDirectionality = correlationDirectionality, 
                                 seedValue = seedValue) 


for(i in 1:1000){
  set.seed(i)
  simulatedData<- getSimulatedData(noOfObservations = noOfObservations, 
                                   noOfIrrelevantFeatures = noOfIrrelevantFeatures,
                                   noOfRelevantFeatures = noOfRelevantFeatures,  
                                   simulatedDataLabelName = simulatedDataLabelName,
                                   probability_of_label_1 = probability_of_label_1,
                                   noiseAmplitudes = noiseAmplitudes, 
                                   correlationDirectionality = correlationDirectionality, 
                                   seedValue = seedValue) 
  
  simulatedData$class_label = as.integer(as.logical(simulatedData$class_label))
  
  
  data_train <- create_train_test(data = simulatedData, size = 0.8, train = TRUE)
  data_test  <- create_train_test(data = simulatedData, size = 0.8, train = FALSE)
  
  
  
  DT.fit = rpart(class_label~., data_train, method = "class")
  logistic.fit= glm(formula = class_label~., family = binomial("logit"),data_train )  
    
    
  
  
  
}




