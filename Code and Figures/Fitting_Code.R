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
  
  DT.fit = 
  
  
  
}




