##############################
#### Code for simulation #####
##############################


seedValue<-1000
noOfObservations <-200  
noOfIrrelevantFeatures<-2  
noOfRelevantFeatures<-7  
simulatedDataLabelName <- "class_label"
probability_of_label_1 <-70
noiseAmplitudes<- "0 , .1,10"
correlationDirectionality<- "1,-1,2,-56"  

simulatedData<- getSimulatedData(noOfObservations = noOfObservations, 
                                 noOfIrrelevantFeatures = noOfIrrelevantFeatures,
                                 noOfRelevantFeatures = noOfRelevantFeatures,  
                                 simulatedDataLabelName = simulatedDataLabelName,
                                 probability_of_label_1 = desiredNumberOfPositives,
                                 noiseAmplitudes = noiseAmplitudes, 
                                 correlationDirectionality = correlationDirectionality, 
                                 seedValue = seedValue) 
