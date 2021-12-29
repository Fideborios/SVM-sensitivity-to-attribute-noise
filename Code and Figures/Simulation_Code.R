##############################
#### Code for simulation #####
##############################




getSimulatedData <- function(  
  noOfObservations = 100,  
  noOfIrrelevantFeatures=1,  
  noOfRelevantFeatures=1,  
  simulatedDataLabelName = "simulatedDataLabel",
  probability_of_label_1=20,
  noiseAmplitudes="0.03, 0.4",
  correlationDirectionality="1",
  seedValue=1000) 
{
  library("data.table")  
  
  ### Set the seeding mechanism
  set.seed(seedValue)
  
  # probability_of_label_1 is either a probability (if <1) or the absolute number of 1 in the Label variable
  desiredNumberOfPositives<-probability_of_label_1
  if (probability_of_label_1<1){desiredNumberOfPositives<-round(abs(probability_of_label_1)*noOfObservations, 0)}else if (probability_of_label_1>=noOfObservations){desiredNumberOfPositives<-noOfObservations-1}
  
  
  
  
  #Transform noiseAmplitudes  and correlationDirectionality into vectors
  
  noiseAmplitude<-as.numeric(unlist(strsplit(gsub(" ", "", noiseAmplitudes),split=",")))
  correlationDirectionality<-as.numeric(unlist(strsplit(gsub(" ", "", correlationDirectionality),split=",")))
  
  
  #expand noiseAmplitude and correlationDirectionality array to a noOfRelevantFeatures long array
  noiseAmplitudes<-noiseAmplitude*rep.int(1, noOfRelevantFeatures)
  correlationDirectionality<-sign(correlationDirectionality)*rep.int(1, noOfRelevantFeatures)
  correlationDirectionality[correlationDirectionality==0]<-1
  
  #generate numeric array that will be used to generate correlated relevant features, and label column through thresholding
  labelColumnNumeric<-runif(noOfObservations)  
  
  
  trainingData<-cbind(  
    #labels column, correlated with labelColumnNumeric that will be discarded
    setnames(data.frame(labelcol=
                          as.factor(labelColumnNumeric>
                                      quantile(labelColumnNumeric,
                                               probs = c((desiredNumberOfPositives/noOfObservations))))),
             eval(simulatedDataLabelName)),
    #Relevant Features
    as.data.frame(lapply(seq(from=1,to=noOfRelevantFeatures,by=1), function(crtRelevantColumnCounter){
      crtNoise<-runif(noOfObservations)
      desiredRange<-range(crtNoise)#will reuse this range for crt feature
      
      #add noise to labelColumnNumeric
      crtRelevantColValues<-
        correlationDirectionality[crtRelevantColumnCounter]*labelColumnNumeric+  
        (noiseAmplitudes[crtRelevantColumnCounter]*(crtNoise-.5))
      
      #reuse the (0,1)-ish range of original noise for the simulated values
      oldRange<-range(crtRelevantColValues)
      crtRelevantColValues<-(crtRelevantColValues-oldRange[1])/(oldRange[2]-oldRange[1])*
        (desiredRange[2]-desiredRange[1])+desiredRange[1]
      
      #return relevant columns
      setnames(data.frame(someCol=  crtRelevantColValues), paste0("RelevCol", crtRelevantColumnCounter))
    })),
    #Irrelevant Features
    as.data.frame(  
      matrix(rnorm(noOfIrrelevantFeatures*noOfObservations),  
             noOfObservations,noOfIrrelevantFeatures))   
  )
  
  cat(summary(trainingData))
  
  trainingData = as.data.frame(trainingData)
}




source("DataSimulatorForMachineLearning.R")
seedValue<-1000
noOfObservations <-200  
noOfIrrelevantFeatures<-2  
noOfRelevantFeatures<-7  
simulatedDataLabelName <- "class_label"
desiredNumberOfPositives <-70
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


corrResults<-cor(as.matrix(as.numeric(simulatedData[,c(simulatedDataLabelName)])), as.matrix(simulatedData[,!names(simulatedData) %in% c(simulatedDataLabelName)]))
corrResults
