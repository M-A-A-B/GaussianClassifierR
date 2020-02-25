trainclassifier <- function(trainX, labels) {
  
  oneClass = labels==1				#this will indicate true for all values for class = 1
  zeroClass = labels==0			#this will indicate true for all values for class = 0
  
  oneDat = trainX[oneClass,]			#this will give you the data matrix for class = 1
  zeroDat = trainX[zeroClass,]    #this will give you the data matrix for class = 0
  
  oneMean = colMeans(oneDat)	#this will give the mean of class = 1
  oneCov    = cov(oneDat)		#this will give the covariance matrix of class = 1	
  
  zeroMean = colMeans(zeroDat)	#this will give the mean of class = 0
  zeroCov    = cov(zeroDat)		#this will give the covariance matrix of class = 0	
  
  fullCov = cov(trainX)   #For case 2, common cov matrix 
  
  prior1=mean(labels)
  prior0= 1-prior1
  
  trainingParamsList<-list(
    "oneMean" = oneMean,
    "oneCov" = oneCov,
    "zeroMean" = zeroMean,
    "zeroCov" = zeroCov,
    "prior1" = prior1,
    "prior0" = prior0,
    "fullCov" = fullCov
  )
  return(trainingParamsList)
}

gauss  <- function(x, meanVe, covM=1,prior) {
  rowCount = nrow(covM)
  if(is.null(rowCount)) {
    rowCount = 0
  }
  
  if (rowCount == 2) {
    var = cbind(covM[1,1],covM[2,2])
    detFullCov = det(covM)
    invFullCov = solve(covM)
  }
  else if (rowCount == 4) {
    oneCov = covM[1:2,1:2]
    zeroCov = covM[3:4,1:2]
    
    var1 = cbind(oneCov[1,1],oneCov[2,2])
    var2 = cbind(zeroCov[1,1],zeroCov[2,2])
    
    detOneCov = det(covM[1:2,1:2])
    invOneCov = solve(covM[1:2,1:2])
    detZeroCov = det(covM[3:4,1:2])
    invZeroCov = solve(covM[3:4,1:2])
  }
  
  #Case1: For C =1
  xMinMew = x - meanVe[,1]
  if (rowCount==2) {
    sumsquarexMinMew=t(xMinMew)%*%invFullCov%*%xMinMew
  } else if(rowCount==4) {
    sumsquarexMinMew=t(xMinMew)%*%invOneCov%*%xMinMew
  } else {
    squarexMinMew = xMinMew^2
    sumsquarexMinMew = sum(squarexMinMew)
  }
  
  exponentVal = exp((-1/2)*sumsquarexMinMew)
  
  if (rowCount==2) {
    pdf1 = (1/(2*pi*(detFullCov^(1/2)))*exponentVal)
  } else if (rowCount==4) {
    pdf1 = (1/(2*pi*(detOneCov^(1/2)))*exponentVal)
  } else {
    pdf1 = ((1/2*pi)*exponentVal)
  }
  
  #Case1: For C = 0
  xMinMew = x - meanVe[,2]
  
  if (rowCount==2) {
    sumsquarexMinMew=t(xMinMew)%*%invFullCov%*%xMinMew
  } else if(rowCount==4) {
    sumsquarexMinMew=t(xMinMew)%*%invZeroCov%*%xMinMew
  } else {
    squarexMinMew = xMinMew^2
    sumsquarexMinMew = sum(squarexMinMew)
  }
  
  exponentVal = exp((-1/2)*sumsquarexMinMew)
  if (rowCount==2) {
    pdf0 = (1/(2*pi*(detFullCov^(1/2)))*exponentVal)
  } else if (rowCount==4) {
    pdf0 = (1/(2*pi*(detZeroCov^(1/2)))*exponentVal)
  } else {
    pdf0 = ((1/2*pi)*exponentVal)
  }
  
  pdf<-list(
    "pdf1"=pdf1,
    "pdf0"=pdf0
  )
  prob=testMap(pdf,prior[,1],prior[,2])
  return(prob)
}


testMap <- function(pdf,prior1,prior0) {
  mapevidence = pdf$pdf1*prior1+pdf$pdf0*prior0
  
  mapclass1 = pdf$pdf1*prior1/mapevidence
  mapclass0 = pdf$pdf0*prior0/mapevidence

  predictedLabel = 0
  if (mapclass1 > mapclass0) {
    predictedLabel = 1
  }
  
  return(c("MAPClass0"=mapclass0,"MAPClass1"= mapclass1,
           "predictedLabel"=predictedLabel))
}
driver<-function(){
  datAll = read.table("GaussTrain.txt")  #Training data
  datAll = data.matrix(datAll)	#this will convert to a matrix data structure
  labels = datAll[,ncol(datAll)]  	#this will store last column of datAll in labels
  trainX = datAll[,-ncol(datAll)]		#this will store all features except label in dat
  oneClass = labels==1				#this will indicate true for all values for class = 1
  zeroClass = labels==0			#this will indicate true for all values for class = 0
  
  oneDat = trainX[oneClass,]			#this will give you the data matrix for class = 1
  zeroDat = trainX[zeroClass,]
  par(mfcol = c(3,3))
  
  
  training = trainclassifier(trainX, labels)
  
  testDatAll = read.table("GaussTest.txt")  #Training data
  testDatAll = data.matrix(testDatAll)	#this will convert to a matrix data structure
  
  #------------------------------------------Case1 Covariance Martix = Identity------------------------------------------
  probabilitiesCase1Test <- apply(testDatAll,1,gauss,cbind(training$oneMean,training$zeroMean),1,cbind(training$prior1,training$prior0))
  probabilitiesCase1Train <- apply(trainX,1,gauss,cbind(training$oneMean,training$zeroMean),1,cbind(training$prior1,training$prior0))
  
  #Case1:Mistake Plotting
  plot(oneDat[,1],oneDat[,2],col='yellow')  	#this will plot training data points for class one in yellow color
  points(zeroDat[,1],zeroDat[,2],col='green') #this will plot training data points for class zero in green color
  probabilitiesCase1Train = t(data.matrix(probabilitiesCase1Train))
  predict = (probabilitiesCase1Train[,ncol((probabilitiesCase1Train))])
  mistake = predict!=labels
  points(trainX[mistake,1],trainX[mistake,2],col='red')
  
  #Case1:Decision Boundry Plotting
  probabilitiesCase1Test = t(data.matrix(probabilitiesCase1Test))
  testLabelsCase1 = probabilitiesCase1Test[,ncol(probabilitiesCase1Test)]  	#this will store last column of datAll in labels
  tOClassCase1 = testLabelsCase1==1				#this will indicate true for all values for class = 1
  tzClassCase1 = testLabelsCase1==0			#this will indicate true for all values for class = 0
  testOneDatCase1 = testDatAll[tOClassCase1,]			#this will give you the data matrix for class = 1
  testZeroDatCase1 = testDatAll[tzClassCase1,]
  plot(testOneDatCase1[,1],testOneDatCase1[,2],col='red')  	#this will plot data points for class one in yellow color
  points(testZeroDatCase1[,1],testZeroDatCase1[,2],col='blue') #this will plot data points for class zero in green color

  #------------------------------------------Case2 Covaraince Matrix = Common------------------------------------------
  probabilitiesCase2Test <- apply(testDatAll,1,gauss,cbind(training$oneMean,training$zeroMean),training$fullCov,cbind(training$prior1,training$prior0))
  probabilitiesCase2Train <- apply(trainX,1,gauss,cbind(training$oneMean,training$zeroMean),training$fullCov,cbind(training$prior1,training$prior0))
  
  #Case2:Mistake Plotting
  plot(oneDat[,1],oneDat[,2],col='yellow')  	#this will plot data points for class one in yellow color
  points(zeroDat[,1],zeroDat[,2],col='green') #this will plot data points for class zero in green color
  probabilitiesCase2Train = t(data.matrix(probabilitiesCase2Train))
  predict = (probabilitiesCase2Train[,ncol((probabilitiesCase2Train))])
  mistake = predict!=labels
  points(trainX[mistake,1],trainX[mistake,2],col='red')
  
  #Case2:Decision Boundry Plotting
  probabilitiesCase2Test = t(data.matrix(probabilitiesCase2Test))
  testLabelsCase2 = probabilitiesCase2Test[,ncol(probabilitiesCase2Test)]  	#this will store last column of datAll in labels
  tOClassCase2 = testLabelsCase2 == 1				#this will indicate true for all values for class = 1
  tzClassCase2 = testLabelsCase2== 0			#this will indicate true for all values for class = 0
  testOneDatCase2 = testDatAll[tOClassCase2,]			#this will give you the data matrix for class = 1
  testZeroDatCase2 = testDatAll[tzClassCase2,]
  plot(testOneDatCase2[,1],testOneDatCase2[,2],col='red')  	#this will plot data points for class one in yellow color
  points(testZeroDatCase2[,1],testZeroDatCase2[,2],col='blue') #this will plot data points for class zero in green color

  #------------------------------------------Case3 Covariance Matrix = Full(Separate for Class 0 and Class 1)------------------------------------------
  probabilitiesCase3Test <- apply(testDatAll,1,gauss,cbind(training$oneMean,training$zeroMean),rbind(training$oneCov,training$zeroCov),cbind(training$prior1,training$prior0))
  probabilitiesCase3Train <- apply(trainX,1,gauss,cbind(training$oneMean,training$zeroMean),rbind(training$oneCov,training$zeroCov),cbind(training$prior1,training$prior0))
  
  #Case3:Mistake Plotting
  plot(oneDat[,1],oneDat[,2],col='yellow')  	#this will plot data points for class one in yellow color
  points(zeroDat[,1],zeroDat[,2],col='green') #this will plot data points for class zero in green color
  probabilitiesCase3Train = t(data.matrix(probabilitiesCase3Train))
  predict = (probabilitiesCase3Train[,ncol((probabilitiesCase3Train))])
  mistake = predict!=labels
  points(trainX[mistake,1],trainX[mistake,2],col='red')
  
  #Case3:Decision Boundry Plotting
  probabilitiesCase3Test = t(data.matrix(probabilitiesCase3Test))
  testLabelsCase3 = probabilitiesCase3Test[,ncol(probabilitiesCase3Test)]  	#this will store last column of datAll in labels
  tOClassCase3 = testLabelsCase3 == 1				#this will indicate true for all values for class = 1
  tzClassCase3 = testLabelsCase3 == 0			#this will indicate true for all values for class = 0
  testOneDatCase3 = testDatAll[tOClassCase3,]			#this will give you the data matrix for class = 1
  testZeroDatCase3 = testDatAll[tzClassCase3,]
  plot(testOneDatCase3[,1],testOneDatCase3[,2],col='red')  	#this will plot data points for class one in yellow color
  points(testZeroDatCase3[,1],testZeroDatCase3[,2],col='blue') #this will plot data points for class zero in green color
  
  print("-------------------------Displaying First Five MAP Probabilities (Case 0, Case 1, Case 2)-------------------------")
  print(cbind((probabilitiesCase1Test[1:5,1:3]),(probabilitiesCase2Test[1:5,1:3]),(probabilitiesCase3Test[1:5,1:3])))
}

driver()
