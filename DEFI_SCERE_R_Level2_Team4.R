#### Defi-R-SCERE_Team4_DUO2021
#
# Level 2: The aim is to write a script that allow the clustering of gene expression in N clusters 
# with the choice of different algorithm methods (K-means or Hierarchical Clustering) 
# and the choice of the distance calculation methods (euclidean or correlation)
# with the creation of the associated graphical representation of the gene expression of the different clusters.


##### FOR THE USER :
#### Choose the Number of Clusters = N
N <- 5

#### Choose the Distance calculation method between euclidean "EUCLIDEAN" or correlation "CORREL" (copy and paste your choice)
DISTANCE <- "CORREL"

#### Choose the Clustering Algorithm method between K-means "KMEANS" or Hierarchical Clustering "HCL" (copy and paste your choice)
ALGORITHM <- "KMEANS"

#### Choose the graphical representation between gene expression profile "EXPRESSION" or heatmap "HEATMAP" or both of these representations "BOTH" (copy and paste your choice)
GRAPH <- "HEATMAP"



#_______________________________________________________________________________#
##### Import the data
## Carefully select the right directory where the data table is located with the command line setwd or the otpion in "More", and "Set As Working Directory"
expData <- read.table("cell-cycle_SCERE_DUO.txt", row.names = 1, sep = "\t", header = TRUE)

#### Check for missing values (NA) in the data (poser la question à Gaelle, comment transformer les NA en valeurs, et définir le nombre de NA par échantillon % à déterminer)
sum(is.na(expData))


#### Graphical reprensentation of all the data
plotGenes2 <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    yMax = max(expData)
    
  }
  
  # Representation of the first expression profile
  plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
       ylim = c(floor(yMin), ceiling(yMax)),
       xlab = "Time point", ylab = "Gene expression level",
       main = title)
  
  # Add expression profile for other genes
  for(i in 2:nrow(expData)){
    
    lines(1:ncol(expData), expData[i,], col = "grey")
    
    # end of for()  
  }
  
  # Average expression profile
  if(meanProfile == TRUE){
    expMean = apply(expData, 2, mean)
    lines(1:ncol(expData), expMean, col = "red", 
          lwd = 1.5, lty = "dashed")
  }
  
  # end of function plotGenes()  
}
plotGenes2(expData)


#### DISTANCE Calculation
if(DISTANCE == "EUCLIDEAN"){
  MatDis <- as.dist(dist(expData, method = "euclidean"))
}else{
  if(DISTANCE == "CORREL")
  {MatDis <- as.dist(1-cor(t(expData)))}
  else{
    print("Miss the choice of the distance method")}
}


#### Clustering 
if(ALGORITHM == "KMEANS"){
  ResultsData <- kmeans(MatDis, centers = N)
}else{
  if(ALGORITHM == "HCL"){
    ResultsData <- hclust(MatDis)
  }else{
    print("Miss the choice of the algorithm method")
  }
}


#### Graphical representation
if(GRAPH == "EXPRESSION"){
  if(ALGORITHM == "KMEANS"){par(mfrow=c(2,2))
    for(i in 1:N){cluster <- expData[which(ResultsData$cluster == i),]
    assign(paste("cluster", i),cluster, envir = parent.frame())
    T <- paste ("Cluster", i)
    plotGenes2(cluster, title = T)}
  }else{par(mfrow=c(2,2))
    for(i in 1:N){cluster <- expData[which(cutree(ResultsData, k = N) == i),]
    assign(paste("cluster", i),cluster, envir = parent.frame())
    T <- paste ("Cluster", i)
    plotGenes2(cluster, title = T)}
  }
}else if(GRAPH == "HEATMAP"){
    if(ALGORITHM == "KMEANS"){par(mfrow=c(1,1))
      for(i in 1:N){
        cluster <- expData[which(ResultsData$cluster == i),]
        assign(paste("cluster", i),cluster, envir = parent.frame())
        T <- paste ("Cluster", i)
        if(nrow(cluster) > 1){
          heatmap(as.matrix(cluster), main = T, col=hcl.colors(12, "viridis"))
        }else{print(paste("no heatmap for this cluster", i))}}
    }else{par(mfrow=c(1,1))
      for(i in 1:N){
        cluster <- expData[which(cutree(ResultsData, k = N) == i),]
        assign(paste("cluster", i),cluster, envir = parent.frame())
        T <- paste ("Cluster", i)
        if(nrow(cluster) > 1){
          heatmap(as.matrix(cluster), main = T, col=hcl.colors(12, "viridis"))
        }else{print(paste("no heatmap for this cluster", i))}
      }
    }
  }else{
    if(ALGORITHM == "KMEANS"){par(mfrow=c(1,1))
      for(i in 1:N){
        cluster <- expData[which(ResultsData$cluster == i),]
        assign(paste("cluster", i),cluster, envir = parent.frame())
        T <- paste ("Cluster", i)
        plotGenes2(cluster, title = T)
        if(nrow(cluster) > 1){
          heatmap(as.matrix(cluster), main = T, col=hcl.colors(12, "viridis"))
        }else{print(paste("no heatmap for this cluster", i))}}
    }else{par(mfrow=c(1,1))
      for(i in 1:N){
        cluster <- expData[which(cutree(ResultsData, k = N) == i),]
        assign(paste("cluster", i),cluster, envir = parent.frame())
        T <- paste ("Cluster", i)
        plotGenes2(cluster, title = T)
        if(nrow(cluster) > 1){
          heatmap(as.matrix(cluster), main = T, col=hcl.colors(12, "viridis"))
        }else{print(paste("no heatmap for this cluster", i))}
      }
    }
  }
