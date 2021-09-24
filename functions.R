
sumStructure <- function(inputFolder, outputFolder) {
  
opar <- graphics::par(no.readonly = TRUE)
  
files = list.files(path = inputFolder, full.names = TRUE)

Kmatrix = matrix(data = NA, nrow = length(files), ncol = 1)
LnPDmatrix = matrix(data = NA, nrow = length(files), ncol = 1)
ancestryList = list()

for (f in 1:length(files)) {
  file = readLines(files[f], warn = FALSE)
  K = grep("populations assumed", file, value = TRUE)
  Kmatrix[f] = K = as.numeric(  gsub("[^0-9.-]", "", K) )

  indivs = grep("individuals", file, value = TRUE)[1]
  indivs = as.numeric(  gsub("[^0-9.-]", "", indivs) )
    
  LnPD = grep("Estimated Ln Prob of Data", file, value = TRUE)
  LnPDmatrix[f] = LnPD = as.numeric(  gsub("[^0-9.-]", "", LnPD) )

  
  indivAncestryLine = grep("Inferred ancestry of individuals", file, value = FALSE)
  indivAncestry = file[(indivAncestryLine+2):(indivAncestryLine+1+indivs)]
  indivAncestry = unlist(strsplit(indivAncestry, split = ":"))[c(FALSE, TRUE)]
  ancestryList[[f]] = matrix(as.numeric(unlist(strsplit(indivAncestry, split = " "))), ncol = K+2, byrow = T)[,- (1:2)]
}


### plot of LnP(D) for a series of structure output files

pdf(file = paste(outputFolder, "/LnP(D).pdf", sep = "")) 
plot(Kmatrix, LnPDmatrix, xlab="K", ylab="Ln P(D)", cex=1.8, frame=TRUE)
dev.off()
plot(Kmatrix, LnPDmatrix, xlab="K", ylab="Ln P(D)", cex=1.8, frame=TRUE)


# Calculation and plot similarity coefficients

simil.structure <- function(rundata) {
  
  runnb <- length(rundata)  
  indnb <- dim(rundata[[1]])[1]
  Ko <- dim(rundata[[1]])[2]
  
  simil <- matrix(NA, runnb*(runnb-1)/2, 3)
  no <- matrix(1/Ko, indnb, Ko)
  s <- 1
  frob <- function(x) {
    return(sqrt(sum(x^2))) }
  
  for (m in 1:(runnb-1)) {
    for (n in (m+1):(runnb)) {
      pair <- ordermat(rundata[[m]], rundata[[n]], indnb, Ko)
      simil[s, 1] <- m
      simil[s, 2] <- n
      simil[s, 3] <- 1 - (frob(pair[ , , 1] - pair[ , , 2])/ sqrt( (frob(pair[ , , 1] - no)) * (frob(pair[ , , 2] - no))))
      s <- s+1
    }
  }
  
  nam <- paste(outputFolder, "/simil", Ko, ".txt", sep="", collapse=NULL)
  write.table(simil, file=nam, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		
  return(simil)
}

ordermat <- function(mat1, mat2, indnr, Ka) {
  
  sr <- rank(apply(mat2, 2, sum))
  mat2 <- rbind(sr, mat2)
  mat2 <- mat2[ ,sort.list(mat2[1, ])]
  mat2 <- mat2[2:(indnr+1), ]
  mat1.sort <- matrix(NA, indnr, Ka)
  
  for (i in 0:(Ka-2)) {
    difference <- matrix(NA, indnr, (Ka-i))
    for (j in 1:(Ka-i)) {
      difference[ ,j] <- abs(mat2[ ,(Ka-i)]-mat1[ ,j])
    }
    dsr <- rank(apply(difference, 2, sum))
    mat1 <- rbind(dsr, mat1)
    mat1 <- mat1[ ,sort.list(mat1[1, ])]
    mat1.sort[ ,(Ka-i)] <- mat1[2:(indnr+1), 1]
    mat1 <- mat1[2:(indnr+1), 2:(Ka-i)]
    rm(difference, dsr)
  }
  mat1.sort[ ,1] <- mat1
  res <- array(data=NA, c(indnr, Ka, 2))
  res[ , , 1] <- mat1.sort
  res[ , , 2] <- mat2
  return(res)
}

simils <- matrix(NA, max(Kmatrix), 4)

for (K in 1:max(Kmatrix)) {

  if (K != 1) {	
    rundata = ancestryList[Kmatrix == K]
      
    similarity <- simil.structure(rundata)
    simils[K, 1] <- K
    simils[K, 2] <- length(rundata)
    simils[K, 3] <- mean(similarity[ , 3])
    simils[K, 4] <- sd(similarity[ , 3]) 
    rm(rundata)
  }
}

plotSimilCoef <- function() {
  plot(simils[ ,1], (simils[ ,3]+ simils[ ,4]), type="p", pch=24, xlab="K", ylab="similarity coefficient", ylim=c(0, 1.25), axes=FALSE, frame=TRUE)
  lines(simils[ ,1], simils[ ,3], type="b")
  points(simils[ ,1], (simils[ ,3]- simils[ ,4]), pch=25)
  axis(1, at=c(1:K))
  axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2))
}
pdf(file = paste(outputFolder, "/similarity.pdf", sep = ""))
plotSimilCoef()
dev.off()
plotSimilCoef()

simils <- rbind(c("K", "nb of runs", "mean similarity coeff.", "standard dev."), simils)
write.table(simils, file=paste(outputFolder, "/simil_coefficients.txt", sep = ""), quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		


# deltaK

averages <- matrix(NA, max(Kmatrix), 7)
for (K in 1:max(Kmatrix)){
  
  lnprobab = LnPDmatrix[Kmatrix == K]
  
  averages[K, 1] <- K
  averages[K, 2] <- length(lnprobab)
  averages[K, 3] <- mean(lnprobab)
  averages[K, 4] <- sd(lnprobab)
  
  
}


for (i in 2:max(Kmatrix)) {
  averages[i, 5]<- averages[i, 3] - averages[(i-1), 3]
}

for (i in 2:(max(Kmatrix)-1)) {
  averages[i, 6]<- abs(averages[(i+1), 5] - averages[i, 5])
}

for (i in 2:(max(Kmatrix)-1)) {
  averages[i, 7]<- averages[i, 6] / averages[i, 4]
}

averages.t <- rbind(c("K", "nb of runs", "mean similarity coeff.", "standard dev.", "L'(K)", "L''(K)", "DeltaK"), averages) 
write.table(averages.t, file=paste(outputFolder, "/DeltaK.txt",sep = ""), quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		

plotDeltaK <- function() {
  par(mfrow=c(2, 2)) 
  plot(averages[ ,1], averages[ ,3], type="p", xlab="K", ylab="Mean L(K)")
  #error.bar(averages[ ,1], averages[ ,3], averages[3]-averages[4], averages[3] + averages[4])
  plot(averages[ ,1], averages[ ,5], type="p", xlab="K", ylab="Mean L'(K)")
  plot(averages[ ,1], averages[ ,6], type="p", xlab="K", ylab="Mean L''(K)")
  plot(averages[ ,1], averages[ ,7], type="b", xlab="K", ylab="Mean DeltaK")
  
}
pdf(file = paste(outputFolder, "/deltaK.pdf", sep = ""))
plotDeltaK()
dev.off()
plotDeltaK()

graphics::par(opar)
}


















