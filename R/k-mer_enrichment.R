

#' Take directly determined motifs where available
#'
#' Given a data.frame of motif enrichments, finds instances where the same TF has directly determined motifs as well as indirectly determined motifs (e.g. inferred by homology), and retains only those directly determined, unless no other motif exists (see \url{http://cisbp.ccbr.utoronto.ca/index.php} for further explanation). 
#' 
#' @param motifEnrichments data.frame as made by getKMerTFEnrichment, and merged with cisbp$TFTable, such that it has column names including  c("Motif_ID", "direction","p","TF_Name","PC")) 
#' @return Returns a filtered version of motifEnrichments with a new column "indirect" containing all the Motif_IDs that represented the same TF, but were inferred indirectly.
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' tfEnrichmentsPBM = preferDirect(tfEnrichmentsPBM, cisbp$TFTable)
#' tfEnrichmentsPBM$Bon.P = tfEnrichmentsPBM$p + log(nrow(tfEnrichmentsPBM)) +  log(3000) #approximate Bonferroni MHT correction; multiply by 3000 for n_max
#' tfEnrichmentsPBM = tfEnrichmentsPBM[tfEnrichmentsPBM$Bon.P < log(0.01),] # cutoff of P<0.01
#' tfEnrichmentsPBM_NR = dropSimilarMotifs(tfEnrichmentsPBM, cisbp$similarMotifs);
#' tfEnrichmentsPBM_NR_onePerTF = bestMotifPerTF(tfEnrichmentsPBM_NR);

preferDirect = function(motifEnrichments, TFInfo){
  motifEnrichments$indirect = "";
  keepTFStatus = T
  if (!"TF_Status" %in% colnames(motifEnrichments)){
    motifEnrichments = merge(motifEnrichments, unique(TFInfo[c("Motif_ID","TF_Status", "TF_Name")]), by=c("Motif_ID","TF_Name"));
    keepTFStatus=F;
  }
  directWhereAvailable = data.frame()
  for (t in unique(as.character(motifEnrichments$TF_Name))){
    if (all(c("D","I") %in% unique(motifEnrichments$TF_Status[motifEnrichments$TF_Name==t]))){
      
      motifEnrichments$indirect[motifEnrichments$TF_Name==t & motifEnrichments$TF_Status=="D"] = paste(unique(motifEnrichments$Motif_ID[motifEnrichments$TF_Name==t & motifEnrichments$TF_Status=="I"]), collapse = ",")
      directWhereAvailable = rbind(directWhereAvailable, motifEnrichments[motifEnrichments$TF_Name==t & motifEnrichments$TF_Status=="D",])
    }else{
      directWhereAvailable = rbind(directWhereAvailable, motifEnrichments[motifEnrichments$TF_Name==t,])
    }
  }
  if (!keepTFStatus){
    directWhereAvailable$TF_Status=NULL
  }
  return(directWhereAvailable)
}

#' Take best motif per TF
#'
#' Given a table of motif enrichments in PCs (as created with getKMerTFEnrichment), selects only the best motif per TF (treating each PC and direction separately) where TFs have more than one motif. Thus, different motifs could be selected for different PCs/directions.
#' 
#' @param motifEnrichments data.frame as made by getKMerTFEnrichment, and merged with cisbp$TFTable, such that it has column names including  c("Motif_ID", "direction","p","TF_Name","PC")) 
#' @return Returns a filtered version of motifEnrichments with a new column "sameTF" containing all the Motif_IDs that represented the same TF, but had a lower enrichment
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' tfEnrichmentsPBM = preferDirect(tfEnrichmentsPBM, cisbp$TFTable)
#' tfEnrichmentsPBM$Bon.P = tfEnrichmentsPBM$p + log(nrow(tfEnrichmentsPBM)) +  log(3000) #approximate Bonferroni MHT correction; multiply by 3000 for n_max
#' tfEnrichmentsPBM = tfEnrichmentsPBM[tfEnrichmentsPBM$Bon.P < log(0.01),] # cutoff of P<0.01
#' tfEnrichmentsPBM_NR = dropSimilarMotifs(tfEnrichmentsPBM, cisbp$similarMotifs);
#' tfEnrichmentsPBM_NR_onePerTF = bestMotifPerTF(tfEnrichmentsPBM_NR);

bestMotifPerTF = function(motifEnrichments){
  motifEnrichments$sameTF = ""
  bestPerTF = data.frame()
  for (p in unique(as.character(motifEnrichments$PC))){
    for (d in unique(as.character(motifEnrichments$direction))){
      for (t in unique(as.character(motifEnrichments$TF_Name))){
        if (sum(motifEnrichments$TF==t & motifEnrichments$PC==p & motifEnrichments$direction==d)>1){ #there are more than one for this TF
          bestP = min(motifEnrichments$p[motifEnrichments$TF==t & motifEnrichments$PC==p & motifEnrichments$direction==d])
          motifEnrichments$sameTF[motifEnrichments$TF==t & motifEnrichments$PC==p & motifEnrichments$direction==d & motifEnrichments$p == bestP] = paste(unique(motifEnrichments$Motif_ID[motifEnrichments$TF==t & motifEnrichments$PC==p & motifEnrichments$direction==d & motifEnrichments$p != bestP]), collapse=",")
          bestPerTF = rbind(bestPerTF,motifEnrichments[motifEnrichments$TF==t & motifEnrichments$PC==p & motifEnrichments$direction==d & motifEnrichments$p == bestP,])
        }else{
          bestPerTF = rbind(bestPerTF,motifEnrichments[motifEnrichments$TF==t & motifEnrichments$PC==p & motifEnrichments$direction==d,])
        }
      }
    }
  }
  return (bestPerTF);
}

#' Drop any motifs that appear to be redundant
#'
#' Uses a table of correlations between motifs' k-mer z scores (as in cisbp$similarMotifs) to remove those that appear to be redundant, taking the best one in cases where motifs are very similar.  Takes the best one independently for each PC/enrichment direction, so different PCs might get different motifs as the best representative motif.
#' 
#' @param motifEnrichments data.frame as made by getKMerTFEnrichment, and merged with cisbp$TFTable, such that it has column names including  c("Motif_ID", "direction","p","TF_Name","PC")) 
#' @param motifSimilarity three column data frame (cisbp$similarMotifs), with columns c("Motif_ID1", "Motif_ID2", "R") describing pairwise motif similarity.
#' @param similarityCutoff a cutoff above which motifs are considered to be redundant.
#' @return Returns a filtered version of motifEnrichments with a new column "wereRedundant" containing all the Motif_IDs that were redundant with the remaining Motif_ID
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' tfEnrichmentsPBM = preferDirect(tfEnrichmentsPBM, cisbp$TFTable)
#' tfEnrichmentsPBM$Bon.P = tfEnrichmentsPBM$p + log(nrow(tfEnrichmentsPBM)) +  log(3000) #approximate Bonferroni MHT correction; multiply by 3000 for n_max
#' tfEnrichmentsPBM = tfEnrichmentsPBM[tfEnrichmentsPBM$Bon.P < log(0.01),] # cutoff of P<0.01
#' tfEnrichmentsPBM_NR = dropSimilarMotifs(tfEnrichmentsPBM, cisbp$similarMotifs);

dropSimilarMotifs = function(motifEnrichments, motifSimilarity, similarityCutoff=0.5){
  motifSimilarity = motifSimilarity[motifSimilarity$R>=similarityCutoff,]
  motifEnrichments = motifEnrichments[order(motifEnrichments$p),]
  motifEnrichments$wereRedundant = "";
  nonRedundantMotifs = data.frame()
  for (p in unique(as.character(motifEnrichments$PC))){
    for (d in unique(as.character(motifEnrichments$direction))){
      curSubset = motifEnrichments[motifEnrichments$PC==p & motifEnrichments$direction==d,]
      z=1;
      keep = c();
      while(z <= nrow(curSubset)){
        keep = c(keep,curSubset$Motif_ID[z])
        chuck = unique(motifSimilarity$Motif_ID1[motifSimilarity$Motif_ID2==curSubset$Motif_ID[z]])
        chuck = chuck[!chuck %in% keep]
        curSubset$wereRedundant[z] = paste(chuck, collapse = ",")
        curSubset = curSubset[!curSubset$Motif_ID %in% chuck,]
        z=z+1;
      }
      nonRedundantMotifs = rbind(nonRedundantMotifs,curSubset)
    }
  }
  return(nonRedundantMotifs)
}



#' Get TF motif enrichment for given gkm-PCs 
#'
#' Applies the minimum hypergeometric test to all PC*motif combinations, for highly-weighted k-mers and lowly-weighted k mers.
#' 
#' @param rotated A numeric matrix containing the k-mer loadings for however many PCs you want analyzed (rows are k-mers, columns are PCs).
#' @param binaryKMerMatchesToTFs A matrix of k-mers by TF motifs, where each row is a k-mer and each column is a TF motif. Each entry in the matrix is 1 if the k-mer matches the TF's motif, and 0 otherwise.
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to 3000. Should be no more than about 10\% of the number of k-mers. 
#' @param verbose Print verbose output, describing progression.
#' @return Returns a data.frame containing the following columns: Motif_ID (the motif ID from binaryKMerMatchesToTFs), PC (the PC from rotated), p (minimum hyper geometric ln(p-values)), k (the number of top k-mers that yielded maximal enrichment), log2OR (the log2 odds ratio (observed/expected) of k-mers for the point of maximal enrichment), i (a unique integer for each PC-motif-direction combination).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' head(tfEnrichmentsPBM[tfEnrichmentsPBM$PC==treatmentPCs$PC[1],], n=20) # show top 20 motifs for the first treatment-distinguishing PC

getKMerTFEnrichment = function(rotated, binaryKMerMatchesToTFs,n_max=3000, verbose=F){
  z=1;
  rotated = rotated[row.names(binaryKMerMatchesToTFs),]; #filter out and sort rows so that they're in the same order
  tfKMerEnrichments = data.frame(Motif_ID = NA, PC=NA, p = NA, k = NA, log2OR = NA, direction = NA, i = 1:(ncol(rotated)*ncol(binaryKMerMatchesToTFs)*2), stringsAsFactors = F );
  for(pci in 1:ncol(rotated)){
    curOrder = order(rotated[,pci]); #increasing order
    if (verbose){
      message(sprintf("PC = %i/%i",pci,ncol(rotated)));
    }
    for (tfi in 1:ncol(binaryKMerMatchesToTFs)){
      if (tfi %% 20==0 && verbose){
        message(sprintf("\tTF = %i/%i",tfi,ncol(binaryKMerMatchesToTFs)));
      }
      tfKMerEnrichments$Motif_ID[z:(z+1)] = colnames(binaryKMerMatchesToTFs)[tfi];
      tfKMerEnrichments$PC[z:(z+1)] = colnames(rotated)[pci];
      tfKMerEnrichments$direction[z:(z+1)] = c("low","high")
      curTest = minHG(binaryKMerMatchesToTFs[curOrder,tfi], n_max=n_max); #check for enrichment among low PC weights
      curTestRev = minHG(rev(binaryKMerMatchesToTFs[curOrder,tfi]), n_max=n_max); #check for enrichment among high PC weights
      
      #makeEnrichmentGraph(binaryKMerMatchesToTFs[curOrder,tfi], n_max=3000)
      #makeEnrichmentGraph(rev(binaryKMerMatchesToTFs[curOrder,tfi]), n_max=3000)
      tfKMerEnrichments$p[z] = curTest$lnMinP;
      tfKMerEnrichments$p[z+1] = curTestRev$lnMinP;
      tfKMerEnrichments$k[z] = curTest$k;
      tfKMerEnrichments$k[z+1] = curTestRev$k;
      tfKMerEnrichments$log2OR[z] = curTest$log2OR;
      tfKMerEnrichments$log2OR[z+1] = curTestRev$log2OR;
      z=z+2;
    }
  }
  return(tfKMerEnrichments);
}


#' Do a minimum hyper-geometric test
#'
#' Applies the minimum hypergeometric test to given data, returning the maximal enrichment achieved.
#' 
#' @param x A sorted vector of binary values, where 1 is a "hit" (e.g. cognate k-mer) and 0 is a "miss" (e.g. non-cognate k-mer).
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to testing all possible entries in the vector for maximal enrichment.
#' @return Returns a list containing the lnMinP (natural log(min(p))), k (the number of vector entries included for maximal enrichment), log2OR (log2 odds ratio (observed hits/expected hits)), and n_max (as given to method).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' testResultsLow = minHG(as.logical(cisbp$binaryPBMZScores[order(pcs$rotation[row.names(cisbp$binaryPBMZScores),treatmentPCs$PC[1]]),"M0312_1.02.PBM"])); # lowly-weighted k-mers
#' testResultsHigh = minHG(as.logical(cisbp$binaryPBMZScores[order(pcs$rotation[row.names(cisbp$binaryPBMZScores),treatmentPCs$PC[1]],decreasing = T),"M0312_1.02.PBM"])); # highly-weighted k-mers

minHG = function(x, n_max = length(x)-1){
  minP = 0;
  mink = 0;
  logOR = 0;
  n=sum(x); # number of black balls in urn (1s in x)
  m=length(x)-n; #number of white balls in urn (0s in x)
  x=1-x; # make TFBSs=0, so now white=1, black=0
  numX = cumsum(x[1:n_max])
  testI = (1:n_max)[x[1:(n_max-1)]!=x[2:(n_max)]]
  for(i in testI){
    curP = phyper(numX[i], m, n, i, log.p=T) # natural log
    #print(curP)
    if (!is.nan(curP) && curP<minP){
      mink=i;
      minP=curP;
      logOR = log2(((i-numX[i])/i)/(n/length(x)));
    }
  }
  return(list(lnMinP = minP, k=mink, log2OR = logOR, n_max=n_max));
}


#' Makes an enrichment graph for a sorted binary vector
#'
#' Applies the minimum hypergeometric test to given data and plots the enrichment of observed over expected across the data, including the ln(P) value (title), log2(observed/expected) (y-axis; e.g. observed/expected cognate k-mers), and element rank (x-axis, e.g. k-mer). The blue line indicates the point of minimum hypergeometric p-value.
#' 
#' @param x A sorted vector of binary values, where 1 is a "hit" (e.g. cognate k-mer) and 0 is a "miss" (e.g. non-cognate k-mer).
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to testing all possible entries in the vector for maximal enrichment.
#' @param sortedBy An alternate x-axis (instead of vector index) to use for the x-axis of the graph. Generally, what the vector x was sorted by.
#' @return Returns a list containing the plot (plot), the data.frame used to make the plot (rawData), and the minHG test (minHGTest).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' plot = makeEnrichmentGraph(as.logical(cisbp$binaryPBMZScores[order(pcs$rotation[row.names(cisbp$binaryPBMZScores),treatmentPCs$PC[1]],decreasing = T),"M0312_1.02.PBM"])); # highly weighted k-mers

makeEnrichmentGraph = function(x, n_max = length(x)-1, sortedBy=NULL){
  testData = data.frame(hits = x, rank =1:length(x), cumu=0, expected = (1:length(x)) * mean(x), OR=1);
  testData$cumu = cumsum(testData$hits)
  testData$cumu[1] = testData$cumu[1]+(testData$expected[1]/2); # add a pseudocount so that it doesn't start at -infinity
  testData$OR= testData$cumu/testData$expected
  testData = rbind(data.frame(hits=0,rank=0,cumu=0,expected=0,OR=1), testData);
  testData$log2OR =log2(testData$OR);
  test = minHG(x, n_max = n_max);
  if (!is.null(sortedBy)){
    testData$sortedBy = sortedBy[1:nrow(testData)];
    p = ggplot(testData, aes(x=sortedBy, y=log2OR)) + geom_line() + theme_bw() + geom_hline(yintercept=0,colour="red")+ggtitle(sprintf("ln(p) = %g",test$lnMinP))+geom_vline(xintercept=testData$sortedBy[testData$rank==test$k], colour="blue") +ylab("log2(O/E)")+xlab("variable"); print(p);
  }else{
    p = ggplot(testData, aes(x=rank, y=log2OR)) + geom_line() + theme_bw() + geom_hline(yintercept=0,colour="red")+ggtitle(sprintf("ln(p) = %g",test$lnMinP))+geom_vline(xintercept=test$k, colour="blue") +ylab("log2(O/E)")+xlab("rank"); print(p);
  }
  
  return (list(plot = p, rawData=testData, minHGTest=test));
}


#' Make a graph showing enrichment of a motif within a PC
#'
#' Applies the minimum hypergeometric test to the current data, and makes a plot where the
#' 
#' @param PC A named numeric vector representing the k-mer loadings of a PC
#' @param binaryData A named logical vector representing k-mer cognate (T) or non-cognate (F) status for a motif.
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to 3000. Should be no more than about 10\% of the number of k-mers. 
#' @return Returns a list containing the plot (plot), the data.frame used to make the plot (rawData), and the minHG test (minHGTest).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$rotation[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' p = makeEnrichmentGraphForPC(pcs$rotation[,treatmentPCs$PC[1]],cisbp$binaryPBMZScores[,head(tfEnrichmentsPBM$Motif_ID[tfEnrichmentsPBM$PC==treatmentPCs$PC[1] & tfEnrichmentsPBM$direction=="low"],n=1)]) #top motif for top treatment-distinguishing PC for lowly-weighted k-mers

makeEnrichmentGraphForPC = function(PC, binaryData, n_max=3000, decreasing=F){
  data = merge(as.data.frame(PC),as.data.frame(binaryData), by=c("row.names"))
  data = data[order(data[[2]],decreasing=decreasing),];
  p = makeEnrichmentGraph(data[[3]], n_max=n_max, sortedBy = data[[2]])
  return(p)
}



