

na.rm = function(x){ x[!is.na(x)]}

ranksumROC = function(x,y,na.rm=T,...){
  if (na.rm){
    x=na.rm(x);
    y=na.rm(y);
  }
  curTest = wilcox.test(x,y,...);
  curTest$AUROC = (curTest$statistic/length(x))/length(y)
  return(curTest)
}

#Compares the original eigen values to permuted matrix values (done ntimes) for each PC, each time removing the last nPCs significant PCs, stopping when one is not significant. 
findSigPCs = function(myPCA, nfolds=10, alpha=0.01){
  permPCAStats = data.frame();
  stillSig=T;
  nPCs = 0
  rePCA = prcomp(t(myPCA$rotation) * matrix(rep(myPCA$sdev,times=nrow(myPCA$rotation)),nrow=ncol(myPCA$rotation), ncol=nrow(myPCA$rotation)),center = F, scale. = F);
  while(stillSig){
    curPCAStats = data.frame();
    for (i in 1:nfolds){
      remainder = ncol(myPCA$rotation) - nPCs;
      x1 = as.vector(t(myPCA$rotation[,(nPCs+1):ncol(myPCA$rotation)])* matrix(rep(myPCA$sdev[(nPCs+1):ncol(myPCA$rotation)],times=nrow(myPCA$rotation)),nrow=remainder,ncol=nrow(myPCA$rotation)))
      x1 = x1[sample(length(x1),length(x1))]
      x1 = matrix(x1,nrow=nrow(myPCA$rotation),ncol=remainder)
      curPCA = prcomp(t(x1),center = F, scale. = F);
      curPCAStats = rbind(curPCAStats, data.frame(sdev = curPCA$sdev, perm=i, eig = curPCA$sdev^2, ranks = 1:length(curPCA$sdev),nPCs=nPCs));
    }
    curPCAStats=curPCAStats[order(curPCAStats$sdev,decreasing = T)[1:nfolds],]
    message(sprintf("PC=%i; P=%g; sdev_act=%g, sdev_max_perm=%g",nPCs+1, mean(rePCA$sdev[nPCs+1] < curPCAStats$sdev),rePCA$sdev[nPCs+1], max(curPCAStats$sdev)))
    
    if (mean(rePCA$sdev[nPCs+1] > curPCAStats$sdev) < (1-alpha)){
      stillSig=F
    }else{
      nPCs = nPCs+1;
    }
    if(nPCs >= ncol(myPCA$rotation)){
      stillSig=F;
    }
    permPCAStats = rbind(permPCAStats, curPCAStats)
  }
  myPCA$permPCA = permPCAStats;
  myPCA$rePCA_sdev = rePCA$sdev;
  myPCA$nPCs = nPCs
  return(myPCA)
}





#' Inputs k-mer frequencies
#'
#' Inputs k-mer frequency data from a given set of file names. Each file is tab delimited and has the header: "kMer\\tPct_FG"
#' Followed by, for example "TTGATTG\\t0.012" ("<kMer>\\t<frequency>").
#' Assumes that all inputs have the same k-mer order.
#' 
#' @param fileNames List of file names to input.
#' @param IDs List of sample IDs, same order as the file names. Defaults to fileNames.
#' @return Returns a numeric matrix of kmer frequencies, kmers (rows) by samples (columns).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)

inputKMerFreqs = function(fileNames, IDs=fileNames){
	stopifnot(length(fileNames)==length(IDs))
	curData = read.table(fileNames[1],sep="\t",header=T,stringsAsFactors=FALSE,quote="",row.names=NULL);

	allData=matrix(nrow = nrow(curData),ncol = nrow(sampleDesc));
	row.names(allData)=curData$kMer;
	colnames(allData) = IDs;
	for (i in 1:length(fileNames)){
	  if (i %% 50==1){
	    message(sprintf("Inputting i=%i/%i",i,length(fileNames)))
	  }
	  curData = read.table(fileNames[i],sep="\t",header=T,stringsAsFactors=FALSE,quote="",row.names=NULL);
	  allData[,i] = curData$Pct_FG;
	}
	return(allData)
}


#' k-mer PCA+tSNE
#'
#' Performs PCA+tSNE of a k-mer frequency matrix. First, the k-mer matrix is scaled, then PCA (prcomp) is performed, then the number of 
#' significant PCs is determined, then tSNE is performed using the significant PCs.
#' 
#' 
#' @param x A matrix of k-mer frequencies, kmers (rows) by samples (columns), as would be produced by inputKMerFreqs.
#' @param nPCs Method for identifying significant PCs, or the number of PCs to use. Either one of ("jackstraw", or "permutation"), or an integer number of PCs. defaults to "jackstraw"
#' @param scale scale the data? defaults to True.
#' @return Returns a PCA object (as created by prcomp), also including $nPCs (the number of significant PCs), $nPCMethod (the method used to get the number of PCs), 
#' $tSNE, the tsne object (as created by tsne), and $tSNEProj, the projection of the samples onto the two tSNE components. 
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' p = ggplot(pcs$tSNEProj, aes(tSNE1, tSNE2)) + geom_point(); print(p) # plot tSNE projection
#' pcs$tSNEProj = merge(pcs$tSNEProj,sampleDesc, by.x="ID",by.y="goodID") # add sample information to tSNE projection
#' p = ggplot(pcs$tSNEProj, aes(tSNE1, tSNE2, colour=celltype)) + geom_point() + theme_classic(); print(p)

doKMerPCA = function(x, nPCs="jackstraw", scale=T){
	if(scale){
		message("Scaling data")
		sd = apply(x,1,sd)
		mean = apply(x,1,mean)
		scaledData = t(scale(t(x),scale=T, center=T));
  }else{
		scaledData=x;
	}
	rm('x');
	message("Doing PCA")
  pcs = prcomp(t(scaledData), scale.=F, center=F);
	if(scale){
		pcs$means = mean;
		pcs$sds = sd;
	}
	if (nPCs=="jackstraw"){
		message("Finding Significant PCs by Jackstraw")
		require(jackstraw)
  	nPCs = permutationPA(scaledData); 
		pcs$nPCs = nPCs$r
		pcs$nPCMethod = "jackstraw";
	}else if (nPCs=="permutation"){
		message("Finding Significant PCs by partial randomization")
		pcs = findSigPCs(pcs, nfolds=1)
		pcs$nPCMethod = "permutation";
	}else if (is.numeric(nPCs)){
		message(sprintf("Using %i PCs",nPCs))
		pcs$nPCs=nPCs;
		pcs$nPCMethod = "fixed";
	}else{
		stop("Unrecognized 'nPCs' parameter.  Must be numeric or one of 'jackstraw' or 'permutation'");
	}
  rm('scaledData')
  
	message("Doing tSNE")
  ##tSNE only top PCs, including PC1
  pcs$tSNE = tsne(pcs$x[,1:pcs$nPCs]);
  pcs$tSNEProj = data.frame(pcs$tSNE);
  names(pcs$tSNEProj) = c('tSNE1','tSNE2');
  pcs$tSNEProj$ID=row.names(pcs$x);
	return(pcs)
}

#' Find distinguishing PCs
#'
#' Identifies which PCs distinguish a feature of interest. If the feature of interest contains two classes, performs a single rank sum test for each PC.
#' If the feature of interest has >2 classes, performs a rank sum test to see if the PC separates the class from all others. 
#' In either case, returns a data frame containing the PC, class, AUROC, and rank sum P-value for how well the PC identifies the class.
#' 
#' @param rotatedData The samples' projections on PC space. Include here as many PCs as you want to inspect
#' @param classLabels Either a data frame (first column containing sample labels, second containing classes), or a list where the names of the list are the sample labels and the entries in the list are the classes.
#' @return Returns a data.frame containing the following columns: PC (which PC was used for each comparison), AUROC (the area under the ROC curve for how well each PC distinguishes each class), P (rank sum p-value for how well the PC distinguishes the class), class (the class).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#` treatmentPCs = findDistinguishingPCs(myPCA$x[,1:myPCA$nPCs], sampleDesc[c("id","treated")])

findDistinguishingPCs = function(rotatedData, classLabels){
  if (is.data.frame(classLabels)){
    classLabels2 = as.list(classLabels[[2]]);
    names(classLabels2)=classLabels[[1]];
    classLabels=classLabels2;
  }
  stopifnot(is.list(classLabels))
  classes = unique(classLabels);
  allResults = data.frame()
  for(c in classes){
    if (length(classes)>2 || classes[1]==c){
      results = data.frame(PC = colnames(rotatedData), AUROC = NA, P=NA,class=c, stringsAsFactors=F)
      for (i in 1:ncol(rotatedData)){
        curTest = ranksumROC(rotatedData[names(classLabels)[classLabels==c],i],rotatedData[names(classLabels)[classLabels!=c],i])
        results$AUROC[i] = curTest$AUROC
        results$P[i] = curTest$p.value
      }
      allResults = rbind(allResults, results)
    }
  }
  return(allResults)
}

