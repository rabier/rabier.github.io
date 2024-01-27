
calcProbSamePara = function (logLam, ti, mStart, nEnd) {
	# calculate probability from mStart genes in a parent to nEnd genes in a child in time=ti
	# here lamda = birth rate = death rate
	# formula Matthew Hahn ("Estimating the tempo and mode of gene family ...") pg1158
	
	if (mStart==0 && nEnd==0) {return(1)}
	
	alpha = (exp(logLam) *ti) / (1+ exp(logLam)*ti)
	minMN = min(mStart, nEnd)
	value = choose(mStart+nEnd-1, mStart-1) *alpha^(mStart+nEnd)	# j=0 in the formula
	Pmn = value
	ab = (1-2*alpha) / alpha^2
	
	if (minMN >= 1) {
		for (j in 1:minMN){
			value = value *(mStart-j+1)/j *(nEnd-j+1)/(mStart+nEnd-j) *ab
			Pmn = Pmn + value
		}
	}
	return (Pmn)
}

calcProb = function (logLamlogMu, ti, mStart, nEnd) {
	# calculate probability from mStart genes in a parent to nEnd genes in a child in time=ti
	# formula Forrest Crawford ("Transition probabilities for general birth-death processes ...") pg566
	
	if (mStart==0 && nEnd==0) {return(1)}
	
	if ((length(logLamlogMu) ==1) || (logLamlogMu[1] == logLamlogMu[2]) ) {
		return( calcProbSamePara(logLamlogMu[1], ti, mStart, nEnd) )
	}
	
	logLam = logLamlogMu[1]; logMu = logLamlogMu[2]
	
	lmdiff = exp(logLam) - exp(logMu)	# lam - mu
	mlratio = exp(logMu) / exp(logLam)	# mu/lam
	
	beta = (exp(lmdiff *ti) -1) / (exp(lmdiff *ti) -mlratio)
	alpha = beta * mlratio
	
	minMN = min(mStart,nEnd)
	value = choose(mStart+nEnd-1, mStart-1) *alpha^mStart *beta^nEnd 	# j=0 in the formula
	Pmn = value
	ab = (1-alpha-beta)/(alpha*beta)
	
	if (minMN >=1) {
		for (j in 1:minMN) {
			value = value *(mStart-j+1)/j *(nEnd-j+1)/(mStart+nEnd-j) *ab
			Pmn = Pmn + value
		}
	}
	return (Pmn)
}

getPt = function (logLamlogMu, ti, nPos, nFamily=NULL, isChildLeaf=F, nLeafCount=NULL, isWgdNode=F, wgdLR=NULL) {
	# nLeafCount is 1 column in the geneCountData, it is a vector of genecounts of 1 species in all families
	# nPos = posible number of genes at each internal node, usually nPos=101 (0:100) 
		
	# if isChildLeaf=T, Pt is a nPos x nFamily matrix in which each entry Pt[i,j]=prob from i-1 genes to nLeafCount[j] genes in family j 
	# if isChildLeaf=F, Pt is nPos x nPos matrix, each entry Pt[i,j] = prob from i-1 genes to j-1 genes
	# isWgdNode = whether or not it is a node before WGD event
	# wgdLR = wgd loss rate
		
	if (isWgdNode) {
		Pt = matrix (0, nrow=nPos, ncol=nPos) 
		for (i in 0:(nPos-1)) {
			N = ifelse( 2*i > (nPos-1), nPos-1, 2*i)
			for (j in i:N) {
				Pt[i+1,j+1] = choose(i, j-i) *(1-wgdLR)^(j-i) *wgdLR^(2*i -j)
				# binomial dist of getting j genes from i genes after wgd, j is between i:2i or i:(nPos-1)
			}
		}
		return (Pt)
	}	
	
	if (isChildLeaf) {
		Pt = matrix(0, nrow=nPos, ncol=nFamily)
		for (i in 1:nPos) {
			for (j in 1:nFamily) {
				Pt[i,j] = calcProb (logLamlogMu, ti, i-1, nLeafCount[j])
			}
		}
	} else {
		Pt = matrix (0, nPos, nPos)
		for (i in 1:nPos) {
			for (j in 1:nPos) {
				Pt[i,j] = calcProb (logLamlogMu, ti, i-1, j-1)
			}
		}
	}
	# calcProb = function (logLamlogMu, ti, mStart, nEnd) {}
	return (Pt)
} 


getMatAndMatDoomed = function(logLamlogMuWgdLR, nLeaf, nFamily, phyloMat, geneCountData, nPos) {
	# phyloMat is matrix with 4 columns: 1st column is Parent, 2nd is Child, 3rd is Time, 4th is Species 
	
	nNode = max(phyloMat$Parent)
	logLamlogMu = logLamlogMuWgdLR[-length(logLamlogMuWgdLR)]
	wgdLR = logLamlogMuWgdLR[length(logLamlogMuWgdLR)]
	
	Mat = array (0, c(nPos, nNode-nLeaf, nFamily))
	# Mat is a 3D matrix of size nPos x (nNode-nLeaf) x nFamily
	# Mat has nNode-nLeaf columns, each column j corresponds to an internal node j+nLeaf
	# each entry Mat[i,j,k] is prob of tree below node (j+nLeaf) in family k given starting with i genes
	

	MatDoomed = array (0, c(nNode,3))
	# MatDoomed is a 2D matrix of size nPos x (nNode-nLeaf) x nFamily
	# MatDoomed has nNode rows and 3 columns
	# each column  corresponds to the probability to be doommed when a lineage starts at
	# the Node considered
	# column 1 probability to be doomed 
	# column 2 probability to be doomed only on the left
	# column 3 probability to be doomed only on the right

	inteNode = nNode : (nLeaf+1) 
	# this is sequence of internal nodes in decreasing order
	# larger nodes are at the bottom of the tree and we calculate from these nodes up

	for (parent in inteNode) {
		#print (paste("node " , parent))
		parentIndex = which(phyloMat$Parent == parent)
		
		child1 = phyloMat$Child[ parentIndex[1]]	
		time1 = phyloMat$Time[ parentIndex[1]]
		prob1 = matrix(0, nrow=nPos, ncol=nFamily)
	
		# print( paste("currentNode is ", parent))
		# print( paste("child1 is ", child1))
		# print( paste("child2 is ", child2))
		
		if (child1 <= nLeaf) { # this child is a leaf
			species1 = phyloMat$Species [phyloMat$Child == child1]
			prob1 = getPt (logLamlogMu, time1, nPos, nFamily=nFamily, isChildLeaf=T,
					nLeafCount=geneCountData[ ,species1])	

			MatDoomed[parent,2]=prob1[2,nFamily]

							
		} else { # this child is an internal node
			if (time1 == 0) { # this parent is the node before wgd
				Pt1 = getPt (logLamlogMu, time1, nPos, isWgdNode=T, wgdLR=wgdLR)
				MatDoomed[parent,3]=1
				MatDoomed[parent,2]= wgdLR*MatDoomed[child1,1] +(1-wgdLR)*MatDoomed[child1,1]^2

			} else {
				Pt1 = getPt (logLamlogMu, time1, nPos)
				
				MatDoomed[parent,2]=Pt1[2,1]
				for (j in 2:nPos) {
				MatDoomed[parent,2]=MatDoomed[parent,2]+ Pt1[2,j] * MatDoomed[child1,1]^(j-1)
				}	

			}
			
			for (i in 1:nPos) {
				prob1[i, ] = Pt1[i, ] %*% Mat[, child1-nLeaf, ]
			}
		}
		# getPt = function (logLamlogMu, ti, nPos, nFamily=NULL, isChildLeaf=F, 
		# 		nLeafCount=NULL, isWgdNode=F, wgdLR=NULL) {}

		if (length(parentIndex) ==2) { # this parent has 2 childs
			child2 = phyloMat$Child[ parentIndex[2]]
			time2 = phyloMat$Time[ parentIndex[2]]
			prob2 = matrix(0, nrow=nPos, ncol=nFamily)

			if (child2 <= nLeaf) { # this child is a leaf
				species2 = phyloMat$Species [phyloMat$Child == child2]
				prob2 = getPt (logLamlogMu, time2, nPos, nFamily=nFamily, isChildLeaf=T,
							 nLeafCount=geneCountData[ ,species2])	
								
				MatDoomed[parent,3]=prob2[2,nFamily]


			} else { # this child is an internal node
				Pt2 = getPt (logLamlogMu, time2, nPos)
				for (i in 1:nPos) {
					prob2[i, ] = Pt2[i, ] %*% Mat[, child2-nLeaf, ]
				}

				MatDoomed[parent,3]=Pt2[2,1]
				for (j in 2:nPos) {
				MatDoomed[parent,3]=MatDoomed[parent,3]+ Pt2[2,j] * MatDoomed[child2,1]^(j-1)
				}	


			}
		}
		
		if (length(parentIndex) ==2) {
			Mat[, parent-nLeaf, ] = prob1 * prob2
			MatDoomed[parent,1]= MatDoomed[parent,2]*MatDoomed[parent,3]

		} else {

			Mat[, parent-nLeaf, ] = prob1	            
			MatDoomed[parent,3]=1
			MatDoomed[parent,1]=MatDoomed[parent,2]

		}				
	}


	return(list(first=Mat,second=MatDoomed))

}






logLikGeneCount = function (logLamlogMuWgdLR, tr, geneCountData, nPos=101, geomMean=NULL, diracMean=NULL, isBestStartNum=F, nCondData=NULL) {
	
	# geneCountData is data table which each row corresponds to a gene family and each column corresponds to a species
	# tr is a tree object from read.tree() method
	# geomMean is the mean of prior geometric dist, diracMean is the mean of prior dirac dist
	# isBestStartNum=T means choosing the best number of gene to start with (MLE approach)
	# nCondData= conditioning (1 for conditioning on observing something, 2 for conditioning on observing ar least two genes, 
	#3 for conditioning on observing at least one gene on each sides of the root )
	
	nLeaf = length(tr$tip.label)
	
	phyloMat = cbind (tr$edge, tr$edge.length)
	phyloMat = data.frame (phyloMat)
	colnames(phyloMat) = c("Parent", "Child", "Time")
	
	phyloMat$Species = rep('_', nrow(phyloMat))
	phyloMat$Species[phyloMat$Child <= nLeaf] = tr$tip.label
	
	if (!is.null(nCondData)) {
		if (nCondData==1) {
			geneCountData = rbind(geneCountData, rep(0, nLeaf))
			# add 1 row of zero's to geneCountData to calculate prob of observe nothing
		
		} else if (nCondData==2) {
			addiFamily = data.frame( rbind( rep(0, nLeaf), diag(1, nrow=nLeaf, ncol=nLeaf) ))
			names(addiFamily) = names(geneCountData)
			# add nLeaf+1 additional families to the geneCountData
			
			geneCountData = rbind(geneCountData, addiFamily)
		} else if (nCondData==3) {
	            geneCountData = rbind(geneCountData, rep(0, nLeaf))
			# add 1 row of zero's to geneCountData to calculate prob of observe nothing
		} else (stop("can't condition on this number of lineages")) }
		
	nFamily = nrow(geneCountData)
		
	myResults = getMatAndMatDoomed(logLamlogMuWgdLR, nLeaf, nFamily, phyloMat, geneCountData, nPos)
	# getMatAndMatDoomed = function(logLamlogMuWgdLR, nLeaf, nFamily, phyloMat, geneCountData, nPos) {}
	
	Mat=myResults$first
	
	MatDoomed=myResults$second

	
	if (isBestStartNum) {
		lkli = apply( Mat[ 2:nPos,1, ], 2, max)
		#print (apply( Mat[ 2:nPos,1, ], 2, function(x) { which( x==max(x) ) })) 
		
	} else if (!is.null(diracMean)) {
		lkli = Mat[diracMean+1, 1,]	
	
	} else if (!is.null(geomMean)){
		rootPriorPmf = dgeom(0:(nPos-2), ifelse(geomMean >1, 1/geomMean, .99))	
		# prior prob mass at the root
		
		lkli = rootPriorPmf %*% Mat[ 2:nPos, 1, ]
		# column 1 is the root node of the tree
		# starting from row 2 since at the root, number of starting genes should be nonzero
	}
	# lkli is a vector whose length equals the number of gene families, 
	# each entry lkli[i]=likelihood of the the tree in family i
	
	if (is.null(nCondData)) {
		logLkli = sum(log( lkli ))
	
	} else if (nCondData==1) {
		logLkli = sum(log( lkli[1:(nFamily-1)] )) - (nFamily-1)*log(1-lkli[nFamily])
		# lkli[nFamily] is prob of observing nothing
	
	} else if (nCondData==2){
		logLkli = sum(log( lkli[1:(nFamily-nLeaf-1)] )) - (nFamily-nLeaf-1)*log(1 - sum(lkli[ (nFamily-nLeaf):nFamily ]))
		# lkli[ nFamily-nLeaf : nFamily] are prob of added families that have only 0 or 1 gene in total
	
	} else if (nCondData==3){		
		#denom = 1-(MatDoomed[nLeaf+1,2])/(1-MatDoomed[nLeaf+1,2])-(MatDoomed[nLeaf+1,3])/(1-MatDoomed[nLeaf+1,3])	+(MatDoomed[nLeaf+1,1])/(1-MatDoomed[nLeaf+1,1])
	 	denom = 1-MatDoomed[nLeaf+1,2]-MatDoomed[nLeaf+1,3]+MatDoomed[nLeaf+1,1]
		logLkli = sum(log( lkli[1:(nFamily-1)] ))  - (nFamily-1)*log(denom)
	}
		
	print (-logLkli)
	return (-logLkli) 
}




MLEGeneCount<- function(tr, geneCountData, nCondData=1, geomMean=NULL, isBestStartNum=F, diracMean=NULL, nPos = 101 , logLamlogMuWgdLR=c(log(0.01), log(0.02), 0.5)){

# tr is a tree object from read.tree() method

result = optim (logLamlogMuWgdLR, logLikGeneCount, method="L-BFGS-B", lower=c(-Inf,-Inf,0), upper=c(Inf,Inf,1), tr=tr, geneCountData=geneCountData, nPos=nPos, diracMean=diracMean, geomMean=geomMean, isBestStartNum=isBestStartNum, nCondData=nCondData)
para = result$par
lambda = exp(para[1])
mu = exp(para[2])
wgdLR = para[3]
loglik=-result$value

list(birthrate=lambda, deathrate=mu, WGDlossrate=wgdLR, loglikelihood=loglik )

}




