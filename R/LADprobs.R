#' LADprobs
#' 2016 Sean Thomas

#' @param binFile a bin file in HiTC format
#' @param countFile a counts file in HiTC format
#' @param ladFile a 2 column file, first column comntains all host genome bin IDs from bin file, second column is an indicator(0/1) showing whether or not each bin overlaps a LAD.
#' @param viralChr a viral genome name: 'chrebv' for example

#' @return vector of LAD state prediction probabilities for all bins in bin file

#' @export
#' @examples
#' ladStateProbs <- LADprobs(binFile,countFile,ladFile,viralChr)

LADprobs <- function(binFile,countFile,ladFile,viralChr) {
	# first load the data
	require(HiTC)
	pcLAD <- read.delim(ladFile) 
	hic <- importC(countFile,binFile)
	
	# next build giant interaction matrix for all autosomal chromosomes #########
	chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
	
	msr <- matrix(0,length(chr),length(chr))
	msc <- matrix(0,length(chr),length(chr))
	for (i in 1:length(chr)) {
		for (j in i:length(chr)) {
			ygiName <- chr[i]
			xgiName <- chr[j]
			chrT <- hic[[paste(ygiName,xgiName,sep="")]]
			chrTm <- as.matrix(chrT@intdata)
			colnames(chrTm) <- chrT@xgi@ranges@NAMES
			rownames(chrTm) <- chrT@ygi@ranges@NAMES
		
			msr[i,j] <- nrow(chrTm)
			msr[j,i] <- ncol(chrTm)
			msc[i,j] <- ncol(chrTm)
			msc[j,i] <- nrow(chrTm)
		}
	}

	msa <- msr[,1]
	for (i in 2:length(msa)) {msa[i] <- msa[i-1]+msa[i]}
	msb <- c(1,msa[1:(length(msa)-1)]+1)
	allm <- matrix(0,sum(msr[,1]),sum(msr[,1]))
	binNames <- rep(".",nrow(allm))
	for (i in 1:length(chr)) {
		for (j in i) {
			ygiName <- chr[i]
			xgiName <- chr[j]
			chrT <- hic[[paste(ygiName,xgiName,sep="")]]
			chrTm <- as.matrix(chrT@intdata)
			colnames(chrTm) <- chrT@xgi@ranges@NAMES
			rownames(chrTm) <- chrT@ygi@ranges@NAMES
		
			binNames[msb[i]:msa[i]] <- rownames(chrTm)
		}
	}
	binState <- 1:length(binNames)*0
	for (i in 1:length(binNames)) {
		binState[i] <- 0;
		wb <- which(pcLAD[,1]==binNames[i])
		if (length(wb)==1) { binState[i] <- pcLAD[wb,2] }
	}

	for (i in 1:length(chr)) {
		for (j in i:length(chr)) {
			ygiName <- chr[i]
			xgiName <- chr[j]
			chrT <- hic[[paste(ygiName,xgiName,sep="")]]
			chrTm <- as.matrix(chrT@intdata)
			colnames(chrTm) <- chrT@xgi@ranges@NAMES
			rownames(chrTm) <- chrT@ygi@ranges@NAMES
		
			allm[msb[i]:msa[i],msb[j]:msa[j]] <- chrTm
		
			tcTm <- t(chrTm)
		 	allm[msb[j]:msa[j],msb[i]:msa[i]] <- tcTm
		}
	}
	###### end of matrix construction step
	
	# now get correlations
	wbs0 <- which(binState==0)
	wbs1 <- which(binState==1)
	s0 <- apply(allm[wbs0,],2,mean)
	s1 <- apply(allm[wbs1,],2,mean)
	cs0 <- apply(allm,1,function(a) cor(a,s0))
	cs1 <- apply(allm,1,function(a) cor(a,s1))

	### now learn ...
	g1 <- glm(binState ~ cs0 + cs1, family="binomial")
	
	### now build viral chr interaction matrix
	ygiName <- "chr1"
	xgiName <- viralChr
	nEbvBins <- 0
	chrT <- hic[[paste(ygiName,xgiName,sep="")]]
	chrTm <- as.matrix(chrT@intdata)
	colnames(chrTm) <- chrT@xgi@ranges@NAMES
	rownames(chrTm) <- chrT@ygi@ranges@NAMES
	ebvBinState <- 1:ncol(chrTm)*0
	ebvM <- matrix(0,nrow(allm),ncol(chrTm))
	colnames(ebvM) <- colnames(chrTm)

	for (i in 1:length(chr)) {
		ygiName <- chr[i]
		xgiName <- "chrebv"
		chrT <- hic[[paste(ygiName,xgiName,sep="")]]
		chrTm <- as.matrix(chrT@intdata)
		colnames(chrTm) <- chrT@xgi@ranges@NAMES
		rownames(chrTm) <- chrT@ygi@ranges@NAMES
		
		ebvM[msb[i]:msa[i],] <- chrTm
	}

	es0 <- apply(ebvM,2,function(a) cor(a,s0))
	es1 <- apply(ebvM,2,function(a) cor(a,s1))

	## next predict LAD state across viralChr using glm
	g1.pred <- predict.glm(g1,newdata=data.frame(cs0=es0,cs1=es1),type="response")
	return(g1.pred)
}

#################### example usage for a directory of result files

# setwd("data")
# baseNames <- unlist(strsplit(dir(pattern="_allBins.txt"),split="_allBins.txt"))
# viralChr <- "chrebv"
# ladFile <- "pcLAD_h1000000.txt"

# first get predictions for each dataset independently
# viralLADState <- c()
# for (i in 1:length(baseNames)) {
# 	binFile <- paste(baseNames[i],"_allBins.txt",sep="")
# 	countFile <- paste(baseNames[i],"_allCounts.txt",sep="")
# 	viralLADState[[baseNames[i]]] <- LADprobs(binFile,countFile,ladFile,viralChr)
# }

# next collect all bin names and initialize a matrix
# vBins <- c();
# for (i in 1:length(viralLADState)) {vBins <- c(vBins,names(viralLADState[[i]]))}
# vBins <- unique(vBins)
# vData <- matrix(0.5,length(baseNames),length(vBins))
# for (i in 1:length(baseNames)) {
# 	tmp <- viralLADState[[baseNames[i]]]
# 	for (j in 1:length(vBins)) {
# 		tmp2 <- which(names(tmp)==vBins[j])
# 		if (length(tmp2)==1) {
# 			vData[i,j] <- tmp[tmp2]
# 		}
# 	}
# }
# colnames(vData) <- vBins
# rownames(vData) <- baseNames


##################################################
# function make a barplot of one dataset

# generate a color scheme
# xc <- 1:1000/1000
# bc <- c(0.1000,0.1000,0.3)
# wc <- c(1,1,1)
# oc <- c(1,.5,0)
# getColGrad <- function(a,b,x) {
# 	rgb(
# 		round((b[1]-a[1])*c(1:x)/x+a[1],2),
# 		round((b[2]-a[2])*c(1:x)/x+a[2],2),
# 		round((b[3]-a[3])*c(1:x)/x+a[3],2)
# 	)
# }
# cvBW <- getColGrad(bc,wc,100)
# cvWO <- getColGrad(wc,oc,100)
# cvBWO <- c(cvBW,cvWO)

# getLADpolygon <- function(ladprobs,filename="LADpolygon.pdf",col1="midnightblue",col2="darkorange2") {
#     wp <- ladprobs; wp[which(ladprobs<=0.5)] <- 0.5; wp[which(is.na(wp))] <- 0.5
#     wn <- ladprobs; wn[which(ladprobs>=0.5)] <- 0.5; wn[which(is.na(wn))] <- 0.5

#    pdf(filename)
#    plot(c(0.5,wp,0.5),type="l",col="white",ylim=c(0,1),axes=F,ylab="LAD probability",xlab="",main=""); lines(c(0.5,wn,0.5),col="white");axis(2)
#    polygon(c(0.5,wp,0.5),col=col2)
#    polygon(c(0.5,wn,0.5),col=col1)
#    dev.off()
#}
# example of how to call the function (colors are set by default in the function but can be changed)
# e.g.: getBarPlot(vData[1,],"testBarPlot.pdf",col1="purple",col2="green")
#getLADpolygon(vData[1,],"testBarPlot.pdf")
##################################################
