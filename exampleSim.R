library(simplePopSimulatR)
set.seed(1234)

nM <- rep(30, 5)
m <- sum(nM)
nQTL <- c(3, 0, 1, 1, 2)
ncM <- 100

pal <- colorRampPalette(c("#441557", "#38728B", "#88CD5A", "#BCD949", "#EFE12C"))

qtl <- sampleLoci(nMarkers = nQTL, cM = ncM, method = "random") # need omit sites
markers <- sampleLoci(nMarkers = nM, cM = ncM, exclude = qtl, method = "random") # this makes sure markers are not on qtl, could do in the other direction too. 

qtlVar <- 4
qtlEff <- list()
for(i in 1:length(qtl)) qtlEff[[i]] <- runif(length(qtl[[i]]), 1, 4) * sample(c(-1, 1), nQTL[i], replace = TRUE)

u <- unlist(qtlEff)

# make founder population (e.g. screened some germplasm for the first time)
f1 <- makeF1(markers, cM = 100)
burnIn <- 3
nFounder <- 100
founder <- makePop(f1, nFounder, type = "outcross")
while(burnIn > 0){
	founder <- mate(founder, nFounder)
	burnIn <- burnIn - 1	
} 


QTL <- genotype(founder, loci = qtl)
colMeans(QTL) / 2
hist(colMeans(do.call(cbind, getSeqMatrix(founder))) / 2)

g <- QTL %*% u

Vg <- mean(g^2) - mean(g)^2
h2 <- 0.5
Ve <- (1-h2) / h2 * Vg

y <- g + rnorm(nFounder, sd = sqrt(Ve))

selThresh <- quantile(y, 0.70) # selection intensity of 30%

selected <- founder[y >= selThresh]



ngen <- 15
popsize <- 100
selectPerc <- 0.30
VgL <- list()
gL <- list()
for(i in 1:ngen){
	# make new gen
	newPop <- mate(selected, popsize)
	#get QTL scores
	QTLnew <- genotype(newPop, loci = qtl)
	# calc true genotypic values
	gnew <- QTLnew %*% u
	#calc and save g, Vg
	gL[[i]] <- gnew
	VgL[[i]] <- mean(gnew^2) - mean(gnew)^2
	# phenotype 
	ynew <- gnew + rnorm(popsize, sd = sqrt(Ve))
	# get new selection threshold
	selThreshNew <- quantile(ynew, 1-selectPerc) # selection intensity of 30%
	# select
	selected <- newPop[ynew >= selThreshNew]
}

meanG <- lapply(gL, mean)
plot(1:ngen, meanG)

plot(1:ngen, unlist(VgL))

colMeans(genotype(newPop, qtl))/ 2
# genotype(newPop, markers)

# plot heritability
plot(1:ngen, unlist(lapply(VgL, function(x) x/(x+Ve))))