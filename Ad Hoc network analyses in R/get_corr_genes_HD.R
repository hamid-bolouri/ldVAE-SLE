library(SummarizedExperiment)
library(Seurat)
library(coop)
library(edgeR)
library(limma)


#-------------------------------------------------------------------------------#
#								  Nakano et al									#
#-------------------------------------------------------------------------------#

# setwd("/mnt/d/SLP_3rdParty/Nakano_136SLE_89HD_2022")
# "sampleSheet.csv" # healthy vs disease status derived from SDRF
# "_count.txt" # have per cell type counts
# We normalize, calculate corrs, & make corrsList.

samplesMap <- read.csv("sampleSheet.csv")

countFiles <- dir(pattern = "_count.txt")

pops <- gsub("_count.txt", "", countFiles)
names(pops) <- countFiles

iHD <- samplesMap$diseaseStatus == "healthy control"
selSamples <- samplesMap$sampleID[iHD]


corrsNakano = list()
for (aFile in countFiles) {
	pop <- pops[aFile]
	print(paste0(pop, "                    ", Sys.time()))
	
	counts <- read.delim(aFile, header = TRUE, sep = "\t")
	
	counts <- avereps(counts[ , 3:ncol(counts)], 
					  ID = counts$Gene_name)
	
	iSel <- match(selSamples, colnames(counts))
	counts <- counts[ , iSel[!is.na(iSel)]] 
	
	
	y <- DGEList(counts)
	y <- calcNormFactors(y)
	log2TMM <- cpm(y, log = FALSE)
	log2TMM <- log2(log2TMM + 1)
	
	maxs <- apply(log2TMM, 1, max)
	log2TMM <- log2TMM[maxs > 1, ]
	
	fracZeros <- apply(log2TMM, 1, function(x) sum(x < 0.1) / length(x))
	log2TMM <- log2TMM[fracZeros < 0.5, ]
	
	corrMat <- coop::pcor(t(log2TMM), 
							use = "pairwise.complete.obs")
	
	corrMat[upper.tri(corrMat, diag = TRUE)] <- 0
	
	indx <- which(abs(corrMat) > 0.5, arr.ind = TRUE)
	r <- corrMat[indx]
	
	indx <- as.data.frame(indx)
	indx[ , 1] <- rownames(indx)
	indx[ , 2] <- colnames(corrMat)[indx[ , 2]]
	
	indx <- t(apply(indx, 1, function(x) sort(x)))
	indx <- as.data.frame(indx)
	
	indx$r <- r
	rownames(indx) <- NULL
	
	corrsNakano[[pop]] <- indx
}

saveRDS(corrsNakano, file = "corrsNakano.RDS")

rm(corrMat, indx, r, counts, y, log2TMM, maxs, fracZeros)
gc()


#-------------------------------------------------------------------------------#
#									Ota et al.									#
#-------------------------------------------------------------------------------#
# Convenience function to make SummarizedExperiment object:
mkSE <- function(filesList = files, 
				 meta = CDEs, grp, phase) {
	
	j <- meta$id[grepl(grp, meta$disease) & 
				 grepl(phase, meta$Phase)]
	
	counts = lbls = list()
	for (aFile in filesList) {
		print(aFile)
		
		counts[[aFile]] <- read.delim(aFile, header = TRUE, 
									  sep = "\t", row.names = 1)
		
		selCols <- colnames(counts[[aFile]]) %in% j
		counts[[aFile]] <- counts[[aFile]][ , selCols]
		
		lbls[[aFile]] <- rep(gsub("_count.txt", "", aFile), 
							 sum(selCols))
	}
	
	counts <- do.call(cbind, counts)
	
	counts <- counts[!duplicated(geneIDs), ]
	rownames(counts) <- geneIDs[!duplicated(geneIDs)]
	colnames(counts) <- gsub("count.txt.", "", colnames(counts))
	
	cData <- data.frame(lbls = do.call(c , lbls))
	
	
	dge <- DGEList(counts = counts, group = cData$lbls, 
				   remove.zeros = TRUE)
	dge <- calcNormFactors(dge, method = "TMM")
	CPM <- edgeR::cpm(dge, log = FALSE)
	
	logCPM <- log2(CPM + 1)
	
	rownames(cData) <- colnames(logCPM)
	
	SummarizedExperiment(assays = 
			 list(logcounts = logCPM), colData = cData)
}

CDEs <- read.delim("clinical_diagnosis_age_sex_v2.txt",
					header = TRUE, sep = "\t", 
					stringsAsFactors = FALSE)

files <- dir(pattern = "_count.txt")
# By inspection dup geneIDs are NOT real genes.
geneIDs <- read.delim(files[1], header = TRUE, 
					  sep = "\t")$Gene_name

# Using only "Phase2" HC data:
# hd1SE <- mkSE(grp = "HC", phase = 1)
hd2SE <- mkSE(grp = "HC", phase = 2)

donorIDs <- unlist(lapply(strsplit(colnames(hd2SE), "_"), 
						  function(x) rev(x)[1]))
colData(hd2SE)$donorID <- donorIDs

saveRDS(hd2SE, file = "hd2SE_phase2HD.RDS")


library(coop)
hd2SE <- readRDS("hd2SE_phase2HD.RDS")
popsOta <- unique(hd2SE$lbls)

corrsOta = list()
for (pop in popsOta) {
	print(paste0(Sys.time(), "  ", pop))
	
	tmpMat <- assays(hd2SE)$logcounts[ , hd2SE$lbls == pop]
	
	maxs <- apply(tmpMat, 1, max)
	tmpMat <- tmpMat[maxs > 1, ]
	
	fracZeros <- apply(tmpMat, 1, function(x) sum(x < 0.1) / length(x))
	tmpMat <- tmpMat[fracZeros < 0.5, ]
	
	corrMat <- coop::pcor(t(tmpMat), use = "pairwise.complete.obs")
	corrMat[upper.tri(corrMat, diag = TRUE)] <- 0
	
	indx <- which(abs(corrMat) > 0.5, arr.ind = TRUE)
	r <- corrMat[indx]
	
	indx <- as.data.frame(indx)
	indx[ , 1] <- rownames(indx)
	indx[ , 2] <- colnames(corrMat)[indx[ , 2]]
	
	indx <- t(apply(indx, 1, function(x) sort(x)))
	indx <- as.data.frame(indx)
	
	indx$r <- r
	rownames(indx) <- NULL
	
	corrsOta[[pop]] <- indx
}

rm(indx, r); gc()

# full corrsTablesOta ~= 60 GB! Only saving r > 0.5
saveRDS(corrsOta, file = "corrsOta.RDS")


#-------------------------------------------------------------------------------#
#						AIFI_GSE214546											#
#-------------------------------------------------------------------------------#

filesAgeAIFI <- dir(pattern = "_corrsPerPop.RDS")

metaData <- read.csv("metaData.CSV", header = TRUE, row.names = 1)

fileGSMs <- strtrim(filesAgeAIFI, 10)
iMeta <- match(fileGSMs, metaData["geo_accession", ])
DoBs <- as.numeric(unlist(metaData["DoB", iMeta]))
names(DoBs) <- filesAgeAIFI

filesAgeAIFI <- filesAgeAIFI[DoBs < 2000] # 32. AOK.

popsAIFI <- names(readRDS(filesAgeAIFI[1]))

sharedGenes = list()
for (pop in popsAIFI) {
	print(paste0("Pop:   ", pop, 
		  "     ", Sys.time()))	
	
	allGenes = NULL
	for (aFile in filesAgeAIFI) {
		print(aFile)
		corrsAgeAIFI <- readRDS(aFile)
		
		allGenes <- c(allGenes, rownames(corrsAgeAIFI[[pop]]))
	}
	
	freqs <- table(allGenes)
	sharedGenes[[pop]] <- names(freqs)[freqs > 
									   round(0.9 * length(filesAgeAIFI))]
	print(paste0(grep("^IRAK1$", sharedGenes[[pop]]),
		  ",          ", length(sharedGenes[[pop]])))
}

saveRDS(sharedGenes, file = "sharedGenes.RDS")

# Based on above 'print' results, run:
noIRAK1 = c("GSM6929104",
		"GSM6929107",
		"GSM6929113",
		"GSM6929114",
		"GSM6929115",
		"GSM6929116",
		"GSM6929118",
		"GSM6929119",
		"GSM7498084")


#-------------------------------------------------------------------------------#
#						AIFI_GSE214546 adults-only portion						#
#-------------------------------------------------------------------------------#

setwd("D:/AIFI_GSE214546")
filesAgeAIFI <- dir(pattern = "_corrsPerPop.RDS")

metaData <- read.csv("metaData.CSV", header = TRUE, row.names = 1)

fileGSMs <- strtrim(filesAgeAIFI, 10)
iMeta <- match(fileGSMs, metaData["geo_accession", ])

DoBs <- as.numeric(unlist(metaData["DoB", iMeta]))
names(DoBs) <- filesAgeAIFI

corrFilesAgeAIFI <- filesAgeAIFI[DoBs < 2000] # 32. AOK.

popsAIFI <- names(readRDS(corrFilesAgeAIFI[1]))

####---------- Shuffle -> get corr thresholds for empiricalP < 0.01 ---------####
library(Seurat)
library(coop)

magicDir = "D:/AIFI_GSE214546/"
magicFiles <- dir(magicDir, pattern = "_labeledMAGIC.RDS")

magicIDs <- gsub("_labeledMAGIC.RDS", "", magicFiles)
corrIDs <- gsub("_corrsPerPop.RDS", "", corrFilesAgeAIFI)

magicFiles <- magicFiles[match(corrIDs, magicIDs)]

setwd("D:/AIFI_GSE214546")
for (aFile in magicFiles) {
	print( paste0(aFile, "           ", Sys.time()) )
	
	seuratObj <- readRDS(paste0(magicDir, aFile))
	
	corrsPvalThreshL <- shuffleCorrs(seuratObj)
	
	saveRDS(corrsPvalThreshL, file = 
			paste0("corrsPvalThreshL_", 
			strtrim(aFile, 26), ".RDS"))
}


####------------------ Select sig corrs per donor per pop -------------------####
setwd("D:/AIFI_GSE214546")

corrFiles <- dir(pattern = "corrsPerPop")
threshFiles <- dir(pattern = "^corrsPvalThreshL_")

iMatch <- match(strtrim(gsub("_corrsPerPop", "", corrFiles), 26), 
				gsub("corrsPvalThreshL_|.RDS", "", threshFiles))
threshFiles <- threshFiles[ iMatch[!is.na(iMatch)] ]
corrFiles <- corrFiles[ which(!is.na(iMatch)) ]


for (i in 16:length(corrFiles)) { 
	print(paste0(i, "      ", corrFiles[i], 
				 "         ", Sys.time()))
	corrsL <- readRDS(corrFiles[i])
	corrTs <- readRDS(threshFiles[i])
	
	edgesL = list()
	for (pop in names(corrsL)) {
		print(pop)
		
		if (!is.na(match("IRAK1", rownames(corrsL[[pop]])))) {
			print(summary(corrsL[[pop]]["IRAK1", ]))
		}
		
		if (!is.null(dim(corrTs[[pop]]))) {
			threshP <- corrTs[[pop]][1, 1]
			threshN <- corrTs[[pop]][2, 1]
			
			corrsL[[pop]][corrsL[[pop]] < threshP &
						  corrsL[[pop]] > threshN] <- NA
			
			indx = which(!is.na(corrsL[[pop]]), arr.ind = TRUE)
			
			edgesTbl <- data.frame(Gene1 = rownames(indx), 
								   Gene2 = colnames(corrsL[[pop]])[indx[ , 2]])
			edgesTbl <- t(apply(edgesTbl, 1, function(x) {
								sort(x)
								}))
			edgesTbl <- as.data.frame(edgesTbl)
			edgesTbl$r <- corrsL[[pop]][indx]
			
			colnames(edgesTbl) <- c("Gene1", "Gene2", "r")
			
			edgesL[[pop]] <- edgesTbl
		} else edgesL[[pop]] <- NA
	}
	
	saveRDS(edgesL, file = paste0("signifEdges_", 
			gsub("corrsPerPop_", "", corrFiles[i])))
}


####------------------ Get edges occuring in > 1 donor ----------------------####
setwd("D:/AIFI_GSE214546")
corrFiles <- dir(pattern = "signifEdges_")

cumulativeFreqs = list()
for (aFile in corrFiles) { 
	print(paste0(gsub("_corrsPerPop.RDS", "", aFile), 
		  "         ", Sys.time()))
	
	corrsL <- readRDS(aFile)
	
	for (pop in names(corrsL)) {			
		if ( !is.null(dim(corrsL[[pop]])) ) {
			edges <- paste(corrsL[[pop]][ , 1], 
						   corrsL[[pop]][ , 2],
						   sep = ":")
			
			edgeTbl <- rep(1, length(edges))
			names(edgeTbl) <- edges
			
			if (length(cumulativeFreqs[[pop]]) == 0 |
				is.null(cumulativeFreqs[[pop]])) {
				cumulativeFreqs[[pop]] <- edgeTbl
			} else {
				iShared <- names(cumulativeFreqs[[pop]]) %in% 
						   names(edgeTbl)
				cumulativeFreqs[[pop]][iShared] <- 
						cumulativeFreqs[[pop]][iShared] + 1
				
				iExtras <- !(names(edgeTbl) %in% 
							 names(cumulativeFreqs[[pop]]))
				cumulativeFreqs[[pop]] <- c(cumulativeFreqs[[pop]],
						edgeTbl[iExtras])
			}
		}
	}
}

saveRDS(cumulativeFreqs, file = "cumulativeFreqsAifiAge.RDS")


#-------------------------------------------------------------------------------#
#	  Longitudinal AIFI_GSE190992: 4 donors, 6 samples each, over 10 weeks		#
#-------------------------------------------------------------------------------#

#### 4 donors, 6 samples each over 10-week period.
#### Filenames contain patient ID and sample collection week as '_PTIDxWx_'.

magicFiles <- dir(magicDir, pattern = "MAGIC")

for (aFile in magicFiles) {
	print(paste0("File:   ", aFile, 
		  "     ", Sys.time()))
	
	sObj <- readRDS(paste0(magicDir, aFile))
	
	allPops <- sort(unique(sObj$popLbls))

	magicRNA <- GetAssayData(sObj, slot = "data", 
							 assay = "MAGIC_SCT")
	countsRNA <- GetAssayData(sObj, slot = "counts", 
							  assay = "RNA")
	
	corrs6WkAIFI = list()
	for (pop in allPops) {
		print(pop)
		
		maxs <- apply(countsRNA[ , grepl(pop, sObj$popLbls)], 1, max)
		genes2Kp <- names(maxs)[maxs > 2]
		
		magicPopRNA <- magicRNA[rownames(magicRNA) %in% genes2Kp, 
								 grepl(pop, sObj$popLbls)]
		
		corrMatMagic <- pcor(t(magicPopRNA))
		corrMatMagic[lower.tri(corrMatMagic, diag = TRUE)] <- NA
		
		corrs6WkAIFI[[pop]] <- corrMatMagic
	}
	
	saveRDS(corrs6WkAIFI, file = paste0("corrsPerPop_", 
			strtrim(aFile, 26), ".RDS"))
}


#-------------------------------------------------------------------------------#
#				Label-shuffled corrs. -> correlation vs p-value table			#
#-------------------------------------------------------------------------------#

shuffleCorrs <- function(seuratObj) { #-----------------------------------------#
	
	magicRNA <- GetAssayData(seuratObj, slot = "data", 
							 assay = "MAGIC_SCT")
	countsRNA <- GetAssayData(seuratObj, slot = "counts", 
							  assay = "RNA")
	
	allPops <- sort(unique(seuratObj$popLbls))
	res = list()
	for (pop in allPops) { 
		print(pop)
		
		j <- seuratObj$popLbls %in% pop
		if (sum(j) > 10) { #### Skip pops smaller than 10 cells
			
			maxs <- apply(countsRNA[ , j], 1, max)
			genes2Kp <- names(maxs)[maxs > 2]
			
			magicPopRNA <- magicRNA[rownames(magicRNA) %in% genes2Kp, 
									 grepl(pop, seuratObj$popLbls)]
			
			N <- ncol(magicPopRNA)
			
			allCorrs = NULL
			for (nShuffles in 1:5) { # 5 is a speed-precision compromise.
				
				# Shuffle cells per gene:
				shuffledMAGIC <- t( apply(magicPopRNA, 1, function(x) {
									  x[sample(1:N, N)] }) )
				
				corrShuffMagic <- pcor(t(shuffledMAGIC))
				corrShuffMagic[lower.tri(corrShuffMagic, diag = TRUE)] <- NA
				
				corrs <- corrShuffMagic[!is.na(corrShuffMagic)]
				if (nShuffles == 1) {
					allCorrs <- corrs
				} else {
					allCorrs <- c(allCorrs, corrs)
				}
			}
			
			maxCorr <- max(allCorrs)
			minCorr <- min(allCorrs)
			
			if (maxCorr < 0.5) {
				res[[pop]] <- data.frame( corrThreshold = maxCorr, 
							  p = (1 / length(allCorrs)) )
			} else {
				corrThreshByP <- quantile(allCorrs, probs = 0.995)
				pByCorr0.5 <- sum(allCorrs > 0.5) / length(allCorrs)
				
				if (pByCorr0.5 < 0.01) {
					res[[pop]] <- data.frame( corrThreshold = 0.5,
											  p = pByCorr0.5 )
				} else {
					res[[pop]] <- data.frame( corrThreshold = corrThreshByP,
											 p = 0.01 )
				}
			}
			if (minCorr > -0.5) {
				res[[pop]] <- rbind(res[[pop]],
									data.frame( corrThreshold = minCorr, 
									p = (1 / length(allCorrs)) )
									)
			} else {
				corrThreshByP <- quantile(allCorrs, probs = 0.005)
				pByCorr0.5 <- sum(allCorrs < -0.5) / length(allCorrs)
				
				if (pByCorr0.5 < 0.01) {
					res[[pop]] <- rbind(res[[pop]],
										data.frame( corrThreshold = 0.5,
										p = pByCorr0.5 )
										)
				} else {
					res[[pop]] <- rbind(res[[pop]],
										data.frame( corrThreshold = corrThreshByP,
										p = 0.01 )
										)
				}
			}
		} else res[[pop]] <- NA
	}
	return(res)
}


#--------------- Calc corr thresholds for p < 0.01 -----------------------------#
library(Seurat)
library(coop)

magicDir = "/nfs/hbolouri/AIFI_GSE190992_Pure/"
magicFiles <- dir(magicDir, pattern = "MAGIC")

setwd("/mnt/d/AIFI_GSE190992")
for (aFile in magicFiles) {
	print( paste0(aFile, "           ", Sys.time()) )
	
	seuratObj <- readRDS(paste0(magicDir, aFile))
	
	corrsPvalThreshL <- shuffleCorrs(seuratObj)
	
	saveRDS(corrsPvalThreshL, file = 
			paste0("corrsPvalThreshL_", 
			strtrim(aFile, 26), ".RDS"))
}


#------------------ Get per sample SIGNIFICANT edges ---------------------------#
corrFiles <- dir(pattern = "corrsPerPop_")
threshFiles <- dir(pattern = "^corrsPvalThreshL_")

threshFiles <- threshFiles[
					match(gsub("corrsPerPop_", "", corrFiles), 
						  gsub("corrsPvalThreshL_", "", threshFiles))
						  ]

setwd("/mnt/d/AIFI_GSE190992")
for (i in 1:length(corrFiles)) {
	print(paste0(i, "      ", corrFiles[i], 
				 "         ", Sys.time()))
	corrsL <- readRDS(corrFiles[i])
	corrTs <- readRDS(threshFiles[i])
	
	edgesL = list()
	for (pop in names(corrsL)) {
		print(pop)
		
		if (!is.na(match("IRAK1", rownames(corrsL[[pop]])))) {
			print(summary(corrsL[[pop]]["IRAK1", ]))
		}
		
		threshP <- corrTs[[pop]][1, 1]
		threshN <- corrTs[[pop]][2, 1]
		
		corrsL[[pop]][corrsL[[pop]] < threshP &
					  corrsL[[pop]] > threshN] <- NA
		
		indx = which(!is.na(corrsL[[pop]]), arr.ind = TRUE)
		
		edgesTbl <- data.frame(Gene1 = rownames(indx), 
							   Gene2 = colnames(corrsL[[pop]])[indx[ , 2]])
		edgesTbl <- t(apply(edgesTbl, 1, function(x) {
							sort(x)
							}))
		edgesTbl <- as.data.frame(edgesTbl)
		edgesTbl$r <- corrsL[[pop]][indx]
		
		colnames(edgesTbl) <- c("Gene1", "Gene2", "r")
		
		edgesL[[pop]] <- edgesTbl
	}
	
	saveRDS(edgesL, file = paste0("signifEdges_", 
			gsub("corrsPerPop_", "", corrFiles[i])))
}


#-------------------------------------------------------------------------------#
#					AIFI longitudinal - get multi-donor edges					#
#-------------------------------------------------------------------------------#

# STEP 1: Get conservedCorrEdgesPerDonor
corrFiles <- dir(pattern = "signifEdges_")

sampleMap <- strsplit(corrFiles, "_|\\.")
sampleMap <- unlist(lapply(sampleMap, function(x) x[4]))
donorIDs <- strtrim(sampleMap, 5)
names(donorIDs) <- corrFiles

perDonorFiles <- split(x = donorIDs, f = donorIDs)
perDonorFiles <- lapply(perDonorFiles, function(x) names(x))

for (donor in 1:length(perDonorFiles)) {
	print(paste0(donor, "          ", Sys.time()))
	
	dFiles = perDonorFiles[[donor]]
	
	cumulativeFreqs = list()
	for (aFile in dFiles) {
		corrsL <- readRDS(aFile)
		
		for (pop in names(corrsL)) {			
			edges <- paste(corrsL[[pop]][ , 1], 
						   corrsL[[pop]][ , 2],
						   sep = ":")
			
			edgeTbl <- rep(1, length(edges))
			names(edgeTbl) <- edges
			
			if (length(cumulativeFreqs[[pop]]) == 0 |
				is.null(cumulativeFreqs[[pop]])) {
				cumulativeFreqs[[pop]] <- edgeTbl
			} else {
				iShared <- names(cumulativeFreqs[[pop]]) %in% names(edgeTbl)
				cumulativeFreqs[[pop]][iShared] <- cumulativeFreqs[[pop]][iShared] + 1
				
				iExtras <- !(names(edgeTbl) %in% names(cumulativeFreqs[[pop]]))
				cumulativeFreqs[[pop]] <- c(cumulativeFreqs[[pop]],
									 edgeTbl[iExtras])
			}
		}
	}
	
	# Remove edges seen in only 1 sample:
	for (pop in names(cumulativeFreqs)) {
		cumulativeFreqs[[pop]] <- cumulativeFreqs[[pop]][cumulativeFreqs[[pop]] > 1]
	}
	
	saveRDS(cumulativeFreqs, file = paste0("recurrentEdges_", 
			names(perDonorFiles)[donor], ".RDS"))
}


# STEP 2: get edges conserved across donors 

freqFiles <- dir(pattern = "recurrentEdges_")

cumulativeFreqs = list()
for (aFile in freqFiles) {
	edges1 <- readRDS(aFile)
	
	for (pop in names(edges1)) {
		print(pop)
		
		if (length(cumulativeFreqs[[pop]]) == 0 |
			is.null(cumulativeFreqs[[pop]])) {
			cumulativeFreqs[[pop]] <- rep(1, length(edges1[[pop]]))
			names(cumulativeFreqs[[pop]]) <- names(edges1[[pop]])
		} else {
			iShared <- names(cumulativeFreqs[[pop]]) %in% names(edges1[[pop]])
			cumulativeFreqs[[pop]][iShared] <- cumulativeFreqs[[pop]][iShared] + 1
			
			iExtras <- !(names(edgeTbl) %in% names(cumulativeFreqs[[pop]]))
			cumulativeFreqs[[pop]] <- c(cumulativeFreqs[[pop]],
								 edgeTbl[iExtras])
		}
	}
}

cumulativeFreqs <- lapply(cumulativeFreqs, function(x) x[x > 1])
saveRDS(cumulativeFreqs, file = "multiDonorRecurrentEdges.RDS")


#-------------------------------------------------------------------------------#
#	Get edges in AIFI sets AND at least 1 of Perez, Yazar, Ped, Nakano, Ota		#
#-------------------------------------------------------------------------------#

#### From above:
setwd("D:/AIFI_GSE190992")
edgesAIFI.6Wk <- readRDS("multiDonorRecurrentEdges.RDS")

setwd("D:/AIFI_GSE214546")
edgesAIFI.age <- readRDS("cumulativeFreqsAifiAge.RDS")

edgesAIFI.age <- lapply(edgesAIFI.age, function(x) x[x > 1])


# Rm dnT in 'age' & 'CD4 CTL' in '6Wk':
sharedPopsAIFI <- intersect(names(edgesAIFI.age), 
							names(edgesAIFI.6Wk)) # 18 pops

edgesAIFI.age <- edgesAIFI.age[sharedPopsAIFI]
edgesAIFI.6Wk <- edgesAIFI.6Wk[sharedPopsAIFI]


# All edges in at least 1 of the 2 AIFI datasets:
allEdgesAIFI = list()
for (pop in sharedPopsAIFI) {
	print(pop)
	allEdgesAIFI[[pop]] <- unique(c(names(edgesAIFI.age[[pop]]), 
									names(edgesAIFI.6Wk[[pop]])))
}


edgesOta <- readRDS("corrsOta.RDS") # DF with cols Edge1, Edge2, r
edgesNakano <- readRDS("corrsNakano.RDS") # As above.


#### From file '/codeBRI/SLE2/mergeCorrNets.R':
corrEdgesMyelPerezHD <- readRDS("corrEdgesGT.5.MyelPerezHD.RDS")
corrEdgesLymphPerezHD <- readRDS("corrEdgesGT.5.LymphPerezHD.RDS")

corrEdgesPerezHD <- c(corrEdgesLymphPerezHD, 
					  corrEdgesMyelPerezHD)


corrEdgesYazar <- readRDS("corrEdgesGt.5.Yazar.RDS")
corrEdgesPedHD <- readRDS("corrEdgesGT.5.PedHD.RDS")


#-------------------------------------------------------------------------------#
#				 allEdgesAIFI %in% Ota/Nakano/Perez/Yazar/Ped 					#
#-------------------------------------------------------------------------------#

popMap <- read.csv("cellTypeMapHD2.csv", header = TRUE)

iSel <- which(popMap$good == "y")

valSets = c("corrEdgesYazar", "corrEdgesPerezHD",
			"edgesNakano", "corrEdgesPedHD", "edgesOta")
names(valSets) <- colnames(popMap)[2:6]

validatedEdges = list()
for (dataset in valSets) {
	print(paste0(dataset, "      ", Sys.time()))
	
	j <- names(valSets)[valSets %in% dataset]
	
	selPops <- popMap[iSel, j]
	selPops <- selPops[!is.na(selPops)]
	
	for (pop in selPops) {
		print(pop)
		
		popAIFI <- popMap$AIFI[match(pop, popMap[ , j])]
		
		qEdges <- paste(get(dataset)[[pop]][ , 1], 
						get(dataset)[[pop]][ , 2], sep = ":")
		
		validE <- intersect(allEdgesAIFI[[popAIFI]], qEdges)
		validatedEdges[[popAIFI]] <- c(validatedEdges[[popAIFI]], 
									   validE)
	}
}

saveRDS(validatedEdges, file = 
	   "validatedEdgesAIFI.plus1Other.RDS")





















# #-------------------------------------------------------------------------------#
# #			Get edges recurrent-AIFI + recurrent-in-other-datasets				#
# #-------------------------------------------------------------------------------#

# #### Step1: Find recurrent edges among Ota, Nakano, Perez, Yazar, & ped datasets:

# recurrPubEdges = list()
# for (i in which(!is.na(popMap$AIFI))) {
	# print(paste0(popAIFI, "      ", Sys.time()))
	
	# tstPop <- popMap$corrsYazar[i]
	# if (!is.na(tstPop)) {
		# edgesYazar <- paste(corrEdgesYazar[[tstPop]][ , 1], 
							# corrEdgesYazar[[tstPop]][ , 2], 
							# sep = ":")
	# } else edgesYazar <- NULL
	
	
	# corrEdgesPerezHD
	# edgesNakano
	# edgesOta
	# corrEdgesPedHD
# }


# valEdgesV2 = list()
# for (dataset in valSets) {
	# print(paste0(dataset, "      ", Sys.time()))
	
	# j <- names(valSets)[valSets %in% dataset]
	
	# selPops <- popMap[iSel, j]
	# selPops <- selPops[!is.na(selPops)]
	
	# for (pop in selPops) {
		# print(pop)
		
		# popAIFI <- popMap$AIFI[match(pop, popMap[ , j])]
		
		# qEdges <- paste(get(dataset)[[pop]][ , 1], 
						# get(dataset)[[pop]][ , 2], sep = ":")
		
		# validE <- intersect(allEdgesAIFI[[popAIFI]], qEdges)
		# valEdgesV2[[popAIFI]] <- c(valEdgesV2[[popAIFI]], 
									   # validE)
	# }
# }

# setwd("D:/SLP_3rdParty")
# saveRDS(valEdgesV2, file = 
	   # "valEdgesV2AIFI.plus1Other.RDS")
















