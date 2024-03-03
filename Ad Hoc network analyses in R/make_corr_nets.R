library(SummarizedExperiment)
library(SingleCellExperiment)
library(Seurat)
library(coop)


#### For Perez, Yazar, and Ped HD scRNAseq datasets:

#-------------------------------------------------------------------------------#
#		scRNA-seq -> label-shuffled corrs. -> correlation vs p-value table		#
#-------------------------------------------------------------------------------#

shuffleCorrs <- function(seuratObj, bulk = FALSE) { #---------------------------#
	if (bulk) {
		countsRNA <- GetAssayData(seuratObj, slot = "counts")
		magicRNA <- GetAssayData(seuratObj, slot = "counts")
	} else {
		magicRNA <- GetAssayData(seuratObj, slot = "data", 
								 assay = "MAGIC_SCT")
		countsRNA <- GetAssayData(seuratObj, slot = "counts", 
								  assay = "RNA")
		}
	
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
					res[[pop]] <- data.frame( corrThreshold = corrThreshByP,
											  p = pByCorr0.5 )
				} else {
					res[[pop]] <- data.frame( corrThreshold = 0.5,
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
										data.frame( corrThreshold = corrThreshByP,
										p = pByCorr0.5 )
										)
				} else {
					res[[pop]] <- rbind(res[[pop]],
										data.frame( corrThreshold = -0.5,
										p = 0.01 )
										)
				}
			}
		} else res[[pop]] <- NA
	}
	return(res)
}


#-------------------------------------------------------------------------------#
signifEdgesBulk <- function(corrTblsL, threshTblsL) { #-------------------------#
	
	edgesL = list()
	for (pop in names(corrTblsL)) {
		print(pop)
		
		if (length(dim(corrTblsL[[pop]])) != 0) {
			threshP <- threshTblsL[[pop]][1, 1]
			threshN <- threshTblsL[[pop]][2, 1]
			
			iBad <- corrTblsL[[pop]]$r < threshP &
					corrTblsL[[pop]]$r > threshN
			
			edgesTbl <- corrTblsL[[pop]][!iBad, ]
			
			# corrTbls are already sorted!
			
			colnames(edgesTbl) <- c("Gene1", "Gene2", "r")
			
			edgesL[[pop]] <- edgesTbl
		} else edgesL[[pop]] <- NA
	}
	return(edgesL)
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


#### For Nakano and Ota sorted-bulk RNA-seq datasets:

#-------------------------------------------------------------------------------#
#		For each pop, shuffled-corrs -> p < 0.01 corr thresholds -> selCorrs	#
#-------------------------------------------------------------------------------#

#---------------------------------- Ota ----------------------------------------#
corrsOta <- readRDS("corrsOta.RDS")

OtaSE <- readRDS("hd2SE_phase2HD.RDS")
popLbls <- unlist(lapply(strsplit(colnames(OtaSE), "_NW"), 
				  function(x) x[1]))

OtaSeurat <- CreateSeuratObject(counts = assays(OtaSE)$logcounts)
OtaSeurat$popLbls <- popLbls

OtaThreshTbl <- shuffleCorrs(OtaSeurat, bulk = TRUE)
saveRDS(OtaThreshTbl, file = "OtaThreshTbl.RDS")

rm(OtaSE, OtaSeurat); gc()


signifEdgesOta <- signifEdgesBulk(corrsOta, OtaThreshTbl)
saveRDS(signifEdgesOta, file = "signifEdgesOta.RDS")


#--------------------------------- Nakano --------------------------------------#
library(limma)
library(edgeR)
library(Seurat)

corrsNakano <- readRDS("corrsNakano.RDS")

samplesMap <- read.csv("sampleSheet.csv")
countFiles <- dir(pattern = "_count.txt")

pops <- gsub("_count.txt", "", countFiles)
names(pops) <- countFiles

iHD <- samplesMap$diseaseStatus == "healthy control" # 79
selSamples <- samplesMap$sampleID[iHD]


NakanoThreshTbl = list()
for (aFile in countFiles) {
	pop <- pops[aFile]
	print(paste0(pop, "                    ", Sys.time()))
	
	counts <- read.delim(aFile, header = TRUE, sep = "\t")
	
	counts <- avereps(counts[ , 3:ncol(counts)], 
					  ID = counts$Gene_name)
	
	iSel <- match(selSamples, colnames(counts))
	counts <- counts[ , iSel[!is.na(iSel)]] # 76 samples
	
	
	y <- DGEList(counts)
	y <- calcNormFactors(y)
	log2TMM <- cpm(y, log = FALSE)
	log2TMM <- log2(log2TMM + 1)
	
	maxs <- apply(log2TMM, 1, max)
	log2TMM <- log2TMM[maxs > 1, ]
	
	fracZeros <- apply(log2TMM, 1, function(x) sum(x < 0.1) / length(x))
	log2TMM <- log2TMM[fracZeros < 0.5, ]
	
	sObj <- CreateSeuratObject(counts = log2TMM)
	sObj$popLbls <- pop
	
	NakanoThreshTbl[[pop]] <- (shuffleCorrs(sObj, bulk = TRUE))[[1]]
	
}

saveRDS(NakanoThreshTbl, file =  "NakanoThreshTbl.RDS")

signifEdgesNakano <- signifEdgesBulk(corrsNakano, 
									 NakanoThreshTbl)

saveRDS(signifEdgesNakano, file = 
		"signifEdgesNakano.RDS")


#------------------------------------ Perez ------------------------------------#
corrEdgesMyelPerezHD <- readRDS("corrEdgesGT.5.MyelPerezHD.RDS")
corrEdgesLymphPerezHD <- readRDS("corrEdgesGT.5.LymphPerezHD.RDS")

corrsPerez <- c(corrEdgesLymphPerezHD, corrEdgesMyelPerezHD)

magicMatPerez <- readRDS("magicPerezHD.RDS")[[1]]
colDataPerez <- readRDS("colData_sceSLE_Perez.RDS")
colDataPerez <- colDataPerez[colDataPerez$SLE_status == "Healthy", ]

iBT <- colDataPerez$cg_cov %in% c("B", "T4", "T8")
popLbls <- as.character(colDataPerez$cg_cov)
popLbls[iBT] <- as.character(colDataPerez$ct_cov[iBT])


sObjPerez <- CreateSeuratObject(counts = t(magicMatPerez))
sObjPerez$popLbls <- popLbls


perezThreshTbl <- shuffleCorrs(sObjPerez, bulk = TRUE) 
				  # 'bulk = TRUE' to use 'counts' (== 'MAGIC')
signifEdgesPerez <- signifEdgesBulk(corrsPerez, 
									perezThreshTbl) 
					# Using corrEdgesTbl


saveRDS(perezThreshTbl, file = "perezThreshTbl.RDS")
saveRDS(signifEdgesPerez, file = "signifEdgesPerez.RDS")


#----------------------------------- Yazar -------------------------------------#
corrEdgesYazar <- readRDS("corrEdgesGt.5.Yazar.RDS")

subYazar <- readRDS("subYazar_MAGIC.RDS") # Sub-sampled to 1500 cells/pop
subYazar$popLbls <- subYazar$PopLbls # For compatibility


yazarThreshTbl <- shuffleCorrs(subYazar, bulk = FALSE) # Using 'MAGIC_SCT' 

signifEdgesYazar <- signifEdgesBulk(corrEdgesYazar, 
									yazarThreshTbl) # Using corrEdgesTbl


saveRDS(yazarThreshTbl, file = "yazarThreshTbl.RDS")
saveRDS(signifEdgesYazar, file = "signifEdgesYazar.RDS")


#----------------------------- Ped (Nehar-Belaid) ------------------------------#
corrEdgesPedHD <- readRDS("corrEdgesGT.5.PedHD.RDS")

refPedHD <- readRDS( "refPedHD16_MAGIC_v2.RDS") # sub-sampled. 393 == all pDCs!

refPedHD$popLbls <- refPedHD$PopLbls # For compatibility


pedThreshTbl <- shuffleCorrs(refPedHD, bulk = FALSE) # Using 'MAGIC_SCT' 

signifEdgesPed <- signifEdgesBulk(corrEdgesPedHD, 
								  pedThreshTbl) # Using corrEdgesTbl


saveRDS(pedThreshTbl, file = "pedThreshTbl.RDS")
saveRDS(signifEdgesPed, file = "signifEdgesPed.RDS")


#-------------------------------------------------------------------------------#
#					Make multi-assay weighted-edges table						#
#-------------------------------------------------------------------------------#

mkEdgePairs <- function(edgesTbl) {
	lapply(edgesTbl, function(x) {
		if (length(dim(x)) > 0) {
			paste(x[ , 1], x[ , 2], sep = ":")
		} else return(NA)
	})
}


edgesAIFI.6Wk <- readRDS("multiDonorRecurrentEdges.RDS")
edgesAIFI.age <- readRDS("cumulativeFreqsAifiAge.RDS")

edgesAIFI.age <- lapply(edgesAIFI.age, function(x) x[x > 1])


signifEdgesOta <- readRDS("signifEdgesOta.RDS")
signifEdgesNakano <- readRDS("signifEdgesNakano.RDS")
signifEdgesPerez <- readRDS("signifEdgesPerez.RDS")
signifEdgesYazar <- readRDS("signifEdgesYazar.RDS")
signifEdgesPed <- readRDS("signifEdgesPed.RDS")


edgePairsAIFI6Wk <- lapply(edgesAIFI.6Wk, function(x) names(x))
edgePairsAIFIAge <- lapply(edgesAIFI.age, function(x) names(x))
edgePairsOta <- mkEdgePairs(signifEdgesOta)
edgePairsNakano <- mkEdgePairs(signifEdgesNakano)
edgePairsPerez <- mkEdgePairs(signifEdgesPerez)
edgePairsYazar <- mkEdgePairs(signifEdgesYazar)
edgePairsPed <- mkEdgePairs(signifEdgesPed)

# grep("IRAK1:|:IRAK1", edgePairsAIFIAge[["CD16 Mono"]]) # OK
# lapply(edgePairsPed, head) # OK


popMap <- read.csv("cellTypeMapHD2.csv", header = TRUE)

iSel <- grep("y", popMap[ , "good"])

pubSets = c("edgePairsYazar", "edgePairsPerezHD",
			"edgePairsNakano", "edgePairsPedHD", 
			"edgePairsOta")
names(pubSets) <- colnames(popMap)[2:6]


for (i in iSel) {
	
	popYazar <- popMap[i, "corrsYazar"]
	popPerez <- popMap[i, "corrsPerez"]
	popPed <- popMap[i, "corrsPed"]
	popNakano <- popMap[i, "corrsNakano"]
	popOta <- popMap[i, "corrsOta"]
	popAIFI <- popMap[i, "AIFI"]
	
	print(paste0(popAIFI, "     ",Sys.time()))
	
	weightedEdges <- 
		table(c(edgePairsYazar[[popYazar]],
				edgePairsPerez[[popPerez]],
				edgePairsPed[[popPed]],
				edgePairsNakano[[popNakano]],
				edgePairsOta[[popOta]],
				edgePairsAIFI6Wk[[popAIFI]],
				edgePairsAIFIAge[[popAIFI]]
			  ))
	
	saveRDS(weightedEdges, file = paste0(
			popAIFI, "_weightedEdges7Datasets.RDS"))
}


#-------------------------------------------------------------------------------#
#		Keep weighted recurrent edges that include DEGs/GWAS SLE genes			#
#-------------------------------------------------------------------------------#
edgeFiles <- dir(pattern = "_weightedEdges7Datasets.RDS")

pops <- unlist(lapply(strsplit(edgeFiles, "_weightedEdges7Datasets"),
			   function(x) x[1]))
names(pops) <- edgeFiles

pops <- gsub("Naïve", "Naive", pops)

#### DEGs and DEG-edges (from 'DoD.R'):
upDEGsL <- readRDS("upDEGsL.RDS")
dnDEGsL <- readRDS("dnDEGsL.RDS")


#### GWAS genes & edges:
GWAS <- read.csv("SLE_GWAS_genes_from_PMC6829439.CSV",
				 header = FALSE, as.is = TRUE)$V1


Q1 <- paste0("^", paste(GWAS, collapse = ":|^"), ":")
Q2 <- paste0(":", paste(GWAS, collapse = "$|:"), "$")
for (pop in pops) {	
	wEdges <- readRDS(names(pops[pops %in% pop]))
	
	i1 <- grepl(Q1, names(wEdges))
	i2 <- grepl(Q2, names(wEdges))
	indx <- i1 | i2
	
	wEdges <- wEdges[indx]
	
	qGenes <- c(upDEGsL[[pop]], dnDEGsL[[pop]])
	DEG1 <- paste0("^", paste(qGenes, collapse = ":|^"), ":")
	DEG2 <- paste0(":", paste(qGenes, collapse = "$|:"), "$")
	
	i1 <- grepl(DEG1, names(wEdges))
	i2 <- grepl(DEG2, names(wEdges))
	indx <- i1 | i2
	
	wEdges <- wEdges[indx]
	
	print(paste0(pop, "    ", length(wEdges),  
				 "    ", Sys.time()))
	
	saveRDS(wEdges, file = 
			paste0(gsub(" ", "_", pop), "_wEdges_DEG_GWAS.RDS"))
}


#-------------------------------------------------------------------------------#
#						Make annotated edge tables								#
#-------------------------------------------------------------------------------#
library(msigdbr)

##### MY gene sets:
LAP	= c("MAP1LC3A", "CASP3", "CASP6", "CASP7", "CASP8", "CASP9")

dyingCells = c("HMGB1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", 
			   "TLR7", "TLR8", "TLR9", "TLR10", "IRAK1", 
			   "IRF3", "NFKB1", "IKBKG", "TRIM5")


DNases = c("DNASEIL3", "DNASE2")

autophagy = c("ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101", "BECN1", 
"ATG14", "PIK3C3", "PIK3R4TG7", "ATG10", "ATG12", "ATG16L1", "ATG5", 
"ATG4", "ATG7", "ATG3", "WIPI2", "ATG9", "ATG2A", "ATG2B", "WIPI4", "ATG9A")

antiAutophagy = c("mTOR", "MLST8", "PRAS40", "DEPTOR")

dontEatMe <- c("SIRPA", "CD47", "PECAM1", "CD300A", "CD300LF", "CD33", 
"CD22", "SIGLEC10", "PDCD1", "LILRB1")

deathReceptors <-c("AIFM1", "FAS", "TNFSF10", "TNFRSF10A", "TNFRSF10B",
				   "TNFRSF1A", "FASLG", "TNF")

#### Combined cell-death gene set:
deathSets <- unique(c(LAP, dyingCells, DNases, autophagy, antiAutophagy, 
					dontEatMe, deathReceptors)) # 67 


# https://www.biorxiv.org/content/10.1101/2023.07.03.547491v1
IFN1response <- c("IFI44L", "ISG15", "IFIT3", "XAF1", "MX1", "IFI6", 
				  "IFIT1", "TRIM22", "MX2", "RSAD2")

# imsig (PMID: 30266715) IFN1 response signature genes:
sig66IFN1 <- c("APOL1","APOL6","BATF2","BST2","C19orf66","C5orf56","CMPK2",
"DDX58","DDX60","DHX58","DTX3L","EPSTI1","FBXO6","GBP1","GBP4","HELZ2","HERC5",
"HERC6","HSH2D","IFI16","IFI35","IFI44","IFI44L","IFI6","IFIH1","IFIT1","IFIT2",
"IFIT3","IFIT5","IFITM1","IRF7","IRF9","ISG15","LAMP3","LAP3","MX1","MX2","OAS2",
"OAS3","OASL","PARP10","PARP12","PARP14","PARP9","PHF11","PML","PSMB9","RNF213",
"RSAD2","RTP4","SAMD9","SAMD9L","SHISA5","SIGLEC1","SP110","STAT1","STAT2","TAP1",
"TRAFD1","TRIM21","TRIM22","TRIM5","UBE2L6","USP18","XAF1","ZNFX1")

# 5-gene IFN1 sig from PMC5751729 
sig5IFN1 <- c("MX1", "IFI6", "OAS1", "ISG15", "IFI44L")

#### IFN1 signature
sig16 <- c("DDX58","DHX58","EIF2AK2", "CCL5",
"IFIH1","IRF7","ISG20","MX1","OAS1","PLSCR1","RSAD2",
"IFI6","IFI16","STAT1","TLR3","XAF1")

#### Combined IFN1 signature genes:
ifn1Sets <- unique(c(IFN1response, sig66IFN1, 
					 sig5IFN1, sig16)) # 72

#### immuneSigDB genes:
immSigDB<- msigdbr(species = "Homo sapiens", 
				   category = "C7", 
				   subcategory = "IMMUNESIGDB") %>% 
				   dplyr::select(gs_name, gene_symbol)



immSigDB <- immSigDB[grep("GSE10325", immSigDB$gs_name), ]

Q <- "_MYELOID_VS_LUPUS_MYELOID_|_BCELL_VS_LUPUS_BCELL_|CD4_TCELL__VS_LUPUS_CD4_TCELL_"
immSigDB <- immSigDB[grep(Q, immSigDB$gs_name), ] # 726 genes


#### GWAS genes & edges:
setwd("D:/SLP_3rdParty/")
GWAS <- read.csv("SLE_GWAS_genes_from_PMC6829439.CSV",
				 header = FALSE, as.is = TRUE)$V1


#### All genes of interest:
myGeneSets <- unique(c(deathSets, ifn1Sets, 
					 immSigDB$gene_symbol)) # 819 genes
# myGeneSets <- unique(c(myGeneSets, GWAS))	# 940 genes


#----------------------- Filter wEdges & make edgeTbls -------------------------#
upDEGsL <- readRDS("upDEGsL.RDS")
dnDEGsL <- readRDS("dnDEGsL.RDS")


edgeFiles <- dir(pattern = "_wEdges_DEG_GWAS.RDS")

pops <- lapply(strsplit(edgeFiles, "_wEdges_DEG_GWAS"),
			   function(x) gsub("_", " ", x[1]))
pops <- gsub(" CD56", "_CD56", unlist(pops))
pops <- gsub("Naïve", "Naive", pops)

names(edgeFiles) <- pops

iKp <- pops %in% names(upDEGsL)
pops <- pops[iKp]
edgeFiles <- edgeFiles[iKp]


for (iFile in 1:length(edgeFiles)) {
	aFile <- edgeFiles[iFile]
	wEdge <- readRDS(aFile)
	
	if (!is.null(dim(wEdge))) {
		edgeTbl <- do.call(rbind, strsplit(names(wEdge), ":"))
		
		edgeTbl <- as.data.frame(edgeTbl)
		colnames(edgeTbl) <-c("Gene1", "Gene2")
		
		edgeTbl$nPubs <- as.vector(wEdge)
		
		# Confirmed that every edge includes a GWAS gene
		edgeTbl$GWAS <- 0
		edgeTbl$GWAS[edgeTbl$Gene1 %in% GWAS] <- 1
		edgeTbl$GWAS[edgeTbl$Gene2 %in% GWAS] <- 2
		
		i1 <- edgeTbl$Gene1 %in% myGeneSets
		i2 <- edgeTbl$Gene2 %in% myGeneSets
		
		edgeTbl <- edgeTbl[i1 | i2, ]
		
		# Mark up DEGs & down DEGs:
		edgeTbl$upInSLE <- 0
		i <- edgeTbl$Gene1 %in% upDEGsL[[names(aFile)]]
		edgeTbl$upInSLE[i] <- 1
		
		i <- edgeTbl$Gene2 %in% upDEGsL[[names(aFile)]]
		edgeTbl$upInSLE[i] <- 2
		
		
		edgeTbl$dnInSLE <- 0
		i <- edgeTbl$Gene1 %in% dnDEGsL[[names(aFile)]]
		edgeTbl$dnInSLE[i] <- 1
		
		i <- edgeTbl$Gene2 %in% dnDEGsL[[names(aFile)]]
		edgeTbl$dnInSLE[i] <- 2
		
		# Remove edges where only the GWAS gene is DE:
		i <- ( edgeTbl$GWAS != edgeTbl$upInSLE &
			   edgeTbl$GWAS != edgeTbl$dnInSLE )
		
		edgeTbl <- edgeTbl[i, ]
		
		# Add gene set annotation
		edgeTbl$g1GeneSet <- "-"
		edgeTbl$g1GeneSet[edgeTbl$Gene1 %in% 
						  immSigDB$gene_symbol] <- "immSigDB"
		edgeTbl$g1GeneSet[edgeTbl$Gene1 %in% 
						  ifn1Sets] <- "IFN1"
		edgeTbl$g1GeneSet[edgeTbl$Gene1 %in% 
						  deathSets] <- "death"
		
		edgeTbl$g2GeneSet <- "-"
		edgeTbl$g2GeneSet[edgeTbl$Gene2 %in% 
						  immSigDB$gene_symbol] <- "immSigDB"
		edgeTbl$g2GeneSet[edgeTbl$Gene2 %in% 
						  ifn1Sets] <- "IFN1"
		edgeTbl$g2GeneSet[edgeTbl$Gene2 %in% 
						  deathSets] <- "death"
		
		# There should be NO edges with un-annotated genes:
		edgeTbl <- edgeTbl[ edgeTbl$g1GeneSet != "-" |
							edgeTbl$g2GeneSet != "-",  ]
		
		saveRDS(edgeTbl, file = 
				paste0(names(aFile), "_selEdgTbl.RDS"))
	}
}


#-------------------------------------------------------------------------------#
#								Plot graphs										#
#-------------------------------------------------------------------------------#
library(igraph)

setwd("D:/SLP_3rdParty/")
eTblFiles <- dir(pattern = "_selEdgTbl.RDS")

pdf("selEdgePerPopGraphs_kkLayout_v2Colors.pdf", height = 25, width = 25)
for (aFile in eTblFiles) {
	
	edgeTbl <- readRDS(aFile)
	pop = unlist(lapply(strsplit(aFile, "_"), function(x) x[1]))
	
	G <- graph_from_edgelist(as.matrix(edgeTbl[ , 1:2]), 
							 directed = FALSE)
	
	# E(G)$color <- "#BBBBBB33" # "#00990011"
	
	iUpInSLE <- names(V(G)) %in% upDEGsL[[pop]]
	iDnInSLE <- names(V(G)) %in% dnDEGsL[[pop]]
	
	V(G)$color <- rgb(0,0,0,0.05)
	V(G)$color[iUpInSLE] <- rgb(1,0.5,0,0.2)
	V(G)$color[iDnInSLE] <- rgb(0,0.5,1,0.2)
	
	iGWAS <-  names(V(G)) %in% GWAS
	
	V(G)$size <- 3 
	V(G)$size[iGWAS] <- 10
	
	V(G)$label.color <- "black"
	V(G)$label <- names(V(G))
	V(G)$label[!iGWAS] <- NA
	
	iIFIH1 <- names(V(G)) %in% "IFIH1"
	V(G)$label[!iIFIH1] <- NA
	
	iDeath <- names(V(G)) %in% deathSets
	iIFN1 <- names(V(G)) %in% ifn1Sets
	
	V(G)$frame.width <- 2
	V(G)$frame.color <- NA
	V(G)$frame.color[iDeath] <- "black"
	V(G)$frame.color[iIFN1] <- "magenta"
	
	
	myLayout <- layout_with_kk(G)
	# myLayout <- layout_in_circle(G)
	# myLayout <- layout.fruchterman.reingold(G, niter = 1E2)
	# myLayout <- layout_with_lgl(G)	
	# myLayout <- layout_with_mds(G)
	# myLayout = layout_with_graphopt(G, charge = 0.5)
	# myLayout <- layout_as_tree(G)
	# myLayout <- layout_as_star(G, center = V(G)[grep("IRAK1", names(V(G)))])
	
	# angle = (-1 * pi / 10)
	# RotMat = matrix(c(cos(angle), sin(angle), 
					 # -sin(angle), cos(angle)), 
					 # ncol = 2)
	# myLayout <- myLayout %*% RotMat
	
	# myLayout <- -1 * myLayout[ , 2:1]
	
	set.seed(42)
	plot(G, vertex.label.family = "Helvetica",
		 layout = myLayout)
	title(pop, cex.main = 2.5, col.main = "black")
}
dev.off()


#-------------------------------------------------------------------------------#
#		Are IRAK1 & UBE2L3 differently connected in CD16+ Monos vs pDCs?		#
#-------------------------------------------------------------------------------#
setwd("D:/SLP_3rdParty/")
eTblFiles <- dir(pattern = "_selEdgTbl.RDS")

# For pDC:
aFile <- eTblFiles[grep("pDC", eTblFiles)]

edgeTbl <- readRDS(aFile)

irak1Edges <- edgeTbl[edgeTbl$Gene1 %in% "IRAK1", ]
irak1Edges <- rbind(irak1Edges, 
			  edgeTbl[edgeTbl$Gene2 %in% "IRAK1", ]) # 356

ubeEdges <- edgeTbl[edgeTbl$Gene1 %in% "UBE2L3", ]
ubeEdges <- rbind(ubeEdges, 
			edgeTbl[edgeTbl$Gene2 %in% "UBE2L3", ]) # 391


length(intersect(c(irak1Edges[ , 1], irak1Edges[ , 2]),
				 c(ubeEdges[ , 1], ubeEdges[ , 2]))) # 355

# 355/391 == 91% shared interactors!

atgEdges <- edgeTbl[edgeTbl$Gene1 %in% "ATG5", ]
atgEdges <- rbind(atgEdges, 
			edgeTbl[edgeTbl$Gene2 %in% "ATG5", ]) # 314


length(intersect(c(irak1Edges[ , 1], irak1Edges[ , 2]),
				 c(atgEdges[ , 1], atgEdges[ , 2]))) # 304

# 304/314 == 97% shared interactions

edgesPDC <- unique(c(unlist(irak1Edges[ , 1:2]), 
					 unlist(ubeEdges[ , 1:2]),
					 unlist(atgEdges[ , 1:2]))) # 395


# For CD14 Monos:
aFile <- eTblFiles[grep("CD14", eTblFiles)]

edgeTbl <- readRDS(aFile)

irak1Edges <- edgeTbl[edgeTbl$Gene1 %in% "IRAK1", ]
irak1Edges <- rbind(irak1Edges, 
			  edgeTbl[edgeTbl$Gene2 %in% "IRAK1", ]) # 856

ubeEdges <- edgeTbl[edgeTbl$Gene1 %in% "UBE2L3", ]
ubeEdges <- rbind(ubeEdges, 
			edgeTbl[edgeTbl$Gene2 %in% "UBE2L3", ]) # 845


length(intersect(c(irak1Edges[ , 1], irak1Edges[ , 2]),
				 c(ubeEdges[ , 1], ubeEdges[ , 2]))) # 833 

# 833/845 == 99% shared interactors!

atgEdges <- edgeTbl[edgeTbl$Gene1 %in% "ATG5", ]
atgEdges <- rbind(atgEdges, 
			edgeTbl[edgeTbl$Gene2 %in% "ATG5", ]) # 854


length(intersect(c(irak1Edges[ , 1], irak1Edges[ , 2]),
				 c(atgEdges[ , 1], atgEdges[ , 2]))) # 849

# 849/854 == 99.4% shared interactions

edgesCD14 <- unique(c(unlist(irak1Edges[ , 1:2]), 
					 unlist(ubeEdges[ , 1:2]),
					 unlist(atgEdges[ , 1:2]))) # 875


#### IFIH1 -partner genes:
edgeTbl <- readRDS(eTblFiles[grep("CD14", eTblFiles)])
iIFIH1 <- grepl("IFIH1", edgeTbl[ , 1]) |
		  grepl("IFIH1", edgeTbl[ , 2])
cd14IFIH1 <- setdiff(sort(unique(unlist(edgeTbl[iIFIH1, 1:2]))),
					 "IFIH1") # 338

edgeTbl <- readRDS(eTblFiles[grep("CD16", eTblFiles)])
iIFIH1 <- grepl("IFIH1", edgeTbl[ , 1]) |
		  grepl("IFIH1", edgeTbl[ , 2])
cd16IFIH1 <- setdiff(sort(unique(unlist(edgeTbl[iIFIH1, 1:2]))),
					 "IFIH1") # 261

edgeTbl <- readRDS(eTblFiles[grep("pDC", eTblFiles)])
iIFIH1 <- grepl("IFIH1", edgeTbl[ , 1]) |
		  grepl("IFIH1", edgeTbl[ , 2])
pdcIFIH1 <- setdiff(sort(unique(unlist(edgeTbl[iIFIH1, 1:2]))),
					 "IFIH1") # 112

library(nVennR)
# https://cran.r-hub.io/web/packages/nVennR/vignettes/nVennR.html
venn <- plotVenn(list(CD14 = cd14IFIH1, 
					  CD16 = cd16IFIH1, 
					  pDC = pdcIFIH1),
					  outFile = "IFIH1_partnerGenesVenn.svg",
					  setColors =  c('coral', 'gold', 'skyblue'), 
					  labelRegions = FALSE, fontScale = 3, 
					  opacity = 0.1,  borderWidth = 3)


#### Plot overlaps:
# For CD16 Monos:
aFile <- eTblFiles[grep("CD16", eTblFiles)]
# aFile <- eTblFiles[grep("CD14", eTblFiles)]
# aFile <- eTblFiles[grep("pDC", eTblFiles)]

edgeTbl <- readRDS(aFile)

irak1Edges <- edgeTbl[edgeTbl$Gene1 %in% "IRAK1", ]
irak1Edges <- rbind(irak1Edges, 
			  edgeTbl[edgeTbl$Gene2 %in% "IRAK1", ]) # 524
irak1Edges <- unique(unlist(irak1Edges[ , 1:2]))

ubeEdges <- edgeTbl[edgeTbl$Gene1 %in% "UBE2L3", ]
ubeEdges <- rbind(ubeEdges, 
			edgeTbl[edgeTbl$Gene2 %in% "UBE2L3", ]) # 481
ubeEdges <- unique(unlist(ubeEdges[ , 1:2]))

atgEdges <- edgeTbl[edgeTbl$Gene1 %in% "ATG5", ]
atgEdges <- rbind(atgEdges, 
			edgeTbl[edgeTbl$Gene2 %in% "ATG5", ]) # 516
atgEdges <- unique(unlist(atgEdges[ , 1:2]))

edgesCD16 <- unique(c(irak1Edges, ubeEdges, atgEdges)) # 559

library(nVennR)
# https://cran.r-hub.io/web/packages/nVennR/vignettes/nVennR.html
venn <- plotVenn(list(ATG5 = atgEdges, 
					  IRAK1 = irak1Edges, 
					  UBE2L3 = ubeEdges),
					  outFile = "pDC_IRAK1_module_Venn.svg",
					  setColors = c('coral', 'gold', 'skyblue'), 
					  labelRegions = FALSE, fontScale = 3, 
					  opacity = 0.1,  borderWidth = 3)


#### Overlap across cell types:
venn <- plotVenn(list('CD14' = edgesCD14, 
					  'CD16' = edgesCD16, 
					  'pDC' = edgesPDC), nCycles = 1E4,
					  outFile = "tmp.svg", # "IRAK1_moduleVenn.svg",
					  setColors = c('coral', 'gold', 'skyblue'), 
					  labelRegions = FALSE, fontScale = 3, 
					  opacity = 0.1,  borderWidth = 3)

