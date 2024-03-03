# Data processing example code using Yazar data for illustration.

library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)


# Colors
qualitativeColorPalette = 
	brewer.pal.info[brewer.pal.info$category == 'qual',]
myColors <- unlist(mapply(brewer.pal, 
				   qualitativeColorPalette$maxcolors, 
				   rownames(qualitativeColorPalette)))


#### Yazar
expYazar <- readRDS("Yazar_Seurat.rds") # 36,571 x 1,248,980
table(expYazar$development_stage)
# 78 donors, ages 19-97, PBMCs


#-------------------------------------------------------------------------------#
#					Filter out questionable "cells"								#
#-------------------------------------------------------------------------------#

selPops = c("B intermediate", "B naive", "B memory", 
			"CD4 Naive","CD4 TCM","CD4 TEM","Treg",
			"CD8 Naive","CD8 TCM","CD8 TEM",
			"gdT","MAIT","NK","NK_CD56bright",
			"CD14 Mono", "CD16 Mono", "cDC","pDC")


iKeep <- 	expYazar$predicted.celltype.l2.score > 0.5 & 
			expYazar$nFeature_RNA < 2000 &
			expYazar$nFeature_RNA > 500 &
			expYazar$nCount_RNA < 10E3 &
			expYazar$nCount_RNA > 1E3 &
			expYazar$percent.mt < 7

expYazar <- expYazar[ , iKeep] # 1,063,635 cells


expYazar$PopLbls <- expYazar$predicted.celltype.l2 
expYazar$PopLbls <- gsub("^cDC1$|^cDC2$", "cDC", expYazar$PopLbls)


expYazar <- expYazar[ , expYazar$PopLbls %in% selPops]
expYazar$PopLbls <- as.factor(expYazar$PopLbls) # 1,043,893

iSel = NULL
for (pop in selPops) {
	iPop <- which(expYazar$PopLbls == pop)
	iSel <- c(iSel, sample(iPop, size = 1500))
}

subYazar <- expYazar[ , iSel]


#### scTransform ####
library(sctransform)
library(future)
plan("multisession", workers = 4) 
# Assign 16GB per core to satify scTransform!
options(future.globals.maxSize = 16 * 1024^3)

subYazar@assays$RNA@counts <- subYazar@assays$RNA@data # Required for SCT

subYazar <- SCTransform(object = subYazar, 
						 vst.flavor = "v2", 
						 verbose = FALSE,
						 variable.features.n = 5000,
						 vars.to.regress = 
							c("pool_number", "age", "percent.mt"),
						 seed.use = 1234) 


subYazar <- RunPCA(subYazar, verbose = FALSE, seed.use = 1234)
subYazar <- FindNeighbors(subYazar, dims = 1:30)
subYazar <- FindClusters(subYazar)
subYazar <- RunUMAP(subYazar, dims = 1:50, verbose = FALSE, seed.use = 12)
subYazar <- RunTSNE(subYazar, dims.use = 1:30, k.seed = 42, seed.use = 1234)

saveRDS(subYazar, file =  "subYazar.RDS")

p1 <- DimPlot(subYazar, cols = rev(myColors), 
			  group.by = "PopLbls", 
			  raster = TRUE, raster.dpi = c(2048, 2048), 
			  reduction = "umap", pt.size = 6, 
			  shuffle = TRUE, seed = 42) # + NoLegend()
p2 <- LabelClusters(p1, id = "PopLbls", 
					repel = TRUE, max.overlaps = 100) 

p3 <- DimPlot(subYazar, cols = c(rev(myColors), "black"), 
			  group.by = "pool_number", 
			  raster = TRUE, raster.dpi = c(2048, 2048), 
			  reduction = "umap", pt.size = 6, 
			  shuffle = TRUE, seed = 42) # + NoLegend()

p4 <- DimPlot(subYazar, 
			  group.by = "age", 
			  raster = TRUE, raster.dpi = c(2048, 2048), 
			  reduction = "umap", pt.size = 6, 
			  shuffle = TRUE, seed = 42) # + NoLegend()

p4b <- DimPlot(subYazar, cols = c(rev(myColors), "black"), 
			  group.by = "seurat_clusters", 
			  raster = TRUE, raster.dpi = c(2048, 2048), 
			  reduction = "umap", pt.size = 6, 
			  shuffle = TRUE, seed = 42) # + NoLegend()

p5a <- DimPlot(subYazar, cols = rev(myColors), 
			  group.by = "PopLbls", 
			  raster = TRUE, raster.dpi = c(2048, 2048), 
			  reduction = "tsne", pt.size = 6, 
			  shuffle = TRUE, seed = 42) # + NoLegend()
p5 <- LabelClusters(p5a, id = "PopLbls", 
					repel = TRUE, max.overlaps = 100) 

pdf("subYazar_labeledUMAP.pdf", width = 9)
p2
p3
p4
p4b
p5
dev.off()


#-------------------------------------------------------------------------------#
#								MAGIC transform 								#
#-------------------------------------------------------------------------------#
library(Rmagic)
library(Seurat)
library(dplyr)
library(purrr)


rootDir = "/Volumes/hbolouri/SLP3_pureWorkspace/" 
setwd(rootDir)

Sys.setenv('R_MAX_VSIZE' = 8000000000)

print( Sys.time() )
suppressWarnings( RNA <- readRDS("subYazar.RDS") )
print( Sys.time() )

DefaultAssay(RNA) <- "SCT"

maxs <- apply(GetAssayData(RNA, assay = "SCT"), 1, max)
sel <- rownames(RNA)[maxs > 0]
RNA <- subset(RNA, features = sel) # 13,863 of 19,964

print(paste0("Started MAGIC   ", Sys.time()))

suppressWarnings(
RNA <- magic( RNA,
			assay = "SCT",
			n_jobs = as.integer(6),
			verbose = FALSE )
)

print(paste0("Finished MAGIC ", Sys.time()))

saveRDS(RNA, file = "subYazar_MAGIC.RDS")
