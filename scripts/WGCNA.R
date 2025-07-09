
install.packages("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)

# prepare expression data ------------------------------------------------------
# load data
dgeDispersion <- readRDS("dgeDispersion.rds")

logCPM <- cpm(dgeDispersion, log=TRUE)

# transpose columns as rows and vice versa
# WGCNA needs genes as columns and samples as rows
datExpr <- t(logCPM)

# checking for missing values (no missing values allowed in WGCNA)
# removes samples with too many missing values 
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# sample clustering ------------------------------------------------------------
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "") # saved in visuals


# choosing threshold -----------------------------------------------------------
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")

# 5 power is right at 0.8 but 6 is above 0.8 
# so choosing 6 for soft threshold power 

# setting chosen power as 6 
chosen_power = 6

# build network and identify modules -------------------------------------------
net <- blockwiseModules(datExpr, power = chosen_power,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, saveTOMFileBase = "TOM",
                        verbose = 3)

names(net)

# save net as rds
saveRDS(net, "WGCNAnet.rds")


# plot modules 
moduleColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03)

# save module assignments
geneModules <- data.frame(Gene = colnames(datExpr), Module = moduleColors)
write.csv(geneModules, "WGCNAgenemodules.csv", row.names = FALSE)

# one-hot encoding -------------------------------------------------------------

# construct trait matrix 
# load 'groups' variable from DEanalysis_prep script
traitData <- data.frame(
  MW = ifelse(groups == "MW", 1, 0),
  FW = ifelse(groups == "FW", 1, 0),
  MC = ifelse(groups == "MC", 1, 0),
  FC = ifelse(groups == "FC", 1, 0)
)

rownames(traitData) <- rownames(datExpr)

moduleTraitCor <- cor(net$MEs, traitData, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# checking to make sure there are no NA values in traitData
dim(moduleTraitCor)
summary(as.vector(moduleTraitCor))
all(rownames(traitData) == rownames(datExpr))
sum(moduleTraitPval < 0.05)

# create heatmap of module-trait relationships
textMatrix <- ifelse(moduleTraitPval < 0.05, signif(moduleTraitCor, 2), "")
pdf("moduleTraitHeatmap.pdf", width = 10, height = 8)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traitData),
  yLabels = names(net$MEs),
  ySymbols = names(net$MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = signif(moduleTraitCor, 2),  # Or textMatrix if you made a custom one
  cex.text = 0.8,
  xLabelsAngle = 45,
  main = "Module-Trait Relationships"
)
dev.off()

# plotting only significant modules on heat map

# identifying modules with at least 1 significant trait cor
sigModules <- apply(moduleTraitPval, 1, function(p_row) any(p_row < 0.05))

# subsetting matrices with significant modules 
sigCor <- moduleTraitCor[sigModules, ]
sigPval <- moduleTraitPval[sigModules, ]

#build new text matrix
textMatrix <- ifelse(sigPval < 0.05, signif(sigCor, 2), "")

# create heatmap with significant modules 
pdf("moduleTraitHeatmapSignificantOnly.pdf", width = 10, height = 8)
labeledHeatmap(
  Matrix = sigCor,
  xLabels = colnames(sigCor),
  yLabels = rownames(sigCor),
  ySymbols = rownames(sigCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  cex.text = 0.8,
  xLabelsAngle = 45,
  main = "Significant Module-Trait Relationships (p < 0.05)"
)
dev.off()



# binary traits ----------------------------------------------------------------

# extract traits 
sex <- substr(groups, 1, 1) 
treatment <- substr(groups, 2, 2)

# create binary columns 
traitdata_binary <- data.frame(
  male = ifelse(sex == "M", 1, 0),
  female = ifelse(sex == "F", 1, 0),
  water = ifelse(treatment == "W", 1, 0),
  control = ifelse(treatment == "C", 1, 0)
)

rownames(traitdata_binary) <- rownames(datExpr)

# correlation between module eigengenes 
moduleTraitCor_bin <- cor(net$MEs, traitdata_binary, use = "p")
moduleTraitPval_bin <- corPvalueStudent(moduleTraitCor_bin, nSamples = nrow(datExpr))

# checking to make sure there are no NA values in traitdata_binary
dim(moduleTraitCor_bin)
summary(as.vector(moduleTraitCor_bin))
all(rownames(traitdata_binary) == rownames(datExpr))
sum(moduleTraitPval_bin < 0.05)

# create heatmap of module-trait relationships for binary traits 
pdf("moduleTraitHeatmap_Binary.pdf", width = 10, height = 8)
labeledHeatmap(
  Matrix = moduleTraitCor_bin,
  xLabels = names(traitdata_binary),
  yLabels = names(net$MEs),
  ySymbols = names(net$MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = signif(moduleTraitCor_bin, 2),
  cex.text = 0.8,
  xLabelsAngle = 45,
  main = "Module-Trait Relationships (Binary Traits)"
)
dev.off()

# plotting only significant modules on heat maps 

# filter significant modules 
sigModules_bin <- apply(moduleTraitPval_bin, 1, function(p_row) any(p_row < 0.05))

sigCor_bin <- moduleTraitCor_bin[sigModules_bin, ]
sigPval_bin <- moduleTraitPval_bin[sigModules_bin, ]

# create matching text matrix
textMatrix_bin <- ifelse(sigPval_bin < 0.05, signif(sigCor_bin, 2), "")

# create heatmap with significant modules for binary traits 
pdf("moduleTraitHeatmap_Binary_SignificantOnly.pdf", width = 10, height = 8)
labeledHeatmap(
  Matrix = sigCor_bin,
  xLabels = colnames(sigCor_bin),
  yLabels = rownames(sigCor_bin),
  ySymbols = rownames(sigCor_bin),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_bin,
  cex.text = 0.8,
  xLabelsAngle = 45,
  main = "Significant Module-Trait Relationships (Binary Traits, p < 0.05)"
)
dev.off()
