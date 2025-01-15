# run_wgcna step3: construct network

###############################################################
##################### define input variable
###############################################################

suppressMessages(library(argparse))
suppressMessages(library(circlize))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))
suppressMessages(library(ComplexHeatmap))

options(stringsAsFactors = FALSE)
enableWGCNAThreads();

parser = ArgumentParser()
parser$add_argument("--input_power", help = 'design this value based step2 result', type = "integer", required = TRUE)
parser$add_argument("--network_type", help = 'network type, unsigned/signed/signed hybrid', default = "unsigned", required = TRUE)
parser$add_argument("--cor_type", help = 'corr type, pearson/bicor', default = "pearson", required = TRUE)
parser$add_argument("--maxBlockSize", help = 'maxBlockSize argument for blockwiseModules function. default: 100000', type = "integer", default = 100000, required = TRUE)
parser$add_argument("--interest_geneset", help = 'gene list file, per gene per line, line 1 must GeneID', required = TRUE)
parser$add_argument("--outdir", help = 'outdir folder', required = TRUE)
args <- parser$parse_args()

print(paste0('###############################################################'))
print(paste0('Input parameters: '))
print(paste0(str(args)))
print(paste0('###############################################################'))

input_power = args$input_power
network_type = args$network_type
cor_type = args$cor_type # default: pearson
maxBlockSize = args$maxBlockSize # 该值与可用内存有关，默认5000。4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可以处理3万个。计算资源允许的情况下最好放在一个block里面。
interest_geneset = args$interest_geneset
outdir = args$outdir

sampleorder = NA # 可通过该参数指定绘图顺序

# input_power = 12
# network_type = "unsigned"
# cor_type = "pearson"
# maxBlockSize = 80000
# interest_geneset = "test.list"
# trait_data = "trait_data.tsv"
# outdir = "run_WGCNA"

###############################################################
##################### load step2 result
###############################################################

load(paste0(outdir, "/02_wgcna_step2_test_power.Rdata")) # datExpr0

interest_geneset = read.table(interest_geneset, header = T, sep = '\t', check.names = F, comment.char = '')$"GeneID"

maxPOutliers = ifelse(cor_type == "pearson", 1, 0.05)
robustY = ifelse(cor_type == "pearson", T, F)

###############################################################
##################### blockwiseModules：自动网络构建和模块检测
###############################################################

net <- blockwiseModules(datExpr0, power = input_power, networkType = network_type, numericLabels = TRUE, 
    TOMType = "signed", mergeCutHeight = 0.25, minModuleSize = 50, 
    maxBlockSize = maxBlockSize, saveTOMs = TRUE, corType = cor_type, maxPOutliers = maxPOutliers, 
    saveTOMFileBase = paste0(outdir, "/03.WGCNA.tom"), verbose = 3)

moduleLabels <- net$colors
moduleLabelsNewName <- paste("ME", moduleLabels, sep = "")  
moduleColors <- labels2colors(moduleLabels)

print(paste0('###############################################################'))
table(moduleColors)
table(moduleLabelsNewName)
print(paste0('###############################################################'))

###############################################################
##################### 统计图形绘制
###############################################################

# 层级聚类树展示各模块。灰色的为**未分类**到模块的基因。
pdf(file = paste0(outdir, "/03.auto_modules_color.pdf"), width = 11, height = 6)
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],  groupLabels = "Module colors", rowText = moduleLabelsNewName, dendroLabels = FALSE, addGuide = TRUE, addTextGuide = TRUE)
dev.off()

# module eigengene, 用于绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
MEs_col = MEs
# colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs), "ME", ""))))
# MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图. marDendro/marHeatmap 设置下、左、上、右的边距
pdf(file = paste0(outdir, "/03.Eigengene_adjacency_heatmap.pdf"))
    plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
    marDendro = c(3, 3, 2, 4), marHeatmap = c(3, 4, 2, 2), plotDendrograms = T, xLabelsAngle = 90)
dev.off()

###############################################################
##################### 04. gene module membership
###############################################################

modNames <- names(MEs)
geneModuleMembership <- cor(datExpr0[, net$goodGenes], MEs, use = "p")
nSamples <- nrow(datExpr0[, net$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership, nSamples)
colnames(MMPvalue) <- paste("p.", modNames, sep = "")

text <- paste("cor=", round(geneModuleMembership, 4), "; p-value=", round(MMPvalue, 4), sep = "")
dim(text) <- dim(geneModuleMembership)
rownames(text) <- rownames(geneModuleMembership)
colnames(text) <- colnames(geneModuleMembership)
text <- cbind(rownames(text), text)
colnames(text)[1] <- "GeneID"
write.table(text, file = paste0(outdir, "/04.genes_module_membership.tsv"), quote = F, sep = "\t", row.names = F)


###############################################################
##################### 05. Module eigengenes Calculation：绘制模块之间的相关性图
###############################################################

# 计算模块特征基因（第一主成分）
# Calculates module eigengenes (1st principal component) of modules in a given single dataset.
MEs0 <- moduleEigengenes(datExpr0[, net$goodGenes], moduleLabels[net$goodGenes])$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames( datExpr0[, net$goodGenes])
text <- cbind(rownames(MEs), MEs)
colnames(text)[1] <- "samples"
write.table(text, file = paste0(outdir, "/05.module_eigengenes.tsv"), quote = F, sep = "\t", row.names = F)

# names(MEs) <- substring(names(MEs), 3)
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
pdf(file = paste0(outdir, "/05.modules_cluster_tree.pdf"), width = 7, height = 5)
    plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# 计算模块间相似度并绘图
moduleCor <- corAndPvalue(MEs, use = "p")
rowLabels <- paste("", names(MEs), sep = "")
textMatrix <- paste(signif(moduleCor$cor, 2), "\n(", signif(moduleCor$p, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleCor$cor)

pdf(file = paste0(outdir, "/05.modules_relationships.pdf"), 12, 8)
    par(mar = c(10, 10, 1, 2))
    labeledHeatmap(Matrix = moduleCor$cor, textMatrix = textMatrix, cex.text = 0.5, xLabels = rowLabels, yLabels = rowLabels, 
    xSymbols = names(MEs), ySymbols = names(MEs), colorLabels = TRUE, colors = blueWhiteRed(50), setStdMargins = FALSE, xLabelsAngle = 90, zlim = c(-1, 1))
dev.off()

text <- paste("cor=", round(moduleCor$cor, 4), "; p-value=", round(moduleCor$p, 4), sep = "")
dim(text) <- dim(moduleCor$cor)
rownames(text) <- rowLabels
colnames(text) <- rowLabels
text <- cbind(rownames(text), text)
colnames(text)[1] <- "modules"
write.table(text, file = paste0(outdir, "/05.modules_relationships_plot_data.tsv"), quote = F, sep = "\t", row.names = F)

# 绘制模块内基因表达趋势
dir.create(paste0(outdir, "/05.expression_ME"))
for(i in 1:(ncol(MEs)-1)) {
    which.module <- paste0("ME", i)
    which.module2 <- labels2colors(i)
    dir <- paste0(outdir, "/05.expression_ME")
    pdf(file = paste(dir, "/05.expression_ME_", which.module, ".pdf", sep = ""), 25, 10)
        ME <- MEs[, which.module]
        ME <- t(as.matrix(MEs[, which.module]))
        colnames(ME) <- rownames(datExpr0[, net$goodGenes])
        layout(matrix(c(1, 2)), heights = c(1.5, 3))
        par(mar = c(0.3, 9, 3, 5))
        plotMat(t(scale(datExpr0[, net$goodGenes][, moduleLabelsNewName[net$goodGenes] == which.module])), nrgcols = 30, rlabels = F, rcols = which.module2, main = paste0(which.module, ": ", which.module2), cex.main = 1)
        par(mar = c(8, 4, 0, 1))
        barplot(ME, col = which.module2, main = "", cex.names = 1, cex.axis = 1, ylab = "module eigengene", las = 3)
    dev.off()
}


###############################################################
##################### 06. Relating modules to external information (samples/traits)
###############################################################

sample_cor <- cor(t(datExpr0[, net$goodGenes]), t(datExpr0[, net$goodGenes]), use = 'pairwise.complete.obs')
moduleSampleCor <- cor(MEs, sample_cor, use = "p")
nSamples <- nrow( datExpr0[, net$goodGenes])
rowLabels <- names(MEs)
# rowLabels <- paste("ME", names(MEs), sep = "")

# 绘制模块和样品的相关性
pdf(file = paste0(outdir, "/06.modules_samples_relationships.pdf"), 0.6*nSamples, 0.6*length(rowLabels))
    par(mar = c(7, 12, 1, 1))

    colorstat = as.matrix(table(moduleLabelsNewName))
    colorstat = colorstat[row.names(moduleSampleCor), ]
    genemodule = data.frame(gene = colnames(datExpr0), module = moduleLabelsNewName)

    test.overlap = lapply(unique(genemodule$module), function(x){
        modulegene = genemodule[genemodule$module == x, ]$gene
        return(length(intersect(modulegene, interest_geneset)))
    })%>%unlist()

    names(test.overlap) = unique(genemodule$module)
    row_anno = HeatmapAnnotation(genenum = anno_barplot(colorstat, add_numbers = TRUE, width = unit(2, 'cm')), interest_genenum = anno_barplot(test.overlap[names(colorstat)], add_numbers = TRUE, width = unit(2, 'cm')), which = 'row')
    col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

    if (!is.na(sampleorder)){moduleSampleCor = moduleSampleCor[, sampleorder]}

    moduleSamplePvalue <- corPvalueStudent(moduleSampleCor, nSamples)
    textMatrix <- paste(signif(moduleSampleCor, 2), "\n(", signif(moduleSamplePvalue, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleSampleCor)
    print(Heatmap(moduleSampleCor, right_annotation = row_anno, cluster_columns = F, name = 'cor', col = col_fun, 
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%s", textMatrix[i, j]), x, y, gp = gpar(fontsize = 10))
    }))
dev.off()

text <- paste("cor=", round(moduleSampleCor, 4), "; p-value=", round(moduleSamplePvalue, 4), sep = "")
dim(text) <- dim(moduleSampleCor)
rownames(text) <- rownames(moduleSampleCor)
colnames(text) <- colnames(moduleSampleCor)
text <- cbind(rownames(text), text)
colnames(text)[1] <- "modules"
colnames(moduleSampleCor)[1] <- "modules"

write.table(text, file = paste0(outdir, "/06.modules_samples_relationships.cor_pvalue.tsv"), quote = F, sep = "\t", row.names = F)
write.table(moduleSampleCor, file = paste0(outdir, "/06.modules_samples_relationships.cor.tsv"), quote = F, sep = "\t", row.names = T)
write.table(genemodule, file = paste0(outdir, "/06.genemodule_information.tsv"), quote = F, sep = "\t", row.names = F)

save(datExpr0, MEs, net, moduleLabelsNewName, moduleColors, geneModuleMembership, file = paste0(outdir, "/03_wgcna_step3_construct_network.Rdata"))

