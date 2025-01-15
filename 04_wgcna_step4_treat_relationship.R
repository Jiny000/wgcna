# run_wgcna step4: module & trait relationship

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
parser$add_argument("--interest_geneset", help = 'gene list file, per gene per line, line 1 must GeneID', required = TRUE)
parser$add_argument("--trait_data", help = 'trait data, can be omitted', default = "")
parser$add_argument("--outdir", help = 'outdir folder', required = TRUE)
args <- parser$parse_args()

print(paste0('###############################################################'))
print(paste0('Input parameters: '))
print(paste0(str(args)))
print(paste0('###############################################################'))

interest_geneset = args$interest_geneset
trait_data = args$trait_data
outdir = args$outdir

sampleorder = NA # 可通过该参数指定绘图顺序

# interest_geneset = "test.list"
# trait_data = "trait_data.tsv"
# outdir = "run_WGCNA"

###############################################################
##################### load step3 result
###############################################################

load(paste0(outdir, "/03_wgcna_step3_construct_network.Rdata")) # datExpr0, MEs, net, moduleLabelsNewName, moduleColors, geneModuleMembership

interest_geneset = read.table(interest_geneset, header = T, sep = '\t', check.names = F, comment.char = '')$"GeneID"
input_trait <- read.table(trait_data, header = T, row.names = 1, sep = '\t', check.names = F, comment.char = '')

###############################################################
##################### 07. module trait relationship
###############################################################

datTraits <- input_trait[intersect(row.names(MEs), row.names(input_trait)), ]
datTraits[is.na(datTraits)] <- 0
MEs1 <- MEs[intersect(row.names(MEs), row.names(datTraits)), ] # MEs1：选取有对应表型性状的样本

moduleTraitCor <- cor(MEs1, datTraits, use = "p")
if (!is.na(sampleorder)){moduleTraitCor = moduleTraitCor[, sampleorder]}
nSamples <- nrow(datExpr0[, net$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file = paste0(outdir, "/06.modules_trait_relationships.pdf"), 15, 8)
    par(mar = c(7, 12, 1, 2))
    colorstat = as.matrix(table(moduleLabelsNewName))
    colorstat = colorstat[row.names(moduleTraitCor), ]
    genemodule = data.frame(gene = colnames(datExpr0), module = moduleLabelsNewName)

    test.overlap = lapply(unique(genemodule$module), function(x){
        modulegene = genemodule[genemodule$module == x, ]$gene
        return(length(intersect(modulegene, interest_geneset)))
    })%>%unlist()
    names(test.overlap) = unique(genemodule$module)
    row_anno = HeatmapAnnotation(genenum = anno_barplot(colorstat, add_numbers = TRUE, width = unit(2, 'cm')), interest_genenum = anno_barplot(test.overlap[names(colorstat)], add_numbers = TRUE, width = unit(2, 'cm')), which = 'row')
    col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

    print(Heatmap(moduleTraitCor, right_annotation = row_anno, cluster_columns = F, name = 'cor', col = col_fun, 
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%s", textMatrix[i, j]), x, y, gp = gpar(fontsize = 10))})
        )
dev.off()

text <- paste("cor=", round(moduleTraitCor, 4), "; p-value=", round(moduleTraitPvalue, 4), sep = "")
dim(text) <- dim(moduleTraitCor)
rownames(text) <- rownames(moduleTraitCor)
colnames(text) <- colnames(moduleTraitCor)
text <- cbind(rownames(text), text)
colnames(text)[1] <- "modules"
write.table(text, file = paste0(outdir, "/06.modules_trait_relationships.tsv"), quote = F, sep = "\t", row.names = F)


# each MM_vs_trait and genes_vs_trait
dir.create(paste0(outdir, '/07.MM_vs_trait'))
dir.create(paste0(outdir, '/07.genes_vs_trait'))

for (t in 1:length(colnames(datTraits))) {
    modNames <- names(MEs1)
    trigly <- as.data.frame(datTraits[, t])
    names(trigly) <- colnames(datTraits)[t]
    geneTraitSignificance <- as.data.frame(cor(datExpr0[, net$goodGenes], trigly, use = "p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    names(geneTraitSignificance) <- paste("GS.", colnames(datTraits)[t], sep = "")
    names(GSPvalue) <- paste("p.GS.", colnames(datTraits)[t], sep = "")
    # head(geneTraitSignificance)
    # head(GSPvalue)
    m = paste0(outdir, "/07.MM_vs_trait/07.mm_vs_", colnames(datTraits)[t])
    dir.create(m)
    for(i in 1:(ncol(MEs)-1)) {
        which.module <- paste0("ME", i)
        which.module2 <- labels2colors(i)
        column <- match(which.module, modNames)
        moduleGenes <- moduleColors[net$goodGenes] == which.module2
        pdf(file = paste(m, "/07.", which.module, "_mm_vs_", colnames(datTraits)[t], ".pdf", sep = ""), 6, 6)
            trigly <- colnames(datTraits)[t]
            verboseScatterplot(geneModuleMembership[moduleGenes, column], geneTraitSignificance[moduleGenes, 1], 
            xlab = paste("Module membership (MM) in", which.module, "module"), 
            ylab = paste("Gene significance for ", colnames(datTraits)[t]), 
            main = paste("Module membership vs. gene significance\n"), col = which.module2)
        dev.off()
    }

    text <- cbind(geneTraitSignificance, GSPvalue)
    text <- cbind(rownames(text), text)
    colnames(text)[1] <- "GeneID"
    o = paste0(outdir, "/07.genes_vs_trait/07.genes_trait_significance_", colnames(datTraits)[t], ".tsv")
    write.table(text, file = o, quote = F, sep = "\t", row.names = F)
}

