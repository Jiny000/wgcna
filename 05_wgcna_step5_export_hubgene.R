# run_wgcna step5: export hubgene

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
parser$add_argument("--outdir", help = 'outdir folder', required = TRUE)
args <- parser$parse_args()

print(paste0('###############################################################'))
print(paste0('Input parameters: '))
print(paste0(str(args)))
print(paste0('###############################################################'))

outdir = args$outdir

sampleorder = NA # 可通过该参数指定绘图顺序

# outdir = "run_WGCNA"

###############################################################
##################### load step4 result
###############################################################

load(paste0(outdir, "/03_wgcna_step3_construct_network.Rdata")) # datExpr0, MEs, net, moduleLabelsNewName, moduleColors, geneModuleMembership
load(file = paste0(outdir, "/03.WGCNA.tom-block.1.RData"))

###############################################################
##################### 08. exporthubgene
###############################################################

ATOM <- as.matrix(TOM)
TOM1 <- ATOM[1:round((nrow(ATOM)/2)), 1:round((nrow(ATOM)/2))]
TOM2 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM), 1:round((nrow(ATOM)/2))]
TOM3 <- ATOM[1:round((nrow(ATOM)/2)), (round(nrow(ATOM)/2)+1):nrow(ATOM)]
TOM4 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM), (round(nrow(ATOM)/2)+1):nrow(ATOM)]

dir.create(paste0(outdir, "/08.module_result"))

for(i in 1:(ncol(MEs)-1)) {
    module <- labels2colors(i)
    module_num <- paste0("ME", i)

    inModule <- moduleColors[net$goodGenes] == module
    genename <- colnames(datExpr0[, net$goodGenes])
    modGenes <- genename[inModule]

    modTOM1 <- TOM1[inModule[1:round((nrow(ATOM)/2))], inModule[1:round((nrow(ATOM)/2))]]
    modTOM2 <- TOM2[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)], inModule[1:round((nrow(ATOM)/2))]]
    modTOM3 <- TOM3[inModule[1:round((nrow(ATOM)/2))], inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
    modTOM4 <- TOM4[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)], inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]

    modTOM <- rbind(cbind(modTOM1, modTOM3), cbind(modTOM2, modTOM4))
    IMConn <- softConnectivity(datExpr0[, net$goodGenes][, modGenes])

    cyt1 <- exportNetworkToCytoscape(modTOM, 
        edgeFile = paste0(outdir, "/08.module_result/CytoscapeInput_edges_", module_num, ".tsv"), 
        nodeFile = paste0(outdir, "/08.module_result/CytoscapeInput_nodes_", module_num, ".tsv"), 
        weighted = TRUE, threshold = 0.02, nodeNames = modGenes, altNodeNames = modGenes, 
        nodeAttr = moduleLabelsNewName[net$goodGenes][inModule])

    out <- cbind(modGenes, IMConn)
    colnames(out) <- c("gene", "connectivity")
    out <- out[order(as.numeric(out[, 2]), decreasing = T), ]
    write.table(out, paste0(outdir, "/08.module_result/", module_num, "_module_gene.tsv"), sep = "\t", quote = F, row.names = F)

    nTop <- 0.05*length(modGenes)
    top <- (rank(-IMConn) <= nTop)
    out <- cbind(modGenes[top], IMConn[top])
    colnames(out) <- c("gene", "connectivity")
    out <- out[order(as.numeric(out[, 2]), decreasing = T), ]
    write.table(out, paste(outdir, "/08.module_result/", module_num, "_percent5_hubgene.tsv", sep = ""), sep = "\t", quote = F, row.names = F)

    nTop <- 0.1*length(modGenes)
    top <- (rank(-IMConn) <= nTop)
    out <- cbind(modGenes[top], IMConn[top])
    colnames(out) <- c("gene", "connectivity")
    out <- out[order(as.numeric(out[, 2]), decreasing = T), ]
    write.table(out, paste(outdir, "/08.module_result/", module_num, "_percent10_hubgene.tsv", sep = ""), sep = "\t", quote = F, row.names = F)

    nTop <- 25
    top <- (rank(-IMConn) <= nTop)
    out <- cbind(modGenes[top], IMConn[top])
    colnames(out) <- c("gene", "connectivity")
    out <- out[order(as.numeric(out[, 2]), decreasing = T), ]
    write.table(out, paste(outdir, "/08.module_result/", module_num, "_top25_hubgene.tsv", sep = ""), sep = "\t", quote = F, row.names = F)
}
