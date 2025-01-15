# run_wgcna step2: test power

###############################################################
##################### define input variable
###############################################################

suppressMessages(library(argparse))
suppressMessages(library(WGCNA))
suppressMessages(library(flashClust))
options(stringsAsFactors = FALSE)
#enableWGCNAThreads();

parser = ArgumentParser()
parser$add_argument("--cut_height", help = 'sampleTree cut_height, used for remove low quality sample. default: 80000. 可根据第一轮结果判断是否有显著离群样本，根据结果调整该值', default = "80000", required = TRUE)
parser$add_argument("--network_type", help = 'network type, unsigned/signed/signed hybrid, default: unsigned', default = "unsigned", required = TRUE)
parser$add_argument("--outdir", help = 'outdir folder', required = TRUE)
args <- parser$parse_args()

print(paste0('###############################################################'))
print(paste0('Input parameters: '))
str(args)
print(paste0('###############################################################'))

cut_height = as.numeric(args$cut_height)
network_type = args$network_type # default: unsigned
outdir = args$outdir

###############################################################
##################### load step1 result
###############################################################

load(paste0(outdir, "/01_wgcna_step1_filter_exp_matrix.Rdata")) # load datExpr0

###############################################################
##################### check sample quality
###############################################################

gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###############################################################
##################### plot sampleTree
###############################################################

sampleTree <- hclust(dist(datExpr0), method = "average")
pdf(paste0(outdir, "/01.samples_cluster_tree.pdf"), width = 25, height = 8)
    par(mar = c(0, 4, 2, 0))
    plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
    #plot中cex = 0.8 调整字体大小
    cutHeight <- cut_height 
    abline(h = cutHeight, col = "red")
dev.off()

#找到离群样本过滤掉
#clust <- cutreeStatic(sampleTree, cutHeight = cut_height)
#keepSamples <- (clust == 1)
#print(clust)
#datExpr0 = datExpr0[keepSamples, ]

###############################################################
##################### 软阈值筛选：构建的网络更符合无标度网络特征
###############################################################
print(dim(datExpr0))
#seqseq=$(echo $x | awk '{print $2}')
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft <- pickSoftThreshold(datExpr0, powerVector = powers, networkType = network_type, verbose = 5)
pdf(file = paste0(outdir, "/02.soft_threshold.pdf"), width = 8, height = 4.5)
    par(mfrow = c(1, 2))
    cex1 <- 0.9
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], xlab = "Soft threshold (power)", ylab = "Scale free topology model fit, signed R^2", type = "n", main = paste("Scale independence"))
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
    abline(h = 0.8, col = "red")
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft threshold (power)", ylab = "Mean connectivity", type = "n", main = paste("Mean connectivity"))
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# power值也可自己设定，但对应至少 >= 0.8, `power = 10`##########
power = sft$powerEstimate
print(paste0('###############################################################'))
print(paste0('powerEstimate: ', power))
print(paste0('###############################################################'))

save(datExpr0, file = paste0(outdir, "/02_wgcna_step2_test_power.Rdata"))


###############################################################
##################### 经验power (无满足条件的power时选用)
###############################################################

# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
# if (is.na(power)){
#     power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
#             ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
#             ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
#             ifelse(type == "unsigned", 6, 12))      
#             )
#             )
# }
