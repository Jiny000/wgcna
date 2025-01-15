# run_wgcna step1: select gene

###############################################################
##################### define input variable
###############################################################

suppressMessages(library(argparse))

parser = ArgumentParser()
parser$add_argument("--input_exp", help = "gene expression file with sep = '\t', genes as rows and samples as cols; gene_id in 1st col. recommend: TPM/FPKM", required = TRUE)
parser$add_argument("--large_exp", help = 'filter low expression gene, remove max expression in all sample less than `large_exp`', default = "3", required = TRUE)
parser$add_argument("--outdir", help = 'outdir folder', required = TRUE)
args <- parser$parse_args()

print(paste0('###############################################################'))
print(paste0('Input parameters: '))
str(args)
print(paste0('###############################################################'))

# exp_file = "fpkm.tsv"
exp_file = args$input_exp
large_exp = as.numeric(args$large_exp)
outdir = args$outdir

if (!file.exists(outdir)) dir.create(outdir) # 创建输出文件夹

###############################################################
##################### Checking data for excessive missing values
###############################################################

print(paste0(exp_file, ' Readign data......'))
dataExpr = read.table(exp_file, header = T, row.names = 1, sep = '\t', check.names = F, comment.char = '')

m.max = apply(dataExpr, 1, max)
dataExpr1 = dataExpr[which(m.max > large_exp), ] # remove max expression in all sample less than `large_exp`

m.mad = apply(dataExpr1, 1, mad)
dataExpr2 = dataExpr1[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01)), ] # 筛选中位绝对偏差前75%的基因，至少MAD大于0.01

print('Raw expression matrix dim:')
dim(dataExpr)
print('Filter expression matrix dim, round1:')
dim(dataExpr1)
print('Filter expression matrix dim, round2:')
dim(dataExpr2)

###############################################################
##################### log2 translation expression, used for wgcna analysis
###############################################################

datExpr0 = as.data.frame(t(dataExpr2))
datExpr0 = log2(datExpr0 + 1)

save(datExpr0, file = paste0(outdir, "/01_wgcna_step1_filter_exp_matrix.Rdata"))
# load(paste0(outdir, "/01_wgcna_step1_filter_exp_matrix.Rdata"))

# Tips
# Ref: https://rstudio-pubs-static.s3.amazonaws.com/554864_6e0f44cdc1474b8ba2c5f8326559bbe4.html
# 共表达网络：定义为加权基因网络。点代表基因，边代表基因表达相关性。加权是指对相关性值进行冥次运算 (冥次的值也就是软阈值 (power, pickSoftThreshold这个函数所做的就是确定合适的power))。无向网络的边属性计算方式为 abs(cor(genex, geney)) ^ power；有向网络的边属性计算方式为 (1+cor(genex, geney)/2) ^ power; sign hybrid 的边属性计算方式为cor(genex, geney)^power if cor > 0 else 0。这种处理方式强化了强相关，弱化了弱相关或负相关，使得相关性数值更符合无标度网络特征，更具有生物意义。如果没有合适的power，一般是由于部分样品与其它样品因为某种原因差别太大导致的，可根据具体问题移除部分 样品或查看后面的经验值。
# Module(模块)：高度內连的基因集。在无向网络中，模块内是高度相关的基因。在有向网络中，模块内是高度正相关的基因。把基因聚类成模块后，可以对每个模块进行三个层次的分析：1. 功能富集分析查看其功能特征是否与研究目的相符；2. 模块与性状进行关联分析，找出与关注性状相关度最高的模块；3. 模块 与样本进行关联分析，找到样品特异高表达的模块。
# Connectivity (连接度)：类似于网络中 “度” (degree)的概念。每个基因的连接度是与其相连的基因的边属性之和。
# Module eigengene E: 给定模型的第一主成分，代表整个模型的基因表达谱。
# Module membership: 给定基因表达谱与给定模型的eigengene的相关性。
# Intramodular connectivity: 给定基因与给定模型内其他基因的关联度，判断基因所属关系。
# Adjacency matrix (邻接矩阵)：基因和基因之间的加权相关性值构成的矩阵。
# TOM (Topological overlap matrix)：把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得的新距离矩阵，这个信息可拿来构建网络或绘制TOM图。其定义依 据是任何两个基因的相关性不近由它们

# WGCNA本质是基于相关系数的网络分析方法，适用于多样品数据模式，一般要求样本数多于15个。样本数多于20时效果更好，样本越多，结果越稳定。
# 基因表达矩阵: 常规表达矩阵即可，即基因在行，样品在列，进入分析前做一个转置。RPKM、FPKM或其它标准化方法影响不大，推荐使用Deseq2的varianceStabilizingTransformation或log2(x+1)对标准化后的数据做个转换。如果数据来自不同的批次，需要先移除批次效应。如果数据存在系统偏移，需要做下quantile normalization。
# 性状矩阵：用于关联分析的性状必须是数值型特征 (如下面示例中的Height, Weight, Diameter)。如果是区域或分类变量，需要转换为0-1矩阵的形 式(1表示属于此组或有此属性，0表示不属于此组或无此属性，如样品分组信息WT, KO, OE)。
# 推荐使用Signed network和Robust correlation (bicor)。(这个根据自己的需要，看看上面写的每个网络怎么计算的，更知道怎么选择)
# 无向网络在 power 小于 15 或有向网络 power 小于 30 内，没有一个power值可以使无标度网络图谱结构R^2达到0.8或平均连接度降到100以下，可能是由于部分样品与其他样品差别太大造成的。这可能由批次效应、样品异质性或实验条件对表达影响太大等造成, 可以通过绘制样品聚类查看分组信息、关联批次信息、处理信息和有无异常样品 。如果这确实是由有意义的生物变化引起的，也可以使用后面程序中的经验power值。
# 网络类型： 官方推荐 "signed" 或 "signed hybrid", 用signed获得的模块包含的基因会少
# 相关性计算, 官方推荐 biweight mid-correlation & bicor

# SAMPLE: args$trait_data
# ID WT KO OE Height Weight Diameter
# samp1 1 0 0 1 2 3
# samp2 1 0 0 2 4 6
# samp3 0 1 0 10 20 50
# samp4 0 1 0 15 30 80
# samp5 0 0 1 NA 9 8
# samp6 0 0 1 4 8 7