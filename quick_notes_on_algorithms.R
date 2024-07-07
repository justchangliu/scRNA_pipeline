
# Basic Pipeline Go Through
# https://wutaoblog.com.cn/archives/seurat_pipeline
# 

#### ------------ Normalize method ------------ ####
# LogNormalize calculation: log2((count/total_counts)*10000+1)
# just a regular method to avoid the depth of the library to impact the variances from barcode to barcode
# This is A MUST DO otherwise there is not much sense to work with the data at all.
# correct barcode direction technical errors.

#### ------------- Find Variable Features ------------- ####
# just try to find a pseudo-linear correlation between mean(expression) and variance
# and get the standard deviation 
# different non-supervised method to estimate the method to pick up features
# 

#### ScaleData() ####
# ScaleData 是用来standardization。 
# Scale.Data is z scores, 或者说how many sigma away from standard deviation.
# 那这个值本质上所表达的是什么？
# 这是因为下一步这个PCA是在提取特征向量，说白了就是看这些gene内部的variance对
# 说白了就是每个feature都有自己的distribution。现在就是根据这个distribution的characteristics来对数据进行降纬
# https://www.jianshu.com/p/a2de65ebb069
# This explanation is decent 
# this blog shows what is stored in each slot in seurat object
# https://blog.csdn.net/Nh_code/article/details/125566701
# regress out genes就是用把 vars.regress.genes用作base vectors来进行投影，剩下的变量进行intercept 
# 


### --------   Data Integration   ------------------   ####
# https://blog.csdn.net/weixin_60734652/article/details/133022518 
# seurat data integration 算法解析
# seurat data integration就是对所有的combined count matrix进行校正，得到一个新的matrix。


###  --------   Find Markers and Find all Markers ------------- ###
# how this algorithm works
# do repeated t tests and rank p values between clusters and see which molecules are highly expressed
# 

###  --------   Findneighbors/clustering/UMAP and tSNE plot ----------- ###

# https://blog.csdn.net/weixin_53637133/article/details/138020905
# Jaccard distance: 比较两个样品之间的相似度度来判断neighboring distance
# https://www.jianshu.com/p/72313d70d9ab
# FindNeighbors是KNN+SNN聚类  KNN计算最近邻，SNN计算共享最近邻-均是计算的过程，可以认为是将细胞进行连线的过程

# FindClusters是Louvain算法分群 将计算好的细胞划分成不同的clusters，可以认为这是一个切割的过程，即将高维空间中连好线的细胞群切割成群
# FindNeighbors耗内存，时间慢，需要先经过PCA降维，去除差异小的基因
# https://juejin.cn/post/7171032224890880007
# 模块度的定义就是描述社区内紧密程度的值, louvian score就是来描述社区内紧密程度的值。
# louvian算法就是根据相邻关系的权重进行划分，是一种clustering的算法
# Findneigbors 是产生距离列表，对每个cell 之间的距离根据gene值进行计算。
# PCA只是用来降维，降低Findneighbors的算法的complexity.

### -------------------  ARACNe algorithm ---------------------   ###
# Use ARACNe 来create a net。
# ARACNe算法的核心是什么：把给定的gene matrix，通过randomly permute这个gene matrix，
# 来算一个MI，然后通过估算MI threshold来剔除不好的tf相连接的target。
# this net就是类似louvain算法里产生的这个graph。代表某个东西以及相关联的连线。
# 其实就是计算出了提供的列表中的gene,以及与其高相关度的gene。
# 根据这个pseudo protein activity来进行clustering，这样看到的东西更加清晰。
# 那viper algorithm说白了就是在做手动的enrichment。挑选出genes来看东西。
# MI threshold estimation. This preprocessing step identifies a significance threshold of MI values from the GEPs provided. The
# threshold depends on the number of samples provided in the
# input.
# 2. Bootstrapping/MI network reconstruction. In this phase MI networks are reconstructed for randomly sampled GEP. For N such
# bootstraps of the data N MI networks are generated. The calculation of the networks involves three steps: (a) Compute MI for
# every TF/Target pair after rank-transformation of the GEPs. (b)
# Removal of non-statistically significant connections using the
# MI threshold. (c) Removal of indirect interactions by applying a
# Data Processing Inequality tolerance filter (DPI, Margolin et al.,
# 2006).
# 3. Building consensus network. A consensus network is computed
# by estimating the statistical significance of the number of times a
# specific edge is detected across all bootstrap runs, based on a
# Poisson distribution. Only significant pairs are kept (P < 0.05,
# Bonferroni corrected).
# doi: 10.1093/bioinformatics/btw216 this work is about

###  -------------------  VIPER Algorithm ------------------------ ###

# Finally, we describe the use of VIPER to evaluate all non-silent
# somatic mutations in TCGA samples and report the aberrant activity
# of all oncogenes listed in the Catalogue Of Somatic Mutations In
# Cancer (COSMIC)22 on an individual sample basis. VIPER can be
# used to systematically assess aberrant activity of oncoproteins for
# which high-affinity inhibitors are available, independent of their
# mutational state, thus establishing them as valuable therapeutic targets on an individual patient basis.
#  The VIPER algorithm tests for regulon enrichment on gene expression signatures. The gene expression signature is first obtained by comparing
# two groups of samples representing distinctive phenotypes or treatments. Any
# method that generates a quantitative measurement of difference between the
# groups can be used (fold change, Student’s t-test, Mann-Whitney U test, etc.).
# Alternatively, single-sample-based gene expression signatures can be obtained
# by comparing the expression levels of each feature in each sample against a set
# of reference samples by any suitable method, including for example Student’s
# t-test, Z-score transformation or fold change; or relative to the average expression level across all samples when clear reference samples are not available.
# Then we compute the enrichment of each regulon on the gene expression
# signature using different implementations of aREA (see below). Finally, we
# estimate the significance, including P value and normalized enrichment score,
# by comparing each regulon enrichment score to a null model generated by
# randomly and uniformly permuting the samples 1,000 times. Alternatively,
# when the number of samples is not enough to support permutation with reposition (at least five samples per group is required), permutation of the genes
# in the gene expression signature or its analytic approximation can be used
# (see below).
# 

### -------------------- The VIPER algorithm -------------------- ###
# https://bioconductor.org/packages/release/bioc/vignettes/viper/inst/doc/viper.pdf
# GES: Gene Expression Signature
# signature <- rowTtest(dset, "description", c("CB", "CC"), "N")
# First, input gene expressiond ata, and in VIPER package, a function will efficiently performs Student's t-test for each row of
# a dataset. 
# > nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000,
# + repos = TRUE, verbose = FALSE)
# and then get a nullmodel by permutating the sample randomly. 
# A regulon object from the ARACNe network
# 我就这么理解吧，这个msVIPER
# regulons include the list of transcriptors and the up/down regulated targets.
# VIPER is the extension of msVIPER to single sample-based analysis. It effectively transforms a gene
# expression matrix to a regulatory protein activity matrix. The simplest implementation of VIPER is based
# on single-sample gene expression signatures obtained by scaling the probes or genes  subtracting the mean
# and dividing by the standard deviation of each row
# msVIPER是在估算master regular gene的enrichment根据这些target genes的activities。
# 换言之，这个regulon network 提供的是什么？
# regulon network，是在一种cell type中找出来的transcriptor factors以及与其强相关target genes.
# 根据这些基因，我们在gene network里面用msVIPER来推算master regulator‘s activities
# viper matrix就是算出来的normalized enrichment score matrix.
# this matrix 是一种pseudo protein activities. 
# 现在的问题是viper是怎么从regulon这个network来倒退这个normalized enrichment score的。
# 这个部分我还没有理解清楚。







