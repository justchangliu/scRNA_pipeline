#### ------------ Normalize method ------------ ####
# LogNormalize calculation: log2((count/total_counts)*10000+1)
# just a regular method to avoid the depth of the library to impact the variances from barcode to barcode
# This is A MUST DO otherwise there is not much sense to work with the data at all.
# correct barcode direction technical errors.

#### ------------- Find Variable Features ------------- ####
# just try to find a pseudo-linear correlation between mean(expression) and variance
# and get the standard deviation 
# different non-supervised method to estimate the method to pick up features


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


