#### Normalize method ####
# LogNormalize calculation: log2((count/total_counts)*10000+1)
# just a regular method to avoid the depth of the library to impact the variances from barcode to barcode
# This is A MUST DO otherwise there is not much sense to work with the data at all.
# correct barcode direction technical errors.

#### Find Variable Features ####
# just try to find a pseudo-linear correlation between mean(expression) and variance
# and get the standard deviation 
# different non-supervised method to estimate the method to pick up features


#### ScaleData ####

