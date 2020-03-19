## try http:// if https:// URLs are not supported
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2")
BiocManager::install("vsn")
BiocManager::install("phyloseq")
library(DESeq2)#bioconductor package
library(ggplot2)
library(vsn)#bioconductor package

#nomalizing using DESEQ2
new.deseq.obj=phyloseq_to_deseq2(phyloseq_obj, ~SampleID)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
vsd<-varianceStabilizingTransformation(new.deseq.obj)#may result in the below error
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means
rld<-rlog(new.deseq.obj)#may result in the below error
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

#Fix for above common error from the support website https://support.bioconductor.org/p/63229/


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x>0]),na.rm = na.rm)/length(x))
}
geoMeans1=apply(counts(new.deseq.obj), 1, gm_mean)
dseq_f1=estimateSizeFactors(new.deseq.obj, geoMeans=geoMeans1)
head(dseq_f1)
plotSparsity(dseq_f1)
dseq_f2=varianceStabilizingTransformation(dseq_f1, fitType = "local")
head(assay(dseq_f2))
dseq_f3=varianceStabilizingTransformation(dseq_f1, blind = TRUE)
head(assay(dseq_f3))
dseq_f4=rlog(dseq_f1, blind = TRUE)#this is much slower than VST 
#rlog requires fitting a shrinkage term for each sample and each gene
#if your dataset has a lot of zeros, you will get the below message

#Warning message:
#In sparseTest(counts(object, normalized = TRUE), 0.9, 100, 0.1) :
#  the rlog assumes that data is close to a negative binomial distribution, an assumption
#which is sometimes not compatible with datasets where many genes have many zero counts
#despite a few very large counts.
#In this data, for 14.1% of genes with a sum of normalized counts above 100, it was the case 
#that a single sample's normalized count made up more than 90% of the sum over all samples.
#the threshold for this warning is 10% of genes. See plotSparsity(dds) for a visualization of this.
#We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).
head(assay(dseq_f4))

#below are visualizations of the sparsity of the data across the normalization strategies
plotSparsity(assay(dseq_f4))
plotSparsity(assay(dseq_f3))
plotSparsity(assay(dseq_f2))

#below are the visualizations of the standard deviation as read counts increase 
#https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
meanSdPlot(log2(counts(dseq_f1,normalized=TRUE) + 1))
meanSdPlot(assay(dseq_f4))
meanSdPlot(assay(dseq_f3))
meanSdPlot(assay(dseq_f2))
#"Note that the vertical axis in such plots is the square root of the variance over all samples, 
#so including the variance due to the experimental conditions. While a flat curve of the square root 
#of variance over the mean may seem like the goal of such transformations, this may be unreasonable 
#in the case of datasets with many true differences due to the experimental conditions."


matrix_obj_vst=assay(dseq_f3)#variance stabilized transformation OTU table in matrix form
matrix_obj_vst=phyloseq_obj_vst+abs(min(phyloseq_obj_vst))#add a constant to remove any zeros
#the normalization results in negatives
#I have not found a consensus on how to handle them 

#now we need to turn it back into a phyloseq object
OTU.table.vst = otu_table(matrix_obj_vst, taxa_are_rows = TRUE)
phyloseq_obj_vst=phyloseq(OTU.table.vst,TAX.file, sample_data(map.file))