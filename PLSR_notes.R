
######################################################
#PLSR notes on code
# Jennifer Jones 2020-03-31
######################################################
# this code is for explaining a single continuous variable with a large otu x sample matrix. 

library(pls)
#link to package instructions: https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf

#first make your otu x sample data frame into a "matrix"
fungi.matrix <-as.matrix(fungi)

#Data is the data file with metadata like percent mass loss. 

#this code combines the metadata file with the otu matrix. 
# you need to have the data in this format in order to run the plsr command. 
mydata <- data.frame(Data, fungi.I = I(fungi.matrix))

#this is code for fungal community composition explaining mass loss. 
#ncomp is the number of axes you want in the output. 
##I chose 5, but you can choose any number
#validation - you might want to try other options for this. I don't really know what this means. 
## "LOO" or "leave one out" is what's in the basic tutorial. 
pls1 <-plsr(ln.massloss~ fungi.I, ncomp = 5, data = mydata, validation = "LOO")
summary(pls1)
# In the output the row for your dependent variable gives the amount of variation explained by the plsr axis. 
# you can check this with a linear model (I do that below.)

#You can use the following graph to figure out how many axes you need to explain your dependent variable. 
## I think it's similar to the graphs for PCAs. 
plot(RMSEP(pls1), legendpos = "topright")

#This explained variance by each axis
explvar(pls1)

#This plots the fit of the plsr to the actual data. 
plot(pls1, ncomp = 1, asp = 1, line = TRUE)

#This is essentially the same thing without the prediction line. 
plot(pls1$scores[,1] ~ mydata$ln.massloss)
plot(pls1$scores[,2] ~ mydata$ln.massloss)

#If you want to check the statistical realtionship between the plsr axis and your dependent variable. . 
lm1<-lm(mydata$ln.massloss ~ pls1$scores[,1])
summary(lm1)#adj R = 0.4127 

#The coefficients are the 
plot(pls1$coefficients)
plot(pls1, "loadings", comps = 1)


### Selecting the top 25 OTUs that are highest weighted in the plsr axis. 
# in order to do this, you need to use the coefficients and you need to use absolute value (because the coefficients have directions). 
# usually do this with one pls axis, but you can do it with more axes too. 
plsr.ln.massloss<-data.frame(pls1$coefficients)

#make the otu names the first column in the dataset. Double check that this looks right! 
plsr.ln.massloss<-data.frame(rownames(plsr.ln.massloss),plsr.ln.massloss)
plsr.ln.massloss<-data.frame(plsr.ln.massloss[,1:2])#select the first plsr axis only. 
colnames(plsr.ln.massloss)<-c("OTU","coefficients")

## The list of top 25 OTUs that contribute to mass loss either positively or negatively. 
plsr.ln.massloss.out<-head(plsr.ln.massloss[order(abs(plsr.ln.massloss$coefficients),decreasing=T),], n=25)
plsr.ln.massloss.out

#make a vector out of the 25 otus that are most important in the pls axis. 
otu.use<-as.character(plsr.ln.massloss.out$OTU)
otu.use

#from the original dataset, select the otus that were chosen via plsr. 
fungi.subset <- fungi[, otu.use]
fungi.subset
head(fungi.subset[,1:5])

#from here you can kind of do whatever you want.
#You can figure out their taxonomy, use them in an nmds, etc. 