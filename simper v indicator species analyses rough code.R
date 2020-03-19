library(vegan)
library(phyloseq)
library(indicspecies)
library(car)
#Simpering this all soil samples
STL.fung.filtered.otab <- t(as(otu_table(STL.fung.filtered), "matrix"))


#First let's try simper with my categorical treatments
#Simper is run on the raw OTU table

#First let's test what OTUs are responding to my soil treatment 
#Clay V Sand
Fstl.soil.simp <- with(STL.fung.shared.map, simper(STL.fung.filtered.otab, SoilType,permutations=999))
#STL.bact.shared.map has my categorical variables still (soil type in this case) 
#and STL.bact.filtered.otab is my OTU table
Fstl.soil.simp
summary(Fstl.soil.simp)#stats summary

#First let's test what OTUs are responding to my soil treatment 
#Sterile V NonSterile
Fstl.raincom.simp <- with(STL.fung.shared.map, simper(STL.fung.filtered.otab, RainCom,permutations=999))
#STL.bact.shared.map has my categorical variables still (soil type in this case) 
#and STL.bact.filtered.otab is my OTU table
Fstl.raincom.simp
summary(Fstl.raincom.simp)#stats summary


#How does this compare to the indicator species analyses

Fstl.raincom.ind = multipatt(as.data.frame(STL.fung.filtered.otab), STL.fung.shared.map$RainCom,duleg=TRUE, 
                             control = how(nperm=9999), func="IndVal.g")
summary(Fstl.raincom.ind, indvalcomp=TRUE)