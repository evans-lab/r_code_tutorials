
library(otuSummary)
library(vegan)
library(ggplot2)
#Here is my randomize community matrix
mock_com=data.frame("samples"=c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10"),
                    "trt"=c(rep("t1",5),rep("t2",5)),
                    "x_corr"=runif(10, min = 1, max = 10),"y_corr"=runif(10, min = 1, max = 10),
                    "sp1"=runif(10, min = 1, max = 100),"sp2"=runif(10, min = 1, max = 100),
                    "sp3"=runif(10, min = 1, max = 100),"sp4"=runif(10, min = 1, max = 100),
                    "sp5"=runif(10, min = 1, max = 100),"sp6"=runif(10, min = 1, max = 100),
                    "sp7"=runif(10, min = 1, max = 100),"sp8"=runif(10, min = 1, max = 100),
                    "sp9"=runif(10, min = 1, max = 100),"sp10"=runif(10, min = 1, max = 100))

#I would make my samples names row names just so I can keep them in my distance matrix

row.names(mock_com)=mock_com$samples


#You may want to subset your data into treatment to reduce the number of pairwise comparisons
#I am not going to do this since this is just a toy datset

#Let's calculate the physical distance between points

mock_com_phy_dis=vegdist(mock_com[,c("x_corr","y_corr")],method = "euclidean")

#Now make this this triangle into a three column matrix using this handy command

mock_com_phy_dis_M=matrixConvert(mock_com_phy_dis, 
              colname = c("sample1", "sample2", "sample_dist"))


#Now calculate the community distance using Bray-Curtis

mock_com_sp_dis=vegdist(mock_com[,5:14],method = "bray")

#Now make this this triangle into a three column matrix using this handy command

mock_com_sp_dis_M=matrixConvert(mock_com_sp_dis, 
                                 colname = c("sample1", "sample2", "bray"))

#Now you want to merge the two matrices you have created
#I am unnecessarily worried about my data getting out of order during merging
#So I like to merge with a column that is made from the combination of samples
#you can probably just merge using the sample1 column but I am paranoid

mock_com_phy_dis_M$s1_s2=with(mock_com_phy_dis_M,interaction(sample1,sample2))
#You will need to remove the sample1 and sample2 column so you will not have a repeat when you merge
mock_com_phy_dis_M[,c("sample1","sample2")]=NULL


mock_com_sp_dis_M$s1_s2=with(mock_com_sp_dis_M,interaction(sample1,sample2))

mock_com_phy_com_dis_M=merge(mock_com_phy_dis_M,mock_com_sp_dis_M, by="s1_s2")

#Now you can put your treatments in the matrix

mock_com_phy_com_dis_M_trt=merge(mock_com_phy_com_dis_M,mock_com[,c("samples","trt")], 
                                 by.x="sample1",by.y="samples")

#Now you need to rename the treatment column to reflex this is the sample1 treatment
colnames(mock_com_phy_com_dis_M_trt)[names(mock_com_phy_com_dis_M_trt)=="trt"]="trt_sample1"

#Do the same thing for sample2 so you have treatments for both sides of the comparsion

mock_com_phy_com_dis_M_trt=merge(mock_com_phy_com_dis_M_trt,mock_com[,c("samples","trt")], 
                                 by.x="sample2",by.y="samples")

#Now you need to rename the treatment column to reflect this is the sample2 treatment
colnames(mock_com_phy_com_dis_M_trt)[names(mock_com_phy_com_dis_M_trt)=="trt"]="trt_sample2"

#Now is when you need to start making decisions about what you want to compare 

#Here is a graph with the treatment in sample1 as the color of the points

ggplot(mock_com_phy_com_dis_M_trt, aes(x=sample_dist,y=bray, color=trt_sample1))+geom_point()

#Notice that t2 is pretty rare in the graph.
#this is because t2 is most of the time in the sample2 portion of the comparison

ggplot(mock_com_phy_com_dis_M_trt, aes(x=sample_dist,y=bray, color=trt_sample2))+geom_point()

#Here are all of the comparisons the dataset
mock_com_phy_com_dis_M_trt$trt_s1_s2=with(mock_com_phy_com_dis_M_trt,
                                          interaction(trt_sample1,trt_sample2))


ggplot(mock_com_phy_com_dis_M_trt, aes(x=sample_dist,y=bray, color=trt_s1_s2))+geom_point()


#I have learned the order of your dataset when you calculated the distance matrix is the driver of 
#the order of your treatment comparisons. Just be careful. 

