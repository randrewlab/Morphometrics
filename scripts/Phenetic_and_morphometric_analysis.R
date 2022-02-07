#Pipeline to explore phenetic and morphometric data


#Check working directory
getwd()

# load packages

library(xlsx) #gives you access to the package contents
library(cluster) #for hierarchical clustering
library(vegan) # NMDS and other statistics
library(ggplot2)
library(tidyverse)

# read in the data set and assign it the name 'dataset0'. Update endrow as you collect more data. 
dataset0 <- read.xlsx(file = "data/rawdata.xlsx", sheetIndex = 1, endRow = 71, row.names = 6)
#remove metadata columns
dataset1 = subset(dataset0, select = -c(1:5) )
# #Remove incomplete data or any other rows you dont want
# dataset1 <- dataset1[-c(8, 11, 19),]

# List column names and pick which variables to use in analysis
colnames(dataset1)
dataset = dataset1[c("Petiole..length.max","Petiole..width.max","Lamina..length.max","Lamina..width.max", "Lamina..length.to.widest.point","Lamina..apex..mucro...length" ,"Bud..length."
                     ,"Bud..width","Bud..length_to_widest_.point","Bract..maximum.length.of.last","Bract..maximum.width.of.last", "Bract..maximum.length.to.widest.point.last",
                     "Bract..apex.angle.of.last" 
                     , "Bract.abaxial.surface..indumentum","Calyx..sepal.length.maximum.",  "Calyx..sepal.width.maximum.", "Calyx..sepal.length.to.widest.point.maximum","Calyx..sepal.apex.angle",
                     "Calyx..sepal.abaxial.surface..indumentum")
                     ]
  
# subset(dataset1, select = c(13:15, 24:26, 33,38,54,65:68, 78:80, 83:84, 92, 99:101, 105, 107) )

#Just select the data rows you want in the analysis (ask Rose for ideas on how to do this more neatly)
row.names(dataset)
datasetisla = dataset[c(26:47, 52:55,69:70), 1:19]
# datasetdatasetisla1 = dataset[35:36, 1:25]
# dataset1 <- rbind(datasetdatasetisla,datasetdatasetisla1)
#remove Gibbergaee 10c, too much missing data
row.names(datasetisla)
datasetisla <- datasetisla[-c(11,22),]




# # Log transform continuous data (or skip)
# cvs <- c(24:26, 33)
# log.dataset <- dataset
# log.dataset[,cvs] <- log10(log.dataset[,cvs])
# dataset <- log.dataset
# View(dataset)

# Association 
str(datasetisla)
#daisy treats factors as nominal
dist <- daisy(datasetisla, metric = "gower", type = list(symm = c(14,19))) #generates the distance matrix
colnames(datasetisla
         )

#dist #show distance matrix
clust <- hclust(dist, method = "average") #runs UPGMA cluster analysis
plot(clust)#plots dendrogram

# clustering analysis
nmds <- metaMDS(dist) #run NMDS analysis
plot(nmds, type = 't') #plots NMDS 1 vs 2

# General Importance
# This analysis provides a relative importance of each character and is carried out by projecting the original data set onto the NMDS ordination.
#The function envfit in vegan runs such an analysis and provides you with a table in which the columns identify 
# 1/2) how each variable relates to the NMDS axes (NMDS1 and 2),
#3) how important that variable was to the entire analysis (r2) 
#4) If the importance of each variable is significantly higher than would be expected by random (i.e., the p-value) <0.05 = significant
#The envfit table can also then be used to generate a biplot, whereby both the NMDS ordination (black text) and the variable importance (red arrows) are included on the same plot.
importance <- envfit(nmds, datasetisla, na.rm = TRUE) #generates table of character importance
#i_short <- importance$vectors$arrows:[-c("Petiole..length.max")]
plot(nmds, type = "t") #plots NMDS 1 vs 2
plot(importance)
#export plot to results

#Extract p, r and NMDS values to inspect in a table that is sortable by magnitude
scores(importance$vectors$pvals)
scores(importance$vectors$arrows)
rvalues <- data.frame(importance$vectors$r)
pvalues <- data.frame(importance$vectors$pvals)
arrows <- data.frame(importance$vectors$arrows)
randpandarrow <- cbind(pvalues, rvalues,arrows)
write.csv(randpandarrow, "C:/Users/hkenned6/Documents/Melichrus/Morphology (morphometrics, PATN, DELTA)/Morphometrics_2021_2022/results/importancedatasetislav1.csv")
#Sort characters for importance to understand which characters are seperating the specimens. 


# cluster Importance
#You can create groups like this
# groups <- as.factor(c(rep("Boonoo", 7), rep("Herberton", 9), rep("Yuraygir", 6), rep("erubescens", 6),
#                       rep("Gardens", 3), rep("procumbens", 1), rep("Gilgandra", 2), rep("Gurulmundi", 2),
#                       rep("Inglewood", 4), rep("Isla", 4), rep("Silent", 4), rep("Tara", 4), rep("gibberagee", 4), rep("hirsutus", 4))) #creates a vector of the groups, same order as the rows in the data set

groups <- as.factor(c (rep("erubescens", 6), rep("Gardens", 3), rep("Isla", 4), rep("Tara", 4),
                    rep("gibberagee", 3), rep("hirsutus", 4), rep ("mareeba", 4))) #creates a vector of the groups, same order as the rows in the data set

datasetisla <- data.frame(datasetisla, groups = groups) #adds the group vector to the data set

## Create boxplots for variables of interest (ggplot can order by mean)
colnames(datasetisla)
# boxplot(Bract..apex.angle.of.last~groups, datasetisla) #generates boxplots for the four pre-defined groups 

datasetisla %>%
  ggplot(aes(x= reorder(groups,Petiole..length.max), y=Petiole..length.max)) +
  geom_boxplot() +
  labs(y="Petiole..length.max", x="OTU")

#maybe write this as a loop
kruskal.test(Lamina..length.max~groups, datasetisla) #runs a Kruskal-Wallis Test to see if the groups are significantly different



colnames(datasetisla)

