#Objective: To subset morphological data by including or excluding samples based on filters
# (e.g. list of characters, groupings, etc.)
# Description of script:* This is a two-step filtering process for the downstream analysis with filtering process 
# (1) based on the characters used and (2) based on the samples to include.
# Exploration of data:* Using data explorer to examine subset of data for missing data, data type, 
# distribution type, correlation of characters and distribution of characters within OTU groups. 

# Install Packages
# install.packages("knitr")
# install.packages("vegan")
# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("naniar")
# install.packages("ggcorrplot")
# install.packages("DataExplorer")
#install.packages("corrplot")

# Check WD
getwd()

#Load packages
library(cluster) #for hierarchical clustering
library(vegan) # NMDS and other statistics
library(ggplot2)
library(tidyverse)
library(GGally)
library(naniar)
library(ggcorrplot)
library(DataExplorer)
library(dplyr)
library(corrplot)

#Define names for input files, directories and subset name
# "indir" is the directory containing the input files

indir <- "data"
rawdata <- "rawdata.csv"
column_metadata <- "column_metadata.csv"
row_metadata <- "row_metadata.csv"
subset_suffix = "group"
variableset_suffix = "set"
id_col="collection_number"

variableset <- "erubgalam"

rawdatafile <- file.path(indir,rawdata)
column_metadatafile <- file.path(indir,column_metadata)
row_metadatafile <- file.path(indir,row_metadata)

## Documentation of data filtering and manipulation for morphological analysis

# The analysis is based on three data files:
# `data/rawdata.csv`
# `data/column_metadata.csv`
# `data/row_metadata.csv`

#Load data and metadata

dataset0<-read_csv(rawdatafile) 
col_meta<-read_csv(column_metadatafile) 
row_meta<-read_csv(row_metadatafile) 

# Check integrity of rawdata headings with column metadata values
# Output is two lists; 
# 1. values which are in metadata but not in rawdata
# 2. values which are in rawdata but not in metadata

# If no values are returned, lists match perfectly.

## Characters in metadata without matches in dataset
col_meta$character[!col_meta$character %in% names(dataset0)]

## Characters in dataset without matches in metadata
names(dataset0)[!names(dataset0) %in% col_meta$character]

# Filtering of variables (characters) to use in analysis based on the column_metadata file:
# The selection of characters for downstream analyses is based on the filtering of groups in the `data/column_metadata.csv` file.
# Characters are included in a group if marked with a yes (`Y`).Then these characters are pulled in the rawdata file 
# based on the list of selected characters in the `data/column_metadata.csv`.

include_strings=c("1","y","Y","yes","Yes")
subsetcol=paste0(variableset,"_",variableset_suffix)
include_columns<-filter(col_meta,.data[[subsetcol]] %in% include_strings)

## Selecting samples to include in analysis (filtering on rows) based on the row_metadata file:
#'include_rows' includes all samples with a "y" in the "isla_group" column.
# Use semi_join to keep all rows in dataset0 that match "include_rows" based on collection number

subsetrow=paste0(variableset,"_",subset_suffix)
include_rows<-filter(row_meta,.data[[subsetrow]]=="y")

## Combine to filter rawdata

dataset1 <- dataset0 %>% 
  select(all_of(c(id_col,include_columns$character))) %>% 
  semi_join(include_rows, by = id_col)

## Inspect filtered dataset for variable type and for missing data. 

# ## Check variable types
# str(dataset1)
# 
# ## Change variable types as needed
# #dataset1$`Petiole; length max/width max` <- as.numeric(dataset1$`Petiole; length max/width max`)
# #dataset1$`Lamina; margin width (of hyaline section)` <- as.numeric(dataset1$`Lamina; margin width (of hyaline section)`)
# 
# ## Calculate proportion of missing data
# missing.prop=sapply(dataset1,function(x)sum(is.na(x))/length(x)) 
# missing.prop
# 
# ## Set and apply missing data threshold 
# max.miss.prop=.5
# too.much.missing=which(missing.prop>max.miss.prop)
# 
# if(any(too.much.missing)){
#   dataset2=dataset1[,-c(too.much.missing)]
#   dataset1 <- dataset2 
# }

## Visualise missing data
#  vis_miss(dataset1)
#  plot_missing(dataset1)
## To remove any row with any amount of missing data
#  dataset1 <- dataset1 %>%
#  na.omit()

## Use DataExplorer to produce a pdf with; data structure, type, missing data, Pearsons correlation heatmap

## See vignette for more detailed explanations: 
##  https://cran.r-project.org/web/packages/DataExplorer/vignettes/dataexplorer-intro.html

boxplot1 <- plot_boxplot(
  dataset1,
  "OTU",
  binary_as_factor = TRUE,
  geom_boxplot_args = list(),
  scale_y = "continuous",
  title = NULL,
  ggtheme = theme_gray(),
  theme_config = list(),
  nrow = 3L,
  ncol = 4L,
  parallel = FALSE)


config <- configure_report(
  add_plot_str = FALSE,
  add_plot_prcomp = FALSE,
  add_plot_boxplot = TRUE,
  add_plot_scatterplot = FALSE,
  plot_boxplot_args = list(
   boxplot1),
  global_ggtheme = quote(theme_minimal(base_size = 14))
)
create_report(dataset1, config = config)
create_report(dataset1, output_file = "report.pdf", output_format = "pdf_document")


##Generate a dataframe with a summary of the pearsons correlation values
#Remove character values
dataset2 <- select(dataset1, -collection_number)
dataset2 <- select(dataset2, -OTU)
str(dataset2)
## USe rquery to produce correlation matrix, then extract as a data frame
source("http://www.sthda.com/upload/rquery_cormat.r")
rquery.cormat(dataset2)
corvalues0 = rquery.cormat(dataset2, type="flatten", graph=FALSE)
corvalues <- corvalues0$r


# ## Broad summary (data type and missing data)
# plot_intro(dataset1) 
# ## Barplots - should show only discrete variables
# plot_bar(dataset1[,-1]) 
# ## Histograms - should show only continuous variables
# ##  - may not vary much in height with few samples, but pattern of bars
# ##    should allow any worrying distributions (non normal, may require transformation, bimodal could be interesting) to be spotted
# plot_histogram(dataset1)
# ## Quantile-quantile plots are a way to visualize the deviation from a specific
# ##  probability distribution - normal by default
# plot_qq(dataset1)
# ## Plot correlation heatmap - probably needs finessing
# corrM <- plot_correlation(na.omit(dataset1))
# ## Plot all vs the first variable - doesn't seem to like name with spaces
# ##  or semicolons
# dataset1_rename=dataset1%>%
#   select(where(is.double))
# names(dataset1_rename)=gsub(" ","_",names(dataset1_rename),fixed = T)
# names(dataset1_rename)=gsub(";","",names(dataset1_rename),fixed = T)
# plot_scatterplot(dataset1_rename,by = names(dataset1_rename[,2]))


