library(MultiDataSet)

## Loading the transcriptomics data
expr_data <- read.csv("data/expression_data.csv", row.names = 1)
expr_data <- as.matrix(expr_data)
expr_samples_info <- read.csv("data/samples_info_expression.csv", row.names = 1)
expr_annot <- read.csv("data/gene_annotation.csv", row.names = 1)

## Loading the snps data
snps_data <- read.csv("data/snps_data.csv", row.names = 1)
snps_data <- as.matrix(snps_data)
snps_samples_info <- read.csv("data/samples_info_snps.csv", row.names = 1)
snps_annot <- read.csv("data/snps_annotation.csv", row.names = 1)

## ----------------------------------------------------- ##
## Step 1: Gather the datasets into standardised classes ##
## ----------------------------------------------------- ##

## Create a SnpSet for the snps data
snp_set <- new("SnpSet", 
               call = snps_data,
               phenoData = Biobase::AnnotatedDataFrame(snps_samples_info),
               featureData = Biobase::AnnotatedDataFrame(snps_annot))


snp_set


## Create an ExpressionSet for the transcriptomics data
expr_set <- ExpressionSet(assayData = expr_data,
                          phenoData = Biobase::AnnotatedDataFrame(expr_samples_info),
                          featureData = Biobase::AnnotatedDataFrame(expr_annot))

expr_set

## MultiDataSet offers SnpSet, ExpressionSet, MethylationSet
## Also accepts eSet and SummarizedExperiment, as well as classes derived from eSet
## https://www.bioconductor.org/packages/release/bioc/vignettes/MultiDataSet/inst/doc/MultiDataSet_Extending_Proteome.html

## -------------------------------------- ##
## Step 2: create the MultiDataSet object ##
## -------------------------------------- ##

multiomics_set <- createMultiDataSet()
multiomics_set <- add_snps(multiomics_set, snp_set)
multiomics_set <- add_genexp(multiomics_set, expr_set)


multiomics_set

## ----------------------------- ##
## Query the MultiDataSet object ##
## ----------------------------- ##

names(multiomics_set) ## name of datasets

sampleNames(multiomics_set) ## list of samples in each dataset
commonIds(multiomics_set) ## list of samples present in all datasets

## --------------------------- ##
## Subset the multi-omics data ##
## --------------------------- ##

## retain only common samples
multiomics_set2 <- commonSamples(multiomics_set) 
multiomics_set2

## retain only some samples, even if not present in all datasets
multiomics_set[c("x0001", "x0002"), ]

## retain only some samples based on phenotype information
subset(multiomics_set, , sex == "Female")

## retain only some features based on annotation
subset(multiomics_set, chromosome == "chr1")

subset(multiomics_set, locus.type == "Coding")

## Extract the datasets
datasets <- as.list(multiomics_set)
str(datasets)

## Extract the features annotation
featureData(multiomics_set)
