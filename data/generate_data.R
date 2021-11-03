library(Biobase)
library(brgedata)
library(dplyr)
library(stringr)
library(tibble)
data("brge_gexp")

## These datasets are used in the MultiDataSet tutorial at:
## https://www.bioconductor.org/packages/devel/bioc/vignettes/MultiDataSet/inst/doc/MultiDataSet.html

brge_gexp ## the expression dataset we're going to copy


expr_df <- Biobase::exprs(brge_gexp)
write.csv(expr_df, file = "data/expression_data.csv", row.names = TRUE)

expr_samples_df <- Biobase::pData(brge_gexp)
write.csv(expr_samples_df, file = "data/samples_info_expression.csv", row.names = TRUE)

expr_annot_df <- Biobase::fData(brge_gexp)
write.csv(expr_annot_df, file = "data/gene_annotation.csv", row.names = TRUE)


## Creating completely fake snp data
set.seed(456)

chroms <- unique(expr_annot_df$chromosome)

snps_annot_df <- tibble(chromosome = rep(chroms, each = 100),
                        position = sample(1:1e6, length(chroms) * 100, replace = TRUE),
                        type = sample(c("intron", "exon", "noncoding"), length(chroms) * 100, replace = TRUE),
                        snpID = paste0(chromosome, "_", position)) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "snpID")

write.csv(snps_annot_df, file = "data/snps_annotation.csv", row.names = TRUE)

samples <- rownames(expr_samples_df)
samples_id <- as.numeric(str_remove(samples, "x"))

snps_samples_df <- tibble(samplesID = c(sample(samples_id, 70, replace = FALSE), (max(samples_id) + 1):(max(samples_id) + 51)),
                  samples = paste0("x", str_pad(samplesID, 4, pad = "0"))) %>% 
  arrange(samplesID) %>% 
  select(-samplesID) %>% 
  mutate(age = 4,
         sex = expr_samples_df[samples, "sex"],
         sex = case_when(is.na(sex) ~ sample(c("Male", "Female"), 1, replace = TRUE),
                         TRUE ~ sex),
         group = sample(c("A", "B", "C", "D"), n(), replace = TRUE)) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "samples")
  
write.csv(snps_samples_df, file = "data/samples_info_snps.csv", row.names = TRUE)


snps_df <- matrix(data = sample(0:2, nrow(snps_annot_df) * nrow(snps_samples_df), replace = TRUE),
                  nrow = nrow(snps_annot_df),
                  ncol = nrow(snps_samples_df),
                  dimnames = list(rownames(snps_annot_df), rownames(snps_samples_df)))

write.csv(snps_df, file = "data/snps_data.csv", row.names = TRUE)
