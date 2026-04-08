library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(openxlsx)

res <- read_excel(snakemake@input[[1]])

keep_cols1 <- c("CHROM","POS","TYPE","REF","ALT")
keep_cols2 <- c("Gencode_43_hugoSymbol","Gencode_43_ncbiBuild","Gencode_43_chromosome","Gencode_43_start","Gencode_43_end","Gencode_43_variantClassification","Gencode_43_genomeChange","Gencode_43_annotationTranscript","Gencode_43_cDnaChange","Gencode_43_codonChange","Gencode_43_proteinChange","ClinVar_VCF_AF_ESP","ClinVar_VCF_AF_EXAC","ClinVar_VCF_AF_TGP","ClinVar_VCF_ALLELEID","ClinVar_VCF_CLNDISDB","ClinVar_VCF_CLNDISDBINCL","ClinVar_VCF_CLNDN","ClinVar_VCF_CLNDNINCL","ClinVar_VCF_CLNHGVS","ClinVar_VCF_CLNREVSTAT","ClinVar_VCF_CLNSIG","ClinVar_VCF_CLNSIGCONF","ClinVar_VCF_CLNSIGINCL","ClinVar_VCF_CLNVC","ClinVar_VCF_CLNVCSO")

res %>%
    dplyr::select(all_of(keep_cols1), ends_with(".GT"), ends_with(".AD"), ends_with(".VAF"), all_of(keep_cols2)) %>%
    write.xlsx(snakemake@output[["v1"]], overwrite=TRUE)
        
res %>%
    dplyr::select(all_of(keep_cols1), ends_with(".GT"), ends_with(".AD"), ends_with(".VAF"), all_of(keep_cols2)) %>%
    dplyr::filter(!ClinVar_VCF_CLNSIG %in% c("Benign","Benign/Likely_benign","Likely_benign")) %>%
    write.xlsx(snakemake@output[["v2"]], overwrite=TRUE)
        
clinvar_columns <- res %>% dplyr::select(starts_with("ClinVar_VCF_"))
res$clinvar_is_missing <- apply(clinvar_columns, 1, function(x) all(is.na(x)))
        
res %>% 
    dplyr::select(all_of(keep_cols1), ends_with(".GT"), ends_with(".AD"), ends_with(".VAF"), all_of(keep_cols2), clinvar_is_missing) %>%
    dplyr::filter(!ClinVar_VCF_CLNSIG %in% c("Benign","Benign/Likely_benign","Likely_benign")) %>%
    dplyr::filter(!clinvar_is_missing) %>%
    dplyr::select(-clinvar_is_missing) %>%
    write.xlsx(snakemake@output[["v3"]], overwrite=TRUE)
