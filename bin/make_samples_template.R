library(readr)
suppressMessages(library(dplyr))
library(stringr)
library(tibble)
library(optparse)
library(gtools)

option_list <- list(
  make_option(c("-f", "--fq_dir"), type="character", default="../raw_data/", 
              help="Directory containing fastqs. [default = %default]", metavar="character"),
  make_option(c("--se"), action="store_true", default=FALSE,
              help="Set if single-end reads."),
  make_option(c("-s", "--sample_rgx"), type="character", default="^[^_]+", 
              help="Regex to extract sample name. [default = %default]", metavar="character"),
  make_option(c("-r", "--r1_rgx"), type="character", default="[._-]R?1([._-]|$)",
              help="R1 plus leading and trailing characters. [default = %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="samples_template.tsv", 
              help="output file name [default = %default]", metavar="character")
) 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


r1_files <- list.files(opt$fq_dir, pattern=paste0(opt$r1_rgx, ".*f(ast)?q(\\.gz)?$"), recursive=TRUE)
if (length(r1_files) == 0){
    stop("No read files found.")
}
if(!opt$se){
    r2_rgx <- str_replace(opt$r1_rgx, "1", "2")
    r2_files <- list.files(opt$fq_dir, pattern=r2_rgx, recursive=TRUE)

    if (length(r2_files) == 0){
        stop("No R2 files found. Please set --se if single end reads.")
    }
    if (length(r1_files) != length(r2_files)){
        stop("Number of R1 and R2 files not the same. Doublecheck raw data directory.")
    }
    # if paired end data, R1 and R2 files names should be identical after removing R1/R2 portion of the name.
    stopifnot(identical(str_remove(r1_files, opt$r1_rgx), str_remove(r2_files, r2_rgx)))
} else{
    r2_files <- NA
}

df <- tibble(fq1=r1_files) %>%
    dplyr::mutate(sample=str_extract(basename(fq1), opt$sample_rgx),
                  fq2=r2_files) %>%
    dplyr::select(sample, fq1, fq2)

# version sort by the sample column
version_sort <- mixedorder(df$sample)
df <- df[version_sort, ]

df %>%
    write_tsv(., opt$out)
