library(GEOquery)
library(SummarizedExperiment)
library(tidyverse)
library(here)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicRanges)
library(txdbmaker)

# retrieve annotation info
txdb_path <- here("data/txdb_gencode_v49.rds")

if (file.exists(txdb_path)) {
    txdb <- AnnotationDbi::loadDb(txdb_path)
} else {
    gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    gtf_tmp <- tempfile(fileext = ".gtf.gz")
    curl::curl_download(gtf_url, gtf_tmp, quiet = FALSE)
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_tmp, format = "gtf")
    AnnotationDbi::saveDb(txdb, txdb_path)
}

hemoglobin_genes <- c(
    ENSG00000206172 = "HBA1",
    ENSG00000188536 = "HBA2",
    ENSG00000244734 = "HBB",
    ENSG00000229988 = "HBBP1",
    ENSG00000223609 = "HBD",
    ENSG00000213931 = "HBE1",
    ENSG00000213934 = "HBG1",
    ENSG00000196565 = "HBG2",
    ENSG00000206177 = "HBM",
    ENSG00000086506 = "HBQ1",
    ENSG00000130656 = "HBZ",
    ENSG00000206178 = "HBZP1"
)

# retrieve phenotype data
gset <- getGEO("GSE236892", GSEMatrix = TRUE, getGPL = TRUE)
gset <- gset[[1]]

pdata <- pData(gset) |>
    as_tibble(rownames = "geo_accession_row") |>
    mutate(
        sepsis_group  = str_extract(title, "(?<=_)(Hyper|Hypo)$"),
        cnts_colnames = str_remove(title, "_(?<=_)(Hyper|Hypo)$")
    )

GSE236892_meta <- pdata |>
    dplyr::select(geo_accession, cnts_colnames, `age:ch1`, `gender:ch1`, `lca:ch1`) |>
    dplyr::rename(age = `age:ch1`, sex = `gender:ch1`, lca_label = `lca:ch1`)

# retrieve count matrix
supp_files <- getGEOSuppFiles("GSE236892", baseDir = tempdir(), fetch_files = TRUE)
cnt_path <- rownames(supp_files)[str_detect(rownames(supp_files), "cnt_data")]

cnt_mat <- read_csv(cnt_path) |>
    dplyr::rename(ensg = `...1`) |>
    column_to_rownames("ensg") |>
    dplyr::select(all_of(GSE236892_meta$cnts_colnames)) |>
    as.matrix()
storage.mode(cnt_mat) <- "integer"

col_data <- GSE236892_meta |>
    column_to_rownames("cnts_colnames") |>
    mutate(
        lca_label   = factor(lca_label, levels = c("Hypo", "Hyper")),
        sex         = factor(sex),
        # label columns which had non-numeric age
        imputed_age = age == "90+",
        # strip non-numeric characters and cast to integer
        age         = as.integer(str_remove(age, "\\+")),
        # scale age to avoid near collinearity with the intercept. Note
        # that this is only necessary or useful if age is used as a predictor
        age_scaled  = as.numeric(scale(age))
    )

gencode_genes <- suppressMessages(genes(txdb))
names(gencode_genes) <- str_remove(names(gencode_genes), "\\.\\d+$")

ensembl_meta <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = names(gencode_genes),
    columns = c("ENSEMBL", "SYMBOL", "GENENAME", "GENETYPE"),
    keytype = "ENSEMBL"
) |>
    as_tibble() |>
    filter(!is.na(ENSEMBL)) |>
    distinct(ENSEMBL, .keep_all = TRUE)

mcols(gencode_genes) <- ensembl_meta[
    match(names(gencode_genes), ensembl_meta$ENSEMBL),
    c("ENSEMBL", "SYMBOL", "GENENAME", "GENETYPE")
]

cnt_ensg <- rownames(cnt_mat)
mapped <- gencode_genes[names(gencode_genes) %in% cnt_ensg]
unmapped_ids <- setdiff(cnt_ensg, names(gencode_genes))

unmapped <- GRanges(
    seqnames = rep("chrUn", length(unmapped_ids)),
    ranges   = IRanges(start = seq_along(unmapped_ids), width = 1),
    ENSEMBL  = unmapped_ids,
    SYMBOL   = NA_character_,
    GENENAME = NA_character_,
    GENETYPE = NA_character_
)
names(unmapped) <- unmapped_ids

row_gr <- c(mapped, unmapped)[cnt_ensg]

mcols(row_gr)$autosomal_protein_coding <- (
    !is.na(mcols(row_gr)$GENETYPE) &
        mcols(row_gr)$GENETYPE == "protein-coding" &
        as.character(seqnames(row_gr)) %in% paste0("chr", 1:22)
)

mcols(row_gr)$hemoglobin_related <- names(row_gr) %in% names(hemoglobin_genes)

se <- SummarizedExperiment(
    assays    = list(counts = cnt_mat),
    colData   = col_data,
    rowRanges = row_gr
)

saveRDS(se, here("data/earli_summarizedexperiment.rds"))
