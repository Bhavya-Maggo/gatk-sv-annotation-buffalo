library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(txdbmaker)
library(GenomicFeatures)
library(GenomeInfoDbData)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)

vcf <- readRDS("c1.rds")
message("Total variants loaded: ", nrow(vcf))

# Check SV types
print(table(info(vcf)$SVTYPE))
stopifnot("SVTYPE" %in% names(info(vcf)))

# Convert VCF into breakpoint ranges.
sv_gr <- breakpointRanges(
  vcf,
  inferMissingBreakends = TRUE
)

sv_gr <- trim(sv_gr)

message("Total breakpoint ranges: ", length(sv_gr))
print(table(sv_gr$svtype))

# Check for out-of-bound coordinates again
sl <- seqlengths(sv_gr)[as.character(seqnames(sv_gr))]
oob <- sum(end(sv_gr) > sl, na.rm = TRUE)

stopifnot(oob == 0)
message("OOB check passed")

# Create TxDb
#txdb <- makeTxDbFromGFF(
#  file = "genes.gtf",
#  format = "gtf",
#  organism = "Bubalus bubalis",
#  taxonomyId = 89462
#)

#saveDb(txdb, "buffalo_txdb.sqlite")

# Load Gene Annotation
txdb <- loadDb("buffalo_txdb.sqlite")

genes_gr <- genes(txdb)
exons_gr <- exons(txdb)
transcripts_gr <- transcripts(txdb)

message("Genes loaded: ", length(genes_gr))

# Harmonise chromosomes
common <- intersect(seqlevels(sv_gr), seqlevels(genes_gr))
seqlevels(sv_gr, pruning.mode="coarse") <- common
seqlevels(genes_gr, pruning.mode="coarse") <- common
stopifnot(all(seqlevels(sv_gr) == seqlevels(genes_gr)))
message("Seqlevels harmonised")

# Import Gene Metadata From GTF
gtf <- import("genes.gtf")

gene_lookup <- unique(as.data.frame(mcols(gtf))[
  ,c("gene_id","gene","gene_biotype")
])

gene_name_map <- setNames(gene_lookup$gene, gene_lookup$gene_id)
gene_biotype_map <- setNames(gene_lookup$gene_biotype, gene_lookup$gene_id)

genes_gr$gene_name <- gene_name_map[genes_gr$gene_id]
genes_gr$gene_biotype <- gene_biotype_map[genes_gr$gene_id]

# Find SV–Gene Overlaps
hits <- findOverlaps(
  sv_gr,
  genes_gr,
  ignore.strand = TRUE
)

message("Total overlaps: ", length(hits))

# Find Multiple Gene Hits
collapse_hits <- function(sv, genes, hits, column){
  
  vals <- tapply(
    mcols(genes)[[column]][subjectHits(hits)],
    queryHits(hits),
    function(x) paste(unique(x), collapse=";")
  )
  
  out <- rep(NA_character_, length(sv))
  out[as.integer(names(vals))] <- vals
  
  out
}

# Annotate SVs
sv_gr$gene_id <- collapse_hits(sv_gr, genes_gr, hits, "gene_id")
sv_gr$gene_name <- collapse_hits(sv_gr, genes_gr, hits, "gene_name")
sv_gr$gene_biotype <- collapse_hits(sv_gr, genes_gr, hits, "gene_biotype")

# Export Table
sv_df <- as.data.frame(sv_gr)

list_cols <- sapply(sv_df, is.list)

sv_df[list_cols] <- lapply(
  sv_df[list_cols],
  function(x) sapply(x, function(v) paste(unlist(v), collapse=";"))
)

write.csv(
  sv_df,
  "annotated_SVs.csv",
  row.names=FALSE
)

# Summary Statistics
sv_explore <- sv_df %>%
  filter(!is.na(gene_name)) %>%
  separate_rows(gene_name, sep=";") %>%
  separate_rows(gene_biotype, sep=";")

message("Total SVs: ", nrow(sv_df))
message("SVs hitting genes: ", sum(!is.na(sv_df$gene_name)))
message("Intergenic SVs: ", sum(is.na(sv_df$gene_name)))

print(table(sv_df$svtype))
print(sort(table(sv_explore$gene_biotype), decreasing=TRUE))

# Top genes affected by structural variants
top_genes <- sort(table(sv_explore$gene_name), decreasing=TRUE)[1:20]

df <- data.frame(
  gene=names(top_genes),
  count=as.numeric(top_genes)
)

p1 <- ggplot(df, aes(reorder(gene,count),count)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(
    x="Gene",
    y="SV count",
    title="Top genes affected by structural variants"
  )
p1

# Gene Biotype Distribution plot
# Flatten multiple biotypes per SV
sv_explore <- sv_df %>%
  filter(!is.na(gene_biotype)) %>%
  separate_rows(gene_biotype, sep=";")

# Count per biotype
biotype_df <- sv_explore %>%
  count(gene_biotype) %>%
  arrange(desc(n))

# Plot
p2 <- ggplot(biotype_df, aes(x = reorder(gene_biotype, n), y = n)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=n), hjust=-0.1, size=3) +
  coord_flip() +
  labs(
    x = "Gene Biotype",
    y = "Number of SVs",
    title = "Gene Biotype Distribution of Structural Variants",
    subtitle = paste0("n = ", nrow(sv_df), " SVs")
  ) +
  theme_minimal(base_size = 12)

p2

# SV Type Distribution
# Count per SV type
svtype_df <- sv_df %>%
  count(svtype) %>%
  arrange(desc(n))

# Plot
p3 <- ggplot(svtype_df, aes(x = reorder(svtype, -n), y = n, fill=svtype)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=n), vjust=-0.5, size=3.5) +
  labs(
    x = "SV Type",
    y = "Count",
    title = "Structural Variant Type Distribution",
    subtitle = paste0("n = ", nrow(sv_df), " SVs")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

p3








