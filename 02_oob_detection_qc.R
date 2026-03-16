# =============================================================================
# 00_qc_oob_diagnosis.R — UPDATED
# chrMT excluded — no RefSeq for MT in GCF_019923935.1
# =============================================================================

library(VariantAnnotation)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomeInfoDb)

# ── CHANGE ONLY THESE 2 LINES ─────────────────────────────────────────────
vcf_file  <- "high_quality_4998.vcf"  
rds_out   <- "c2.rds"               
# ──────────────────────────────────────────────────────────────────────────

# Seqinfo — 25 nuclear chromosomes only (no MT)
new_seqinfo <- Seqinfo(
  seqnames = c(
    "NC_059157.1", "NC_059158.1", "NC_059159.1", "NC_059160.1",
    "NC_059161.1", "NC_059162.1", "NC_059163.1", "NC_059164.1",
    "NC_059165.1", "NC_059166.1", "NC_059167.1", "NC_059168.1",
    "NC_059169.1", "NC_059170.1", "NC_059171.1", "NC_059172.1",
    "NC_059173.1", "NC_059174.1", "NC_059175.1", "NC_059176.1",
    "NC_059177.1", "NC_059178.1", "NC_059179.1", "NC_059180.1",
    "NC_059181.1"
  ),
  seqlengths = c(
    "NC_059157.1" = 202348575, "NC_059158.1" = 188164321,
    "NC_059159.1" = 174872328, "NC_059160.1" = 164971372,
    "NC_059161.1" = 132500020, "NC_059162.1" = 120418900,
    "NC_059163.1" = 116997125, "NC_059164.1" = 119318788,
    "NC_059165.1" = 110262714, "NC_059166.1" = 104551495,
    "NC_059167.1" = 102416932, "NC_059168.1" = 106386409,
    "NC_059169.1" = 89364132,  "NC_059170.1" = 83479656,
    "NC_059171.1" = 81832314,  "NC_059172.1" = 85115218,
    "NC_059173.1" = 72603365,  "NC_059174.1" = 65791760,
    "NC_059175.1" = 71582194,  "NC_059176.1" = 69611172,
    "NC_059177.1" = 60624307,  "NC_059178.1" = 61771150,
    "NC_059179.1" = 52095120,  "NC_059180.1" = 42172481,
    "NC_059181.1" = 143192433
  ),
  isCircular = rep(FALSE, 25),
  genome     = "GCF_019923935.1"
)

chr_to_nc <- c(
  "chr1"  = "NC_059157.1", "chr2"  = "NC_059158.1",
  "chr3"  = "NC_059159.1", "chr4"  = "NC_059160.1",
  "chr5"  = "NC_059161.1", "chr6"  = "NC_059162.1",
  "chr7"  = "NC_059163.1", "chr8"  = "NC_059164.1",
  "chr9"  = "NC_059165.1", "chr10" = "NC_059166.1",
  "chr11" = "NC_059167.1", "chr12" = "NC_059168.1",
  "chr13" = "NC_059169.1", "chr14" = "NC_059170.1",
  "chr15" = "NC_059171.1", "chr16" = "NC_059172.1",
  "chr17" = "NC_059173.1", "chr18" = "NC_059174.1",
  "chr19" = "NC_059175.1", "chr20" = "NC_059176.1",
  "chr21" = "NC_059177.1", "chr22" = "NC_059178.1",
  "chr23" = "NC_059179.1", "chr24" = "NC_059180.1",
  "chrX"  = "NC_059181.1"
  # chrMT excluded — no RefSeq ID in GCF_019923935.1
)

# Load VCF
vcf <- readVcf(vcf_file, genome = "GCF_019923935.1")
message("Raw SVs loaded: ", nrow(vcf))

# Step 1 — MT hatao
seqlevels(vcf, pruning.mode = "coarse") <- grep(
  "chrMT|NC_049568.1", seqlevels(vcf), value = TRUE, invert = TRUE
)
message("After MT removal: ", nrow(vcf), " SVs")

# Step 2 — rename if chr style
if (any(grepl("^chr", seqlevels(vcf)))) {
  seqlevels(vcf) <- chr_to_nc[seqlevels(vcf)]
  message("Renamed chr → NC_")
}

print(seqlevels(vcf))

# Step 3 — seqinfo assign
seqinfo(vcf) <- new_seqinfo

# Step 4 — OOB diagnosis
vcf_ranges <- rowRanges(vcf)
sl_raw     <- seqlengths(vcf_ranges)[as.character(seqnames(vcf_ranges))]
oob_raw    <- end(vcf_ranges) > sl_raw
message("Total OOB: ", sum(oob_raw, na.rm = TRUE))

# Agar sirf BND hain ya OOB = 0 → seedha save karo
bnd_only <- sum(oob_raw, na.rm = TRUE) == 0 ||
  all(info(vcf[oob_raw])$SVTYPE == "BND")

if (bnd_only) {
  message("OOB = 0 or BND only → no cleaning needed → saving directly")
  saveRDS(vcf, rds_out)
  message("Saved: ", rds_out)
} else {
  # OOB by SVTYPE
  message("OOB by SVTYPE:")
  print(table(info(vcf[oob_raw])$SVTYPE))
  message("OOB by chromosome:")
  print(sort(table(as.character(seqnames(vcf_ranges[oob_raw]))), decreasing = TRUE))
  message("→ Investigate further before saving")
}

