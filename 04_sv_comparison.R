# Load libraries--------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggvenn) 

# Load annotated SV tables ---------------------------------------------------
sv1 <- read.csv("annotated_SVs_d1.csv", stringsAsFactors = FALSE)
sv2 <- read.csv("annotated_SVs_c1.csv", stringsAsFactors = FALSE)

sv1$sample <- "Donor"
sv2$sample <- "Clone"

message("Donor SVs: ", nrow(sv1))
message("Clone SVs: ", nrow(sv2))

# Generate unique SV identifiers ---------------------------------------------
sv1 <- sv1 %>%
  mutate(SV_ID = paste(seqnames, start, end, svtype, sep = "_"))

sv2 <- sv2 %>%
  mutate(SV_ID = paste(seqnames, start, end, svtype, sep = "_"))

# Identify shared and unique structural variants -----------------------------
shared_ids <- intersect(sv1$SV_ID, sv2$SV_ID)
unique_sv1 <- setdiff(sv1$SV_ID, sv2$SV_ID)
unique_sv2 <- setdiff(sv2$SV_ID, sv1$SV_ID)

message("Shared SVs:     ", length(shared_ids))
message("Donor only SVs: ", length(unique_sv1))
message("Clone only SVs: ", length(unique_sv2))

# Subset SV rows -------------------------------------------------------------
shared_sv  <- sv1 %>% dplyr::filter(SV_ID %in% shared_ids)
donor_only <- sv1 %>% dplyr::filter(SV_ID %in% unique_sv1)
clone_only <- sv2 %>% dplyr::filter(SV_ID %in% unique_sv2)

# Optional: export SV comparison tables --------------------------------------
#write.csv(shared_sv,  "comparison_outputs/shared_SVs.csv", row.names = FALSE)
#write.csv(donor_only, "comparison_outputs/donor_only_SVs.csv", row.names = FALSE)
#write.csv(clone_only, "comparison_outputs/clone_only_SVs.csv", row.names = FALSE)

# Summary statistics ---------------------------------------------------------
summary_table <- data.frame(
  Metric       = c("Total SVs", "Genic SVs", "Intergenic SVs",
                   "Shared SVs", "Unique SVs"),
  Donor        = c(nrow(sv1),
                   sum(!is.na(sv1$gene_name)),
                   sum(is.na(sv1$gene_name)),
                   length(shared_ids),
                   length(unique_sv1)),
  Clone        = c(nrow(sv2),
                   sum(!is.na(sv2$gene_name)),
                   sum(is.na(sv2$gene_name)),
                   length(shared_ids),
                   length(unique_sv2))
)

print(summary_table)

# write.csv(summary_table, "comparison_outputs/summary_table.csv", row.names = FALSE)

# SV type distribution -------------------------------------------------------
message("Donor SV types:")
print(table(sv1$svtype))

message("Clone SV types:")
print(table(sv2$svtype))

message("Shared SV types:")
print(table(shared_sv$svtype))

# Gene level comparison-------------------------------------------------------
# Expand rows where multiple genes are annotated per SV
explode_genes <- function(df) {
  df %>%
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    tidyr::separate_rows(gene_name,    sep = ";") %>%
    tidyr::separate_rows(gene_biotype, sep = ";") %>%
    tidyr::separate_rows(gene_id,      sep = ";")
}

sv1_exp <- explode_genes(sv1)
sv2_exp <- explode_genes(sv2)

# Identify shared and unique genes -------------------------------------------
shared_genes     <- intersect(unique(sv1_exp$gene_name), unique(sv2_exp$gene_name))

unique_gene_sv1  <- setdiff(unique(sv1_exp$gene_name),   unique(sv2_exp$gene_name))

unique_gene_sv2  <- setdiff(unique(sv2_exp$gene_name),   unique(sv1_exp$gene_name))

message("Shared genes:          ", length(shared_genes))
message("Donor only genes:      ", length(unique_gene_sv1))
message("Clone only genes:      ", length(unique_gene_sv2))
message("Donor unique genes:    ", length(unique(sv1_exp$gene_name)))
message("Clone unique genes:    ", length(unique(sv2_exp$gene_name)))

# Optional: export gene comparison tables ------------------------------------
#write.csv(data.frame(gene = shared_genes),"shared_genes.csv", row.names = FALSE)

#write.csv(data.frame(gene = unique_gene_sv1), "donor_only_genes.csv", row.names = FALSE)

#write.csv(data.frame(gene = unique_gene_sv2), "clone_only_genes.csv", row.names = FALSE)

# Plot 1: SV type comparison -------------------------------------------------
combined <- bind_rows(sv1, sv2)

svtype_df <- combined %>%
  dplyr::count(sample, svtype) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(pct = round(n / sum(n) * 100, 1))

p1 <- ggplot(svtype_df, aes(x = svtype, y = n, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
            position = position_dodge(width = 0.9),
            vjust = -0.2, size = 2.8) +
  scale_fill_manual(values = c("Donor" = "#2196F3", "Clone" = "#FF5722")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top") +
  labs(x = "SV Type", y = "Count", fill = "",
       title    = "SV Type Distribution — Donor vs Clone",
       subtitle = paste0("Donor n=", nrow(sv1), " | Clone n=", nrow(sv2)))
p1

# Plot 2: Shared vs unique SV counts -----------------------------------------
overlap_df <- data.frame(
  Category = factor(c("Donor only", "Shared", "Clone only"),
                    levels = c("Donor only", "Shared", "Clone only")),
  Count    = c(length(unique_sv1), length(shared_ids), length(unique_sv2))
)

p2 <- ggplot(overlap_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.4, size = 4.5) +
  scale_fill_manual(values = c(
    "Donor only" = "#2196F3",
    "Shared"     = "#9C27B0",
    "Clone only" = "#FF5722")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "", y = "Number of SVs",
       title = "Shared vs Unique SVs — Donor vs Clone")

p2

# Plot 3: Venn diagram of SV overlap -----------------------------------------
p3 <- ggvenn(
  list(Donor = sv1$SV_ID, Clone = sv2$SV_ID),
  fill_color  = c("#2196F3", "#FF5722"),
  stroke_size = 0.5,
  set_name_size = 5
) + labs(title = "SV Overlap — Donor vs Clone")

p3

# Plot 4: Venn diagram of gene overlap ---------------------------------------
p4 <- ggvenn(
  list(Donor = unique(sv1_exp$gene_name), Clone = unique(sv2_exp$gene_name)),
  fill_color  = c("#2196F3", "#FF5722"),
  stroke_size = 0.5,
  set_name_size = 5
) + labs(title = "Gene Overlap — Donor vs Clone")

p4

# Plot 5: Chromosome distribution -------------------------------------------- chr_df <- combined %>% dplyr::count(sample, seqnames)

p5 <- ggplot(chr_df, aes(x = seqnames, y = n, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Donor" = "#2196F3", "Clone" = "#FF5722")) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "top") +
  labs(x = "Chromosome", y = "Number of SVs", fill = "",
       title = "Chromosome Distribution — Donor vs Clone")

p5

# Plot 6: Gene biotype comparison --------------------------------------------
biotype_df <- bind_rows(sv1_exp, sv2_exp) %>%
  dplyr::count(sample, gene_biotype)

p6 <- ggplot(biotype_df, aes(x = reorder(gene_biotype, n), y = n, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(0.9), hjust = -0.2, size = 2.8) +
  coord_flip() +
  scale_fill_manual(values = c("Donor" = "#2196F3", "Clone" = "#FF5722")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top") +
  labs(x = "Gene Biotype", y = "Number of SVs", fill = "",
       title = "Gene Biotype Distribution — Donor vs Clone")

p6

# Plot 7: Top genes affected by SVs ------------------------------------------
top_combined <- bind_rows(
  sv1_exp %>% dplyr::count(gene_name, sort = TRUE) %>% head(20) %>% mutate(sample = "Donor"),
  sv2_exp %>% dplyr::count(gene_name, sort = TRUE) %>% head(20) %>% mutate(sample = "Clone")
)

p7 <- ggplot(top_combined, aes(x = reorder(gene_name, n), y = n, fill = sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), hjust = -0.2, size = 2.8) +
  coord_flip() +
  facet_wrap(~sample, scales = "free") +
  scale_fill_manual(values = c("Donor" = "#2196F3", "Clone" = "#FF5722")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = "Gene", y = "Number of SVs",
       title = "Top 20 Most-Hit Genes — Donor vs Clone")

p7
