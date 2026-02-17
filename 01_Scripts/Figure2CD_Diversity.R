# Microbiome_Analysis_Pipeline.R

# 1. Load Libraries
packages <- c("dplyr", "tidyr", "tibble", "vegan", "car")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
lapply(packages, library, character.only = TRUE)

# 2. Data Loading
# Ensure input files are in the working directory
feature_file <- "exported-table/feature-table.tsv"
taxonomy_file <- "exported-taxonomy/taxonomy.tsv"

if (!file.exists(feature_file) || !file.exists(taxonomy_file)) stop("Input files missing!")

feature_table <- read.delim(feature_file, sep = "\t", skip = 1, header = TRUE, row.names = 1, check.names = FALSE)
taxonomy_raw <- read.delim(taxonomy_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 3. Data Processing
# Parse Taxonomy
taxonomy_clean <- taxonomy_raw %>%
  separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  mutate(across(Kingdom:Species, ~ sub("^.__", "", .x))) %>%
  column_to_rownames("Feature ID")

# Calculate Relative Abundance & Dominant Taxa
otu_mat <- as.matrix(feature_table)
otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), "/")

df_top <- as.data.frame(otu_rel) %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "Sample", values_to = "RelAbundance") %>%
  group_by(Sample) %>%
  slice_max(RelAbundance, n = 1) %>%
  ungroup()

# Annotate Taxa
df_annotated <- df_top %>%
  left_join(rownames_to_column(taxonomy_clean, "ASV"), by = "ASV") %>%
  mutate(Best_Taxon = coalesce(na_if(Genus, ""), na_if(Family, ""), na_if(Order, ""), "Unclassified")) %>%
  filter(Best_Taxon != "Unclassified")

# 4. Grouping & Matrix Construction
# Parse Sample Name (Format: R10306 -> Treatment:1, Plate:03)
df_parsed <- df_annotated %>%
  mutate(
    Sample_Clean = sub("^R", "", Sample),
    Group = ifelse(substr(Sample_Clean, 1, 1) == "1", "Co_Callus", "Control"),
    Plate = substr(Sample_Clean, 2, 3),
    Sample_ID = paste0(Group, "_", as.numeric(Plate))
  )

# Aggregate Count Table
freq_mat <- df_parsed %>%
  count(Sample_ID, Best_Taxon) %>%
  pivot_wider(names_from = Best_Taxon, values_from = n, values_fill = 0) %>%
  column_to_rownames("Sample_ID") %>%
  as.matrix()

# 5. Alpha Diversity (Shannon)
alpha_df <- data.frame(Sample_ID = rownames(freq_mat), Shannon = diversity(freq_mat, "shannon")) %>%
  mutate(Group = sub("_[0-9]+$", "", Sample_ID))

# Normality & Homogeneity Test
shapiro_p <- min(tapply(alpha_df$Shannon, alpha_df$Group, function(x) shapiro.test(x)$p.value))
levene_p <- car::leveneTest(Shannon ~ Group, data = alpha_df)$`Pr(>F)`[1]

# Statistical Test Selection
if (shapiro_p > 0.05) {
  test_res <- t.test(Shannon ~ Group, data = alpha_df, var.equal = (levene_p > 0.05))
  test_name <- ifelse(levene_p > 0.05, "Student's t-test", "Welch's t-test")
} else {
  test_res <- wilcox.test(Shannon ~ Group, data = alpha_df)
  test_name <- "Wilcoxon Rank Sum test"
}

# Export Alpha Stats
write.csv(data.frame(Metric="Shannon", Test=test_name, P_value=test_res$p.value), "Alpha_Stats.csv", row.names = FALSE)

# 6. Beta Diversity (Bray-Curtis)
dist_bc <- vegdist(freq_mat, method = "bray")
metadata <- data.frame(Sample_ID = rownames(freq_mat)) %>%
  mutate(Group = sub("_[0-9]+$", "", Sample_ID))

# PERMANOVA & PERMDISP
permanova <- adonis2(dist_bc ~ Group, data = metadata, permutations = 999)
permdisp <- anova(betadisper(dist_bc, metadata$Group))

# Export Beta Stats
beta_stats <- data.frame(
  Test = c("PERMANOVA", "PERMDISP"),
  F_Statistic = c(permanova$F[1], permdisp$`F value`[1]),
  R2 = c(permanova$R2[1], NA),
  P_value = c(permanova$`Pr(>F)`[1], permdisp$`Pr(>F)`[1])
)
