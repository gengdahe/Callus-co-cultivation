# Wilcoxon_Pathway_Analysis.R

# 1. Libraries
packages <- c("readr", "dplyr", "tidyr", "stringr", "KEGGREST")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "KEGGREST") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("KEGGREST")
    } else install.packages(pkg)
  }
}
library(readr); library(dplyr); library(tidyr); library(stringr); library(KEGGREST)

# 2. Parameters
file_in <- "C:/Users/gengd/Desktop/研一实验数据/Metabolome/Callus植物非靶向代谢组数据n=3.csv"
out_dir <- "C:/Users/gengd/Desktop/Metabolome/Optimized_Results"
MIN_PATHWAY_SIZE <- 15 

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 3. Data Processing
read_data_robust <- function(path) {
  tryCatch({ read_csv(path, col_types = cols(.default = "c"), skip = 1) }, 
           error = function(e) stop("Error reading file"))
}

dat_raw <- read_data_robust(file_in)
colnames(dat_raw) <- gsub("[[:cntrl:]]|\\s+", "", colnames(dat_raw))

sample_cols <- colnames(dat_raw)[grepl("^Callus", colnames(dat_raw), ignore.case = TRUE)]
kegg_col <- colnames(dat_raw)[grepl("kegg_map|kegg", colnames(dat_raw), ignore.case = TRUE)][1]
name_col <- colnames(dat_raw)[1]

if (length(sample_cols) == 0 || is.na(kegg_col)) stop("Missing sample or KEGG columns.")

dat_processed <- dat_raw %>%
  mutate(across(all_of(sample_cols), ~as.numeric(as.character(.)))) %>%
  rowwise() %>%
  mutate(log_abundance = mean(log10(c_across(all_of(sample_cols)) + 1), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(log_abundance > 0, !is.na(!!sym(kegg_col)), !!sym(kegg_col) != "", !!sym(kegg_col) != "-")

# 4. Pathway Mapping & Background
pathway_map <- dat_processed %>%
  select(Feature = all_of(name_col), Abundance = log_abundance, Raw_KEGG = all_of(kegg_col)) %>%
  mutate(PathwayID = str_split(Raw_KEGG, ";|,")) %>%
  unnest(PathwayID) %>%
  mutate(PathwayID = trimws(PathwayID), MapID = str_replace(PathwayID, "^ko", "map")) %>%
  filter(str_detect(MapID, "^map\\d+")) %>%
  distinct(Feature, MapID, .keep_all = TRUE)

background_abundance <- unique(pathway_map$Abundance)

# 5. Wilcoxon Test
stats_list <- list()
unique_pathways <- unique(pathway_map$MapID)

for (pw in unique_pathways) {
  pw_abundances <- pathway_map %>% filter(MapID == pw) %>% pull(Abundance)
  
  if (length(pw_abundances) >= MIN_PATHWAY_SIZE) {
    wt <- wilcox.test(pw_abundances, background_abundance, alternative = "greater", exact = FALSE)
    stats_list[[pw]] <- data.frame(PathwayID = pw, Detected_Count = length(pw_abundances), P_value = wt$p.value)
  }
}

df_stats <- bind_rows(stats_list)

if (nrow(df_stats) > 0) {
  df_stats$FDR <- p.adjust(df_stats$P_value, method = "BH")
  
  # 6. KEGG Metadata Retrieval
  kegg_names <- keggList("pathway")
  names(kegg_names) <- str_replace(names(kegg_names), "path:", "")
  
  final_meta <- list()
  for (i in 1:nrow(df_stats)) {
    pid <- df_stats$PathwayID[i]
    p_name <- if (pid %in% names(kegg_names)) str_split(kegg_names[[pid]], " - ")[[1]][1] else "Unknown"
    
    total_cpds <- NA
    try({ Sys.sleep(0.1); links <- keggLink("compound", paste0("path:", pid)); total_cpds <- length(unique(links)) }, silent = TRUE)
    
    final_meta[[i]] <- data.frame(PathwayID = pid, PathwayName = p_name, Total_N = total_cpds)
  }
  
  # 7. Export
  final_df <- df_stats %>%
    left_join(bind_rows(final_meta), by = "PathwayID") %>%
    arrange(P_value) %>%
    mutate(Coverage_Pct = round(Detected_Count / Total_N * 100, 2),
           Significance = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "ns")) %>%
    select(PathwayID, PathwayName, Detected_Count, Total_N, Coverage_Pct, P_value, FDR, Significance)
  
  write_csv(final_df, file.path(out_dir, "Supplementary_Table_Strict_Wilcoxon.csv"))
}