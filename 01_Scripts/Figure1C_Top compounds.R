# Figure1C_TopCompounds.R

# 1. Libraries
library(readxl); library(dplyr); library(ggplot2); library(writexl)

# 2. Parameters
file_path <- "C:/Users/gengd/Desktop/研一实验数据/Metabolome/20241225/MWXS-24-7906-a 3个样本植物非靶向代谢组结题报告.xlsx"
out_dir <- dirname(file_path)

# 3. Data Processing & QC
met_data <- read_excel(file_path, skip = 1) %>%
  mutate(across(c(QC01, QC02, QC03, Callus_1, Callus_2, Callus_3), as.numeric)) %>%
  mutate(QC_mean = rowMeans(select(., QC01, QC02, QC03), na.rm = TRUE),
         QC_sd = apply(select(., QC01, QC02, QC03), 1, sd, na.rm = TRUE),
         QC_RSD = (QC_sd / QC_mean) * 100,
         Callus_avg = rowMeans(select(., Callus_1, Callus_2, Callus_3), na.rm = TRUE),
         Callus_sd = apply(select(., Callus_1, Callus_2, Callus_3), 1, sd, na.rm = TRUE)) %>%
  filter(QC_RSD <= 30 | is.na(QC_RSD))

# 4. Data Preparation
total_abundance <- sum(met_data$Callus_avg, na.rm = TRUE)

abundance_table <- met_data %>%
  mutate(Relative_Abundance = (Callus_avg / total_abundance) * 100) %>%
  arrange(desc(Relative_Abundance))

top10_table <- abundance_table %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number(), across(c(Callus_avg, Callus_sd, QC_RSD, Relative_Abundance), ~round(., 2))) %>%
  select(Rank, Compounds, Callus_1, Callus_2, Callus_3, Callus_avg, Callus_sd, QC_RSD, Relative_Abundance)

write_xlsx(top10_table, file.path(out_dir, "Top10_Compounds_Statistics.xlsx"))
write.csv(top10_table, file.path(out_dir, "Top10_Compounds_Statistics.csv"), row.names = FALSE)

plot_data <- abundance_table %>%
  filter(Relative_Abundance > 1) %>%
  mutate(Compounds = factor(Compounds, levels = rev(Compounds)))

# 5. Visualization
p_abundance <- ggplot(plot_data, aes(x = Compounds, y = Relative_Abundance)) +
  geom_bar(stat = "identity", fill = "#5B9BD5", width = 0.7) +
  coord_flip() +
  geom_text(aes(label = sprintf("%.1f%%", Relative_Abundance)), hjust = 1.2, size = 4, color = "black", fontface = "bold") +
  labs(title = "Top compounds (Relative abundance > 1%)", x = NULL, y = "Relative abundance (%)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.line = element_line(linewidth = 1)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 6. Export
ggsave(file.path(out_dir, "Figure1C_Top Compounds.pdf"), plot = p_abundance, width = 8, height = 6)
