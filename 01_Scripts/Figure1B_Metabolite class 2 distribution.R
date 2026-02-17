# ClassII_PieChart.R

# 1. Libraries
library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(scales)

# 2. Parameters
file_path <- "C:/Users/gengd/Desktop/研一实验数据/Metabolome/20241225/MWXS-24-7906-a 3个样本植物非靶向代谢组结题报告.xlsx"
out_dir <- dirname(file_path)
OTHER_THRESHOLD <- 3  # Percentage threshold to group small categories into "Others"

# 3. Data Processing & QC
met_data <- read_excel(file_path, skip = 1)

# QC Calculation: Filter features with RSD <= 30%
met_data_filtered <- met_data %>%
  mutate(across(c(QC01, QC02, QC03), as.numeric)) %>%
  mutate(QC_mean = rowMeans(select(., QC01, QC02, QC03), na.rm = TRUE),
         QC_sd   = apply(select(., QC01, QC02, QC03), 1, sd, na.rm = TRUE),
         QC_RSD  = (QC_sd / QC_mean) * 100) %>%
  filter(QC_RSD <= 30)

# Calculate Average Abundance
met_data_filtered <- met_data_filtered %>%
  mutate(across(c(Callus_1, Callus_2, Callus_3), as.numeric)) %>%
  mutate(Callus_avg = rowMeans(select(., Callus_1, Callus_2, Callus_3), na.rm = TRUE))

# 4. Data Preparation for Plotting
classII_raw <- met_data_filtered %>%
  filter(!is.na(`Class II`)) %>%
  group_by(`Class II`) %>%
  summarise(Total = sum(Callus_avg, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Percent = Total / sum(Total) * 100)

# Group small categories (< 3%) into "Others"
classII_plot <- classII_raw %>%
  mutate(ClassII_mod = ifelse(Percent < OTHER_THRESHOLD, "Others", `Class II`)) %>%
  group_by(ClassII_mod) %>%
  summarise(Total = sum(Total)) %>%
  ungroup() %>%
  mutate(Percent = Total / sum(Total) * 100) %>%
  arrange(desc(Percent))

classII_plot$ClassII_mod <- factor(classII_plot$ClassII_mod, levels = classII_plot$ClassII_mod)

# 5. Visualization (Blue Gradient)
colors_palette <- colorRampPalette(c("#08519C", "#DEEBF7"))(nrow(classII_plot))

p_pie <- ggplot(classII_plot, aes(x = "", y = Percent, fill = ClassII_mod)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Distribution of Metabolite Classes (Class II)", fill = "Class (Percentage)") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(
    values = colors_palette, 
    labels = paste0(classII_plot$ClassII_mod, " (", sprintf("%.1f", classII_plot$Percent), "%)")
  )

# 6. Export
ggsave(file.path(out_dir, "ClassII_PieChart_Blue.png"), p_pie, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "ClassII_PieChart_Blue.pdf"), p_pie, width = 7, height = 6)