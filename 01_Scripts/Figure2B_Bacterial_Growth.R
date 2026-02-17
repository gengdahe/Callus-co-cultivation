# Script: 02_Growth_Analysis_Plot.R

# 1. Libraries
library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

# 2. Parameters & Data Loading
file_path <- "C:/Users/gengd/Desktop/研一实验数据/96孔板共培养/数据/Bacterial_Growth_RawData.csv"
out_dir <- dirname(file_path)

data_all <- read.csv(file_path)

# Set Factor Levels (Crucial for plotting order)
data_all$Rep <- factor(data_all$Rep, levels = c("Rep1", "Rep2", "Rep3"))
data_all$Treatment <- factor(data_all$Treatment, levels = c("Control", "Co_cultivation"))

# Colors: Blue for Control, Red for Co-cultivation
col_ctrl <- "#6699CC"; col_exp <- "#961D4C"
my_colors <- c("Control" = col_ctrl, "Co_cultivation" = col_exp)

# 3. Statistical Analysis
# 3.1 Calculate Mean & Growth Rate
stats_df <- data_all %>%
  group_by(Rep, Treatment) %>%
  summarise(Mean = mean(Value), Max = max(Value), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = c(Mean, Max)) %>%
  mutate(
    Increase_Pct = (Mean_Co_cultivation - Mean_Control) / Mean_Control * 100,
    Label_Pct = paste0("+", sprintf("%.1f", Increase_Pct), "%")
  )

avg_growth_rate <- mean(stats_df$Increase_Pct)

# 3.2 T-test for Significance
p_vals <- data_all %>%
  group_by(Rep) %>%
  summarise(p_val = t.test(Value ~ Treatment)$p.value, .groups = "drop") %>%
  mutate(sig_label = case_when(
    p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ "ns"
  ))

# Merge for plotting annotations
plot_anno <- left_join(stats_df, p_vals, by = "Rep") %>%
  mutate(y_star = pmax(Max_Control, Max_Co_cultivation) + 1.5)

# 4. Visualization
# 4.1 Top Plot: Distribution (Boxplot + Jitter)
p_top <- ggplot(data_all, aes(x = Rep, y = Value, fill = Treatment)) +
  geom_violin(position = position_dodge(0.8), width = 0.7, alpha = 0.3, color = NA, trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.15, outlier.shape = NA, alpha = 0.8, size = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
             shape = 21, color = "white", size = 1.8, stroke = 0.3, show.legend = FALSE) +
  geom_text(data = plot_anno, aes(x = Rep, y = y_star, label = sig_label), 
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Bacterial Growth (Counts)", y = "Counts", x = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = c(0.12, 0.92),
    legend.background = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(), plot.margin = margin(b = 2)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# 4.2 Bottom Plot: Growth Rate (Barplot)
p_bottom <- ggplot(plot_anno, aes(x = Rep, y = Increase_Pct)) +
  geom_col(width = 0.5, fill = "#F39C12", alpha = 0.85) +
  geom_hline(yintercept = avg_growth_rate, linetype = "dashed", color = "#D35400", linewidth = 0.8) +
  annotate("text", x = 0.5, y = avg_growth_rate + 2, 
           label = paste0("Avg Increase: ", sprintf("%.1f", avg_growth_rate), "%"), 
           hjust = 0, color = "#D35400", fontface = "italic", size = 4) +
  geom_text(aes(label = Label_Pct), vjust = -0.5, fontface = "bold", size = 4.5) +
  scale_y_continuous(limits = c(0, max(plot_anno$Increase_Pct) * 1.35), expand = c(0,0)) +
  labs(y = "Growth Ratio (%)", x = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), plot.margin = margin(t = 2),
    axis.title.y = element_text(size = 12), axis.text.x = element_text(face = "bold", size = 12)
  )

# 4.3 Combine Plots
final_plot <- (p_top / p_bottom) + plot_layout(heights = c(2, 1))

# 5. Export
ggsave(file.path(out_dir, "Corrected_Bacterial_Growth.pdf"), final_plot, width = 7, height = 8)
ggsave(file.path(out_dir, "Corrected_Bacterial_Growth.png"), final_plot, width = 7, height = 8, dpi = 600)

print(final_plot)