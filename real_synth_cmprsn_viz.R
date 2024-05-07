library(ggplot2)
library(ggcorrplot)
library(ggsci)
library(patchwork)
library(ggpubr)
library(dplyr)
library(tidyr)

###############################################################################
############### Data prep and summary statistics calculations ################
###############################################################################

df <- read.csv("Project/sampleData_digitalAg.csv")
df <- df[, -c(1:3)]

df2 <- read.csv("Project/synth_data.csv", row.names = 1)

df %>% select(-c(group, sb, afb)) %>% summary()
df2 %>% summary()

###############################################################################
###### Histograms of synthetically generated group, SB, and AFB values ########
###############################################################################

hist_simGroup <- ggplot(df2, aes(x=group)) + 
  geom_histogram(aes(y=after_stat(density)), 
                 alpha=0.5, 
                 position="identity", 
                 bins=75) +
  geom_density(alpha=0.2) + 
  geom_vline(aes(xintercept=0.55),
             color="black", 
             linetype="dashed", 
             size=1) +
  xlab("Simulated Group") +
  ylab("Density") +
  geom_segment(aes(x = 0.6, y = 2, xend = 1.2, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "#E64B35FF") +
  geom_segment(aes(x = 0.5, y = 2, xend = -0.1, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "#4DBBD5FF") +
  annotate("text", x=c(0.9, 0.2), y=2.05, label=c("VR", "VS"), color=c("#E64B35FF", "#4DBBD5FF"), size=3) +
  theme_bw()

hist_simSB <- ggplot(df2, aes(x=sb)) + 
  geom_histogram(aes(y=after_stat(density)), 
                 alpha=0.5, 
                 position="identity", 
                 bins=75) +
  geom_density(alpha=0.2) + 
  geom_vline(aes(xintercept=0.5),
             color="black", 
             linetype="dashed", 
             size=1) +
  xlab("Simulated SB") +
  ylab("Density") +
  geom_segment(aes(x = 0.6, y = 2, xend = 1.2, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "#E64B35FF") +
  geom_segment(aes(x = 0.4, y = 2, xend = -0.2, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "#4DBBD5FF") +
  annotate("text", x=c(0.9, 0.1), y=2.05, label=c("SB+", "SB-"), color=c("#E64B35FF", "#4DBBD5FF"), size=3) +
  theme_bw()

hist_simAFB <- ggplot(df2, aes(x=afb)) + 
  geom_histogram(aes(y=after_stat(density)), 
                 alpha=0.5, 
                 position="identity", 
                 bins=75) +
  geom_density(alpha=0.2) + 
  xlab("Simulated AFB") +
  ylab("Density") +
  geom_vline(aes(xintercept=0.55),
             color="black", 
             linetype="dashed", 
             size=1) +
  geom_segment(aes(x = 0.6, y = 2, xend = 1.2, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "#E64B35FF") +
  geom_segment(aes(x = 0.5, y = 2, xend = -0.1, yend = 2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "#4DBBD5FF") +
  annotate("text", x=c(0.9, 0.2), y=2.05, label=c("AFB+", "AFB-"), color=c("#E64B35FF", "#4DBBD5FF"), size=3) +
  theme_bw()

###############################################################################
### Disease prevalence comparisons between original and synthetic datasets ####
###############################################################################

df3 <- df2 %>%
  mutate(group_label=ifelse(group>0.55, 1, 0)) %>%
  mutate(sb_label=ifelse(sb>0.5, 1, 0)) %>%
  mutate(afb_label=ifelse(afb>0.55, 1, 0))

prev_real <- df %>% 
  summarise(vr=(sum(group)/n())*100,
            vs=(1-(sum(group)/n()))*100,
            sb_pos=(sum(sb)/n())*100,
            sb_neg=(1-(sum(sb)/n()))*100,
            afb_pos=(sum(afb)/n())*100,
            afb_neg=(1-(sum(afb)/n()))*100) %>%
  pivot_longer(cols = everything(),
               names_to = "var", 
               values_to = "value") %>%
  mutate(dataset="original")

prev_sim <- df3 %>% 
  summarise(vr=(sum(group_label)/n())*100,
            vs=(1-(sum(group_label)/n()))*100,
            sb_pos=(sum(sb_label)/n())*100,
            sb_neg=(1-(sum(sb_label)/n()))*100,
            afb_pos=(sum(afb_label)/n())*100,
            afb_neg=(1-(sum(afb_label)/n()))*100) %>%
  pivot_longer(cols = everything(),
               names_to = "var", 
               values_to = "value") %>%
  mutate(dataset="synthetic")

prev_combined <- rbind(prev_real, prev_sim)

prev_barp <- ggplot(prev_combined, aes(y=value, x=var, fill=dataset)) + 
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_npg() +
  ylab("Proprotion of total observations (%)") +
  xlab("") +
  scale_x_discrete(labels=c("vr"="VR", "vs"="VS", "sb_pos"="SB+", "sb_neg"="SB-", "afb_pos"="AFB+","afb_neg"="AFB-")) +
  theme_bw()

png("Project/real_synth_discrete_comparison.png", width = 4000, height = 3000, res = 300)
(hist_simGroup + hist_simSB) / (hist_simAFB + prev_barp) + plot_annotation(tag_levels = "A")
dev.off()

###############################################################################
############################ Correlation matrices #############################
###############################################################################

corr_labels <- c("Varroa", "Nosema", "ABPV", "BQCV", "DWV", "KV", "VDV1")

og_corr_r <- df %>%
  select(varroa_count, nosema, abpv, bqcv, dwv, kv, vdv1) %>%
  cor(method = "pearson")

colnames(og_corr_r) <- corr_labels
rownames(og_corr_r) <- corr_labels

og_corr_p <- df %>%
  select(varroa_count, nosema, abpv, bqcv, dwv, kv, vdv1) %>%
  cor_pmat()

colnames(og_corr_p) <- corr_labels
rownames(og_corr_p) <- corr_labels

synth_corr_r <- df2 %>%
  select(varroa_count, nosema, abpv, bqcv, dwv, kv, vdv1) %>%
  cor(method = "pearson")

colnames(synth_corr_r) <- corr_labels
rownames(synth_corr_r) <- corr_labels

synth_corr_p <- df2 %>%
  select(varroa_count, nosema, abpv, bqcv, dwv, kv, vdv1) %>%
  cor_pmat()

colnames(synth_corr_p) <- corr_labels
rownames(synth_corr_p) <- corr_labels

og_cor_plot <- ggcorrplot(og_corr_r, outline.col = "white",
                          lab = TRUE, p.mat = og_corr_p, insig = "blank",
                          colors = c(pal_npg("nrc")(2)[2], "white", pal_npg("nrc")(2)[1]))

synth_cor_plot <- ggcorrplot(synth_corr_r, outline.col = "white",
                             lab = TRUE, p.mat = synth_corr_p, insig = "blank",
                             colors = c(pal_npg("nrc")(2)[2], "white", pal_npg("nrc")(2)[1]))

png("Project/real_synth_corr_comparison.png", width = 3000, height = 2000, res = 300)
og_cor_plot + synth_cor_plot + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect")
dev.off()
