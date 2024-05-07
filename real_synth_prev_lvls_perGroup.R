library(ggpubr)
library(ggsci)
library(patchwork)
library(dplyr)

###############################################################################
############################### data preparation ##############################
###############################################################################

df <- read.csv("Project/sampleData_digitalAg.csv")
df <- df[, -c(1:3)] %>%
  mutate(dataset="original")

df2 <- read.csv("Project/synth_data.csv", row.names = 1)
df3 <- df2 %>%
  mutate(group_label=ifelse(group>0.55, 1, 0), .keep="unused", .before = "varroa_count") %>%
  mutate(sb_label=ifelse(sb>0.5, 1, 0), .keep="unused", .after = "varroa_count") %>%
  mutate(afb_label=ifelse(afb>0.55, 1, 0), .keep="unused", .after = "varroa_count") %>%
  mutate(dataset="synthetic")
colnames(df3)[1:4] <- c("group", "varroa_count", "sb", "afb")

df_all <- rbind(df, df3)

sb_afb_prev <- df_all %>% 
  group_by(group, dataset) %>%
  summarise(prev_sb=(sum(sb)/n())*100,
            prev_afb=(sum(afb)/n())*100) %>%
  pivot_longer(cols = starts_with("prev_"),
               names_to = "disease", 
               names_prefix = "prev_",
               values_to = "prev")

###############################################################################
############# plotting prevalence or levels per group per dataset #############
###############################################################################

# prevalence of SB and AFB per group per dataset
sb_afb_barp <- ggplot(sb_afb_prev, aes(y=prev, x=disease, fill=as.factor(group))) + 
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_npg(labels=c("0"="VS", "1"="VR")) +
  scale_x_discrete(labels=c("afb"="AFB", "sb"="SB")) +
  facet_grid(~dataset) +
  ylab("Prevalence (%)") +
  xlab("") +
  theme_bw() +
  theme(legend.title = element_blank())

# boxplot of Varroa counts per group per dataset
vr_bp <- ggboxplot(df_all, "group", "varroa_count", 
                   ylab = "Varroa mite count", 
                   color = "group", 
                   facet.by = "dataset",
                   palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, label.y=68, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# boxplot of Nosema levels per group per dataset
nos_bp <- ggboxplot(df_all, "group", "nosema", 
                    ylab = "N. ceranae levels", 
                    color = "group", 
                    facet.by = "dataset",
                    palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, label.y=1.45, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# boxplot of ABPV levels per group per dataset
abpv_bp <- ggboxplot(df_all, "group", "abpv", 
                     ylab = "ABPV levels", 
                     color = "group", 
                     facet.by = "dataset",
                     palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, label.y=0.033, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# boxplot of BQCV levels per group per dataset
bqcv_bp <- ggboxplot(df_all, "group", "bqcv", 
                     ylab = "BQCV levels", 
                     color = "group", 
                     facet.by = "dataset",
                     palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, label.y=0.032, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# boxplot of DWV levels per group per dataset
dwv_bp <- ggboxplot(df_all, "group", "dwv", 
                    ylab = "DWV levels", 
                    color = "group", 
                    facet.by = "dataset",
                    palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, label.y=0.43, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# boxplot of KV levels per group per dataset
kv_bp <- ggboxplot(df_all, "group", "kv", 
                   ylab = "KV levels", 
                   color = "group", 
                   facet.by = "dataset",
                   palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, label.y=0.19, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# boxplot of VDV1 levels per group per dataset
vdv1_bp <- ggboxplot(df_all, "group", "vdv1", 
                     ylab = "VDV1 levels", 
                     color = "group", 
                     facet.by = "dataset",
                     palette = "npg") + 
  xlab("") +
  stat_compare_means(label = "p.signif", label.x.npc=0.5, , label.y=0.47, size = 6) +
  scale_x_discrete(labels=c("0"="VS", "1"="VR")) +
  theme(legend.position="none")

# combining all plots into one and saving it
png("Project/real_synth_bar_boxplot.png", width = 4000, height = 5000, res = 300)
wrap_plots(sb_afb_barp, vr_bp, nos_bp, abpv_bp, bqcv_bp, dwv_bp, kv_bp, vdv1_bp,
           ncol = 2, nrow = 4) + plot_annotation(tag_levels = "A")
dev.off()
