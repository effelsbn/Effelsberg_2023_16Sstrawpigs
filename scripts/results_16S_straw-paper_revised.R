here::i_am("scripts/results_16S_straw-paper.R")

library(here)
library(tidyverse)
library(phyloseq)
library(decontam)
library(scales)
library(ggbump)
library(wesanderson)
library(patchwork)
library(vegan)

# Import files from DADA2-pipeline
readtrack <- read.delim(file = here("input", "readtrack.tsv"), row.names = 1)
metadata <- read.delim(file = here("input", "metadata_16S_Strohstall.tsv"), row.names = 1)
load(here("input", "phyloseq_original.RData"))
ps_original <- ps

# Create filtered objects
all_info <- merge(metadata, readtrack, by = 0)
all_info_kept <- all_info %>%
  filter(type != "control" & ID != "2_394" & type != "na")

ps <- subset_samples(ps_original, type != "control" & ID != "2_394") # continue working without controls & empty sample
OTU <- as.data.frame(otu_table(ps))
taxa <- as.data.frame(tax_table(ps))
melt <- psmelt(ps)

# Remove contamination
sample_data(ps_original)$is.neg <- sample_data(ps_original)$source == "NTC"
contamdf.prev <- isContaminant(ps_original, method = "prevalence", neg = "is.neg", threshold = 0.5)

taxa_contam <- as.data.frame(tax_table(ps_original)[which(contamdf.prev$contaminant)])
genus_contam <- taxa_contam$Genus

ASV_contaminants <- paste0("ASV", which(contamdf.prev$contaminant))
n_samples_contam <- melt %>%
  filter(OTU %in% ASV_contaminants & Abundance > 0) %>%
  summarize(ID) %>%
  nrow()

ps.contam <- ps
ps <- prune_taxa(!contamdf.prev$contaminant, ps)

# Subset phyloseq objects
ps_straw <- psadd::subset_samples_no_zero(ps, type == "straw")
ps_pigs <- psadd::subset_samples_no_zero(ps, type == "pig")

# Create krona plots, commented out because only works in Linux environment with kronatools installed
#plot_krona(ps_straw, output = "krona/straw", variable = "week")
#plot_krona(ps_pigs, output = "krona/pigs", variable = "week")

# Get top taxa
n_reads_straw <- sum(sample_sums(ps_straw))
n_reads_pigs <- sum(sample_sums(ps_pigs))

straw_summary <- psmelt(ps_straw) %>%
  select(OTU, ID, Abundance, Phylum, Class, Order, Family, Genus, Species) %>%
  filter(Abundance > 0) %>%
  group_by(OTU) %>%
  mutate(Sum = sum(Abundance), Sample_count = n_distinct(ID)) %>%
  ungroup() %>%
  select(-ID, -Abundance) %>%
  distinct()

pig_summary <- psmelt(ps_pigs) %>%
  select(OTU, ID, Abundance, Phylum, Class, Order, Family, Genus, Species) %>%
  filter(Abundance > 0) %>%
  group_by(OTU) %>%
  mutate(Sum = sum(Abundance), Sample_count = n_distinct(ID)) %>%
  ungroup() %>%
  select(-ID, -Abundance) %>%
  distinct()

top_phyla_straw <- straw_summary %>%
  group_by(Phylum) %>%
  summarise(Rel_Abun = sum(Sum) / n_reads_straw) %>%
  arrange(desc(Rel_Abun)) %>%
  filter(Rel_Abun > 0.05)

top_phyla_pigs <- pig_summary %>%
  group_by(Phylum) %>%
  summarise(Rel_Abun = sum(Sum) / n_reads_pigs) %>%
  arrange(desc(Rel_Abun)) %>%
  filter(Rel_Abun > 0.05)

top_fam_straw <- straw_summary %>% 
  group_by(Family) %>% 
  summarise(Rel_Abun = sum(Sum)/n_reads_straw) %>% 
  arrange(desc(Rel_Abun)) %>% 
  filter(Rel_Abun > 0.05)

top_fam_pigs <- pig_summary %>% 
  group_by(Family) %>% 
  summarise(Rel_Abun = sum(Sum)/n_reads_pigs) %>% 
  arrange(desc(Rel_Abun)) %>% 
  filter(Rel_Abun > 0.05)

# assess core microbiome
ps_pigs_rel <- microbiome::transform(ps_pigs, "compositional")
core_pigs <- microbiome::core_members(ps_pigs_rel, detection = 0.005, prevalence = 0.5) # get taxa at 0.5% abundance in 50% samples

core_pigs_full <- psmelt(ps_pigs) %>%
  filter(Abundance > 0) %>% 
  group_by(OTU) %>% 
  mutate(Sample_count = n_distinct(ID),
         Sum_Abund = sum(Abundance)) %>% 
  ungroup() %>% 
  filter(OTU %in% core_pigs) %>% 
  select(OTU, Phylum, Class, Order, Family, Genus, Species, Sample_count, Sum_Abund) %>% 
  distinct() %>% 
  mutate(Rel_Abund = Sum_Abund/sum(Sum_Abund))

# Create Phyla bar plot

phyla <- melt %>%
  select(type, week, Phylum, Abundance) %>%
  filter(type != "control" & Abundance > 0) %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla_pigs$Phylum, Phylum, "other")) %>%
  group_by(type, week, Phylum) %>%
  summarize(sum_abund = sum(Abundance), .groups = "drop") %>%
  group_by(type, week) %>%
  mutate(week_total = sum(sum_abund)) %>%
  ungroup() %>%
  # mutate(percentage = percent(sum_abund/week_total, accuracy = 0.1)) %>%
  mutate(rel_abund = round((sum_abund / week_total) * 100, digits = 1)) %>%
  mutate(Phylum = factor(Phylum, levels = c("other", "Actinobacteriota", "Bacteroidota", "Proteobacteria", "Firmicutes")))

bars <- ggplot(phyla, aes(x = type, y = sum_abund, fill = Phylum)) +
  geom_col(position = "fill") +
  geom_text(aes(label = rel_abund), position = position_fill(vjust = 0.5), size = 3.5, show.legend = F) +
  facet_wrap(~week, nrow = 1) +
  scale_fill_manual(
    breaks = c("other", "Actinobacteriota", "Bacteroidota", "Proteobacteria", "Firmicutes"),
    values = c("#888888", "#CC6677", "#88CCEE", "#DDCC77", "#44AA99", "#6699CC")
  ) +
  labs(y = "relative abundance", x = NULL) +
  scale_y_continuous(
    labels = percent,
    expand = expansion(add = c(0, 0.03))
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 11, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  )

ggsave("figures/phyla_barplot.pdf", width = 170, height = 120, units = "mm")

# Create bump chart classes
pig_tax_by_week <- psmelt(ps_pigs) %>%
  select(OTU, ID, week, Abundance, Phylum, Class, Order, Family, Genus, Species) %>%
  filter(Abundance > 0) %>%
  group_by(week, OTU) %>%
  mutate(Sum = sum(Abundance)) %>%
  ungroup() %>%
  select(-ID, -Abundance) %>%
  distinct() %>%
  mutate(Rel_Abund = Sum / n_reads_pigs)

pig_class_by_week_rank <- pig_tax_by_week %>%
  mutate(week = as.numeric(week), .keep = "unused") %>%
  group_by(week, Class) %>%
  filter(is.na(Class) == FALSE) %>%
  mutate(Class_Abund = sum(Sum)) %>%
  ungroup() %>%
  select(week, Class, Class_Abund) %>%
  distinct() %>%
  group_by(week) %>%
  mutate(Rank = rank(desc(Class_Abund))) %>%
  filter(Rank <= 10) %>%
  select(-Class_Abund) %>%
  arrange(week, Rank)

bumps <- ggplot(pig_class_by_week_rank, mapping = (aes(x = week, y = Rank, color = Class))) +
  geom_bump(size = 1.5) +
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Darjeeling2", 13, type = "continuous")) +
  scale_x_continuous(
    limits = c(0, 14),
    breaks = c(0, 1, 5, 10, 11, 14)
  ) +
  scale_y_reverse(breaks = 1:10) +
  labs(x = "week", y = "abundance rank") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  )

ggsave("figures/class_bumps.pdf", width = 170, height = 120, units = "mm")

## Staphylococci

staphs <- melt %>%
  mutate(week = as.numeric(week), .keep = "unused") %>%
  mutate(source = case_when(
    (type == "straw" & cleaning == "clean") ~ "St1",
    (type == "straw" & cleaning == "desinfect") ~ "St2",
    TRUE ~ source
  )) %>%
  filter(Abundance > 0 & type != "control") %>%
  group_by(Sample) %>%
  mutate(
    sample_sum = sum(Abundance),
    rel_abund = Abundance / sample_sum
  ) %>%
  ungroup() %>%
  filter(Genus == "Staphylococcus") %>%
  select(type, week, MRSA, source, Sample, Genus, rel_abund) %>%
  group_by(Sample, Genus) %>%
  mutate(sum_rel_abund = sum(rel_abund), .keep = "unused") %>%
  distinct()

ggplot(staphs, aes(x = week, y = sum_rel_abund, group = source, color = type)) +
  geom_point() +
  geom_line() +
  geom_smooth(aes(group = type), se = FALSE, size = 2) +
  scale_color_manual(
    breaks = c("pig", "straw"),
    # values = c("#E58B8B", "#E9BD64")) +
    values = c("#E58B8B", "#EDA931")
  ) +
  scale_x_continuous(
    limits = c(0, 14),
    breaks = c(0, 1, 5, 10, 11, 14),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_y_continuous(
    labels = percent,
    expand = expansion(add = c(0.003, 0.002))
  ) +
  labs(x = "week", y = "relative abundance of staphylococci") +
  theme_classic() +
  theme(
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.position = c(0.5, 0.95),
    legend.title = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  )

ggsave("figures/staphylococci.pdf", width = 170, height = 100, units = "mm")

## Diversity
ps_pigs_rare <- rarefy_even_depth(ps_pigs, sample.size = min(sample_sums(ps_pigs)), rngseed = 17, replace = FALSE)
alpha_pig_rare <- microbiome::alpha(ps_pigs_rare, index = "all")

alpha_complete <- merge(metadata, alpha_pig_rare, by = "row.names", all = FALSE)

shannon <- ggplot(alpha_complete, aes(x = week, y = diversity_shannon)) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "Shannon diversity"
  ) +
  theme_classic()

observed <- ggplot(alpha_complete, aes(x = week, y = observed)) +
  geom_boxplot() +
  labs(
    x = "week",
    y = "Observed species"
  ) +
  theme_classic()

OTU_pigs <- as.data.frame(otu_table(ps_pigs))
curvedata <- vegan::rarecurve(OTU_pigs, step = 100)

curve_df <- map_dfr(curvedata, bind_rows) %>%
  bind_cols(Sample = rownames(OTU_pigs), .) %>%
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name)

rcurve <- ggplot(curve_df, aes(x = n_seqs, y = value, group = Sample)) +
  geom_line(color = "gray30", size = 0.65) +
  geom_vline(xintercept = min(sample_sums(ps_pigs)), color = "sienna1", size = 0.8) +
  annotate("text", x = 16000, y = 1200, label = min(sample_sums(ps_pigs)), size = 3, color = "sienna1", fontface = 2) +
  labs(x = "Number of sequences", y = "Number of ASVs") +
  scale_x_continuous(expand = c(0.005, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

rcurve + shannon / observed +
  plot_annotation(tag_levels = "A")

ggsave("figures/diversity_combined.pdf", width = 170, units = "mm")




## Bray-Curtis
pig_dist_bray <- vegan::avgdist(OTU_pigs, sample = min(sample_sums(ps_pigs))) # 100 random vegdist, sample = 10358

pig_meta <- all_info %>%
  rename(Sample = Row.names) %>%
  filter(type == "pig" & ID != "2_394") %>%
  select(Sample, source, week, MRSA) %>%
  mutate(week = as.numeric(week), .keep = "unused")

pig_meta_B <- pig_meta %>%
  mutate(
    name = Sample, sourceB = source, weekB = week, MRSAB = MRSA,
    .keep = "unused"
  )

pig_dist_tbl <- pig_dist_bray %>%
  as.matrix() %>%
  as_tibble(rownames = "Sample") %>%
  pivot_longer(-Sample) %>%
  inner_join(pig_meta, by = "Sample") %>%
  distinct() %>%
  inner_join(pig_meta_B, by = "name") %>%
  distinct() %>%
  mutate(sampleA = Sample, sampleB = name, dist = value, .keep = "unused") %>%
  select(dist, sampleA, sampleB, source, sourceB, week, weekB, MRSA, MRSAB)

pig_dist_small <- pig_dist_tbl %>%
  filter(source == sourceB & week == 0 & weekB != 0)

dayzero <- ggplot(pig_dist_small, aes(x = weekB, y = dist, group = source, color = source)) +
  geom_point() +
  geom_line() +
  geom_smooth(aes(group = 1), se = FALSE, color = "sienna1", size = 2) +
  scale_x_continuous(breaks = c(0, 1, 5, 10, 11, 14)) +
  scale_color_grey() +
  labs(
    x = "Weeks after housing in",
    y = "Bray-Curtis distance to day 0"
  ) +
  theme_classic() +
  theme(legend.position = "none")

## NMDS
ord_all <- ordinate(ps, "NMDS", "bray")

NMDS_cols <- feathers::get_pal("cassowary")[1:6]

NMDS_labels <- data.frame(
  x = c(-1.1, -0.6, 0.2, 0.7, 0.1, 1.25),
  y = c(-0.5, 0.4, 0.9, -0.2, -0.7, 0.6),
  label = c("00", "01", "05", "10", "11", "14"),
  week = c("00", "01", "05", "10", "11", "14")
)

NMDS_all <- plot_ordination(ps, ord_all, color = "week", shape = "type") +
  geom_point(size = 4) +
  stat_ellipse() +
  scale_color_manual(values = feathers::get_pal("cassowary")) +
  scale_shape_manual(values = c(16, 18, 17)) +
  geom_label(data = NMDS_labels, aes(x = x, y = y, label = label, color = week), fontface = 2) +
  theme_classic() +
  theme(legend.position = "none")

NMDS_all / dayzero  +
  plot_annotation(tag_levels = "A")

ggsave("figures/beta_combined_revised.pdf", width = 170, height = 200, units = "mm")

## Statistics
dist_bray <- vegan::avgdist(OTU, sample = min(sample_sums(ps))) # 100 random vegdist, sample = 6245
bray_meta <- all_info_kept %>%
  rename(Sample = Row.names) %>%
  select(Sample, source, week, MRSA) %>%
  mutate(week = as.numeric(week), .keep = "unused")

adonis_inkSt <- adonis2(dist_bray~bray_meta$week)

bray_test_inkSt <- pairwiseAdonis::pairwise.adonis2(dist_bray~week, data = bray_meta)

bray_df_inkSt <- data.frame(do.call(rbind.data.frame, bray_test_inkSt)) 

bray_sigs_inkSt <- bray_df_inkSt %>% 
  filter(rownames(bray_df_inkSt) != "parent_call") %>% 
  drop_na() %>% 
  mutate(p = as.numeric(Pr..F.))

## Supplemental figure mock
melt_all <- psmelt(ps_original)
mocks <- melt_all %>% 
  filter(source == "mock") %>%
  filter(Abundance > 0) %>%
  group_by(Genus, ID) %>% 
  mutate(abs_abund = sum(Abundance)) %>%
  ungroup() %>% 
  select(ID, Genus, abs_abund) %>%
  distinct() %>% 
  group_by(ID) %>% 
  mutate(read_total = sum(abs_abund)) %>% 
  ungroup() %>%
  mutate(rel_abund = (abs_abund/read_total)*100, .keep = "unused") %>% 
  pivot_wider(names_from = ID, values_from = rel_abund) %>% 
  arrange(Genus) %>% 
  filter(Mock1 > 0.05 & Mock2 > 0.05)

zymo <- read.delim(file = "input/zymo_mock_genus.tsv")

comptab <- merge(zymo, mocks, by = "Genus", all = TRUE)

cols26 <- read_lines("input/26colors.txt")

comptab_stack <- arrange(comptab, desc(Zymo)) %>%
  pivot_longer(cols = !Genus, names_to = "sample")
comptab_stack$Genus <- factor(comptab_stack$Genus, levels = unique(comptab$Genus))

ggplot(comptab_stack, aes(fill=Genus, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity", color = "black") +
  scale_fill_manual(values = cols26) +
  ylab("relative abundance") +
  xlab(NULL) +
  scale_x_discrete(labels = c("Zymo", "Run1", "Run2")) +
  theme_classic()

ggsave("figures/mock.pdf", width = 170, units = "mm")



