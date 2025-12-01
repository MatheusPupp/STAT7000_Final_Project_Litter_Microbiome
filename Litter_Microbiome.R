# LITTER MICROBIOME ANALYSIS SCRIPT - WAQTER TREATMENTS DAY20 ONLY #

setwd("C:\\Users\\mathe\\Documents") # set your working directory

# Load necessary packages
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(vegan)
library(pairwiseAdonis)
library(picante)
library(ALDEx2)
library(zCompositions)
library(MASS)
library(dendextend)
library(rms)
library(breakaway)
library(dplyr)
library(microViz)
library(ggplot2)
library(ANCOMBC)
library(FSA)
library(DT)

# Creating the phyloseq object from QIIME2 outputs

LitterMicrobiome_Analysis<-qza_to_phyloseq(
  features="table-240-210.qza",
  tree="rooted-tree-240-210.qza",
  taxonomy="taxonomy-240-210.qza",
  metadata = "Day20_metadata.tsv" # metadata has to be .tsv
)

################################################################################
# FILTERING AND SORTING ########################################################

# Sort the samples by number of reads
sort(phyloseq::sample_sums(LitterMicrobiome_Analysis))
     
# Filter out everything not related to bacteria
psJ <- LitterMicrobiome_Analysis %>%
  subset_taxa(Kingdom == "d__Bacteria" &  Family  != "Mitochondria" &    #filter out mitochondria and chloroplast
                Class   != "Chloroplast")

# Filter out everything with less than 5000 reads and empty otus
psJ <- prune_samples(sample_sums(psJ) > 5000, psJ)

# Get unique sampling days
unique_days <- unique(sample_data(psJ)$Day)
print(unique_days)
# Assuming your sample metadata has a column named "Day"
physeq_day20 <- subset_samples(psJ, Day == "20")

# Get the total counts per taxa
taxa_counts <- taxa_sums(physeq_day20)

# Identify taxa with structural zeros across groups
taxa_with_zeros <- taxa_counts[taxa_counts == 0]
print(taxa_with_zeros)

# arrange the treatments for day 20
physeq_day20@sam_data$Water.Treatment <- factor(physeq_day20@sam_data$Water.Treatment, 
                                    levels=c("Control", "Acetic", "Citric"))

################################################################################
# Relative abundance ###########################################################

# Relative abundance by samples
plot_top_taxa <- function(physeq_obj, taxrank = "Genus", topN = 20, plot_title = "") {
  # Check inputs
  if (!inherits(physeq_obj, "phyloseq")) stop("physeq_obj must be a phyloseq object")
  if (!is.character(taxrank) || length(taxrank) != 1) stop("taxrank must be a single string, e.g. 'Genus' or 'Family'")
  
  # Agglomerate to the taxonomic level (phyloseq)
  phy_t <- tax_glom(physeq_obj, taxrank)
  
  # Transform to relative abundance (compositional)
  phy_rel <- microbiome::transform(phy_t, "compositional")
  
  # Melt to long format
  df_rel <- psmelt(phy_rel)
  
  # Ensure taxrank column exists
  if (! taxrank %in% colnames(df_rel)) {
    stop(paste0("Taxonomic column '", taxrank, "' not found in psmelt output. Available columns: ",
                paste(colnames(df_rel)[grepl("Rank|Family|Genus|Phylum|Class|Order", colnames(df_rel))], collapse = ", ")))
  }
  
  # Replace NA or empty names with "Unclassified"
  df_rel[[taxrank]] <- as.character(df_rel[[taxrank]])
  df_rel[[taxrank]][is.na(df_rel[[taxrank]]) | df_rel[[taxrank]] == ""] <- "Unclassified"
  
  # Compute total abundance per taxon across all samples (base R aggregate)
  total_by_taxon <- aggregate(Abundance ~ ., data = df_rel[, c(taxrank, "Abundance")], FUN = sum)
  # aggregate produced two columns: taxrank name and Abundance. Rename
  colnames(total_by_taxon) <- c("TaxonName", "TotalAbundance")
  
  # Select top N taxa by total abundance
  total_by_taxon <- total_by_taxon[order(total_by_taxon$TotalAbundance, decreasing = TRUE), ]
  top_taxa <- head(total_by_taxon$TaxonName, n = topN)
  
  # Create a new column collapsing non-top taxa into "Other"
  df_rel$Taxon2 <- ifelse(df_rel[[taxrank]] %in% top_taxa, df_rel[[taxrank]], "Other")
  
  # Aggregate by Sample x Taxon2 to get relative abundance per sample
  agg_sample_taxon <- aggregate(Abundance ~ Sample + Taxon2, data = df_rel, FUN = sum)
  
  # Turn Taxon2 into a factor with top taxa first, then "Other"
  taxa_levels <- c(top_taxa, setdiff(unique(agg_sample_taxon$Taxon2), top_taxa))
  agg_sample_taxon$Taxon2 <- factor(agg_sample_taxon$Taxon2, levels = taxa_levels)
  
  # OPTIONAL: order samples by treatment (if Water.Treatment exists in sample data)
  sample_order <- rownames(phyloseq::sample_data(physeq_obj))
  if ("Water.Treatment" %in% colnames(phyloseq::sample_data(physeq_obj))) {
    meta_df <- data.frame(sample_data(physeq_obj))
    meta_df$.SampleID <- rownames(meta_df)
    # create order: group samples by treatment then sample name
    sample_order <- meta_df %>%
      arrange(Water.Treatment, rownames(meta_df)) %>%
      pull(.SampleID)
  }
  agg_sample_taxon$Sample <- factor(agg_sample_taxon$Sample, levels = sample_order)
  
  # Plot with ggplot2 (stacked bar)
  p <- ggplot(agg_sample_taxon, aes(x = Sample, y = Abundance, fill = Taxon2)) +
    geom_bar(stat = "identity", width = 0.9) +
    labs(x = "Sample", y = "Relative abundance", fill = taxrank, title = plot_title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7)
    )
  
  return(p)
}

# Generate plots for Genus and Family
# Top 20 Genera - Day 20
rel_abund_genus_day20 <- plot_top_taxa(
  physeq_day20,
  taxrank = "Genus",
  topN = 20,
  plot_title = "Relative Abundance per sample: Top 20 Genera"
)

# Top 20 Families - Day 20
rel_abund_family_day20 <- plot_top_taxa(
  physeq_day20,
  taxrank = "Family",
  topN = 20,
  plot_title = "Relative Abundance per sample: Top 20 Families"
)

# Display plots
rel_abund_genus_day20
rel_abund_family_day20

# Relative abundance by Treatments
plot_top_taxa_grouped <- function(physeq, taxrank, group_var, topN = 20,
                                  plot_title = NULL) {
  
  # Convert phyloseq objects to data frames
  otu_df <- as.data.frame(otu_table(physeq))
  otu_df$OTU <- rownames(otu_df)
  
  tax_df <- as.data.frame(as.matrix(tax_table(physeq)))
  tax_df$OTU <- rownames(tax_df)
  
  meta_df <- as.data.frame(sample_data(physeq))
  meta_df$Sample <- rownames(meta_df)
  
  # Long format
  otu_long <- otu_df %>%
    tidyr::pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance")
  
  # Merge OTU → taxonomy → metadata
  df_rel <- otu_long %>%
    left_join(tax_df, by = "OTU") %>%
    left_join(meta_df, by = "Sample")
  
  # Compute relative abundance per sample
  df_rel <- df_rel %>%
    group_split(Sample) %>%
    purrr::map_df(~ {
      x <- .x
      total <- sum(x$Abundance)
      x$Abundance <- x$Abundance / total
      x
    })
  
  # Total abundance for each taxon (all samples combined)
  tax_totals <- df_rel %>%
    group_split(.data[[taxrank]]) %>%
    purrr::map_df(~ {
      x <- .x
      tibble(
        Taxon = x[[taxrank]][1],
        Total = sum(x$Abundance)
      )
    })
  
  # Identify top N
  top_taxa <- tax_totals %>%
    arrange(desc(Total)) %>%
    slice_head(n = topN) %>%
    pull(Taxon)
  
  # Add "Other" category
  df_rel <- df_rel %>%
    mutate(Taxon2 = ifelse(.data[[taxrank]] %in% top_taxa,
                           .data[[taxrank]],
                           "Other"))
  
  # Summaries per group
  df_group <- df_rel %>%
    group_split(.data[[group_var]], Taxon2) %>%
    purrr::map_df(~ {
      x <- .x
      tibble(
        Group = x[[group_var]][1],
        Taxon = x$Taxon2[1],
        Total = sum(x$Abundance)
      )
    })
  
  # Normalize to 1 per group
  df_group <- df_group %>%
    group_split(Group) %>%
    purrr::map_df(~ {
      x <- .x
      total_group <- sum(x$Total)
      x$RelAbundance <- x$Total / total_group
      x
    })

  palette_colors <- RColorBrewer::brewer.pal(12, "Set3")  
  # Repeat colors if topN > 12
  palette_full <- rep(palette_colors, length.out = topN)
  names(palette_full) <- top_taxa
  
  # Add gray for Other (always last)
  fill_colors <- c(
    palette_full,
    Other = "#BFBFBF"
  )
  
  # Order taxa so Other is last (on top of stack)
  df_group$Taxon <- factor(df_group$Taxon,
                           levels = c(top_taxa, "Other"))
  
  # Plot
  ggplot(df_group,
         aes(x = Group, y = RelAbundance, fill = Taxon)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = group_var,
      y = "Relative Abundance",
      fill = taxrank,
      title = plot_title
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_genus_groups <- plot_top_taxa_grouped(
  physeq_day20,
  taxrank = "Genus",
  group_var = "Water.Treatment",
  topN = 20,
  plot_title = "Water treatments Relative abundance: Top 20 Genera"
)

p_family_groups <- plot_top_taxa_grouped(
  physeq_day20,
  taxrank = "Family",
  group_var = "Water.Treatment",
  topN = 20,
  plot_title = "Water treatments Relative abundance: Top 20 Families"
)

p_genus_groups
p_family_groups

## Pie charts
plot_top_taxa_pie <- function(physeq_obj, taxrank = "Genus", topN = 20,
                              plot_title = "") {
  
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # --- Agglomerate to the taxonomic level ---
  phy_t <- tax_glom(physeq_obj, taxrank)
  
  # --- Relative abundance ---
  phy_rel <- microbiome::transform(phy_t, "compositional")
  
  # --- Long format ---
  df_rel <- psmelt(phy_rel)
  
  # Replace NA/empty names
  df_rel[[taxrank]] <- as.character(df_rel[[taxrank]])
  df_rel[[taxrank]][is.na(df_rel[[taxrank]]) | df_rel[[taxrank]] == ""] <- "Unclassified"
  
  # --- Total abundance per taxon (all samples) ---
  totals <- df_rel %>%
    group_by(.data[[taxrank]]) %>%
    summarise(Total = sum(Abundance)) %>%
    arrange(desc(Total))
  
  # --- Identify top N ---
  top_taxa <- totals %>% slice_head(n = topN) %>% pull(!!taxrank)
  
  # --- Add "Other" ---
  df_rel <- df_rel %>%
    mutate(Taxon2 = ifelse(.data[[taxrank]] %in% top_taxa,
                           .data[[taxrank]], "Other"))
  
  # --- Collapse to one row per taxon ---
  df_pie <- df_rel %>%
    group_by(Taxon2) %>%
    summarise(Total = sum(Abundance)) %>%
    ungroup()
  
  # --- Percentages ---
  df_pie <- df_pie %>%
    mutate(Percent = Total / sum(Total) * 100,
           Label = paste0(Taxon2, " (", sprintf("%.1f", Percent), "%)"))
  
  # --- Color palette (Set3 for taxa, gray for Other) ---
  palette_colors <- RColorBrewer::brewer.pal(12, "Set3")
  palette_full <- rep(palette_colors, length.out = topN)
  names(palette_full) <- top_taxa
  
  fill_colors <- c(
    palette_full,
    Other = "#BFBFBF"
  )
  
  # Factor order: top taxa first, Other last
  df_pie$Taxon2 <- factor(df_pie$Taxon2,
                          levels = c(top_taxa, "Other"))
  
  # --- Plot ---
  ggplot(df_pie, aes(x = "", y = Total, fill = Taxon2)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = grDevices::rainbow(length(unique(df_pie$Taxon2)))) +
    geom_text(
      aes(
        x = 1.6,   
        label = sprintf("%.1f%%", Percent)),
      position = position_stack(vjust = 0.5),
      size = 3) +
    labs(title = plot_title, fill = taxrank) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
}

# Generate pie chart for Genus level
pie_genus_day20 <- plot_top_taxa_pie(
  physeq_day20,
  taxrank = "Genus",
  topN = 20,
  plot_title = "Top 20 Genera relative abundance"
)

pie_genus_day20

# Generate pie chart for Genus level
pie_family_day20 <- plot_top_taxa_pie(
  physeq_day20,
  taxrank = "Family",
  topN = 20,
  plot_title = "Top 20 Families relative abundance"
)

pie_family_day20

################################################################################
# ANCOMBC2 - Differential Abundance ############################################
# Day 20 analysis (only water treatments)
da_result_day20 <- ancombc(
  data = physeq_day20,
  assay.type = NULL,
  assay_name = "counts",
  rank = NULL,
  tax_level = "Genus",
  formula = "Water.Treatment",
  p_adj_method = "BH",             # Benjamini–Hochberg adjust method - to account for sparse features, more useful to use a False Discovery Rate approach ## it takes all the alfas and multiply them by the 
  prv_cut = 0.1,
  lib_cut = 0,
  group = "Water.Treatment",
  struc_zero = TRUE,               # To detect possible structural-zeros and avoid false positive detection of taxa
  neg_lb = FALSE,
  tol = 1e-05,
  max_iter = 200,                  # increase iterations if convergence warnings
  conserve = FALSE,
  alpha = 0.05,
  global = TRUE,                   # consider overall test per taxon
  n_cl = 1,
  verbose = TRUE
)

# Name the columns for the p-values tables
col_name20 = c("taxon", "(Intercept)", "Water.TreatmentAcetic", "Water.TreatmentCitric")

# ANCOMBC p-values
#Day 20
tab_p20 = da_result_day20$res$p_val
colnames(tab_p20) = col_name20
tab_p20 %>% 
  datatable(caption = "P-values from the Primary Result Day 20") %>%
  formatRound(col_name20[-1], digits = 3)

# Filtering only lines with p-values < 0.05 for da_result_day20
tab_p_day20_filtered <- tab_p20 %>%
  filter(`Water.TreatmentAcetic` < 0.05 | `Water.TreatmentCitric` < 0.05)

# ANCOMBC Adjusted p-values (q_val = adjusted p-values)
# Day 20
tab_q20 = da_result_day20$res$q_val
colnames(tab_q20) = col_name20
tab_q20 %>% 
  datatable(caption = "Adjusted p-values from the Primary Result Day 20") %>%
  formatRound(col_name20[-1], digits = 3)
tab_q_day20_filtered <- tab_q20 %>%
  filter(`Water.TreatmentAcetic` < 0.05 | `Water.TreatmentCitric` < 0.05)

# Save the adjusted p-values table of day 20 in TSV format
write.table(tab_q_day20_filtered, file = "q_values_day20_Phylum_BH.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


## Looking for Log Fold Changes (LFC)

# Day 20
# Extracting LFC and preparing table structure
df_lfc = data.frame(da_result_day20$res$lfc[, -1] * da_result_day20$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = da_result_day20$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Renaming columns with treatments
colnames(df_lfc) <- c("taxon_id", "Nothing", "Acetic", "Citric")

# Filtering to keep only Acetic and Citric treatments (if significant and present)
df_fig_water <- df_lfc %>%
  dplyr::select(taxon_id, Acetic, Citric) %>%
  tidyr::pivot_longer(cols = c("Acetic", "Citric"), names_to = "treatment", values_to = "LFC") %>%
  dplyr::filter(!is.na(LFC) & LFC != 0) %>%
  dplyr::mutate(direct = ifelse(LFC > 0, "Positive LFC", "Negative LFC"))

# Check which treatments actually have significant taxa
present_treatments <- unique(df_fig_water$treatment)

# Create a dynamic title that reflects the available treatments
treatments_label <- paste(present_treatments, collapse = " & ")
plot_title <- paste("Log Fold Changes on Genus for", treatments_label, "Treatment on Day 20")

# Adjust factor levels to only include treatments that exist
df_fig_water$treatment <- factor(df_fig_water$treatment, levels = present_treatments)

# Create the plot dynamically
p_water_d20 <- ggplot(data = df_fig_water, aes(y = taxon_id, x = LFC, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  facet_wrap(~treatment, nrow = 1, scales = "free_y") +  # only plots treatments with data
  labs(y = "Taxon", x = "Log Fold Change (LFC)", title = plot_title) + 
  scale_fill_manual(
    name = "Direction",
    values = c("Positive LFC" = "#D9534F", "Negative LFC" = "steelblue")  # customize or invert
  ) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 14),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(hjust = 1, size = 8, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    plot.margin = unit(c(1,1,1,1), "cm"), 
    aspect.ratio = 1.5,
    panel.spacing = unit(1, "lines")
  )

# Displaying the graph
p_water_d20

## Computing impact of treatments from a differential abundance perspective
# Compute impact metrics: Metric 1 - Number of significant taxa (richness of effect)
#                         Metric 2 - Sum of |LFC| across taxa (magnitude of effect)
#                         Metric 3 - Mean |LFC| (effect per taxon)
#                         Metric 4 - Largest single LFC (extreme effect)
impact_metrics <- df_fig_water %>%
  group_by(treatment) %>%
  summarise(
    n_significant_taxa = n(),                           # Metric 1
    total_abs_LFC = sum(abs(LFC)),                     # Metric 2
    mean_abs_LFC = mean(abs(LFC)),                     # Metric 3
    max_abs_LFC = max(abs(LFC))                        # Metric 4
  ) %>%
  arrange(desc(total_abs_LFC))

impact_metrics

# Two proportion Z-test - Does Treatment A produce significantly more differential taxa than Treatment B?
n_total <- nrow(tab_q20)

n_acetic <- sum(tab_q20$Water.TreatmentAcetic < 0.05)
n_citric <- sum(tab_q20$Water.TreatmentCitric < 0.05)

prop.test(
  x = c(n_acetic, n_citric),
  n = c(n_total, n_total),
  alternative = "two.sided",
  correct = FALSE
) # p-value = 1.336e-11 -> one of the treatments significantly have more differential taxa than the other

# Wilcoxon rank-sum test - Does one treatment cause stronger shifts (LFC magnitude) than another?
df_significant <- df_fig_water %>% 
  filter(!is.na(LFC) & LFC != 0)
wilcox.test(
  abs(LFC) ~ treatment,
  data = df_significant,
  alternative = "two.sided"
) # W = 507, p-value = 0.0003119 - LFC magnitude significantly differ between treatments

# global test from ANCOMBC2 which tests overall differential abundance for each taxon across all groups rather than comparing specific pairs of groups.
global_day20 <- da_result_day20$res_global # extract global results
# In the global result we observe that the water treatments are having significant effect in changing the microbiome of the litter

################################################################################
# ALPHA DIVERSITY ##############################################################

# DAY 20
#Generate a data.frame with adiv measures
adiv20 <- data.frame(
  "Observed" = phyloseq::estimate_richness(physeq_day20, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(physeq_day20, measures = "Shannon"),
  "Group" = phyloseq::sample_data(physeq_day20)$Water.Treatment)
#Generate Pielou
"Pielou" = adiv20$Pielou <- adiv20$Shannon / log(adiv20$Observed)
#Plot adiv measures
alpha_plot20 <- adiv20 %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Pielou")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Pielou"))) %>%
  ggplot(aes(group=Group, x = factor(Group, level=c("Nothing", "Acetic", "Citric")), y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Group), height = 0, width = .2, cex = 2.5) +
  scale_color_manual(values=c( "#53B74C", "#B3A033", "#E77D72",
                               "#55BCC2", "#6F9BF8",  "#E46EDD" )) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  ggtitle("Day 20") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(color = "black", size = 20, hjust=1, angle=45),
        axis.text.y = element_text(color = "black", size = 14),
        strip.text = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, color = "black"),
        panel.spacing.x = unit(2, "lines"))
ggsave("alpha diversity OSS 20.png", alpha_plot20, height = 6, width = 13, units = "in", dpi = 300 )

###############################
# NORMALITY AND VARIANCE TESTS
###############################

library(car)  # for Levene test

alpha_vars <- c("Observed", "Shannon", "Pielou")

shapiro_results <- list()

for (metric in alpha_vars) {
  cat("\n-------------------------------\n")
  cat("Shapiro-Wilk Test for:", metric, "\n")
  cat("-------------------------------\n")
  
  # Split dataset by Group
  res_list <- lapply(split(adiv20[[metric]], adiv20$Group), shapiro.test)
  
  # Convert results into a tibble
  shapiro_results[[metric]] <- tibble(
    Group = names(res_list),
    W = sapply(res_list, function(x) x$statistic),
    p.value = sapply(res_list, function(x) x$p.value)
  )
  
  print(shapiro_results[[metric]])
}


# Q-Q PLOTS FOR EACH GROUP (faceted)
for (metric in alpha_vars) {
  p <- ggplot(adiv20, aes(sample = !!sym(metric))) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ Group) +
    ggtitle(paste("Q-Q Plot for", metric, "(Day 20)")) +
    theme_bw()
  
  print(p)
  
  ggsave(
    filename = paste0("QQplot_Day20_", metric, ".png"),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
}

# LEVENE'S TEST FOR EQUALITY OF VARIANCES
levene_results <- list()

for (metric in alpha_vars) {
  cat("\n-------------------------------\n")
  cat("Levene Test for:", metric, "\n")
  cat("-------------------------------\n")
  
  levene_results[[metric]] <- leveneTest(as.formula(paste(metric, "~ Group")), data = adiv20)
  print(levene_results[[metric]])
}

#Kruskal test of location
kruskal.test(Observed ~ Group, data=adiv20)
kruskal.test(Shannon ~ Group, data = adiv20)            
kruskal.test(Pielou ~ Group, data = adiv20)

# Pairwise comparisons (with Bonferroni correction)
dunnTest(Shannon ~ Group,
         data = adiv20,
         method= "bonferroni")

# Pairwise comparisons (with Bonferroni correction)
dunnTest(Observed ~ Group,
         data = adiv20,
         method= "bonferroni")

# Pairwise comparisons (with Bonferroni correction)
dunnTest(Pielou ~ Group,
         data = adiv20,
         method= "bonferroni")

################################################################################
# BETA DIVERSITY ###############################################################

# DAY 20
ps20_clr <- microbiome::transform(physeq_day20, "clr")

#PCA via phyloseq
ord_clr20 <- phyloseq::ordinate(ps20_clr, "RDA")

#Plot scree plot
pcoa20 <-phyloseq::plot_scree(ord_clr20) + 
  geom_bar(stat="identity", fill = "royalblue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle("Day 20") +
  theme_bw() 
theme(axis.text.x = element_text(size = 12, angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size = 12))
ggsave("pcoa_plot_20.png", pcoa20, height = 8, width = 13, units = "in", dpi = 300 )

#Scale axes and plot ordination with first two PCs
clr1 <- ord_clr20$CA$eig[1] / sum(ord_clr20$CA$eig)
clr2 <- ord_clr20$CA$eig[2] / sum(ord_clr20$CA$eig)


ord_plot20 <- phyloseq::plot_ordination(ps20_clr, ord_clr20, type="samples", color="Water.Treatment") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Water.Treatment), linetype = 2) +
  ggtitle("Day 20") +
  theme_bw()
ggsave("ordination plot 20.png", ord_plot20, height = 8, width = 10, units = "in", dpi = 300 )

#Extract the OTU table:
otu_table20 <- otu_table(physeq_day20)

#Convert the OTU table to a matrix:
otu_matrix20 <- as.matrix(otu_table20)

# Plot the ordination
plot_ordination(physeq_day20, ord_clr20, type = "samples", color = "Water.Treatment") +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Water.Treatment), linetype = 2) +
  ggtitle("Day 20") +
  theme_bw()

#Perform PERMANOVA with Bray-Curtis dissimilarity
bray_curtis_distance20 <- phyloseq::distance(physeq_day20, method ="bray") 
adonis2_resultbraycurtis20 <- adonis2(bray_curtis_distance20 ~ phyloseq::sample_data(ps20_clr)$Water.Treatment)
pairwise.adonis(bray_curtis_distance20, phyloseq::sample_data(physeq_day20)$Water.Treatment)


# Perform PERMANOVA with Unweighted UniFrac distance
unweighted_unifrac20 <- UniFrac(physeq=physeq_day20, weighted=FALSE, normalized=T, parallel=T, fast=T)
adonis2_resultunweighted20 <- adonis2(unweighted_unifrac20 ~ phyloseq::sample_data(ps20_clr)$Water.Treatment)
pairwise.adonis(unweighted_unifrac20, phyloseq::sample_data(physeq_day20)$Water.Treatment)


# Perform PERMANOVA with weighted UniFrac distance
weighted_unifrac20 = UniFrac(physeq=physeq_day20, weighted=TRUE, normalized=T, parallel=T, fast=T)
adonis2_resultweighted20 <- adonis2(weighted_unifrac20 ~ phyloseq::sample_data(ps20_clr)$Water.Treatment)
pairwise.adonis(weighted_unifrac20, phyloseq::sample_data(physeq_day20)$Water.Treatment)

#################################################
