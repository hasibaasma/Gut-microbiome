## ---- setup ----
knitr::opts_chunk$set(echo = TRUE)

## ---- Load packages 1 ----
library(dplyr)
library(magrittr)
library(knitr)
library(tidyverse)
library(lattice)
library(latticeExtra)

## ---- Load packages 2 ----
library(phyloseq)
library(microbiome)
library(microViz)
library(ggplot2)
library(microbiomeutilities)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(grid)
library(DirichletMultinomial)
library(ggalluvial)
library(broom)
library(broom.mixed)
library(dotwhisker)
library(tibble)
library(nlme)
library(glmmTMB)
library(splines)
library(lme4)
library(MuMIn)
library(aod)
library(DHARMa)
library(jtools)

## ---- Load Data ----
physeq_mOTU <- readRDS("~/Downloads/mOTUs_phyloseq.rds")

## ---- Preprocessing of the data ----
# Make compositional
physeq_mOTU.rel <- microbiome::transform(physeq_mOTU, "compositional")

# Aggregate to family level
physeq_mOTU.rel <- aggregate_taxa(physeq_mOTU.rel, "family")
physeq_mOTU.rel
get_taxa_unique(physeq_mOTU.rel, "family")

# there was a change in family name
family_to_change <- "Firmicutesfam.incertaesedis"
new_family_name <- "Bacillotafam.incertaededis"

tax_table_physeq <- tax_table(physeq_mOTU.rel)
rows_to_change <- tax_table_physeq[, "family"] == family_to_change
tax_table_physeq[rows_to_change, "family"] <- new_family_name
tax_table(physeq_mOTU.rel) <- tax_table_physeq
get_taxa_unique(physeq_mOTU.rel, "family")

## ---- Subset data and select top 15 patient microbiota ----
temp <- aggregate_top_taxa2(physeq_mOTU.rel, 15, "family")
temp@otu_table %>% rownames()
temp@otu_table <- temp@otu_table[c(11, 1:10, 12:16), ]
temp@tax_table %>% rownames()
temp@tax_table <- temp@tax_table[c(11, 1:10, 12:16), ]
temp

temp.donorA <- temp %>% ps_filter(subject_id == "Donor A", .keep_all_taxa = TRUE)
temp.donorA@sam_data[["timepoint.new"]] <- factor(temp.donorA@sam_data[["timepoint.new"]],
                                                  levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"))

temp.donorB <- temp %>% ps_filter(subject_id == "Donor B", .keep_all_taxa = TRUE)
temp.donorB@sam_data[["timepoint.new"]] <- factor(temp.donorB@sam_data[["timepoint.new"]],
                                                  levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14"))

temp.base <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Baseline", .keep_all_taxa = TRUE)

temp.fmt1 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Pre-FMT", .keep_all_taxa = TRUE)

temp.fmt2 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Post-1", .keep_all_taxa = TRUE)

temp.fmt3 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Post-2", .keep_all_taxa = TRUE)

temp.fmt4 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Post-3", .keep_all_taxa = TRUE)

temp.wk7 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Post-4", .keep_all_taxa = TRUE)

temp.wk8 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Week8", .keep_all_taxa = TRUE)

temp.wk10 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Week10", .keep_all_taxa = TRUE)

temp.wk14 <- temp %>% ps_select(clinical_outcome_wk14, timepoint.new) %>%
  ps_filter(timepoint.new == "Week14", .keep_all_taxa = TRUE)

## ---- Functions for plots ----
make_barplot1_manual <- function (dfm, group_by) {
  dfm <- dfm %>% arrange(Tax)
  dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
  p <- ggplot(dfm, aes(x = Sample, y = Abundance, fill = Tax, colour = "Black")) +
    geom_bar(position = "stack", stat = "identity", colour = "Black") +
    scale_x_discrete(labels = dfm$xlabel, breaks = dfm$Sample)
  p <- p + labs(y = "Abundance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) +
    guides(fill = guide_legend(reverse = FALSE))
  if (!is.null(group_by)) {
    p <- p + facet_grid(. ~ Group, drop = TRUE, space = "free", scales = "free")
  }
  p
}

plot_composition_manual <- function (x, sample.sort = NULL, otu.sort = NULL, x.label = "sample",
                                     plot.type = "barplot", verbose = FALSE, average_by = NULL,
                                     group_by = NULL, ...) {
  if (!is.null(x@phy_tree)) x@phy_tree <- NULL
  xorig <- x
  abu <- abundances(x)
  if (!is.null(average_by)) {
    dff <- as.data.frame(t(abu))
    dff$group <- sample_data(x)[[average_by]]
    if (is.numeric(dff$group) || is.character(dff$group)) {
      dff$group <- factor(dff$group, levels = sort(unique(dff$group)))
    }
    dff <- dff %>% filter(!is.na(group))
    dff$group <- droplevels(dff$group)
    av <- aggregate(. ~ group, data = dff, mean)
    rownames(av) <- as.character(av$group)
    av$group <- NULL
    abu <- t(av)
  }
  if (is.null(sample.sort) || sample.sort == "none" || !is.null(average_by)) {
    sample.sort <- colnames(abu)
  }
  else if (length(sample.sort) == 1 && sample.sort %in% taxa(xorig)) {
    tax <- sample.sort
    sample.sort <- rev(sample_names(x)[order(abundances(x)[tax, ])])
  }
  else if (length(sample.sort) == 1 && sample.sort %in% names(sample_data(x)) &&
           is.null(average_by)) {
    sample.sort <- rownames(sample_data(x))[order(sample_data(x)[[sample.sort]])]
  }
  else if (all(sample.sort %in% sample_names(x)) & is.null(average_by)) {
    sample.sort <- sample.sort
  }
  else if (length(sample.sort) == 1 && sample.sort == "neatmap") {
    sample.sort <- neatsort(x, method = "NMDS", distance = "bray",
                            target = "sites", first = NULL)
  }
  else if (is.vector(sample.sort) && length(sample.sort) > 1) {
    sample.sort <- sample_names(x)[sample.sort]
  }
  else if (!sample.sort %in% names(sample_data(x))) {
    warning(paste("sample.sort argument", sample.sort,
                  "not in sample_data(x). Using original order."))
    sample.sort <- sample_names(x)
  }
  
  if (is.null(otu.sort) || otu.sort == "none") {
    otu.sort <- taxa(x)
  }
  else if (length(otu.sort) == 1 && otu.sort == "abundance2") {
    otu.sort <- rev(c(rev(names(sort(rowSums(abu)))[seq(1, nrow(abu), 2)]),
                      names(sort(rowSums(abu)))[seq(2, nrow(abu), 2)]))
  }
  else if (length(otu.sort) == 1 && otu.sort == "abundance") {
    otu.sort <- rev(names(sort(rowSums(abu))))
  }
  else if (length(otu.sort) == 1 && otu.sort %in% colnames(tax_table(x))) {
    otu.sort <- rownames(sample_data(x))[order(tax_table(x)[[otu.sort]])]
  }
  else if (all(otu.sort %in% taxa(x))) {
    otu.sort <- otu.sort
  }
  else if (length(otu.sort) == 1 && otu.sort == "neatmap") {
    otu.sort <- neatsort(x, method = "NMDS", distance = "bray",
                         target = "species", first = NULL)
  }
  
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("Tax", "Sample", "Abundance")
  dfm$Sample <- factor(dfm$Sample, levels = sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels = otu.sort)
  if (!is.null(group_by)) {
    if (!is.null(average_by)) {
      dfm$Group <- meta(x)[[group_by]][match(as.character(dfm$Sample),
                                             meta(x)[[average_by]])]
    } else {
      dfm$Group <- meta(x)[[group_by]][match(as.character(dfm$Sample),
                                             sample_names(x))]
    }
  }
  if (x.label %in% colnames(sample_data(x)) & is.null(average_by)) {
    meta <- sample_data(x)
    dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample), x.label]))
    if (is.factor(meta[, x.label])) {
      lev <- levels(meta[, x.label])
    } else {
      lev <- unique(as.character(unname(unlist(meta[, x.label]))))
    }
    dfm$xlabel <- factor(dfm$xlabel, levels = lev)
  } else {
    dfm$xlabel <- dfm$Sample
  }
  if (plot.type == "barplot") {
    p <- make_barplot1_manual(dfm, group_by)
  } else if (plot.type == "heatmap") {
    p <- make_heatmap1(x, otu.sort, sample.sort, verbose)
  } else if (plot.type == "lineplot") {
    p <- make_lineplot1(dfm)
  } else stop("plot.type argument not recognized")
  p
}

## ---- Colors ----
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(16)

## ---- Average Plots per timepoint ----
pdonorA <- plot_composition_manual(temp.donorA, average_by = "treated_with_donor") +
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Average", y = "Relative Abundance") +
  scale_fill_manual(values = mycolors) +
  scale_x_discrete(labels = c("Donor A"))

pdonorB <- plot_composition_manual(temp.donorB, average_by = "treated_with_donor") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.6),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line("white"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Average", y = "Relative Abundance") +
  scale_fill_manual(values = mycolors) +
  scale_x_discrete(labels = c("Donor B"))

# (Repeat p1 … p9 code blocks exactly as in Rmd)
# ...
plot.compositions <- (pdonorA + pdonorB + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) +
  plot_layout(nrow = 1)
plot.compositions

## ---- Phyloseq2 object ----
physeq_mOTU <- aggregate_taxa(physeq_mOTU, "family")
physeq_mOTU
df <- abundances(physeq_mOTU)
df.count <- as.matrix(t(df))

## ---- PCA ----
sample_data_df <- data.frame(sample_data(physeq_mOTU), stringsAsFactors = FALSE)
sample_data_df$clinical_outcome_wk14 <- as.character(sample_data_df$clinical_outcome_wk14)
sample_data_df$clinical_outcome_wk14 <- replace_na(sample_data_df$clinical_outcome_wk14, "Donor")
sample_data_df$Timepoint_Category <- factor(sample_data_df$timepoint.new,
                                            levels = c("Pre-FMT", "FMT", "Post-FMT", "Donor"))

sample_data_df$type <- c(rep("Subject", 180), rep("Donor A", 13), rep("Donor B", 14))
sample_data_df$type <- as.factor(sample_data_df$type)

sample_data_df$type_outcome <- ifelse(sample_data_df$type == "Subject" &
                                        sample_data_df$clinical_outcome_wk14 == "Good", "Subject_Good",
                                      ifelse(sample_data_df$type == "Subject" &
                                               sample_data_df$clinical_outcome_wk14 == "None", "Subject_None",
                                             sample_data_df$type))

sample_data_df$type_outcome <- as.factor(sample_data_df$type_outcome)
sample_data(physeq_mOTU) <- sample_data_df

pca.1.2 <- physeq_mOTU %>%
  tax_transform("clr", rank = "unique") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    axes = c(1, 2),
    color = "clinical_outcome_wk14",
    shape = "type_outcome",
    plot_taxa = c("Prevotellaceae","Ruminococcaceae","Lachnospiraceae","Clostridialesfam.incertaesedis"),
    size = 3
  ) +
  scale_color_manual(
    breaks = c("Good", "None", "Donor","Partial"),
    values = c("#85C1E9", "#EC7063", "#999999", "#F4D03F")
  ) +
  scale_shape_manual(values = c("Subject_Good" = 19, "Subject_None" = 1,
                                "Donor A" = 2, "Donor B" = 17)) +
  ggtitle("PCA Based on Aitchison Distance") +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 14)) +
  labs(caption = "180 samples, 93 taxa") +
  theme_bw() +
  ggplot2::stat_ellipse(aes(colour = clinical_outcome_wk14))

pca.1.2

# ------------------------------------------------------------
# 5. Data preprocessing and Bacteroidaceae analysis
# ------------------------------------------------------------

# Load the data
# We use the compositional count abundance dataset
# to avoid convergence problems.

# ------------------------------------------------------------
# Preprocessing of the data
# ------------------------------------------------------------

# Make compositional
physeq_mOTU <- microbiome::transform(physeq_mOTU, "compositional")

# Aggregate to family level
physeq_mOTU <- aggregate_taxa(physeq_mOTU, "family")
physeq_mOTU

# Remove donors from data
physeq_mOTU <- subset_samples(physeq_mOTU, subject_id != "Donor A")
physeq_mOTU <- subset_samples(physeq_mOTU, subject_id != "Donor B")

# make new class - Coriobacteriaceae
physeq_mOTU <- merge_taxa2(physeq_mOTU,
                           pattern = "_Coriobacteriaceae",
                           name = "Coriobacteriaceae")

# ------------------------------------------------------------
# Sample data processing
# ------------------------------------------------------------
sample.data <- as.data.frame(as.matrix(physeq_mOTU@sam_data))
sample.data$timepoint.new <- as.factor(sample.data$timepoint.new)
sample.data$timepoint.new <- factor(
  sample.data$timepoint.new,
  levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2",
             "Post-3", "Post-4", "Week8", "Week10", "Week14")
)
sample.data$timepoint.new.num <- as.numeric(sample.data$timepoint.new)
sample.data$clinical_outcome_wk14 <- factor(
  sample.data$clinical_outcome_wk14,
  levels = c("None", "Good")
)
sample.data$sex <- as.factor(sample.data$sex)
sample.data$age <- as.numeric(sample.data$age)
sample.data$subject_id <- as.factor(sample.data$subject_id)
sample.data$pretreatment <- as.factor(sample.data$pretreatment)
sample.data$treated_with_donor <- as.factor(sample.data$treated_with_donor)

# Abundance data
abund <- as.data.frame(as.matrix(t(physeq_mOTU@otu_table)))

# Combine
df <- data.frame(abund, sample.data)

# ------------------------------------------------------------
# Bacteroidaceae
# ------------------------------------------------------------

# Transform the data
# hist(df$Bacteroidaceae)
df$Bacteroidaceae.new <- asin(sqrt(df$Bacteroidaceae))
# hist(df$Bacteroidaceae.new)

# Proportion of zeros
sum(df$Bacteroidaceae.new == 0) / length(df$Bacteroidaceae.new)

# Split by clinical outcome
df.none <- subset(df, clinical_outcome_wk14 == "None")
df.good <- subset(df, clinical_outcome_wk14 == "Good")

# ------------------------------------------------------------
# Plot Bacteroidaceae
# ------------------------------------------------------------

a <- xyplot(
  Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14,
  data = df.none,
  xlab = "Timepoint",
  ylab = "Relative abundance (transformed)",
  main = "Bacteroidaceae",
  panel = function(x, y) {
    panel.average(x, y, horizontal = FALSE, col = "#EC7063", lwd = 4)
  },
  key = list(
    space = "right",
    lines = list(col = c("#EC7063", "#85C1E9"), lty = c(1, 1), lwd = 2),
    text  = list(c("Non-Responders", "Responders"))
  )
)

b <- xyplot(
  Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14,
  data = df.good,
  xlab = "Timepoint",
  panel = function(x, y) {
    panel.average(x, y, horizontal = FALSE, col = "#85C1E9", lwd = 4,
                  type = "l", lty = 1)
  }
)

c <- xyplot(
  Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14,
  data = df.none, type = "p", col = "#EC7063"
)

d <- xyplot(
  Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14,
  data = df.good, type = "p", col = "#85C1E9"
)

Bacteroidaceae.plot <- a + as.layer(b) + as.layer(c) + as.layer(d)
Bacteroidaceae.plot
