# ========================
#      1. Load Data
# ========================

library(ape)
library(cluster)
library(dplyr)
library(fpc)
library(ggplot2)
library(lubridate)
library(stringr)
library(tidyverse)
library(vegan)


otu_tab <- read.csv("./genus_counts.csv", header = TRUE, row.names = 1)
sample_df <- read.csv("./sample_data.csv", header = TRUE)

# ========================
#   2. Data Description 
# ========================

# ------------------------
#    2.1. Sample data
# ------------------------

# Variables of sample data
colnames(sample_df)

# Type of Sequence Techniques
SeqTech <- unique(sample_df$SeqTech)
data.frame(SeqTech)
table(sample_df$SeqTech)

# Sample ID
ID_sample <- unique(sample_df$SampleID)[c(1, 2, 99, 100, 145, 267)]
ID_sample

# Project
project <- unique(sample_df$Project)
project
table(sample_df$Project)

# Clinical Status
c_status <- unique(sample_df$ClinicalStatus)
c_status
table(sample_df$ClinicalStatus)

# ------------------------
#  2.2. Genus Counts (OTU)
# ------------------------

otu_matrix <- as.matrix(otu_tab)


# ========================
#   3. Client's Approach 
# ========================

# --------------------------
# 3.1. Bray-Curtis Distance
# --------------------------

pyro_indices = which(sample_df$SeqTech == "Pyro454")
otu_tab_pyro = otu_tab[,pyro_indices]
sample_df_pyro = sample_df[pyro_indices,]

d_pyro = vegdist(t(otu_tab_pyro), method = "bray")

as.matrix(d_pyro)[1:5, 1:5]

# --------------------------------
# 3.2. PAMK & Optimal # of cluster
# --------------------------------

pamk_fit_pyro = pamk(d_pyro, krange = 1:6)

# --------------------------------
#     3.3. Silhouette Score
# --------------------------------

# Silhouette Score
plot(pamk_fit_pyro$crit)

# Optimal Number of Clusters
pamk_fit_pyro$nc

# ================================
#           4.PERMANOVA
# ================================

pam_fit_pyro = pam(d_pyro, k = pamk_fit_pyro$nc)
cluster_group_pyro = as.factor(pam_fit_pyro$clustering)
PERM_pyro <- adonis2(d_pyro ~ cluster_group_pyro, permutations = 999)
print(PERM_pyro)

# --------------------------------
#     4.1. Visualization (PCoA)
# --------------------------------

clustering_pyro <- pamk_fit_pyro$pamobject$clustering

cluster_df_pyro <- data.frame(
  SampleID = names(clustering_pyro),
  Cluster_pyro = as.factor(clustering_pyro)
)
rownames(cluster_df_pyro) <- NULL

sample_cluster_df_pyro <- sample_df_pyro |> 
  left_join(cluster_df_pyro, by = "SampleID")

otu_tab_pyro <- as.data.frame(t(otu_tab_pyro))

otu_tab_pyro <- otu_tab_pyro %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sample_cluster_df_pyro %>% select(SampleID, Cluster_pyro), by = "SampleID")

rownames(otu_tab_pyro) <- otu_tab_pyro$SampleID
otu_tab_pyro$SampleID <- NULL

otu_tab_pyro <- otu_tab_pyro %>%
  relocate(Cluster_pyro, .after = last_col())

pcoa_res <- pcoa(d_pyro)
pcoa_scores <- as.data.frame(pcoa_res$vectors[,1:2])

pcoa_scores$SampleID <- rownames(pcoa_scores)

pcoa_scores <- pcoa_scores |>
  left_join(
    otu_tab_pyro |> 
      rownames_to_column("SampleID") |> 
      select(SampleID, Cluster_pyro),
    by = "SampleID"
  )

ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2, color = Cluster_pyro)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis) of Microbiome Composition")


# Figure 3. Boxplot
axis1_df <- pcoa_res$vectors[, 1, drop = FALSE] |>
  as.data.frame() |>
  rownames_to_column("SampleID") |>
  rename(Axis1 = Axis.1) |>
  left_join(
    otu_tab_pyro |>
      rownames_to_column("SampleID") |>
      select(SampleID, Cluster_pyro),
    by = "SampleID"
  )

axis2_df <- pcoa_res$vectors[, 2, drop = FALSE] |>
  as.data.frame() |>
  rownames_to_column("SampleID") |>
  rename(Axis2 = Axis.2) |>
  left_join(
    otu_tab_pyro |>
      rownames_to_column("SampleID") |>
      select(SampleID, Cluster_pyro),
    by = "SampleID"
  ) |>
  filter(!is.na(Cluster_pyro))

# 3. combine two axes
axis_df <- axis1_df |>
  left_join(axis2_df |> select(SampleID, Axis2), by = "SampleID") |>
  filter(!is.na(Cluster_pyro))

# 4. Long format
axis_long <- axis_df |>
  pivot_longer(cols = c(Axis1, Axis2), names_to = "Axis", values_to = "Score")

# 5. Boxplot + Facet
ggplot(axis_long, aes(x = Cluster_pyro, y = Score, fill = Cluster_pyro)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  facet_wrap(~ Axis, scales = "fixed") + 
  theme_minimal() +
  labs(
    title = "Separation Along PCoA Axes 1 and 2",
    x = "Cluster",
    y = "PCoA Score"
  ) +
  theme(legend.position = "none")


# =========================
# 5. Influence of Variables
# =========================

# 5.1. Clustering 
sanger_indices = which(sample_df$SeqTech == "Sanger")
otu_tab_sanger = otu_tab[,sanger_indices]
sample_df_sanger = sample_df[sanger_indices,]

d_sanger = vegdist(t(otu_tab_sanger), method = "bray")

pamk_fit_sanger = pamk(d_sanger, krange = 1:6)

plot(pamk_fit_sanger$crit)

clustering_sanger <- pamk_fit_sanger$pamobject$clustering

cluster_df_sanger <- data.frame(
  SampleID = names(clustering_sanger),
  Cluster_sanger = as.factor(clustering_sanger)
)

rownames(cluster_df_sanger) <- NULL

sample_cluster_df_sanger <- sample_df_sanger |> 
  left_join(cluster_df_sanger, by = "SampleID")

otu_tab_sanger <- as.data.frame(t(otu_tab_sanger))


otu_tab_sanger <- otu_tab_sanger %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sample_cluster_df_sanger %>% select(SampleID, Cluster_sanger), by = "SampleID")

rownames(otu_tab_sanger) <- otu_tab_sanger$SampleID
otu_tab_sanger$SampleID <- NULL

otu_tab_sanger <- otu_tab_sanger |>
  relocate(Cluster_sanger, .after = last_col())

otu_rel_abund_sanger <- otu_tab_sanger |> 
  select(-Cluster_sanger) |> 
  apply(1, function(x) x / sum(x))|> 
  t() |> 
  as.data.frame()


otu_rel_abund_sanger$Cluster_sanger <- otu_tab_sanger$Cluster_sanger


pcoa_res_sanger <- pcoa(d_sanger)
pcoa_scores_sanger <- as.data.frame(pcoa_res_sanger$vectors[,1:2])
pcoa_scores_sanger$Cluster_sanger <- otu_tab_sanger$Cluster_sanger

ggplot(pcoa_scores_sanger, aes(x = Axis.1, y = Axis.2, color = Cluster_sanger)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis) of Microbiome Composition")


# 5.2. Significance test (Sanger)
pam_fit_sanger = pam(d_sanger, k = pamk_fit_sanger$nc)

cluster_group_sanger = as.factor(pam_fit_sanger$clustering)

adonis_res_sanger <- adonis2(d_sanger ~ cluster_group_sanger, permutations = 999)

print(adonis_res_sanger)

# 5.3. chi-square test and ANOVA by variables
sample_df_sanger$Cluster_sanger <- factor(cluster_group_sanger)

# Nationality vs Cluster
table_nat <- table(sample_df_sanger$Nationality, sample_df_sanger$Cluster_sanger)
chisq.test(table_nat)

# Gender vs Cluster
table_gender <- table(sample_df_sanger$Gender, sample_df_sanger$Cluster_sanger)
chisq.test(table_gender)

# ClinicalStatus vs Cluster
table_status <- table(sample_df_sanger$ClinicalStatus, sample_df_sanger$Cluster_sanger)
chisq.test(table_status)

# 5.3.1. Monte Carlo Simulation with categorical variables
chisq.test(table_nat, simulate.p.value = TRUE, B = 10000)
chisq.test(table_gender, simulate.p.value = TRUE, B = 10000)
chisq.test(table_status, simulate.p.value = TRUE, B = 10000)

# 5.3.2. ANOVA test for age
anova_age <- aov(Age ~ Cluster_sanger, data = sample_df_sanger)
summary(anova_age)

# Residuals for ANOVA
residuals_anova <- residuals(anova_age)

qqnorm(residuals_anova)
qqline(residuals_anova, col = "red", lwd = 2)

# 5.4. Relationship between Seqtech and Clustering
d_all = vegdist(t(otu_tab), method = "bray")
pamk_fit_all = pamk(d_all, krange = 1:7)
plot(pamk_fit_all$crit)

# 5.4.1. Visualization
clustering_all <- pamk_fit_all$pamobject$clustering

cluster_df_all <- data.frame(
  SampleID = names(clustering_all),
  Cluster_all = as.factor(clustering_all)
)
rownames(cluster_df_all) <- NULL


sample_cluster_df_all <- sample_df |> 
  left_join(cluster_df_all, by = "SampleID")

otu_tab <- as.data.frame(t(otu_tab))


otu_tab <- otu_tab |> 
  rownames_to_column(var = "SampleID") |> 
  left_join(sample_cluster_df_all |>  select(SampleID, Cluster_all), by = "SampleID")

rownames(otu_tab) <- otu_tab$SampleID
otu_tab$SampleID <- NULL

otu_tab <- otu_tab |> 
  relocate(Cluster_all, .after = last_col())

otu_rel_abund_all <- otu_tab |> 
  select(-Cluster_all) |> 
  apply(1, function(x) x / sum(x))|> 
  t() |> 
  as.data.frame()

otu_rel_abund_all$Cluster_all <- otu_tab$Cluster_all

pcoa_res <- pcoa(d_all)
pcoa_scores <- as.data.frame(pcoa_res$vectors[,1:2])
pcoa_scores$Cluster_all <- otu_tab$Cluster_all

ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2, color = Cluster_all)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis) of Microbiome Composition")

# 5.4.2. Significance test
pam_fit_all = pam(d_all, k = pamk_fit_all$nc)

cluster_group_all = as.factor(pam_fit_all$clustering)

adonis_res_all <- adonis2(d_all ~ cluster_group_all, permutations = 999)

print(adonis_res_all)

# 5.4.2. chi-square test for Seqtech
table_seq <- table(sample_cluster_df_all$SeqTech, sample_cluster_df_all$Cluster_all)
chisq.test(table_seq)

# ==========================
# 6. Checking zero-inflation
# ==========================

zero_prop_by_sample <- colSums(otu_matrix == 0) / nrow(otu_matrix)
summary(zero_prop_by_sample)




