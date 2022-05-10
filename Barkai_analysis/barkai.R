setwd("path/to/directory")

# libraries
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)

# read in evenness functions
source("path/to/Tissue_environment_specificity_X_perc_cis/Evenness_metrics/Evenness_functions.R")

# read in percent cis data
percent_cis <- read.csv("path/to/ASE_tidy_from_metzger_supp2.csv", header = TRUE)

# read in tf regulators
gene_regulators <- read_delim("path/to/Tissue_environment_specificity_X_perc_cis/Yeast_data/yeastract_TFregulator_counts.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

# read in standard deviation data
cpm_counts <- read.csv("path/to/Tissue_environment_specificity_X_perc_cis/Yeast_data/cpm_counts.csv")
cpm_counts %>%
  rename(Name = X)

comparative <- read_excel("./Barkai_analysis/Supplemental_Table_S7.xlsx", sheet = 2)
dataset_control <- read_excel("./Barkai_analysis/Supplemental_Table_S7.xlsx", sheet = 3)

#####################################################################
############################## DATA #################################
#####################################################################

# # read in data for Rc metric
# comparative <- read_excel("./Barkai_analysis/Supplemental_Table_S7.xlsx", sheet = 2)
# dataset_control <- read_excel("./Barkai_analysis/Supplemental_Table_S7.xlsx", sheet = 3)
# 
# Barkai999 <- read_excel("C:/Users/moonl/Documents/Wittkopp Lab/DDivergence/DataFromPapers/bioSample1to999.xlsx")
# BarkaitoEnd <- read_excel("C:/Users/moonl/Documents/Wittkopp Lab/DDivergence/DataFromPapers/bioSample1000toend.xlsx")
# Barkai_cer <- read_excel("C:/Users/moonl/Documents/Wittkopp Lab/DDivergence/DataFromPapers/Supplemental_Table_S11.xlsx", sheet =2)
# Barkai_par <- read_excel("C:/Users/moonl/Documents/Wittkopp Lab/DDivergence/DataFromPapers/Supplemental_Table_S12.xlsx", sheet =2)
# sumExp <- read_excel("C:/Users/moonl/Documents/Wittkopp Lab/Barkai_analysis/DESeq/sumExp.xlsx", sheet=1)
# 
# #create combined metadata file
# Barkai_combined <- rbind(Barkai999, BarkaitoEnd)
# Barkai_combined <- Barkai_combined %>%
#   select(-bioproject_accession, -sample_title, -strain, -isolation_source, -host,
#          -geo_loc_name)

# #cerivisiae data
# Barkai_cer <- t(Barkai_cer)
# Barkai_cer <- as.data.frame(Barkai_cer)
# colnames(Barkai_cer) <- Barkai_cer[1,]
# Barkai_cer <- Barkai_cer[-1,]
# Barkai_cer$sample_name <- c(rownames(Barkai_cer))
# Barkai_cer <- as.numeric(Barkai_cer)
# 
# # join by sample name with combined metadata
# Barkai_cer <- left_join(Barkai_cer, Barkai_combined, by= "sample_name")
# rownames(Barkai_cer) <- Barkai_cer$sample_name
# 
# # remove columns with NAs
# Barkai_cer <- Barkai_cer[, colSums(is.na(Barkai_cer)) != nrow(Barkai_cer)]
# 
# #join with experiment summary
# Barkai_cer <- left_join(Barkai_cer, sumExp, by="sample_name")
# rownames(Barkai_cer) <- Barkai_cer$sample_name
# 
# #remove extra columns
# Barkai_cer <- Barkai_cer %>%
#   select(-well_in_lib, -inner_barcode_index, -library_ID, -glycerol_stock_ID, 
#          -well_in_experiment, -experiment_ID, -condition, -strain, -"repeat",
#          -well_index, -time_point, -experiment_type, -genotype, -collected_by,
#          -sample_type, -collection_date, -isolate, -organism, -sample_name)
# #####################################################################################
# #paradoxus data
# Barkai_par <- as.data.frame(t(Barkai_par))
# 
# colnames(Barkai_par) <- Barkai_par[1,]
# Barkai_par <- Barkai_par[-1,]
# Barkai_par$sample_name <- c(rownames(Barkai_par))
# Barkai_par <- as.numeric(Barkai_par)
# 
# # join by sample name with combined metadata
# Barkai_par <- left_join(Barkai_par, Barkai_combined, by= "sample_name")
# rownames(Barkai_par) <- Barkai_par$sample_name
# 
# # remove columns with NAs
# Barkai_par <- Barkai_par[, colSums(is.na(Barkai_par)) != nrow(Barkai_par)]
# 
# #join with experiment summary
# Barkai_par <- left_join(Barkai_par, sumExp, by="sample_name")
# rownames(Barkai_par) <- Barkai_par$sample_name
# 
# #remove extra columns
# Barkai_par <- Barkai_par %>%
#   select(-well_in_lib, -inner_barcode_index, -library_ID, -glycerol_stock_ID, 
#          -well_in_experiment, -experiment_ID, -condition, -strain, -"repeat",
#          -well_index, -time_point, -experiment_type, -genotype, -collected_by,
#          -sample_type, -collection_date, -isolate, -organism, -sample_name)
# #####################################################################################
# Barkai_par <- Barkai_par %>%
#   filter(species == "par")
# 
# Barkai_cer <- Barkai_cer %>%
#   filter(species == "cer")
# #####################################################################################
# Barkai_par_sum <- Barkai_par[,5371:5375]
# par_rownames <- rownames(t(Barkai_par[,1:5370]))
# Barkai_par <- t(Barkai_par[,1:5370])
# Barkai_par <- as.data.frame(Barkai_par)
# Barkai_par <- as.data.frame(lapply(Barkai_par,as.numeric))
# Barkai_par <- as.data.frame(apply(Barkai_par,2,function(x){(x/sum(x)*1000000)}))
# 
# 
# rownames(Barkai_par) <- par_rownames
# #Barkai_par <- Barkai_par_save 
# #####################################################################################
# cer_rownames <- rownames(t(Barkai_cer[,1:5461]))
# Barkai_cer_sum <- Barkai_cer[,5462:5466]
# Barkai_cer <- t(Barkai_cer[,1:5461])
# Barkai_cer <- as.data.frame(Barkai_cer)
# Barkai_cer <- as.data.frame(lapply(Barkai_cer,as.numeric))
# Barkai_cer <- as.data.frame(apply(Barkai_cer,2,function(x){(x/sum(x)*1000000)}))
# 
# rownames(Barkai_cer) <- cer_rownames
# #####################################################################################
# 
# Barkai_cer <- calculate_tau(Barkai_cer)
# Barkai_par <- calculate_tau(Barkai_par)
# 
# # cer_tau_save <- Barkai_cer
# # par_tau_save <- Barkai_par
# # write.csv(cer_tau_save, "C:/Users/moonl/Documents/Wittkopp Lab/Barkai_analysis/cer_tau.csv", row.names=FALSE)
# # write.csv(par_tau_save, "./Barkai_analysis/par_tau.csv", row.names=TRUE)
# 
# Barkai_cer <- read.csv("C:/Users/moonl/Documents/Wittkopp Lab/Barkai_analysis/cer_tau.csv")
# Barkai_cer <- Barkai_cer[,2:ncol(Barkai_cer)]
# Barkai_cer <- tibble::rownames_to_column(Barkai_cer, "Name")
# Barkai_par <- tibble::rownames_to_column(Barkai_par, "Name")
# 
# # Below if reading in files
# Barkai_cer <- read_csv("./Barkai_analysis/cer_tau.csv")
# Barkai_par <- read_csv("./Barkai_analysis/par_tau.csv")

############################# Tau #############################

# read in Barkai data for S. cerevisiae with tau
Barkai_cer <- read_csv("./Barkai_analysis/cer_tau.csv")

# filtering tau for data of interest: genes with non-even, non-zero expression
tau_filter <- Barkai_cer %>%
  filter(!(tau==0 & PHD1_TP3_cer_C1_E4_H7_GS3==0))

############################# Binning by Tau value #############################

# create tau bin category in dataframe
df_tau_binned <- tau_filter
df_tau_binned$tau_bin <- NA

# bin genes by tau intervals of 0.25 (i.e 0-0.25, 0.25-0.50, etc.)
for (i in 1:nrow(df_tau_binned)){
  if (df_tau_binned$tau[i] < 0.25){
    df_tau_binned$tau_bin[i] <- "First quarter"
  }
  else if (df_tau_binned$tau[i] < 0.5 & df_tau_binned$tau[i] >= 0.25){
    df_tau_binned$tau_bin[i] <- "Second quarter"
  }
  else if (df_tau_binned$tau[i] < 0.75 & df_tau_binned$tau[i] >= 0.5){
    df_tau_binned$tau_bin[i] <- "Third quarter"
  }
  else if (df_tau_binned$tau[i] <= 1 & df_tau_binned$tau[i] >= 0.75){
    df_tau_binned$tau_bin[i] <- "Fourth quarter"
  }
}

# join binned and tf regulator dataframes together
df_bin_gene <- left_join(df_tau_binned, gene_regulators, by = "Name")

# omit NAs
df_bin_gene_clean <- df_bin_gene %>%
  na.omit() %>%
  unique()

# reorder tau_bins
df_bin_gene_clean$tau_bin <- fct_relevel(df_bin_gene_clean$tau_bin , levels = c("First quarter", "Second quarter", "Third quarter", "Fourth quarter"))

############################# % cis and Total difference #############################

# join with tau binned data
df_bin_cis <- left_join(df_bin_gene_clean, percent_cis, by = "Name") %>%
  na.omit() %>%
  unique()

# re-order by divergence time
df_bin_cis$comparison <- fct_relevel(df_bin_cis$comparison , levels = c("cer-cer", "cer-par", "cer-mik", "cer-bay"))

############################# Expression Variability Analysis #############################

# classifying specificity and binning based on 0 threshold
cpm_cer <- Barkai_cer
cpm_cer$min <- as.numeric(apply(cpm_cer, 1, FUN=min))

# creating specificity category: gene is CActive (Constitutively Active)
# if minimum count across environments is > 0. Gene is Non_CActive otherwise.x  
cpm_cer$specificity <- NA

for (i in 1:nrow(cpm_cer)){
  if (cpm_cer$min[i] > 0) {
    cpm_cer$specificity[i] <- "CActive"
  }
  else {
    cpm_cer$specificity[i] <- "Non_CActive"
  }
}

# redo naming
#cpm_cer <- tibble::rownames_to_column(cpm_cer, "Name")
cpm_name_spec <- as.data.frame(cpm_cer$Name)
colnames(cpm_name_spec)[1] <- "Name"

# joining specificity column with comparison data
cpm_name_spec$specificity <- cpm_cer$specificity
df_bin_cis <- left_join(df_bin_cis, cpm_name_spec, by = "Name")
df_bin_cis$class <- NA

# separating into classes by specificity and tau:
# For CActive genes, if tau >0.6, classify as CA_EnvVar (Constitutively
# Active, Environmentally Variable). If tau <0.6, classify as CA_EnvConst
# (Constitutively Active, Environmentally COnstant). All other genes 
# continue to be characterized as nonCA (non-Constitutively Active)
for (i in 1:nrow(df_bin_cis)) {
  if (df_bin_cis$tau[i] > 0.6) {
    if (df_bin_cis$specificity[i] == "CActive") {
      df_bin_cis$class[i] <- "CA_EnvVar"
    }
    else {
      df_bin_cis$class[i] <- "nonCA"
    }
  }
  else {
    if (df_bin_cis$specificity[i] == "CActive") {
      df_bin_cis$class[i] <- "CA_EnvConst"
    }
    else {
      df_bin_cis$class[i] <- "nonCA"
    }
  }
}

############################# Rc Metric #############################

# data for Rc metric
comparative <- comparative %>%
  rename(Name = ORF) %>%
  drop_na()
dataset_control <- dataset_control %>%
  rename(Name = ORF) %>%
  drop_na()

# join tau binned data with standard deviation data; filter genes with
# SD <0.2
df_bin_cis <- left_join(df_bin_cis, cpm_counts, by = "Name")
df_bin_cis <- df_bin_cis %>%
  filter(tau_SD<0.2)

# join all data together
df_all <- left_join(df_bin_cis, comparative, by = "Name") %>%
  na.omit() %>%
  unique()
df_all <- left_join(df_all, dataset_control, by = "Name") %>%
  na.omit() %>%
  unique()

# filter data for looking at "diverged" genes (using thresholds form
# Barkai paper)
df_all_filter <- df_all %>%
  filter(bw_cer_par<0.2 & cer>0.5 & par>0.5)

#####################################################################
############################# PLOTS #################################
#####################################################################

# histogram of genes by unfiltered tau
p_tau_unfiltered <- ggplot(Barkai_cer, aes(x=tau)) + 
  geom_histogram(fill = "#69b3a2", color = "#e9ecef", alpha=0.9) +
  ggtitle("Tau measurement of Yeast genes in different environments") + 
  theme(plot.title = element_text(size=15))
p_tau_unfiltered

# histogram of genes by filtered tau
p_tau_filtered <- ggplot(tau_filter, aes(x=tau)) +
  geom_histogram(fill = "#69b3a2", color = "#e9ecef", alpha=0.9) +
  ggtitle("FIltered Tau measurement of Yeast genes in different environments") +
  theme(plot.title = element_text(size=15))
p_tau_filtered

# boxplot of tau bins vs. number of transcription factor regulators
p_bin_regulator_box <- ggplot(df_bin_gene_clean, aes(x=tau, y=Num_of_regulators, fill=tau_bin)) +
  geom_boxplot(notch=TRUE) 
p_bin_regulator_box

# scatter plot of tau bins vs. number of transcription factor regulators
p_bin_regulator_scatter <- ggplot(df_bin_gene_clean, aes(x=tau, y=Num_of_regulators)) +
  geom_point()
p_bin_regulator_scatter

# boxplot of tau bins vs. % cis, grouped by comparison of yeast species
p_bin_cis <- ggplot(df_bin_cis, aes(x=tau, y=cis_percent, fill=tau_bin)) +
  geom_boxplot(notch=TRUE) + 
  facet_wrap(~comparison)
p_bin_cis

# boxplot of species comparison vs. % cis, grouped by tau bin
p_comparison_cis <- ggplot(df_bin_cis, aes(x=comparison, y=cis_percent)) +
  geom_boxplot(notch=TRUE) + 
  facet_wrap(~tau_bin)
p_comparison_cis

# boxplot of species comparison vs. total difference (trans and cis)
# grouped by tau bin
p_comparison_diff <- ggplot(df_bin_cis, aes(x=comparison, y=abs(total_diff))) +
  geom_boxplot(notch=TRUE) + 
  facet_wrap(~tau_bin)
p_comparison_diff

# boxplot of species compariosn vs. % cis, grouped by class and tau bin
p_expvar_cis <- ggplot(df_bin_cis, aes(x=comparison, y=cis_percent, fill=class,
                                       alpha=comparison)) +
  geom_boxplot(notch=TRUE) + 
  facet_wrap(~class, labeller = as_labeller(c('CA_EnvConst' = "Const. Active, Env. Const.", 
                                              'CA_EnvVar' = "Const. Active, Env. Variable", 'nonCA' = "Non-Const. Active"))) +
  scale_fill_discrete(guide=FALSE) +
  ylab("% cis") + 
  xlab("Species Comparison") + 
  ggtitle("Species Comparison of Variably Expressed Genes vs. % cis")
p_expvar_cis

# boxplot of class vs. % cis, grouped by species comparison?
p_class_cis <- ggplot(df_bin_cis, aes(x=class, y=cis_percent, fill=class)) +
  geom_boxplot(notch=TRUE) + 
  facet_wrap(~comparison) +
  scale_fill_discrete(guide = FALSE) +
  ylab("% cis") + 
  theme(axis.title.x = element_blank()) + 
  ggtitle("Expression Variability vs. % cis by diverging species of Yeast")
p_class_cis

# boxplot of class vs. log of abs. value of total difference
p_ca_diff <- ggplot(df_bin_cis, aes(x=class, y=log(abs(total_diff)), fill=class)) +
  geom_boxplot(notch=TRUE) +
  facet_wrap(~comparison) +
  scale_fill_discrete(guide = FALSE) +
  ylab("abs(total diff)") +
  theme(axis.title.x = element_blank())
p_ca_diff

# boxplot of class vs. RC; unfiltered
p_unfiltered_class_rc <- ggplot(df_all, aes(x=class, y=bw_cer_par, fill=class)) +
  geom_boxplot(notch=TRUE) +
  ggtitle("Barkai Unfiltered: Class vs. Rc")
p_unfiltered_class_rc

# boxplot of tau bins vs. RC; unfiltered
p_unfiltered_bin_rc <- ggplot(df_all, aes(x=tau_bin, y=bw_cer_par, fill=tau_bin)) +
  geom_boxplot(notch=TRUE) +
  ggtitle("Barkai: Tau vs. Rc") +
  ylab("Rc") + 
  xlab("Tau") +
  theme(legend.position = "none")
p_unfiltered_bin_rc

df_all %>%
  group_by(tau_bin) %>%
  tally()

# boxplot of class vs. RC; filtered
p_filter_class_rc <- ggplot(df_all_filter, aes(x=class, y=bw_cer_par, fill=class)) +
  geom_boxplot(notch=TRUE) +
  ggtitle("Barkai Filtered: Class vs. Rc")
p_filter_class_rc

# boxplot of tau bins vs. RC; filtered
p_filter_bin_rc <- ggplot(df_all_filter, aes(x=tau_bin, y=bw_cer_par, fill=tau_bin)) +
  geom_boxplot(notch=TRUE) +
  ggtitle("Barkai Filtered: Tau vs. Rc")
p_filter_bin_rc

