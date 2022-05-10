library(tidyverse)
library(ggplot2)
library(igraph)
library(network3d)

setwd("path/to/directory")
network <- read.table(file = 'elife-51254-code3-v3.tsv', sep = '\t', header=TRUE)
allele.sp <- read.csv("C:/Users/moonl/Documents/Wittkopp Lab/DDivergence/GeneratedFiles/ASE_tidy_from_metzger_supp2.csv", header = TRUE)

allele.sp <- allele.sp[,2:ncol(allele.sp)]
allele.sp <- allele.sp[allele.sp$comparison == "cer-cer",] # choose cer-par or cer-cer

network.sv <- network[,1]
network[,2:ncol(network)] <- abs(network[,2:ncol(network)])

network$sums <- rowSums(network[,2:ncol(network)])
network <- network %>%
  rename(Name = target)
###########################################################################
# binomial testing
# change comparison to either cer-cer or cer-par

pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;


for (i in 1:nrow(allele.sp))
{
  Par1_Bi <- binom.test(allele.sp[[i,3]], (allele.sp[[i,3]]+allele.sp[[i,4]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  Hyb1_Bi <- binom.test(allele.sp[[i,5]], (allele.sp[[i,5]]+allele.sp[[i,6]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  # collect p-values from binomial tests
  pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
  pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);
}
#FDR correct pvalues
pvalsHyb1.adj <- p.adjust(pvalsHyb1,method="fdr");
pvalsPar1.adj <- p.adjust(pvalsPar1,method="fdr");
#append on adjusted pvals
allele.sp <- cbind(allele.sp,pvalsPar1.adj,pvalsHyb1.adj)

# pval distribution plots
ggplot(allele.sp, aes(x=pvalsPar1.adj)) + geom_histogram() + xlab("p-values: Parental") + ggtitle("Distribution of Parental p-values from Binomial Test")
ggplot(allele.sp, aes(x=pvalsHyb1.adj)) + geom_histogram() + xlab("p-values: Hybrid") + ggtitle("Distribution of Hybrid p-values from Binomial Test")
###########################################################################
# in-degree plots

# Hybrid differential expression: Yes or No based on 0.05 significance cutoff
allele.sp$DiffYN <- ""
for (i in 1:nrow(allele.sp)) {
  if (allele.sp$pvalsHyb1.adj[i] < 0.05) {
    allele.sp$DiffYN[i] <- "Yes"
  }
  else {
    allele.sp$DiffYN[i] <- "No"
  }
}

# join Gresham network data with allele specific and filter
network <- left_join(network, allele.sp, by = "Name")
network <- network %>%
  na.omit() %>%
  unique()

# Parental differential expression: Yes or No based on 0.05 significance cutoff
network$DiffYN_par <- ""
for (i in 1:nrow(network)) {
  if (network$pvalsPar1.adj[i] < 0.05) {
    network$DiffYN_par[i] <- "Yes"
  }
  else {
    network$DiffYN_par[i] <- "No"
  }
}


#write.csv(network, "./network.csv", row.names= FALSE)

# plots for distribution of in-degree nodes vs. parental/hybrid diff. expressed genes
ggplot(network, aes(x=DiffYN_par, y=sums)) + geom_boxplot(notch=TRUE) + xlab("Expression Differences") + ylab("In-degree") + ggtitle("cer-cer hybrid expression differences vs. in-degree"
)
ggplot(network, aes(x=DiffYN, y=sums)) + geom_boxplot(notch=TRUE) + xlab("Expression Differences") + ylab("In-degree") + ggtitle("cer-cer parental expression differences vs. in-degree"
)
########################################################################
# out-degree plots

# transpose gresham data to make network.out matrix
network.out <- read.table(file = 'elife-51254-code3-v3.tsv', sep = '\t', header=TRUE)
network.out[,2:ncol(network.out)] <- abs(network.out[,2:ncol(network.out)])
network.out$sums <- rowSums(network.out[,2:ncol(network.out)])
network.out <- network.out %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))
target <- network.out$target
network.out <- as.data.frame(t(network.out[,-1]))
colnames(network.out) <- target
network.out$Name <- factor(row.names(network.out))

# join with allele specific data and filter
network.out <- left_join(network.out, allele.sp, by = "Name")
network.out <- network.out %>%
  na.omit() %>%
  unique() 

# Parental differential expression: Yes or No based on 0.05 significance cutoff
network.out$DiffYN_par <- ""
for (i in 1:nrow(network.out)) {
  if (network.out$pvalsPar1.adj[i] < 0.05) {
    network.out$DiffYN_par[i] <- "Yes"
  }
  else {
    network.out$DiffYN_par[i] <- "No"
  }
}

# plots for distribution of out-degree nodes vs. parental/hybrid diff. expressed genes
ggplot(network.out, aes(x=DiffYN_par, y=Total)) + geom_boxplot(notch=TRUE)+ xlab("Expression Differences") + ylab("Out-degree") + ggtitle("cer-cer parental expression differences vs. out-degree"
)
ggplot(network.out, aes(x=DiffYN, y=Total)) + geom_boxplot(notch=TRUE)+ xlab("Expression Differences") + ylab("Out-degree") + ggtitle("cer-cer hybrid expression differences vs. out-degree"
)

##########################################################################
# recreating figure 1a and 1d from Bing & Wittkopp

# read in proportions data from python script
proportions <- read.csv("./results.csv",)
proportions <- proportions[,1:ncol(proportions)-1]
proportions <- as.data.frame(t(proportions))
proportions <- proportions %>%
  rename(ratio=V1)
proportions$Name <- rownames(proportions)

# join with DiffYN_par (from allele.sp) data and filter
DiffYN_par <- as.data.frame(network$DiffYN_par)
DiffYN_par$Name <- network$Name
proportions <- left_join(proportions, DiffYN_par, by="Name") %>%
  na.omit() %>%
  unique()

proportions <- proportions %>%
  rename(YN='network$DiffYN_par')

# recreation of figure 1D with gresham network
ggplot(proportions, aes(x=YN, y=ratio)) + geom_boxplot(notch=TRUE) + xlab("Expression Differences") +
  ylab("Proportion of regulators changing expression") + ggtitle("cer-cer expression differences vs. proportion of regulators ")

# same as above, but for figure 1A with gresham network
proportions_in <- read.csv("./results_in.csv")
proportions_in <- proportions_in[,1:ncol(proportions_in)-1]
proportions_in <- as.data.frame(t(proportions_in))
proportions_in <- proportions_in %>%
  rename(ratio=V1)
proportions_in$Name <- rownames(proportions_in)

proportions_in <- left_join(proportions_in, DiffYN_par, by="Name") %>%
  na.omit() %>%
  unique()

proportions_in <- proportions_in %>%
  rename(YN='network$DiffYN_par')

# recreation of figure 1A with gresham network
ggplot(proportions_in, aes(x=YN, y=ratio)) + geom_boxplot(notch=TRUE)+ xlab("Expression Differences") +
  ylab("Proportion of targets changing expression") + ggtitle("cer-cer expression differences vs. proportion of targets")

proportions_in %>%
  group_by(YN) %>%
  tally()
#######################################################################
# same objectives as code block above, but using "gold standard" network (Tchourine et al. 2018)
gold.net <- read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
gold.net <- abs(gold.net)
gold.net$sums <- rowSums(gold.net)

ggplot(gold.net, aes(x=sums)) + geom_histogram() + xlab("In-degree") +
  ggtitle("Distribution of targets by in-degree number")



# recreation of figure 1A using gold standard network; transpose data for TFs
gold.tf <- as.data.frame(t(read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
))
gold.tf <- abs(gold.tf)
gold.tf$sums <- rowSums(gold.tf)

ggplot(gold.tf, aes(x=sums)) + geom_histogram(binwidth = 5) + xlab("Out-degree") +
  ggtitle("Distribution of transcription factors by out-degree number")
gold.tf %>%
  group_by(sums) %>%
  tally()
gold.tf$Name <- rownames(gold.tf)
gold.tf <- left_join(gold.tf, allele.sp) %>%
  na.omit() %>%
  unique()

# assign YN for parental differentially expressed genes
gold.tf$DiffYN_par <- ""
for (i in 1:nrow(gold.tf)) {
  if (gold.tf$pvalsPar1.adj[i] < 0.05) {
    gold.tf$DiffYN_par[i] <- "Yes"
  }
  else {
    gold.tf$DiffYN_par[i] <- "No"
  }
}

ggplot(gold.tf, aes(x=DiffYN_par, y=sums)) + geom_boxplot(notch=TRUE)

#write.csv(gold.tf, "./goldtf.csv", row.names=TRUE)

gold.net <- abs(gold.net)
gold.net$Name <- rownames(gold.net)
gold.net <- left_join(gold.net, allele.sp, by="Name") %>%
  na.omit() %>%
  unique()

gold.net$DiffYN_par <- ""
for (i in 1:nrow(gold.net)) {
  if (gold.net$pvalsPar1.adj[i] < 0.05) {
    gold.net$DiffYN_par[i] <- "Yes"
  }
  else {
    gold.net$DiffYN_par[i] <- "No"
  }
}

# produce network matrix for ratio generation using python script
#write.csv(gold.net, "./goldnet.csv", row.names=TRUE)

# read in gold standard network ratio data from python script
gold.ratio <- read.csv("./goldratio.csv")
gold.ratio <- as.data.frame(t(gold.ratio))
gold.ratio$Name <- rownames(gold.ratio)
gold.ratio <- gold.ratio %>%
  rename(ratio=V1)

# join ratio data with tf matrix
gold.tf <- left_join(gold.tf, gold.ratio, by="Name") %>%
  na.omit() %>%
  unique()

# plot YN Diff. expressed genes against ratio of targets changing expression
ggplot(gold.tf, aes(x=DiffYN_par, y=ratio)) + geom_boxplot(notch=TRUE) + xlab("Expression Differences")+
  ylab("Proportion of targets changing expression") + ggtitle("cer-cer expression differences vs. proportion of targets")

# find n for Yes/No categories: 67 for Yes, 37 for No
gold.tf %>%
  group_by(DiffYN_par) %>%
  tally()
#######################################################################
network.graph <- read.table(file = 'elife-51254-code3-v3.tsv', sep = '\t', header=TRUE)
network.graph[,2:ncol(network.graph)] <- abs(network.graph[,2:ncol(network.graph)])

counter <- 1
bip.net <- as.data.frame(matrix())
for (i in 1:nrow(network.graph)) {
  for (j in 2:ncol(network.graph)) {
    if (network.graph[i,j] == 1) {
      bip.net[counter,1] <- network.graph$target[i]
      bip.net[counter,2] <- colnames(network.graph[j])
      counter <- counter+1
    }
  }
}

bip.net <- as_tibble(bip.net)
bip.net <- bip.net %>%
  rename(source=V1) %>%
  rename(target=V2)

edges <- bip.net


vertices <- network.graph$target
vertices <- as_tibble(vertices)
vertices <- vertices %>%
  rename(id=value)
vertices$name <- c(1:2445)
vertices$name <- as.character(vertices$name)


data <- list(edges, vertices)
names(data) <- c("edges", "vertices")
data$vertices %>% head() %>% knitr::kable()
data$edges %>% head() %>% knitr::kable()
network3d(data$vertices, data$edges)

data2 <- collaboration_networks
data2$vertices %>% head() %>% knitr::kable()
data2$edges %>% head() %>% knitr::kable()
network3d(data2$vertices, data2$edges)

#i don't know why this isn't working :(
#########################################################################


