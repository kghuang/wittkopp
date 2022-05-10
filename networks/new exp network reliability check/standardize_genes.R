library(tidyverse)
library(stringr)
setwd("path/to/directory")

cer267 <- read.table("./cer_267.txt", sep="\t")
par.both <- read.table("./par_264_268.txt", sep="\t")
sacgenes <- read.table("./sacgenes.tsv", header=TRUE, sep="\t")

newgenes <- matrix(ncol=1, nrow=7417)
for (i in 1:nrow(cer267)) {
  newgenes[i] <- unlist(str_extract_all(cer267[i,4], "Parent=.+?;"))
}
newgenes <- as.data.frame(newgenes)

newgenes2 <- matrix(ncol=1, nrow=7417)
for (j in 1:nrow(newgenes)) {
  newgenes2[j] <- unlist(str_extract_all(newgenes[j,1], "Scer.+?;"))
}
newgenes2 <- as.data.frame(newgenes2)

str_sub(newgenes2[1,1],1,nchar(newgenes2[1,1])-1)

i <- 0
for (i in 1:nrow(cer267)) {
  cer267[i,4] <- str_sub(newgenes2[i,1],1,nchar(newgenes2[i,1])-1)
}

k <- 0
for (k in 1:nrow(par.both)) {
  par.both[k,4] <- substring(par.both[k,4], 8)
}
rm(newgenes)
rm(newgenes2)

cer267 <- cer267 %>%
  rename(Scer_Name=V4,sc267=V5)
par.both <- par.both %>%
  rename(Spar_Name=V4,sp268=V6)
##########################################################################
sacgenes <- left_join(sacgenes, cer267, by = "Scer_Name", copy=TRUE)
sacgenes <- sacgenes %>%
  select(-c(Scer_Gene, V1, V2, V3))
sacgenes <- left_join(sacgenes, par.both, by="Spar_Name", copy=TRUE)
sacgenes <- sacgenes %>%
  select(-c(V1, V2, V3, V5))
sacgenes <- sacgenes %>%
  select(-Spar_Name)

sacgenes <- na.omit(sacgenes)

sacgenesY <- read.table("./sacgenes.tsv", header=TRUE, sep="\t")
sacgenes <- left_join(sacgenes, sacgenesY, by = "Scer_Name", copy=TRUE)
sacgenes <- sacgenes %>%
  select(-c(Spar_Name, Scer_Name))
##########################################################################
cpm_genes <- as.data.frame(sacgenes$sc267*1000000/44099617)
cpm_genes <- cpm_genes %>%
  rename(sc267='sacgenes$sc267 * 1e+06/44099617')
cpm_genes$sp268 <- sacgenes$sp268*1000000/47740035
cpm_genes$Name <- sacgenes$Scer_Gene
##########################################################################
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;

cpm_genes$sc267 <- round(cpm_genes$sc267, digits = 0)
cpm_genes$sp268 <- round(cpm_genes$sp268, digits = 0)

#write.table(cpm_genes, "./cpm_genes.tsv", sep = "\t")

cpm_genes[cpm_genes==0] <- NA
cpm_genes <- na.omit(cpm_genes)



for (i in 1:nrow(cpm_genes)) {
  Par1_Bi <- binom.test(cpm_genes[[i,1]], (cpm_genes[[i,1]]+cpm_genes[[i,2]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  # collect p-values from binomial tests
  pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);
}
pvalsPar1.adj <- p.adjust(pvalsPar1,method="fdr");
cpm_genes <- cbind(cpm_genes,pvalsPar1.adj)

ggplot(cpm_genes, aes(x=pvalsPar1.adj)) + geom_histogram() + xlab("p-values") + ggtitle("Distribution of p-values from Binomial Test")
##########################################################################
network <- read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
network <- abs(network)
cpm_genes$DiffYN <- ""
for (i in 1:nrow(cpm_genes)) {
  if (cpm_genes$pvalsPar1.adj[i] < 0.05) {
    cpm_genes$DiffYN[i] <- "Yes"
  }
  else {
    cpm_genes$DiffYN[i] <- "No"
  }
}
network$Name <- rownames(network)
network <- left_join(network, cpm_genes, by="Name")
network <- na.omit(network)
network <- unique(network)
write_tsv(network, "./new_network_out.tsv")

net.out.ratio <- read.table(file="./new_out_results.csv", sep=',')
net.out.ratio <- as.data.frame(t(net.out.ratio))
net.out.ratio <- net.out.ratio %>%
  rename(Name=V1, ratio=V2)

net.out.ratio <- left_join(net.out.ratio, cpm_genes, by="Name")
net.out.ratio <- net.out.ratio %>%
  na.omit() %>%
  unique()

net.out.ratio <- as.data.frame(net.out.ratio)
net.out.ratio$ratio <- as.numeric(net.out.ratio$ratio)

ggplot(net.out.ratio, aes(x=DiffYN, y=ratio)) + geom_boxplot()+ xlab("Expression Differences") +
  ylab("Proportion of regulators changing expression") + ggtitle("cer-par expression differences vs. proportion of regulators ")+
  stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red")


net.out.ratio %>%
  group_by(DiffYN) %>%
  tally()

net.yes <- net.out.ratio %>%
  filter(DiffYN=="Yes")
mean(net.yes$ratio)

net.no <- net.out.ratio %>%
  filter(DiffYN=="No")
mean(net.no$ratio)

wilcox.test(net.yes$ratio, net.no$ratio)

net.yes %>%
  filter(ratio==0.5) %>%
  tally()
net.no %>%
  filter(ratio==0.5) %>%
  tally()

sums <- colSums(network[,1:198])
sums <- as.data.frame(t(sums))
sums %>%
  group_by(V1) %>%
  tally()

ggplot(net.out.ratio, aes(x=ratio, y=abs(log2(sc267/sp268)))) + geom_point()
