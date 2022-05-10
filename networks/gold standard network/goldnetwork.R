# replicating network.R plots with Tchourine et al. network
library(tidyverse)
library(ggplot2)
library(igraph)
library(network3d)

setwd("path/to/directory")
network <- read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
allele.sp <- read.csv("C:/Users/moonl/Documents/Wittkopp Lab/DDivergence/GeneratedFiles/ASE_tidy_from_metzger_supp2.csv", header = TRUE)

allele.sp <- allele.sp[,2:ncol(allele.sp)]
allele.sp <- allele.sp[allele.sp$comparison == "cer-cer",] # choose cer-par or cer-cer

network <- abs(network)

network$sums <- rowSums(network)
network$Name <- rownames(network)
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
network.out <- read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
network.out <- abs(network.out)
network.out <- network.out %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))
#target <- rownames(network.out)
network.out <- as.data.frame(t(network.out))
#colnames(network.out) <- target
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

network.out <- network.out %>%
  rename(Total=...1395)
# plots for distribution of out-degree nodes vs. parental/hybrid diff. expressed genes
ggplot(network.out, aes(x=DiffYN_par, y=log(Total))) + geom_boxplot(notch=TRUE)+ xlab("Expression Differences") + ylab("log Out-degree") + ggtitle("cer-cer parental expression differences vs. out-degree"
)
ggplot(network.out, aes(x=DiffYN, y=log(Total))) + geom_boxplot(notch=TRUE)+ xlab("Expression Differences") + ylab("log Out-degree") + ggtitle("cer-cer hybrid expression differences vs. out-degree"
)

network.out %>%
  group_by(DiffYN_par) %>%
  tally()
network.out %>%
  group_by(DiffYN) %>%
  tally()

#########################################################################

network2 <- read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
network2 <- abs(network2)
network2 <- as.data.frame(t(network2))
network2$Name <- factor(row.names(network2))

pval.name <- as.data.frame(matrix(nrow=3133,ncol=2))
pval.name[,1] <- allele.sp$Name
pval.name[,2] <- allele.sp$pvalsPar1.adj
pval.name <- pval.name %>%
  rename(Name=V1)

network2 <- left_join(network2, pval.name, by = "Name")
network2 <- network2 %>%
  na.omit() %>%
  unique() 
rownames(network2) <- network2$Name
network2 <- network2 %>%
  rename(pvalsPar1.adj=V2)

network2$DiffYN_par <- ""
for (i in 1:nrow(network2)) {
  if (network2$pvalsPar1.adj[i] < 0.05) {
    network2$DiffYN_par[i] <- "Yes"
  }
  else {
    network2$DiffYN_par[i] <- "No"
  }
}

network2 <- network2 %>%
  select(-c(Name,pvalsPar1.adj))
network2 <- as.data.frame(t(network2))


tf.change <- c()
for (i in 1:(nrow(network2)-1)) {
  for (j in 1:ncol(network2)) {
    if (network2[i,j]==1 & network2[1395,j]=="Yes") {
      tf.change[i] <- "Yes"
      break
    }
    else if (network2[i,j]==0 & j==ncol(network2)) {
      tf.change[i] <- "No"
    }
  }
}
tf.change <- as.data.frame(tf.change)
rownames(tf.change) <- rownames(network2[1:(nrow(network2)-1),])

tf.change$Name <- rownames(tf.change)

target.pval <- as.data.frame(matrix(nrow=3133, ncol=2))
target.pval[,1] <- allele.sp$Name
target.pval[,2] <- allele.sp$pvalsPar1.adj

target.pval <- target.pval %>%
  rename(Name=V1, pvalsPar1.adj=V2)

target.pval$DiffYN_par <- ""
for (i in 1:nrow(target.pval)) {
  if (target.pval$pvalsPar1.adj[i] < 0.05) {
    target.pval$DiffYN_par[i] <- "Yes"
  }
  else {
    target.pval$DiffYN_par[i] <- "No"
  }
}

tf.change <- left_join(tf.change, target.pval, by="Name") 
tf.change <- na.omit(tf.change)

tf.change <- tf.change %>%
  rename(tfchange=tf.change)

finaldata <- as.data.frame(matrix(0,2,2)) 
rownames(tf.change) <- c(1:632)
for (i in 1:nrow(tf.change)) {
    if (tf.change[i,1]=="Yes" & tf.change[i,4]=="Yes") {
      finaldata[1,1] <- finaldata[1,1] + 1
    } 
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="Yes") {
      finaldata[1,2] <- finaldata[1,2] + 1
    }
    else if (tf.change[i,1]=="Yes" & tf.change[i,4]=="No") {
      finaldata[2,1] <- finaldata[2,1] + 1
    }
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="No") {
      finaldata[2,2] <- finaldata[2,2] + 1
    }
}  

finaldata$genechange <- c("Yes", "No")
finaldata <- finaldata %>%
  rename(TFYes=V1, TFNo=V2)

ggplot(finaldata, aes(x=genechange, y=TFYes/(TFYes+TFNo))) + geom_bar(stat="identity") +
  xlab("Genes Differentially Expressed?") + ylab("Proportion of Genes with >=1 D.E. Regulators")

facet.net <- read.table(file = './SuppNetwork1.tsv', sep = '\t', header=TRUE)
facet.net <- abs(facet.net)
facet.net$sums <- rowSums(facet.net)
facet.net$Name <- rownames(facet.net)
facet.net <- facet.net %>%
  select(Name, sums)

tf.change <- left_join(tf.change, facet.net, by="Name")
tf.change <- tf.change %>%
  na.omit() %>%
  unique()

facet.final <- as.data.frame(matrix(0,8,4))
for (i in 1:nrow(tf.change)) {
  if (tf.change$sums[i]==1) {
    if (tf.change[i,1]=="Yes" & tf.change[i,4]=="Yes") {
      facet.final[1,1] <- facet.final[1,1] + 1
    } 
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="Yes") {
      facet.final[1,2] <- facet.final[1,2] + 1
    }
    else if (tf.change[i,1]=="Yes" & tf.change[i,4]=="No") {
      facet.final[2,1] <- facet.final[2,1] + 1
    }
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="No") {
      facet.final[2,2] <- facet.final[2,2] + 1
    }
  }
  else if (tf.change$sums[i]==2) {
    if (tf.change[i,1]=="Yes" & tf.change[i,4]=="Yes") {
      facet.final[3,1] <- facet.final[3,1] + 1
    } 
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="Yes") {
      facet.final[3,2] <- facet.final[3,2] + 1
    }
    else if (tf.change[i,1]=="Yes" & tf.change[i,4]=="No") {
      facet.final[4,1] <- facet.final[4,1] + 1
    }
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="No") {
      facet.final[4,2] <- facet.final[4,2] + 1
    }
  }
  else if (tf.change$sums[i]==3) {
    if (tf.change[i,1]=="Yes" & tf.change[i,4]=="Yes") {
      facet.final[5,1] <- facet.final[5,1] + 1
    } 
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="Yes") {
      facet.final[5,2] <- facet.final[5,2] + 1
    }
    else if (tf.change[i,1]=="Yes" & tf.change[i,4]=="No") {
      facet.final[6,1] <- facet.final[6,1] + 1
    }
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="No") {
      facet.final[6,2] <- facet.final[6,2] + 1
    }
  }
  else if (tf.change$sums[i]==4) {
    if (tf.change[i,1]=="Yes" & tf.change[i,4]=="Yes") {
      facet.final[7,1] <- facet.final[7,1] + 1
    } 
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="Yes") {
      facet.final[7,2] <- facet.final[7,2] + 1
    }
    else if (tf.change[i,1]=="Yes" & tf.change[i,4]=="No") {
      facet.final[8,1] <- facet.final[8,1] + 1
    }
    else if (tf.change[i,1]=="No" & tf.change[i,4]=="No") {
      facet.final[8,2] <- facet.final[8,2] + 1
    }
  }
}

facet.final <- facet.final %>%
  rename(TFYes=V1, TFNo=V2, genechange=V3, Sums=V4)

facet.final$genechange <- c("Yes","No","Yes","No","Yes","No","Yes","No")
facet.final$Sums <- c(1,1,2,2,3,3,4,4)

ggplot(facet.final, aes(x=genechange, y=TFYes/(TFYes+TFNo))) + geom_bar(stat="identity") +
  xlab("Genes Differentially Expressed?") + ylab("Proportion of Genes with >=1 D.E. Regulators") +
  facet_grid(cols=vars(Sums))

tf.change %>%
  group_by(sums) %>%
  tally()
