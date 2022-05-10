library(tidyverse)
library(ggplot2)
library(sna)
library(igraph)
setwd("path/to/directory")

# read in all the data to make "bigtable" TM
tau <- read_csv("path/to/Tissue_environment_specificity_X_perc_cis/Yeast_data/Gresham/scRNA_tau_env.csv")
dnds <- read.table("./dn_ds_cer_par.txt", sep='\t', header=TRUE)
ppi <- read.table("./yeast_PPI_bipart.txt", sep='\t', header=FALSE)
network <-  read.table("../SuppNetwork1.tsv", sep = '\t', header=TRUE)
exp <- read.table("./cpm_genes.tsv", sep='\t', header=TRUE)

# create big table to fill with all data
bigtable <- as.data.frame(matrix(nrow=198))

# bigtable add name and outdegree: (start with 198)
network <- as.data.frame(t(abs(network)))
network$sums <- rowSums(network)
network$Name <- rownames(network)

bigtable$Name <- network$Name
bigtable$Outdegree <- network$sums
bigtable <- bigtable %>%
  select(!V1)

# bigtable add tau: (198 remain after joining)
tau$sd <- 0
for (i in 1:nrow(tau)) {
  tau$sd[i] <- sd(tau[i,2:12])
}

bigtable <- left_join(bigtable, tau[,c("Name","tau", "sd")], by="Name")
bigtable <- bigtable %>%
  na.omit() %>%
  unique()

# bigtable add exp: (181 remain after joining)
bigtable <- left_join(bigtable, exp[,c("Name","sc267")], by="Name")
bigtable <- left_join(bigtable, exp[,c("Name","sp268")], by="Name")
bigtable$expdiv <- abs(log2(bigtable$sc267/bigtable$sp268))

ggplot(bigtable, aes(x=sc267, y=tau)) + geom_point()
ggplot(bigtable, aes(x=sc267, y=sd)) + geom_point()

# bigtable add ppi: (66 remain after joining)
ppi <- as.data.frame(table(ppi$V1))
ppi <- ppi %>%
  rename(Name=Var1)
bigtable <- left_join(bigtable, ppi, by="Name")

#bigtable add dn/ds: (153 remain after joining)
dnds <- dnds %>%
  rename(Name=ORF.code)
bigtable <- left_join(bigtable, dnds[,c("Name", "dN.dS")], by="Name")

#
bigtable <- bigtable %>%
  na.omit() %>%
  unique()
###########################################################################
# create tau bins by thirds
bigtable$tau.thirds <- 0
for (i in 1:nrow(bigtable)) {
  if (bigtable$tau[i] < 0.33) {
    bigtable$tau.thirds[i] <- 1
  }
  else if (0.33 <= bigtable$tau[i] & bigtable$tau[i] <0.66) {
    bigtable$tau.thirds[i] <- 2
  }
  else {
    bigtable$tau.thirds[i] <- 3
  }
}

bigtable %>%
  group_by(tau.thirds) %>%
  tally()

# order elements of row by ascending order
for (i in 1:nrow(tau)) {
  tau[i,2:12] <- sort(as.data.frame(tau[i,2:12]))
}

tf_tau <- left_join(tau, bigtable[,c("Name", "tau.thirds")], by="Name")
tf_tau <- tf_tau %>%
  na.omit() %>%
  unique()

tf_tau %>%
  pivot_longer(cols=2:12, names_to = "indicators", values_to = "values") %>%
  ggplot(data=.,aes(x=indicators, y=log(values+0.00001))) + geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="red", fill="red")+ facet_grid(rows=vars(tau.thirds))+
  ylim(-2.5,5)

# plots 
## tau vs dnds
ggplot(bigtable, aes(x=tau, y=dN.dS)) + geom_point() + geom_smooth(method=loess)

## tau vs exp
ggplot(bigtable, aes(x=tau, y=expdiv)) + geom_point() + geom_smooth(method=loess)

##  outdegree vs dnds
ggplot(bigtable, aes(x=log(Outdegree+1), y=dN.dS)) + geom_point() + geom_smooth(method=loess)

## outdegree vs exp
ggplot(bigtable, aes(x=log(Outdegree+1), y=expdiv)) + geom_point() + geom_smooth(method=lm)

####################### multivariate regression 
bigtable <- bigtable %>% 
  filter_if(~is.numeric(.), all_vars(!is.infinite(.)))

multi.fit<-lm(cbind(expdiv, dN.dS) ~ Outdegree + tau, data=bigtable)
summary(multi.fit)

############## betweenness

edgelist <- as.character(matrix(data=NA,ncol=2))

for (i in 1:nrow(network)) {
  for (j in 1:(ncol(network)-2)) {
    if (network[i,j]==1) {
      edgelist <- rbind(edgelist, c(rownames(network)[i], colnames(network)[j]))
    }
  }
}

edgelist <- edgelist[2:nrow(edgelist),]
rownames(edgelist) <- c(1:1634)

graph <- graph_from_edgelist(edgelist, directed=TRUE)
btwn <- as.data.frame(betweenness(graph))

colnames(btwn) <- "Betweenness"
btwn$Name <- rownames(btwn)

btwn <- left_join(btwn, network, by="Name") 
btwn <- btwn %>%
  na.omit() %>%
  unique()
btwn <- btwn[,1:2]
rownames(btwn) <- btwn$Name

ggplot(btwn, aes(x=Betweenness)) + geom_bar()
btwn %>%
  group_by(Betweenness) %>%
  tally()

# append btwn to bigtable
bigtable <- left_join(bigtable, btwn, by="Name")
bigtable <- bigtable %>%
  na.omit() %>%
  unique()

#closeness
components <- components(graph)
components$no
