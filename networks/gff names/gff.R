library(tidyverse)
setwd("path/to/directory")
gff <- read_tsv("Spar.gff.txt")
gff.new <- t(colnames(gff))
gff.new <- as.data.frame(gff.new)
gff.new[2:7309,] <- gff


gff.new <- gff.new %>%
  select(V1,V4,V5,V9)

rownames(gff.new) <- c(1:7309)

gff.new %>%
  group_by(V1) %>%
  tally()

gff.new$V1[gff.new$V1==1] <- "chrI"
gff.new$V1[gff.new$V1==2] <- "chrII"
gff.new$V1[gff.new$V1==3] <- "chrIII"  
gff.new$V1[gff.new$V1==4] <- "chrIV"
gff.new$V1[gff.new$V1==5] <- "chrV"
gff.new$V1[gff.new$V1==6] <- "chrVI"
gff.new$V1[gff.new$V1==7] <- "chrVII"
gff.new$V1[gff.new$V1==8] <- "chrVIII"
gff.new$V1[gff.new$V1==9] <- "chrIX"
gff.new$V1[gff.new$V1==10] <- "chrX"
gff.new$V1[gff.new$V1==11] <- "chrXI"
gff.new$V1[gff.new$V1==12] <- "chrXII"
gff.new$V1[gff.new$V1==13] <- "chrXIII"
gff.new$V1[gff.new$V1==14] <- "chrXIV"
gff.new$V1[gff.new$V1==15] <- "chrXV"
gff.new$V1[gff.new$V1==16] <- "chrXVI"

gff.new %>%
  group_by(V1) %>%
  tally()

write_tsv(gff.new, "./Spar_new.gff")
