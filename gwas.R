setwd("/Users/patrickkampmeyer/Dropbox/Bioinformatics")
#install.packages("BiocManager")
#library(BiocManager)
# Mac (terminal): xcode-select --install, to update Rcpp
#install.packages("Rcpp")
#BiocManager::install("snpStats") # Yes -- install from source, update all
library(snpStats)

# too slow to use map/ped
#file_prefix <- "LIBR_1-7_2014-02-03"
#file.ped <- paste(file_prefix,".ped",sep="")
#file.map <- paste(file_prefix,".map",sep="")
#gwas.data <- read.pedfile(file=file.ped, snps=file.map)

# convert to bed: inbix command line
#inbix --file LIBR_1-7_2014-02-03 --make-bed
#data.dir = '/test'
#gwas.fn <- lapply(c(bed='bed',bim='bim',fam='fam',gds='gds'), function(n) sprintf("%s/GWAStutorial.%s", data.dir, n))
# Read in PLINK files
#file_prefix <- "inbix"
#file.bed <- paste(getwd(),"/",file_prefix,".bed",sep="")
#file.bim <- paste(getwd(),"/",file_prefix,".bim",sep="")
#file.fam <- paste(getwd(),"/",file_prefix,".fam",sep="")

gwas.data <- read.plink("inbix.bed", "inbix.bim", "inbix.fam", na.strings = ("-9"))
str(gwas.data)
# contains 3 lists
#1. phenotype, subject ids, not that informative for this dataset
#gwas.data$fam
#phenotype <- ex.data$fam$affected-1  # change pheno from 1/2 to 0/1
#2. SnpMatrix of genotypes (240 subjects by 4,301,332 SNPs )
genotypes <- gwas.data$genotypes  
dim(genotypes)
geno.stats <- col.summary(genotypes)
colnames(geno.stats) # Call.rate, MAF, z.HWE, P.AB...
maf.keep <- geno.stats$MAF>.01  # minor allele frequencies
na.rm <- is.na(geno.stats$MAF)   # will need to remove these na's
sum(as.integer(maf.keep),na.rm=T)  # keep about 1.7 million. 

# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(geno.stats, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotypes)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") 
#1960222 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotypes <- genotypes[,use]
geno.stats <- geno.stats[use,]

print(genotypes) # 2,341,110 SNPs remain
rownames(genotypes)

#3? map of SNP ids
# gwas.data$map$snp.name  # rs numbers and kgp (illumina) snp ids
snp.names <- as.character(colnames(genotypes))

#4. remove samples with NA phenotypes and merge with genotype data
pheno.data <- read.csv("gwas.subjdata.csv")
colnames(pheno.data)

library(dplyr)
pheno.subcols <- pheno.data %>% select(record_id, bioax_kynurenine_1,
                                      bioax_tryptophan_1,
                                      bioax_kynurenicacid_1, 
                                      bioax_quinolinicacid_1)             %>%
  mutate(kyn_trp_ratio = ifelse(bioax_tryptophan_1 != 0,       # make new column and avoid divide by 0
                              bioax_kynurenine_1/bioax_tryptophan_1, NA)) %>%
  mutate(kyna_quina_ratio = ifelse(bioax_quinolinicacid_1 != 0,       # make new column and avoid divide by 0
                                  bioax_kynurenicacid_1/bioax_quinolinicacid_1, NA)) %>% 
  filter(!is.na(kyn_trp_ratio))  # remove rows that have NA in this column

#kyn<-pheno.data$bioax_kynurenine_1
#trp<-pheno.data$bioax_tryptophan_1
#kynacid<-pheno.data$bioax_kynurenicacid_1
#quinacid<-pheno.data$bioax_quinolinicacid_1

common_samples <- intersect(rownames(genotypes),pheno.subcols$record_id)
genotypes <- genotypes[common_samples,]
dim(genotypes)
pheno.subsamples <- pheno.subcols %>% filter(record_id %in% common_samples)
# phenotype vector
KQrat <- pheno.subsamples$kyna_quina_ratio

#5 ld prune with SNPRelate and gds file
ld.thresh <- 0.2    # LD cut-off

# Create gds file, required for SNPRelate functions
#library(BiocManager)
#BiocManager::install("SNPRelate")
library(SNPRelate)
# already ran to create gds file
#snpgdsBED2GDS("inbix.bed", "inbix.fam", "inbix.bim", "inbix.gds")

gdsfile <- openfn.gds("inbix.gds", readonly = FALSE)
gdsfile$root
snp_pos <- read.gdsn(index.gdsn(gdsfile, "snp.position"))
snp_id <- read.gdsn(index.gdsn(gdsfile, "snp.id"))
snp_chr <- read.gdsn(index.gdsn(gdsfile, "snp.chromosome"))
map.df <- data.frame(chr=snp_chr,snp_id=snp_id,cm=0,bp=snp_pos) 
map.df <- map.df %>% mutate_at("snp_id", as.character) 
map.df[10000,]
#map.df <- read.table("LIBR_1-7_2014-02-03.map",header=F,sep="\t")

# Automatically added "-1" sample suffixes are removed
#gds.ids <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
#gds.ids <- sub("-1", "", gds.ids)
#add.gdsn(gdsfile, "sample.id", gds.ids, replace = TRUE)

#LD Prune SNPs for IBD analysis
set.seed(1000)
snpSUB <- snpgdsLDpruning(gdsfile, ld.threshold = ld.thresh,
                          sample.id = rownames(genotypes), # Only analyze the filtered samples
                          snp.id = colnames(genotypes)) # Only analyze the filtered SNPs

# names of SNPs to keep from ld pruning
snpset <- unlist(snpSUB, use.names=FALSE)

ld.usesnps <- colnames(genotypes) %in% snpset 
#length(ld.usesnps)
#sum(ld.usesnps)
genotypes <- genotypes[,ld.usesnps]
dim(genotypes)

#5. gwas qtl with KQ ratio phenotype and ld/maf filtered snps
library(npdr)
genoNum <- as(genotypes, "numeric")

#fit <- summary(lm(KQrat ~ genoNum[,1]))
#coeffs <- fit$coeff
#coeffs[2,]

KQrat.univar <- npdr::uniReg(outcome=KQrat, dataset=genoNum, regression.type="lm")

# check eqtl
colnames(KQrat.univar)
KQrat.univar.addsnpcol <- data.frame(KQrat.univar) %>% 
                    rownames_to_column() %>% rename(snp_id=rowname) %>% 
                    mutate_at("snp_id", as.character)
KQrat.univar.pfilt <- KQrat.univar.addsnpcol %>% 
                      filter(p.adj<.05) # %>% pull(snp_id) 

# Illumina mapping of kgp to rs numbers
rsmap <- read.table("InfiniumOmni5-4v1-2_A1_b144_rsids.txt",header=T,sep="\t")
colnames(rsmap)

# add the RsID column for all SNPs
KQrat.univar.pfilt <- merge(KQrat.univar.pfilt, rsmap, by.x="snp_id", by.y="Name") %>% 
                      mutate_at("RsID", as.character)  # factor to char
colnames(KQrat.univar.pfilt)

KQrat.univar.pfilt$RsID[1]

# map SNPs to genes
# Read in file containing protein coding genes coords
gene_coords <- read.csv("ProCodgene_coords.csv", stringsAsFactors = FALSE)
colnames(gene_coords)
colnames(map.df)

# map snps to genes
#snp_2b_mapped <- KQrat.univar.pfilt$snp_id[1]  # id and rs might be different
#rs_2b_mapped <- KQrat.univar.pfilt$RsID[1]

rs2gene <- function(snp_vector){
  # global variables map.df and gene_cooords
  # map.df uses illumina ids, so need to uses these as snp_vector input (kgp etc.)
  left_genes=rep(NA,length(c))
  right_genes=rep(NA,length(snp_vector))
  gene_chromes=rep(NA,length(snp_vector))

  i<-1
  for (snp_2b_mapped in snp_vector){
    # find the basepair position and chromosome num of snp to be mapped to gene
    # find the row of plink-style map file for the target snp
    # get chr and bp of snp
    map_row <- which(snp_2b_mapped == map.df$snp_id)
    snp_chr <- map.df[map_row,]$chr
    snp_bp <- map.df[map_row,]$bp
    
    # get gene coords for the chr that contains the snp
    chr_coord_rows <- gene_coords[which(snp_chr == gene_coords),]
    starts <- chr_coord_rows$start
    stops <- chr_coord_rows$stop
    genes <- chr_coord_rows$gene
    # check where the snp falls compared to gene coords
    left_idx <- tail(which(snp_bp > starts),1)
    right_idx <- head(which(snp_bp < stops),1)
    gene2left <- genes[left_idx]
    gene2right <- genes[right_idx]
    
    #cat(snp_chr," ",snp_bp,"\nleft: ",gene2left,"\nright: ",gene2right,"\n")
    #cat(str(gene2right),"\n")
    
    if (length(gene2left) != 0){ # sometimes left or right coord not found 
      left_genes[i] <- gene2left
    }
    if (length(gene2right) != 0){ 
      right_genes[i] <- gene2right
    }
    gene_chromes[i] <- snp_chr
    i <- i+1
  }
  return(list(left_genes=left_genes,right_genes=right_genes,gene_chromes=gene_chromes))
}

univ.map <- rs2gene(KQrat.univar.pfilt$snp_id)

KQrat.univar.pfilt.genes <- KQrat.univar.pfilt %>%  # add columns
          mutate(chr=univ.map$gene_chromes,
                 lft_gene=univ.map$left_genes,
                 rt_gene=univ.map$right_genes) %>%
          arrange(pval)  # sort by pval

write.csv(KQrat.univar.pfilt.genes,file="KQ_ratio_results.csv",row.names=F)

# NPDR MDD analysis
start_time <- Sys.time()
KQrat.npdr <- npdr::npdr(outcome=KQrat, dataset=genoNum, regression.type="lm", 
                      attr.diff.type="allele-sharing",
                      nbd.method="multisurf", nbd.metric = "manhattan", 
                      msurf.sd.frac=0.5,
                      dopar.nn = T, dopar.reg=T,
                      padj.method="bonferroni", verbose = T)
end_time <- Sys.time()
end_time - start_time  # about 8 m

colnames(KQrat.npdr)

KQrat.npdr.pfilt <- KQrat.npdr %>% filter(pval.adj<.05) # %>% pull(snp_id) 

# Illumina mapping of kgp to rs numbers
colnames(rsmap)

# add the RsID column for all SNPs
KQrat.npdr.pfilt <- merge(KQrat.npdr.pfilt, rsmap, by.x="att", by.y="Name") %>% 
  mutate_at("RsID", as.character)  # factor to char

# time consuming to get genes
npdr.map <- rs2gene(KQrat.npdr.pfilt$att)  
# use illumina ids (att not rs) to match coords in plink map

KQrat.npdr.pfilt.genes <- KQrat.npdr.pfilt %>%  # add columns
  mutate(chr=npdr.map$gene_chromes,
         lft_gene=npdr.map$left_genes,
         rt_gene=npdr.map$right_genes) %>%
  arrange(pval.att)  # sort by pval

write.csv(KQrat.npdr.pfilt.genes,file="KQ_ratio_npdr_results.csv",row.names=F)


# genotypes is not a dataframe and probably don't want to convert
#rownames_to_column() add names to genotypes for merge
# merge pheno.subset and genotypes and rm na phenotypes
# merge(mdd.cmv.RN,geneExpr.RN,by="rowname")
snp1 <- as.numeric(genotypes[,1])
chisq.test(phenotype,snp1)

# create your own chi-square table
observed.table <- table(phenotype,snp8)
rowmarg <- rowSums(observed.table)    # case and control counts
colmarg <- colSums(observed.table)    # genotype counts (regardless of case/control)
colmarg.probs <- colmarg/sum(colmarg) # marginal probability of each genotype
# geno-prob-vec * case-count, geno-prob-vec * control-count
expected.table <- rbind(colmarg.probs*rowmarg[1],colmarg.probs*rowmarg[2])
# or
matrix(rowmarg,ncol=1)%*%matrix(colmarg.probs,nrow=1)

#Fisher test
ft <- fisher.test(observed.table)
ft$p.value

# logistic regression
# need to be factors
pheno.factor <- factor(phenotype)   
levels(pheno.factor)[levels(pheno.factor)==0] <- 0  # you still need to do this to make 0/1
levels(pheno.factor)[levels(pheno.factor)==1] <- 1
snp8.factor <- factor(snp8)
lr.fit <- glm(pheno.factor~snp8.factor,family=binomial)
summary(lr.fit)

# dominant encoding
snp8.DOM <- factor(snp8==1)
lr.DOM.fit <- glm(pheno.factor~snp8.DOM,family=binomial)
summary(lr.DOM.fit)

library(ggplot2)

#plot(snp8.factor,lr.fit$fitted.values)
#as(genotypes[, "rs7835221"], "numeric")

# plot the additive encoding model
snp8.df <- data.frame(cbind(snp8.factor,phenotype)) # phenotype is numeric
colnames(snp8.df) <- c("geno","pheno")
A1<-ex.data$map$allele.1[8]
A2<-ex.data$map$allele.2[8]
geno.labels <- c(paste(A1,A1,sep=""),paste(A1,A2,sep=""),paste(A2,A2,sep=""))

lr.plot <- ggplot(snp8.df, aes(x=geno, y=pheno)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0.1)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), formula=y~x, se=FALSE) + 
  xlim(geno.labels) +  # add x-tick genotype labels
  geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
  ggtitle(ex.data$map$snp.names[8]) 

print(lr.plot)