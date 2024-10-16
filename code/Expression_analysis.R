##### ENVIRONMENT
# General
library(tidyverse)
# Differential expression
library(DESeq2)
# GO enrichment
library(topGO)
library(GOSim)
# Network analysis
library(WGCNA) 
library(Hmisc) 
library(phyloseq) 
# annotate expression modules
library(GO.db)
# Plotting
library(UpSetR) 
library(gplots) 
options(stringsAsFactors = FALSE)

##### 1. DATA PREP
# load sample table
sampleTable <- read.csv(file = "data/sample_data.csv", row.names = 1, header = T)
sampleTable <- sampleTable[order(rownames(sampleTable)),]
sampleTable$tr_wk <- paste0(sampleTable$Treatment, sampleTable$WeekOfTissueCollection)
# load gene info
em <- read.csv("data/emapper.emapper.annotations", comment.char = "#", sep = "\t")
colnames(em) <- c("query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
tx2gene <- read.table("pipeline_output/star_salmon/tx2gene.tsv")
em$gene <- tx2gene[match(em$query, tx2gene[,1]),"V3"]
# load data and round salmon read counts
dat <- read.table("pipeline_output/star_salmon/salmon.merged.gene_counts_length_scaled.tsv", 
                  header = T, 
                  row.names = 1)[,2:91] %>%
  as.matrix() %>%
  round() 
# datasets
wks <- paste0("Week", c(1,3))
trs <- paste0("P", c(1:4))
dds <- list()

##### 2. DIFFERENTIAL GENE EXPRESSION
# total dataset (for network analysis)
dds$all <- DESeqDataSetFromMatrix(dat, colData = sampleTable, design = ~ WeekOfTissueCollection + Treatment)
# week specific datasets (for DE analysis)
for(wk in wks){
  st <- sampleTable[which(sampleTable$WeekOfTissueCollection == wk),]
  ds <- dat[,rownames(st)]
  dds[[wk]] <- DESeqDataSetFromMatrix(ds, colData = st, design = ~ tr_wk)
}
# filter out rarely expressed genes
for(i in names(dds)){
  my.dds <- dds[[i]]
  smallestGroupSize <- min(table(colData(my.dds)$tr_wk)) ; print(smallestGroupSize)
  keep <- rowSums(counts(my.dds) >= 10) >= smallestGroupSize
  my.dds <- my.dds[keep,]
  dds[[i]] <- my.dds
}
# run deseq
dds <- lapply(dds, DESeq)
save(dds, file = "results/DESeq2/dds.RData")
load("results/DESeq2/dds.RData")
# get results
res <- list()
for(tr in trs){
  for(wk in wks){
    tw <- paste0(tr, wk)
    cw <- paste0("Control", wk)
    # get results
    my.res <- results(dds[[wk]], contrast = c("tr_wk", tw, cw), alpha = 0.05)
    # +ve log2foldchange means higher expression in treatment relative to control
    # -ve log2foldchange means lower expression in treatment relative to control
    # add gene info
    my.res$Gene <- rownames(my.res)
    my.res$Preferred_name <- em[match(rownames(my.res), em$gene), "Preferred_name"]
    my.res$Description <- em[match(rownames(my.res), em$gene), "Description"]
    my.res$GOs <- em[match(rownames(my.res), em$gene), "GOs"]
    # reorder by adjusted P value
    my.res <- my.res[order(my.res$padj),]
    # add to output list
    res[[paste(tr, wk, sep = ".")]] <- my.res
  }
}
# save results
save(res, file = "results/DESeq2/results.RData")
load("results/DESeq2/results.RData")

##### 3. GO TERM ENRICHMENT
GO_en <- list()
GO_en[["dat"]] <- list()
GO_en[["tests"]] <- list()
GO_en[["results"]] <- list()
ont="BP"
alg="elim"
GOanc <- GOSim::getAncestors()
# run for each set of DEGs
for(i in names(res)){
  my.res <- res[[i]]
  # Gene-GO term mapping
  geneID2GO <- as.list(my.res$GOs)
  names(geneID2GO) <- my.res$Gene
  geneID2GO <- lapply(geneID2GO, function(x) str_split(x, pattern = ",", simplify = F) %>% unlist())
  # DE
  geneList <- my.res$padj
  names(geneList) <- my.res$Gene
  # GO datasets
  GOdata <- new("topGOdata",
                ontology = ont,
                geneSelectionFun = function(x) x < 0.05,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO,
                nodeSize = 10)
  GO_en[["dat"]][[i]] <- GOdata
  # run tests
  GOtests <- runTest(GOdata, algorithm = alg, statistic = "fisher")
  # get results
  goRes <- GenTable(GOdata, 
                    P = GOtests,
                    topNodes = length(which(GOtests@score < 0.05)),
                    numChar = 1000)
  # add is immune
  goRes$Immune.related <- unlist(lapply(as.list(goRes$GO.ID), function(x) "GO:0002376" %in% GOanc[[x]]))
  # add to output
  GO_en[["tests"]][[i]] <- GOtests
  GO_en[["results"]][[i]] <- goRes
}
# Remove unnecessary cols and combine all GO enrichment results
all_GO <- GO_en[["results"]]
for(i in names(all_GO)) colnames(all_GO[[i]])[6] <- paste(i, colnames(all_GO[[i]])[6], sep = ".")
for(i in names(all_GO)) all_GO[[i]] <- all_GO[[i]][,c(1,2,6,7)]
all_GO <- lapply(all_GO, as.data.frame)
all_GO <- Reduce(function(x,y) merge(x, y, all = T),  all_GO)
all_GO <- all_GO %>% replace(is.na(.), "n.s")
# put immune related terms first
all_GO <- rbind(all_GO[all_GO$Immune.related,],all_GO[!(all_GO$Immune.related),])
# save
dir.create("results/GO")
write.csv(all_GO, "results/GO/GOenrichment.csv", row.names = F)

##### 4. CREATE GENE EXPRESSION MODULES
dir.create("results/networks")
# variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds$all, blind = T)
# get counts
vsd_counts <- assay(vsd) %>% t()
# pick threshold
sft <- pickSoftThreshold(vsd_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed")
# save
save(sft, file = "results/networks/sft.RData")
load("results/networks/sft.RData")
# get power estimates and plot 
pow <- sft$powerEstimate
plot(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, xlab = "Soft threshold (power)", ylab = "R2")
abline(v = pow, col = "deepskyblue1")
abline(v = pow+1, col = "deepskyblue3")
abline(v = pow+2, col = "deepskyblue4")
# make bwnet
bwnet <- blockwiseModules(vsd_counts,
                          maxBlockSize = 15000, 
                          TOMType = "signed", 
                          power = pow, 
                          numericLabels = TRUE, 
                          randomSeed = 1234)
save(bwnet, file = "results/networks/bwnet.RData")
load("results/networks/bwnet.RData")

##### 5. ANNOTATE GENE EXPRESSION MODULES
# get top level BP GOs
BPendnode <- c("GO:0040011","GO:0043473","GO:0044848","GO:0050789","GO:0002376","GO:0032501","GO:0042592","GO:0051179","GO:0032502","GO:0098754","GO:0048519","GO:0044419","GO:0048518","GO:0022414","GO:0051703","GO:0065007","GO:0048511","GO:0040007","GO:0009987","GO:0050896","GO:0016032")
# get table with GO term names
topBPGO <- select(GO.db, columns=c("GOID","TERM"), keys=BPendnode, keytype="GOID")
topBPGO <- topBPGO[order(topBPGO$TERM),]
# get gene expression module names
modnums <- sort(unique(bwnet$colors))
modnames <- paste0("ME", modnums)
# gene2GO
gene2GO <- as.list(em$GOs)
names(gene2GO) <- em$gene
gene2GO <- lapply(gene2GO, function(x) str_split(x, pattern = ",", simplify = F) %>% unlist())
# reduce to top-level terms
gene2GOslim <- lapply(gene2GO, function(x) x[x %in% BPendnode])
# table
topBPGO <- cbind(topBPGO, data.frame(matrix(NA, ncol = length(modnames), nrow = length(BPendnode), dimnames = list(NULL, modnames))))
# count for all
for(i in modnums){
  my.genes <- names(bwnet$colors[which(bwnet$colors == i)])
  my.allGOs <- unlist(gene2GOslim[my.genes], use.names = F)
  my.table <- table(my.allGOs)[topBPGO$GOID]
  my.table[which(is.na(my.table))] <- 0
  nam <-  paste0("ME", i)
  topBPGO[,nam] <- my.table
}
# remove empty rows
topBPGO <- topBPGO[rowSums(topBPGO[,3:ncol(topBPGO)]) > 0,]
# convert to proportions
topBPGO_prop <- cbind(topBPGO[,1:2], proportions(as.matrix(topBPGO[,3:ncol(topBPGO)]), margin = 2))
# plot
hm_dat <- as.matrix(topBPGO_prop[,3:ncol(topBPGO_prop)])
rownames(hm_dat) <- topBPGO_prop$TERM
heatmap.2(hm_dat,
          col =  hcl.colors(256,palette = "heat2"), 
          density.info= "none", 
          trace="none", 
          key = T,
          key.title = "Proportion of genes in module",
          key.xlab = "",
          dendrogram = "none", 
          Rowv = FALSE, 
          Colv = FALSE, 
          margins=c(6,30), 
          rowsep = c(8,9),
          colsep = c(7,8,16,17,27,28),
          sepcolor = "white", 
          adjCol = c(NA,0.5),
          sepwidth = c(0.01, 0.01))

##### 6. GENE MODULE - MICROBIOME CORRELATION 
# load MB data
# microbiome data files created using code in https://github.com/kvasir7/Smithsonian_Xenopus_Probiotic_Project
# ASV table
featfile <- "data/metabarcoding_final/vJDM_XENGPRO2022_feature_table_comborunFINAL_rarefy_3_7_24.csv"
featureTab <- read.csv(featfile, header = T, row.names = 1)
featureTab <- otu_table(featureTab, taxa_are_rows = TRUE)
# taxonomy
taxofile <- "data/metabarcoding_final/vJDM_XENGPRO2022_taxonomy_comborun.csv"
taxonomy <- as.matrix(read.csv(taxofile, row.names = 1))
taxonomy <- tax_table(taxonomy)
# metadata
metafile <- "data/metabarcoding_final/XENGPro2022_meta_FINAL_3_23_24_with_anti_qpcr.csv"
metadata <- read.csv(metafile, header = T, row.names = 2)
metadata <- sample_data(metadata)
# phyloseq
ps <- merge_phyloseq(featureTab, taxonomy, metadata) 
# keep only samples in RNAseq dataset
ps_filt <- subset_samples(ps, SampleID %in% rownames(bwnet$MEs))
# remove missing samples from bwnet
ME_filt <- bwnet$MEs[sample_names(ps_filt),]
# merge to genus
ps_glom <- tax_glom(ps_filt, taxrank = "Genus")
# filter taxa in less than 50% samples
ps_glom <- filter_taxa(ps_glom, function (x) {sum(x > 0) > nsamples(ps_glom)/2}, prune=TRUE) 
# get matrices of proportional abundance
gtab <- t(as(otu_table(ps_glom), "matrix"))
colnames(gtab) <- tax_table(ps_glom)[,"Genus"]
gtab <- proportions(gtab, margin = 2)
gtab <- gtab[row.names(ME_filt),]
# combine expression modules and bacterial genera
all_dat <- as.matrix(cbind(ME_filt, gtab))
# get correlations
all_cors <- rcorr(all_dat, type = "spearman")
# adjust P vals
all_cors$Q <- p.adjust(all_cors$P[,], method = "fdr") %>% 
  matrix(ncol = ncol(all_cors$P), dimnames = list(rownames(all_cors$P), colnames(all_cors$P)))
# make links and node tables for input into Gephi
# get all pairs
pairs <- combn(rownames(all_cors$Q), m = 2)
# make table
colnames <- c("Source","Target","cor","p","q","AbCor","DirCor")
links <- data.frame(matrix(NA, ncol = length(colnames), nrow = ncol(pairs), dimnames = list(NULL, colnames)))
links$Source <- pairs[1,]
links$Target <- pairs[2,]
# fill table for all pairs
for(n in 1:ncol(pairs)){
  Sou <- pairs[1,n]
  Tar <- pairs[2,n]
  r <- all_cors$r[Sou,Tar]
  links[n,"cor"] <- r
  links[n,"p"] <- all_cors$P[Sou,Tar]
  links[n,"q"] <- all_cors$Q[Sou,Tar]
  links[n,"AbCor"] <- abs(r)
  if(r > 0) dir <- "pos" else dir <- "neg"
  links[n,"DirCor"] <- dir
}
# filter by q value
links <- filter(links, q < 0.05)
# remove spaces from names
links$Source <- gsub(" ", ".", links$Source)
links$Target <- gsub(" ", ".", links$Target)
# make nodes table
nodenames <- sort(unique(c(links$Source, links$Target)))
nodes <- data.frame(Id = nodenames, Label = nodenames, Type = NA)
nodes[grep("ME", nodes$Id),"Type"] <- "Expression.module"
nodes[grep("ME", nodes$Id, invert = T),"Type"] <- "Bacterial.taxon"
# output files
write.table(links, file = "results/networks/network.links.csv", sep = ";", col.names = T, row.names = F, quote = F)
write.table(nodes, file = "results/networks/network.nodes.csv", sep = ";", col.names = T, row.names = F, quote = F)
# heatmap
all_cors_r <- all_cors$r
all_cors_r[which(all_cors$Q > 0.05)] <- 0
# Modules in rows, bacteria in cols
all_cors_r <- all_cors_r[grep("ME", rownames(all_cors_r)), grep("ME", colnames(all_cors_r), invert = T)]
# order by phylum
tax <- as(tax_table(ps_glom), "matrix")
rownames(tax) <- tax[,"Genus"]
gen2phyla <- tax[colnames(all_cors_r),"Phylum"]
# shorten bacterial names
colnames(all_cors_r) <- sub(" .*", "", colnames(all_cors_r)) %>% 
  sub("-.*", "", .)
names(gen2phyla) <- colnames(all_cors_r)
all_cors_r <- all_cors_r[,names(gen2phyla)[order(gen2phyla)]]
all_cors_r <- t(all_cors_r)
# reorder gene modules
all_cors_r <- all_cors_r[,paste0("ME",0:31)]
# colour palettes
palette.breaks <- seq(-1, 1, 0.01)
my.cols <- colorRampPalette(colors = c("darkred","white","royalblue"))(length(palette.breaks) - 1)
phy.cols <- hcl.colors(length(unique(gen2phyla)), palette = "Pastel") 
phy.cols <- colorRampPalette(colors = c("#B6DBC0","#435147"))(length(unique(gen2phyla)))
names(phy.cols) <- sort(unique(gen2phyla))
# plot
pdf("results/networks/heatmap_Mod_Bac_Corr.pdf")
heatmap.2(all_cors_r, 
          col =  my.cols, 
          breaks = palette.breaks, 
          density.info= "none", 
          trace="none", 
          key = F,
          dendrogram = "none", 
          Rowv = FALSE, 
          Colv = FALSE, 
          margins=c(12,8), 
          keysize = 1, 
          RowSideColors = phy.cols[gen2phyla[rownames(all_cors_r)]], 
          ColSideColors = rep("#FFC3AB", ncol(all_cors_r)),
          rowsep = c(0,1,7,8,12,38,39), 
          colsep = c(0,32), 
          adjCol = c(NA,0.5),
          sepcolor = "black", 
          sepwidth = c(0.01, 0.01),
          labRow=as.expression(lapply(rownames(all_cors_r), function(a) bquote(italic(.(a))))))
dev.off()

##### 7. GENE EXPRESSION TABLE
# combine DESeq results
all_DEG <- res
for(i in names(all_DEG)) {
  colnames(all_DEG[[i]])[1:6] <- paste(i, colnames(all_DEG[[i]])[1:6], sep = ".")
}
for(i in names(all_DEG)){
  all_DEG[[i]] <- all_DEG[[i]][,c(2, 5:10)] 
}
all_DEG <- lapply(all_DEG, as.data.frame)
all_DEG <- Reduce(function(x,y) merge(x, y, all = T), all_DEG)
# add column indicating whether gene is annotated with immune GO terms
is_im <- function(x){
  gos <- unlist(gene2GOslim[x])
  return("GO:0002376" %in% gos)
} 
all_DEG$isImmune <-unlist(lapply(as.list(all_DEG$Gene), FUN = is_im))
# add expression module
all_DEG$Expression_module <- paste0("ME",bwnet$colors[all_DEG$Gene])
# add N tests significant
all_DEG$nsig <-  apply(X = all_DEG[,c(7,10,13,16,19,22,25,28)], MARGIN = 1, FUN = function(x) length(which(x < 0.05)), simplify = T)
# reorder columns
all_DEG <- all_DEG[,c(1:4, (ncol(all_DEG) - c(2,1,0)), 5:(ncol(all_DEG)-3))]
write.csv(all_DEG, "results/DESeq2/DESeq.csv", row.names = F)
# filter significant only
sig_DEG <- all_DEG[all_DEG$nsig > 0,]
sig_DEG <- sig_DEG[order(sig_DEG$nsig, decreasing = T),]
# filter immune only
sig_DEG_im <- sig_DEG[sig_DEG$isImmune == TRUE,]
write.csv(all_DEG, "results/DESeq2/DESeq_sig_immune.csv", row.names = F)

##### 8. OVERLAP BETWEEN DIFFERENTIALLY EXPRESSION GENE SETS
# Upset plot for overlap between DEG sets
USdat_all <- list()
USdat_imm <- list()
for(tr in paste0("P", 1:4)){
  for(wk in paste0("Week", c(1,3))){
    col <- paste(tr, wk, "padj", sep = ".")
    P_all <- sig_DEG[which(sig_DEG[,col] <= 0.05), "Gene"]
    P_imm <- sig_DEG_im[which(sig_DEG_im[,col] <= 0.05), "Gene"]
    nam <- paste0(wk, ": ", tr, " vs Control")
    USdat_all[[nam]] <- P_all
    USdat_imm[[nam]] <- P_imm
  }
}
# colours
trcols <- c(P1 = "#433D85", P2 = "#277F8E", P3 = "#95D840", P4 = "#FDE725")
# metadata
USmd <- as.data.frame(cbind(names(USdat_all),str_split_fixed(names(USdat_all), pattern = ' ', n = 4)))
names(USmd)[c(1,3)] <- c("set", "Treatment")
# plot - all genes
pdf("results/DESeq2/upset_DEG_all.pdf")
upset(fromList(USdat_all), 
      sets = rev(names(USdat_all)), 
      keep.order = T, 
      nintersects = 100,
      set_size.show = T,
      order.by = "freq",  
      set.metadata = list(data = USmd, plots = list(list(type = "matrix_rows", column = "Treatment", colors = trcols, alpha = 0.8))), 
      matrix.color = "black",
      set_size.scale_max = 150)
dev.off()
# plot - immune genes only
pdf("results/DESeq2/upset_DEG_imm.pdf")
upset(fromList(USdat_imm), 
      sets = rev(names(USdat_imm)), 
      keep.order = T, 
      nintersects = 100,
      set_size.show = T,
      order.by = "freq",  
      set.metadata = list(data = USmd, plots = list(list(type = "matrix_rows", column = "Treatment", colors = trcols, alpha = 0.8))), 
      matrix.color = "black",
      set_size.scale_max = 22)
dev.off()
