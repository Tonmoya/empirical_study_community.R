# Extraction of differentially expressed genes for each dataset

library(edgeR)
library(e1071)
library(dplyr)
library(igraph) 
library(DESeq2)
library(apeglm)
library(ggplot2)
library(tibble)
library(stringr)
#library(fossil) # for RI, ARI
#library(mclust) # ARI
#library(aricode)
library(bayesbio)
library(NMI)

setwd("DE_genes/")    # set a working directory

#RNA-Seq data
data <- read.csv(file = "DATASETS/GSE138082_SZ_csv.csv", stringsAsFactors = FALSE)    #using dataset GSE138082
dim(data)  # 20255 genes x  78 samples, 1 genename column

outfolder = "SZ_1/"

#pre-processing the data--------
geneNames=data[,1]
rawData<-data[-1] # 20255 x 78

rownames(rawData) <- data[,1]

## NORMALIZATION USING CPM -----------------

normalized_counts <- cpm(rawData)  # normalizaed matrix
thresh1 <- normalized_counts > 0.5
keep1 <- rowSums(thresh1) >= 2
thresData <- normalized_counts[keep1,]
dim(thresData)


##-------------- DESEQ2 DE analysis --------------

colData <- data.frame(condition=ifelse(grepl("CTL", colnames(rawData)), "control", "sz")) #metadata
rownames(colData) <- colnames(thresData)

dds1 <- DESeqDataSetFromMatrix(countData = rawData, colData, formula(~ condition)) # generate object

de_genes <- DESeq(dds1)
res <- results(de_genes)
summary(res)
#res

res <- res[order(res$padj),] # sort summary list acc to log2FC ----------------------
res_logfc <- res[order(res$log2FoldChange),]   # sort summary list acc to log2FC ----------- ***************

sum(res$padj < 0.05, na.rm=TRUE) #98 - number of genes with padj < 0.05
sum(res$padj < 0.1, na.rm=TRUE) # 206

res_matrix <- as.data.frame(res) %>% rownames_to_column("GeneNames")
res_logfc_matrix <- as.data.frame(res_logfc) %>% rownames_to_column("GeneNames")
#res_matrix
plotMA(res)  # DE analysis plot


# moderated log2 fold changes - shrink LFC for visualization and ranking
resultsNames(de_genes)
plotCounts(dds1, gene=which.min(res$padj), intgroup="condition")

resLFC <- lfcShrink(de_genes, coef="condition_sz_vs_control", type="apeglm")
resLFC
plotMA(resLFC)
resLFC_sort <- resLFC[order(resLFC$log2FoldChange),]  # sort acc to log2FC
resLFC[order(resLFC$padj),]  # sort acc to padj


## ----------------- EDGER DE analysis -------------------------------

# use colData for sample names - DataGroups
datagroups <- t(colData)
dge_edger <- DGEList(counts=rawData,group=factor(datagroups))
dim(dge_edger)  # 20255    78
dge_edger_full <- dge_edger  #keep the full list separate

apply(dge_edger$counts, 2, sum) ## total gene counts per sample
keep <- rowSums(cpm(dge_edger)>20) >= 2
dge_edger <- dge_edger[keep,]
dim(dge_edger)  # 9907   78
dge_edger$samples$lib.size <- colSums(dge_edger$counts)  # reset library sizes
dge_edger$samples
dge_edger <- calcNormFactors(dge_edger) # normalize


plotMDS(dge_edger, method="bcv", col=as.numeric(dge_edger$samples$group))
legend("bottomleft", as.character(unique(dge_edger$samples$group)), col=1:3, pch=20) #plotting the samples

#estimate a generalized linear model (GLM) fit using edgeR
design.mat <- model.matrix(~ 0 + dge_edger$samples$group)
colnames(design.mat) <- levels(dge_edger$samples$group)
d2 <- estimateGLMCommonDisp(dge_edger,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

#GLM testing for differential expression
fit <- glmFit(d2, design.mat)
# compare (group 1 - group 2) to 0:
lrt12 <- glmLRT(fit, contrast=c(1,-1))

de_edger_res <- topTags(lrt12, n=10, adjust.method = "BH", sort.by = 'PValue', p.value = 0.1) # most DE genes
#topTags(lrt12, n=10)

de_edger_logfc <- as.data.frame(lrt12) %>% rownames_to_column("GeneNames")    # sort logFC ********************
de_edger_logfc <- de_edger_logfc[order(de_edger_logfc$logFC),]

#plot the log-fold changes of all the genes, and the highlight those that are differentially expressed.
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.1) #significant DE genes (0.05 or 0.1)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-1, 1), col = "blue")

de2_matrix <- as.data.frame(de2) %>% filter(de2 != 0) 




## ---------------------- LIMMA - VOOM DE analysis ----------------------------------------

d0 <- DGEList(rawData)
d0 <- calcNormFactors(d0) # normalization factors

#Filter low-expressed genes
cutoff <- 2.5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d0_genes <- d0[-drop,] 
dim(d0_genes)        
snames <- colnames(rawData) # Sample names
#experiment information from sample names
sample_data <- t(colData)
group <- interaction(sample_data)
#plotMDS(d0_genes, col = as.numeric(group))

# VOOM transformation and calculation of variance weights
mm <- model.matrix(~0 + group)         # specify a model to be fitted -  each coefficient corresponds to a group mean
y_voom <- voom(d0_genes, mm, plot = T)
#tmp <- voom(d0, mm, plot = T)

#Fitting linear models in limma -  using weighted least squares for each gene:
fit <- lmFit(y_voom, mm)
head(coef(fit))

#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models
contr <- makeContrasts(groupcontrol - groupsz, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)  # estimate contrast for each gene
tmp <- eBayes(tmp)  #Empirical Bayes smoothing of standard errors

top.table <- topTable(tmp, sort.by = "P", n = Inf)  # most DE genes ###############
head(top.table, 10)

length(which(top.table$adj.P.Val < 0.05)) 
length(which(top.table$adj.P.Val < 0.1)) # 6

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
#write.table(top.table, file = paste0(outfolder,"gse138082_sz_DE_limma.txt"), row.names = F, sep = "\t", quote = F)

limma_voom_logfc <- top.table[order(top.table$logFC),]  ##*****



##------------------ DE GENES FROM 3 METHODS --------------------------------------

for (i in 1:length(res_matrix$padj)) {
  if(!is.na(res_matrix$padj[i])){
    if(res_matrix$padj[i] < 0.1){
      de_deseq = res_matrix[i,1]
      write.table(de_deseq, paste0(outfolder,"diff_genes_gse138082_SZ_deseq.csv"), 
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}

for(ii in 0:10){
  de_deseq_logfc = res_logfc_matrix[ii,1]  # takr from the begining
  write.table(de_deseq_logfc, paste0(outfolder,"diff_genes_gse138082_SZ_deseq_logfc.csv"), 
              row.names = FALSE, col.names = FALSE, append = TRUE)
  de_deseq_logfc = res_logfc_matrix[(length(res_logfc_matrix$log2FoldChange)-ii),1]  # take from the end
  write.table(de_deseq_logfc, paste0(outfolder,"diff_genes_gse138082_SZ_deseq_logfc.csv"), 
              row.names = FALSE, col.names = FALSE, append = TRUE)
}


de_edger_genes <- as.matrix(de2tags12)
write.table(de_edger_genes, paste0(outfolder,"diff_genes_gse138082_SZ_edger.csv"), row.names = FALSE, col.names = FALSE, 
            append = TRUE)

for(ii in 0:10){
  de_edger_genes = de_edger_logfc[ii,1]  # takr from the begining
  write.table(de_edger_genes, paste0(outfolder,"diff_genes_gse138082_SZ_edger_logfc.csv"), 
              row.names = FALSE, col.names = FALSE, append = TRUE)
  de_edger_genes = de_edger_logfc[(length(de_edger_logfc$logFC)-ii),1]  # take from the end
  write.table(de_edger_genes, paste0(outfolder,"diff_genes_gse138082_SZ_edger_logfc.csv"), 
              row.names = FALSE, col.names = FALSE, append = TRUE)
}


for (i in 1:length(top.table$adj.P.Val)) {
  if(!is.na(top.table$adj.P.Val[i])){
    if(top.table$adj.P.Val[i] < 0.1){
      de_limma = top.table[i,1]
      write.table(de_limma, paste0(outfolder,"diff_genes_gse138082_SZ_limma.csv"), 
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}

for(ii in 0:10){
  de_limma_logfc = limma_voom_logfc[ii,1]  # takr from the begining
  write.table(de_limma_logfc, paste0(outfolder,"diff_genes_gse138082_SZ_limma_logfc.csv"), 
              row.names = FALSE, col.names = FALSE, append = TRUE)
  de_limma_logfc = limma_voom_logfc[(length(limma_voom_logfc$logFC)-ii),1]  # take from the end
  write.table(de_limma_logfc, paste0(outfolder,"diff_genes_gse138082_SZ_limma_logfc.csv"), 
              row.names = FALSE, col.names = FALSE, append = TRUE)
}



##------------------ TAKE DE GENES OF 3 METHODS --------------------------------


de_deseq_list <- t(read.table(paste0(outfolder,"diff_genes_gse138082_SZ_deseq.csv")))
de_deseq_list <- append(t(read.table(paste0(outfolder,"diff_genes_gse138082_SZ_deseq_logfc.csv"))), de_deseq_list)
de_edger_list <- t(read.table(paste0(outfolder,"diff_genes_gse138082_SZ_edger.csv")))
de_edger_list <- append(t(read.table(paste0(outfolder,"diff_genes_gse138082_SZ_edger_logfc.csv"))), de_edger_list)
de_limma_list <- t(read.table(paste0(outfolder,"diff_genes_gse138082_SZ_limma.csv")))
de_limma_list <- append(t(read.table(paste0(outfolder,"diff_genes_gse138082_SZ_limma_logfc.csv"))), de_limma_list)

#int_del <- Reduce(intersect,list(de_deseq_list,de_edger_list,de_limma_list))  

all_del <- append(de_deseq_list,de_edger_list)  
all_del <- append(all_del,de_limma_list)
all_del <- unique(all_del)   # union of the lists




###----------------- get normalized matrix with DE genes --------------------------

norm_diff <- data.frame(matrix(ncol = 0, nrow = dim(colData)[1]))

for (j in 1:length(all_del)) {
  gene_j <- all_del[j]
  de <- normalized_counts[gene_j,]   #get normalized values of samples for the gene
  de <- as.data.frame(de)
  norm_diff[gene_j]<- de
}

#write.table(norm_diff, paste0(outfolder,"diff_genes_gse138082_SZ_normalized.csv"))


##--------------------------- Correlation matrix of normalized DE genes ----------------

# GRAPH FOR DISEASE - 39 samples
norm_D <- norm_diff[c(14,15,16,17,18,19,20,21,22,23,24,25,26,40,41,42,43,44,45,46,47,48,49,50,51,52,66,67,68,69,70,71,72,73,74,75,76,77,78),]  # disease samples

write.csv(norm_D, file = paste0(outfolder,"gse138082_norm_D.csv"), col.names = TRUE)

corr_D <- cor(norm_D, method = c("pearson"))
corr_D1 <- abs(corr_D)  # ## USE ABSOLUTE VALUES
corr_D1[is.na(corr_D1)] <- 0
corr_D1[lower.tri(corr_D1)] <- 0   #lower triangle = 0
corr_graph_D <- graph_from_adjacency_matrix(corr_D1, weighted = TRUE, mode = "undirected")
corr_graph_D <- simplify(corr_graph_D)
plot(corr_graph_D)
comp_corr_graph_D <- components(corr_graph_D)  # membership, csize(cluster sizes), no(no. of clusters)
mod1_D <- modularity(corr_graph_D, membership(comp_corr_graph_D))


# GRAPH FOR CONTROL - 39 samples
norm_C <- norm_diff[c(1,2,3,4,5,6,7,8,9,10,11,12,13,27,28,29,30,31,32,33,34,35,36,37,38,39,53,54,55,56,57,58,59,60,61,62,63,64,65),]   # control samples

write.csv(norm_C, file = paste0(outfolder,"gse138082_norm_C.csv"), col.names = TRUE)

corr_C <- cor(norm_C, method = c("pearson"))
corr_C1 <- abs(corr_C)
corr_C1[is.na(corr_C1)] <- 0
corr_C1[lower.tri(corr_C1)] <- 0 
corr_graph_C <- graph_from_adjacency_matrix(corr_C1, weighted = TRUE, mode = "undirected")
corr_graph_C <- simplify(corr_graph_C)
plot(corr_graph_C)
comp_corr_graph_C <- components(corr_graph_C)  # membership, csize(cluster sizes), no(no. of clusters)
mod1_C <- modularity(corr_graph_C, membership(comp_corr_graph_C))




#nodes - 225 ...... edges - 25200

write.table(corr_D, paste0(outfolder,"corr_D1_gse138082.csv"), row.names = TRUE, col.names = TRUE, append = TRUE)

write.table(corr_C, paste0(outfolder,"corr_C1_gse138082.csv"), row.names = TRUE, col.names = TRUE, append = TRUE)






# COMMUNITY DETECTION USING DIFFERENT METHODS -------------------------------------------------------------------

#--------------------------- COMMUNITIES IN DISEASE ----------------------------------

# LOUVAIN
clv <- cluster_louvain(corr_graph_D)
communities(clv)
sizes(clv)
modularity(clv)


# FAST GREEDY
cfg <- cluster_fast_greedy(corr_graph_D)
communities(cfg)
sizes(cfg)
modularity(cfg)

# SPINGLASS
csg <- cluster_spinglass(corr_graph_D, spins = 4)  
communities(csg)
sizes(csg)
modularity(csg)

# LEADING EIGEN
lec <- cluster_leading_eigen(corr_graph_D)

lec1 <- cluster_leading_eigen(corr_graph_D, start=membership(lec))
communities(lec1)
sizes(lec1)
modularity(lec1)

# INFOMAP
imc <- cluster_infomap(corr_graph_D)
communities(imc)
sizes(imc)
modularity(imc)

# #WALKTRAP
cwt <- cluster_walktrap(corr_graph_D, weights = E(corr_graph_D)$weight, steps = 10)
communities(cwt)
sizes(cwt)
modularity(cwt)


# LABEL PROPAGATION
clp <- cluster_label_prop(corr_graph_D, weights = E(corr_graph_D)$weight)
communities(clp)
sizes(clp)
modularity(clp)

# EDGE_BETWEENNESS - GIRVAN NEWMAN
ceb <- cluster_edge_betweenness(corr_graph_D)
communities(ceb)
sizes(ceb)
modularity(ceb)

## commands for community - modularity, sizes, communities, membership


## ---------------------------- COMMUNITIES IN CONTROL ---------------------------------------------

# LOUVAIN
clv1 <- cluster_louvain(corr_graph_C)
communities(clv1)
sizes(clv1)
modularity(clv1)


# FAST GREEDY
cfg1 <- cluster_fast_greedy(corr_graph_C)
communities(cfg1)
sizes(cfg1)
modularity(cfg1)

# SPINGLASS
csg1 <- cluster_spinglass(corr_graph_C, spins = 4)  
communities(csg1)
sizes(csg1)
modularity(csg1)

# LEADING EIGEN
lec11 <- cluster_leading_eigen(corr_graph_C)

lec12 <- cluster_leading_eigen(corr_graph_C, start=membership(lec))
communities(lec12)
sizes(lec12)
modularity(lec12)

# INFOMAP
imc1 <- cluster_infomap(corr_graph_C)
communities(imc1)
sizes(imc1)
modularity(imc1)

# #WALKTRAP
cwt1 <- cluster_walktrap(corr_graph_C, weights = E(corr_graph_C)$weight, steps = 10)
communities(cwt1)
sizes(cwt1)
modularity(cwt1)


# LABEL PROPAGATION
clp1 <- cluster_label_prop(corr_graph_C, weights = E(corr_graph_C)$weight)
communities(clp1)
sizes(clp1)
modularity(clp1)

# EDGE_BETWEENNESS - GIRVAN NEWMAN
ceb1 <- cluster_edge_betweenness(corr_graph_C)
communities(ceb1)
sizes(ceb1)
modularity(ceb1)


## -------------------------------------------------**---------------------------------------------


# COMPARISON OF COMUUNITIES - JACCARD INDEX, RAND INDEX, ADJUSTED RAND INDEX

# JACCARD INDEX - DISEASE --------------
# 1- lovain, fast greedy
for (u in 1:length(clv)) {
  g1 <- unlist(clv[u])
  for (v in 1:length(cfg)) {
    g2 <- unlist(cfg[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
    
  }
}



# 2 - lovain, spinglass
for (u in 1:length(clv)) {
  g1 <- unlist(clv[u])
  for (v in 1:length(csg)) {
    g2 <- unlist(csg[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
  
}

# 3 - lovain, leading eigen
for (u in 1:length(clv)) {
  g1 <- unlist(clv[u])
  for (v in 1:length(lec1)) {
    g2 <- unlist(lec1[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
  
}

# 4 - lovain, Walktrap
for (u in 1:length(clv)) {
  g1 <- unlist(clv[u])
  for (v in 1:length(cwt)) {
    g2 <- unlist(cwt[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
  
}



# Normalized mutual information - NMI - disease

#as.matrix(membership(clv))

for (y in 1:5) {
  member1 <- read.csv(file = paste0("E:/Tonmoya/R_scripts/community/empirical/SZ_1/GSE138082_communities/membership_",y,".csv"), stringsAsFactors = FALSE)
  for (z in y:5) {
    member2 <- read.csv(file = paste0("E:/Tonmoya/R_scripts/community/empirical/SZ_1/GSE138082_communities/membership_",z,".csv"), stringsAsFactors = FALSE)
    print(paste0("method_",y,"----method_",z,"---- NMI = ",NMI(member1,member2))) 
  }
}






# JACCARD INDEX - CONTROL ---------------
# 1- lovain, fast greedy
for (u in 1:length(clv1)) {
  g1 <- unlist(clv1[u])
  for (v in 1:length(cfg1)) {
    g2 <- unlist(cfg1[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
}



# 2 - lovain, spinglass
for (u in 1:length(clv1)) {
  g1 <- unlist(clv1[u])
  for (v in 1:length(csg1)) {
    g2 <- unlist(csg1[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
  
}

# 3 - lovain, leading eigen
for (u in 1:length(clv1)) {
  g1 <- unlist(clv1[u])
  for (v in 1:length(lec12)) {
    g2 <- unlist(lec12[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
  
}

# 4 - lovain, Walktrap
for (u in 1:length(clv1)) {
  g1 <- unlist(clv1[u])
  for (v in 1:length(cwt1)) {
    g2 <- unlist(cwt1[v])
    print(paste0("g1(",u,")-g2(",v,")-jaccard = ",jaccardSets(g1,g2)))
  }
  
}




# Normalized mutual information - NMI - control

for (y in 1:5) {
  member1 <- read.csv(file = paste0("E:/Tonmoya/R_scripts/community/empirical/SZ_1/GSE138082_communities/membership_",y,y,".csv"), stringsAsFactors = FALSE)
  for (z in y:5) {
    member2 <- read.csv(file = paste0("E:/Tonmoya/R_scripts/community/empirical/SZ_1/GSE138082_communities/membership_",z,z,".csv"), stringsAsFactors = FALSE)
    print(paste0("method_",y,"----method_",z,"---- NMI = ",NMI(member1,member2))) 
  }
}



