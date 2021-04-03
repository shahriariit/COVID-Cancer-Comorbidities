
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GOSemSim", "org.Hs.eg.db"))


# 1. Set your work folder
setwd("D:/result1")

# 2. Load libraries for script two
library(GOSemSim)
library(readtext)
library(stringr)
library(factoextra)
library(dendextend)
library(corrplot)
library(RColorBrewer)

# 3. List the available files 
wkpath <- c("Result")
files <- readtext(wkpath)
split_files <- strsplit(files[,1],"_")
id <- unlist(lapply(split_files, function(z) paste0(z[1],"_",z[2])))

# 4. Create list of genes related to each dataset
term_go <- lapply(files[,2], function(z) str_extract_all(z,"GO:.{7}"))
term_go_sb <- lapply(term_go, function(z) z[[1]][1:30])
names(term_go_sb) <- id


grasp_gene <- function(txt_gene){
  ax_1 <- str_replace_all(txt_gene,"GO:.{7} |\n|genes:", "")
  ax_2 <- strsplit(ax_1, " ")
  ax_3 <- strsplit(ax_2[[1]][2:6],",")
  ax_4 <- unique(unlist(ax_3))
  return(ax_4)
}

term_gene_sb <- lapply(files[,2],grasp_gene)
names(term_gene_sb) <- id


# 5. Select the gene ontology and prepare the annotation
goHs1 <- godata('org.Hs.eg.db', ont="BP")
goHs2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)

# 6. Create semantic similarity matrices (GO terms, genes)
len <- length(id)
mat_go_sem_sim = mat_gene_sem_sim <- matrix(data = 0, nrow = len, ncol = len)
rownames(mat_go_sem_sim) = rownames(mat_gene_sem_sim) <- id
colnames(mat_go_sem_sim) = colnames(mat_gene_sem_sim) <- id

for(x in 1:len){
  for(y in 1:len){
    mat_go_sem_sim[x,y] <- mgoSim(term_go_sub[[x]],term_go_sb[[y]],semData=goHs1, measure="Wang", combine="BMA")
    print(paste("x =",x,"y =",y,"per =", round(((x-1)*len+y)/len^2,digits = 3)*100,"%"))
  }
}

for(x in 1:len){
  for(y in 1:len){
    mat_gene_sem_sim[x,y] <- clusterSim(term_gene_sb[[x]], term_gene_sb[[y]],semData=goHs2, measure="Wang", combine="BMA")
    cat("x =",x,"y =",y,"per =",round(((x-1)*len+y)/len^2,digits = 3)*100,"%\n")
  }
}

load("2nd_code_global_environment.RData")

library(DOSE)

pdf(file = "mat_go.pdf", width = 16, height = 13)
simplot(mat_go_sem_sim, color.low="white", color.high="red", labs=TRUE, digits=2, labs.size=4.5, font.size=16, xlab="", ylab="")
dev.off()

pdf(file = "mat_gene.pdf", width = 16, height = 13)
simplot(mat_gene_sem_sim,color.low="white", color.high="red",labs=TRUE, digits=2, labs.size=4.5, font.size=16, xlab="", ylab="")
dev.off() 

do_id <-  c("DOID:1612","DOID:219","DOID:263","DOID:3571","DOID:1781","DOID:2945")
do_ac <-  c("BC","CC","KD","LC","TC","PMBC")

mat_do_sem_sim <- doSim(do_id, do_id,measure="Wang")
rownames(mat_do_sem_sim) = colnames(mat_do_sem_sim) <- do_ac

pdf(file = "mat_do.pdf", width = 10)
simplot(mat_do_sem_sim, color.low="white", color.high="red", labs=TRUE, digits=2, labs.size=6, font.size=20, xlab="", ylab="")
dev.off()

# 8. Create KEGG Enrichment graph
library(clusterProfiler)

rispcor <- lapply(term_gene_sb,function(z) bitr(z, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
rispcor <- lapply(rispcor,"[","ENTREZID")
rispcor <- lapply(rispcor, unlist)
names(rispcor) = substr(names(rispcor),1,3)

comcluster <- compareCluster(geneCluster = rispcor, fun = "enrichKEGG")
pdf(file = "enrich_KEGG_PMBC2.pdf", width = 18, height = 15)
dotplot(comcluster, font.size = 16)
dev.off()

save.image(".RData")