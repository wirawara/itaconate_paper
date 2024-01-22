

####### -------------------------------

##### Author: Andrea Komljenovic
##### Affiliation: Medical University of Vienna, Vienna, Austria

###### --------------------------------


# ShinyGO v0.77
# R v4.2.1
# fgsea v1.25.1
# GOdb v3.15.0
# enrichplot v1.16.2
# clusterProfiler v4.4.4
# edgeR v3.38.4
# limma v3.52.4
# biomaRt 2.56.1


# 1. volcano plots 5mMvsOM 

library(utils)
library(edgeR)
library(fgsea)
library(openxlsx)
library(org.Hs.eg.db)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(fgsea)
library(GO.db)


path <- "~/Documents/kieler/P198_Markus_Kieler_230228/2023-02-27_198_mRNAreport/"
setwd(path)


####### ---------------------------------
counts <- read.xlsx("genesTable.xlsx") 
counts.orig <- counts
### ensembl gene id
rownames(counts) <- counts$gene
counts <- counts[, 4:ncol(counts)]

#### 
dim(counts)
# 30562    36

### select only protein-coding ones - check the version of the biomart
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "http://www.ensembl.org")
annotation <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
               values = rownames(counts), 
               mart = mart)

annotation.protein.coding <- annotation[annotation$gene_biotype == "protein_coding",]  
dim(annotation.protein.coding)
# 16475     3

counts.protein.coding <- counts[na.omit(match(annotation.protein.coding$ensembl_gene_id, rownames(counts))),]
dim(counts.protein.coding)


colnames(counts.protein.coding) <- c(rep("5mM", 4), rep("1mM", 4), rep("OM", 5), rep("GM", 5), rep("5mM", 5), rep("1mM", 5), rep("OM", 4), rep("GM", 4))

### select only OM, GM, 5mM

counts.protein.coding.osteo.om.gm.5m <- counts.protein.coding[,colnames(counts.protein.coding) %in% c("GM", "5mM", "OM")]
colnames(counts.protein.coding.osteo.om.gm.5m) <- gsub("\\.[0-9]+","", colnames(counts.protein.coding.osteo.om.gm.5m))

### order the columns
counts.protein.coding.osteo.om.gm.5m <- counts.protein.coding.osteo.om.gm.5m[,c(1,  2 , 3 , 4 ,15, 16, 17, 18 ,19, 5 , 6 , 7 , 8 , 9 ,20, 21, 22 ,23,10, 11, 12, 13, 14, 24, 25, 26, 27)]
colnames(counts.protein.coding.osteo.om.gm.5m) <- gsub("\\.[0-9]+","", colnames(counts.protein.coding.osteo.om.gm.5m))
# define the treatment
treatment <- colnames(counts.protein.coding.osteo.om.gm.5m)


#### for differential expression
dge.om.gm.5m<- DGEList(counts=counts.protein.coding.osteo.om.gm.5m, group=treatment)
design.om.gm.5m<- model.matrix(~0+treatment)
rownames(design.om.gm.5m) <- colnames(dge.om.gm.5m)


keep.om.gm.5m <- rowSums(cpm(dge.om.gm.5m)>1) >= 5
dge.om.gm.5m <- dge.om.gm.5m[keep.om.gm.5m,]
dge.om.gm.5m$samples$lib.size <- colSums(dge.om.gm.5m$counts)
dge.om.gm.5m <- calcNormFactors(dge.om.gm.5m)


### PCA
## pca
logCPM.om.gm.5m <- cpm(dge.om.gm.5m ,log=TRUE,prior.count=5) 
# 12607    27

### PCA
pca.om.gm.5m <- prcomp(t(logCPM.om.gm.5m), center = TRUE, scale = TRUE)
score.om.gm.5m <- data.frame(pca.om.gm.5m$x)[,1:2]
eigs.om.gm.5m <- pca.om.gm.5m$sdev^2
perc.vars.om.gm.5m <- round(eigs.om.gm.5m[1:2] / sum(eigs.om.gm.5m),digits=2)*100
# colnames of the conditions
score.om.gm.5m$grouping <- treatment


pdf("./PCA_Treatments_OM_GM_5M.pdf", 5, 5)
g1 <- ggplot(score.om.gm.5m, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = grouping, colour = grouping),show.legend = FALSE) + 
  scale_shape_manual(values=c(15,15,15,15,15,15,15,15,15,15)) +
  geom_label_repel(label = score.om.gm.5m$grouping) +
  # cannot calculate too few points
  stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill =grouping)) +
  xlab(paste("PC1:", perc.vars.om.gm.5m[1], "%")) + 
  ylab(paste("PC2:", perc.vars.om.gm.5m[2], "%")) +
  ggtitle("PCA Treatments") + theme_classic()

g1
dev.off()

####### -------------


## limma-voom
DEG <- voom(dge.om.gm.5m, design.om.gm.5m, plot = T) # now it's better
nrow(DEG)
# [1]  12607
fit.voom <- lmFit(DEG,design.om.gm.5m)
# fit.voom <- eBayes(fit.voom)
head(coef(fit.voom))


### contrasts OM vs GM
contrGMvsOM <- makeContrasts(treatmentGM - treatmentOM, levels = colnames(coef(fit.voom)))
contrGMvsOM
tmpGMvsOM <- contrasts.fit(fit.voom, contrGMvsOM)
efitGMvsOM <- eBayes(tmpGMvsOM)
top.tableGMvsOM <- topTable(efitGMvsOM, sort.by = "P", n = Inf)
head(top.tableGMvsOM, 20)
sum(top.tableGMvsOM$adj.P.Val<0.05)
# 5322

singif.top.GMvsOM.005 <- top.tableGMvsOM[top.tableGMvsOM$adj.P.Val < 0.05,]
# 5322


summary(decideTests(efitGMvsOM))
#        treatmentGM - treatmentOM
# Down                        2690
# NotSig                      7285
# Up                          2632

### add gene names to the top table
topGMvsOM <- top.tableGMvsOM
annotation.protein.coding.genes <- annotation.protein.coding[!duplicated(annotation.protein.coding$external_gene_name),]
top.gene.names.gmvsom <- annotation.protein.coding.genes[na.omit(match(rownames(topGMvsOM), annotation.protein.coding.genes$ensembl_gene_id)),]
top.table.gmvsom  <-topGMvsOM[na.omit(match(top.gene.names.gmvsom$ensembl_gene_id,rownames(topGMvsOM))),]
top.table.gmvsom$external_gene_name <- top.gene.names.gmvsom$external_gene_name


openxlsx::write.xlsx(top.table.gmvsom, "./results_diffexp_GMvsOM_005.xlsx", rowNames=TRUE)


### contrasts OM vs 5M
contr5MvsOM <- makeContrasts(treatment5mM - treatmentOM, levels = colnames(coef(fit.voom)))
contr5MvsOM
tmp5MvsOM <- contrasts.fit(fit.voom, contr5MvsOM)
efit5MvsOM <- eBayes(tmp5MvsOM)
top.table5MvsOM <- topTable(efit5MvsOM, sort.by = "P", n = Inf)
head(top.table5MvsOM, 20)
sum(top.table5MvsOM$adj.P.Val<0.1)
# 126 genes - 0.05
# 260

singif.top.5MvsOM <- top.table5MvsOM[top.table5MvsOM$adj.P.Val < 0.1,]
# 126 genes - 0.05
# 260

summary(decideTests(efit5MvsOM))
#        treatmentOM - treatment5mM
# Down                           69
# NotSig                      12481
# Up                             57


#        treatment5mM - treatmentOM
# Down                           57
# NotSig                      12481
# Up                             69



### add gene names to the top table
top5MvsOM <- top.table5MvsOM
top.gene.names.5mvsom <- annotation.protein.coding.genes[na.omit(match(rownames(top5MvsOM), annotation.protein.coding.genes$ensembl_gene_id)),]
top.table.5mvsom  <-top5MvsOM[na.omit(match(top.gene.names.5mvsom$ensembl_gene_id,rownames(top5MvsOM))),]
top.table.5mvsom$external_gene_name <- top.gene.names.5mvsom$external_gene_name


openxlsx::write.xlsx(top.table.5mvsom, "./results_diffexp_5MvsOM_01.xlsx", rowNames=TRUE)



volcano_plot_summary_update <- function(results, title){

    #results are dataframe with log2FC and adj_pval

    library(ggplot2)
    library(gridExtra)
    library(grid)

    # add a column of NAs
    results$diffexpressed <- "NOT SIG"

    # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP"
    results$diffexpressed[results$logFC > 0.5 & results$adj.P.Val < 0.1] <- "UP"

    # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
    results$diffexpressed[results$logFC < -0.5 & results$adj.P.Val < 0.1] <- "DOWN"

    ### to change the colors of the legend
    print(levels(as.factor(results$diffexpressed)))
    results$delabel <- NA
    results$delabel[results$diffexpressed != "NOT SIG"] <- results$external_gene_name[results$diffexpressed != "NOT SIG"]

    if(length(levels(as.factor(results$diffexpressed)))==3){
    g <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val),label = delabel,
        , col=diffexpressed)) +
            geom_point() +  scale_color_manual(values=c("blue", "gray", "red")) +
            geom_hline(aes(yintercept = -log10(0.1)), linetype = "dashed") +
            labs(x = bquote(~Log[2]~ "fold change"), y = bquote(~-Log[10]~italic(FDR))) +
            theme_bw() + geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed") + ggtitle(title) + geom_text_repel() 
            g2 <- grid.arrange(g,
                bottom = textGrob(title,
                x = 0.40, y = 1, gp = gpar(fontsize = 9)))



        } else {
          g <- ggplot(results, aes(x = logFC, y = -log10(adj_pval), col=diffexpressed)) +
                  geom_point() +  scale_color_manual(values=c( "gray")) +
                  geom_hline(aes(yintercept = -log10(0.1)), linetype = "dashed") +
                  labs(x = bquote(~Log[2]~ "fold change"), y = bquote(~-Log[10]~italic(FDR))) +
                  theme_bw() + geom_vline(xintercept=c(-1, 1), linetype = "dashed") +
                  ggtitle(title)

          g2 <- grid.arrange(g,
                      bottom = textGrob(title,
                      x = 0.45, y = 1, gp = gpar(fontsize = 9)))


        }

    return(g2)

}

##### volcano plot
pdf("./volcano_plot_5mMvsOM.pdf", 10, 10)
volcano_plot_summary_update(top.table.5mvsom, "5mM vs OM")
dev.off()


#### barplots from the ShinyGO
###### ------------------------------------------------
##### go enrichment analysis on the overlap

overlap.topgmvsomsignif.top5mvsom <- top.table5MvsOM[na.omit(match(rownames(singif.top.GMvsOM.005), rownames(singif.top.5MvsOM))),]
nrow(overlap.topgmvsomsignif.top5mvsom)
# 179

# small subset of the genes
write.table(rownames(overlap.topgmvsomsignif.top5mvsom), "./overlap_topgenesGMvsOM_5mVsOM.txt", quote=F, row.names= F)
# background
write.table(rownames(top.table5MvsOM), "./top_table_5MvsOM.txt", quote = F, row.names = F)


### adding the gene names just to print out the overlap
overlap.list.5mvsom <- overlap.topgmvsomsignif.top5mvsom
top.gene.names.5mmvsom <- annotation.protein.coding.genes[na.omit(match(rownames(overlap.list.5mvsom), annotation.protein.coding.genes$ensembl_gene_id)),]
top.table.5mmvsom  <-overlap.list.5mvsom[na.omit(match(top.gene.names.5mmvsom$ensembl_gene_id,rownames(overlap.list.5mvsom))),]
top.table.5mmvsom$external_gene_names <- top.gene.names.5mmvsom$external_gene_name

openxlsx::write.xlsx(top.table.5mmvsom, "./overlap_diffexp_5MvsOM.xlsx", rowNames=TRUE)

#### shinyGO
shinygo <- read.csv("./enrichment_all_shinyGO.csv")

### shinygo xls 
openxlsx::write.xlsx(shinygo, "./enrichment_shinygo_overlap.xlsx", rowNames=TRUE)

### selected FDR < 0.01
selected.pathways.shinygo <- c("Extracellular matrix organization ", "Angiogenesis ", "Response to hypoxia ", "Regulation of cell adhesion ", "Response to lipid ")
# selected pathways
selected.shinygo.fdr <- shinygo[shinygo$Pathway %in% selected.pathways.shinygo,]
selected.shinygo.fdr$Enrichment.FDR.log <- -log10(selected.shinygo.fdr$Enrichment.FDR)
selected.shinygo.fdr <- selected.shinygo.fdr[order(selected.shinygo.fdr$Enrichment.FDR.log, decreasing = TRUE),]
selected.shinygo.fdr$Order <- c(1,2,3,4,5)


## barplots
pdf("./barplots_shinygo_overlap_gmvsom_5mvsom_.pdf", 8,3)
plt <- ggplot(selected.shinygo.fdr) +
  geom_col(aes(x = Enrichment.FDR.log ,y = reorder(Pathway, -Order)), fill = "dodgerblue4", width = 0.8) + theme_cowplot(12) +
  xlab(bquote(Enrichment ~-Log[10]~italic(FDR))) + ylab(" ") + geom_vline(xintercept=-log10(0.05), colour = "red")

plt
dev.off()


###### GSEA 5mM vs OM
###### GSEA erichment analysis


goTerms <- toTable(GOTERM)
# choose only unique
goTerms.unique <- goTerms[!duplicated(goTerms$go_id),]
goTermsBP <- goTerms.unique[which(goTerms.unique$Ontology == "BP"),]
goTermsMF <- goTerms.unique[which(goTerms.unique$Ontology == "MF"),]
goTermsCC <- goTerms.unique[which(goTerms.unique$Ontology == "CC"),]



######
xx <- as.list(org.Hs.egGO2ALLEGS)

xx.gobp <- xx[match(goTermsBP$go_id,names(xx))]
xx.gobp <- xx.gobp[!is.na(names(xx.gobp))]

xx.gomf <- xx[match(goTermsMF$go_id,names(xx))]
xx.gomf <- xx.gomf[!is.na(names(xx.gomf))]


xx.gocc <- xx[match(goTermsCC$go_id,names(xx))]
xx.gocc <- xx.gocc[!is.na(names(xx.gocc))]


genes.5mvsom <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = rownames(top.table5MvsOM), 
               mart = mart)

genes.5mvsom<-  genes.5mvsom[!duplicated(genes.5mvsom$entrezgene_id),]

topgene.list.entrez.ranked.5MvsOM <- top.table5MvsOM[match(genes.5mvsom$ensembl_gene_id, rownames(top.table5MvsOM)),]
topgene.list.entrez.ranked.5MvsOM$entrez <- genes.5mvsom$entrezgene_id


original_gene_list_t_5mvsom <- topgene.list.entrez.ranked.5MvsOM$t

# # name the vector
names(original_gene_list_t_5mvsom) <- topgene.list.entrez.ranked.5MvsOM$entrez

# # omit any NA values 
gene_list_t_5mvsom<-na.omit(original_gene_list_t_5mvsom)

# # sort the list in decreasing order (required for clusterProfiler)
gene_list_t_5mvsom <- sort(gene_list_t_5mvsom, decreasing = TRUE)


fgsea_go_bp_t_5mvsom <- fgsea(pathways = xx.gobp, 
                        stats = gene_list_t_5mvsom,
                         minSize=15,
                         maxSize=500)

sum(fgsea_go_bp_t_5mvsom[, padj < 0.1]) # 27
head(fgsea_go_bp_t_5mvsom[order(padj), ])

ffbp.om5m <- fgsea_go_bp_t_5mvsom

# ###
ffbp.om5m <- ffbp.om5m[match(goTermsBP$go_id, ffbp.om5m$pathway),] 
ffbp.om5m$Term <- goTermsBP$Term
ffbp.om5m <- na.omit(ffbp.om5m)

# order according to NES and adj p values


ffbp.om5m$p.adj <- p.adjust( ffbp.om5m$pval, method="BH") ### new
ffbp_signif_omvs5m <- ffbp.om5m[which(ffbp.om5m$p.adj < 0.1),]
ffbp_signif_omvs5m$Term <- factor(ffbp_signif_omvs5m$Term,                                    # Factor levels in decreasing order
                  levels = ffbp_signif_omvs5m$Term[order(ffbp_signif_omvs5m$NES, decreasing = FALSE)])

#### ------------
ffgobp_tops <- ffbp_signif_omvs5m[order(ffbp_signif_omvs5m$NES, decreasing = FALSE),]
ff_gobp_highlighted <- ffgobp_tops

ff_gobp_highlighted$Term <- factor(ff_gobp_highlighted$Term,                                    # Factor levels in decreasing order
                  levels = ff_gobp_highlighted$Term[order(ff_gobp_highlighted$NES, decreasing = FALSE)])


c4 <- c(rep("gray80", 7),rep("indianred1", 20))
ff_gobp_highlighted <- cbind(ff_gobp_highlighted, c4)

pdf("./5mvsOM_fgsea_GOBP_signif_ordered.pdf", 10,10)

### this works
ggplot(ff_gobp_highlighted, aes(Term, NES)) + theme_classic() +
    # scale_fill_brewer("", palette="Paired") +
    geom_bar(position= "identity", stat="identity", width=0.8, fill = c4) +
    coord_flip() +
    scale_y_continuous(position = "right") +
    labs( y = "GO BP Normalized Enrichment Score") +
    theme(axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    text = element_text(size=14))
dev.off()

openxlsx::write.xlsx(ff_gobp_highlighted, file = "gsea_results_5mMvsOM.xlsx", rowNames = FALSE)


####### leading edge 

#### run this function before plotting the GSEA plot
leading_edge <- function(observed_info) {
    core_enrichment <- lapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            leading_gene <- runningES$gene[1:i]
        } else {
            i <- which.min(runningES$runningScore)
            leading_gene <- runningES$gene[-c(1:(i-1))]
        }
        return(leading_gene)
    })

    rank <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            rr <- which.max(runningES$runningScore)
        } else {
            i <- which.min(runningES$runningScore)
            rr <- nrow(runningES) - i + 1
        }
        return(rr)
    })

    tags <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    ll <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    N <- nrow(observed_info[[1]]$runningES)
    setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
    signal <- tags * (1-ll) * (N / (N - setSize))

    tags <- paste0(round(tags * 100), "%")
    ll <- paste0(round(ll * 100), "%")
    signal <- paste0(round(signal * 100), "%")
    leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)

    res <- list(rank = rank,
                tags = tags,
                list = ll,
                signal = signal,
                leading_edge = leading_edge,
                core_enrichment = core_enrichment)
    return(res)
}


ffbp.om5m <- fgsea_go_bp_t_5mvsom

# ###
ffbp.om5m <- ffbp.om5m[match(goTermsBP$go_id, ffbp.om5m$pathway),] 
ffbp.om5m$Term <- goTermsBP$Term
ffbp.om5m <- na.omit(ffbp.om5m)

p.adj <- p.adjust( ffbp.om5m$pval, method="BH")

Description <- ffbp.om5m$Term
    

res <- data.frame(
        ID = as.character(ffbp.om5m$pathway),
        Description = unname(Description),
        setSize = ffbp.om5m$size,
        enrichmentScore = ffbp.om5m$ES,
        NES = ffbp.om5m$NES,
        pvalue = ffbp.om5m$pval,
        p.adjust = p.adj,
        # qvalue = qvalues,
        stringsAsFactors = FALSE
    )

res <- res[!is.na(res$pvalue),]
res <- res[ res$pvalue <= 0.1, ]
res <- res[ res$p.adjust <= 0.1, ]
idx <- order(res$p.adjust, -abs(res$NES), decreasing = FALSE)
res <- res[idx, ]

row.names(res) <- res$ID
gseaScores <- getFromNamespace("gseaScores", "DOSE")
observed_info <- lapply(xx.gobp[res$ID], function(gs)
        gseaScores(geneSet=gs,
                   geneList=gene_list_t_5mvsom,
                   exponent=1)
        )

ledge <- leading_edge(observed_info)

res$rank <- ledge$rank
res$leading_edge <- ledge$leading_edge
res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')

params <- list(pvalueCutoff = 0.1,
                        pAdjustMethod = "BH",
                        exponent = 1,
                        minGSSize = 15,
                        maxGSSize = 500
                    )

results <- new("gseaResult",
        result     = res,
        geneSets   = xx.gobp,
        geneList   = gene_list_t_5mvsom,
        params     = params,
        readable   = FALSE
        )

results@organism <- "hs"
results@setType <- "GOBP"
results@keytype <- "UNKNOWN"
 

pdf("./insulin_like.pdf",5,5)
gseaplot2(results, geneSetID = 14, title = results$Description[14])
dev.off()


pdf("./ossification.pdf",5,5)
gseaplot2(results, geneSetID = 20, title = results$Description[20])
dev.off()

pdf("./regulation_pi3k.pdf",5,5)
gseaplot2(results, geneSetID = 22, title = results$Description[22])
dev.off()
