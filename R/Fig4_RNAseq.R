# =====================================================================================================
# Figure 4A, 4J and Figure S4G related to
# Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
# Wang AG, Mercado Vasquez JG, Son M, et al. Corresponding author: tays@uchicago.edu
#
# Re-analysis of bulk RNA-seq of BMDMs (GEO: GSE260996; counts in data/rnaseq/merged_FC.csv, 11 conditions x 3
# replicates; column names = phase I stimulus + phase II stimulus: F = media, C = CpG, P = polyI:C, L = LPS, T = TNF).
#   - edgeR: filterByExpr, TMM normalization, logCPM (prior.count = 3)
#   - limma: lmFit / eBayes, contrasts vs FF (media), decideTests (BH < 0.01, |log2 FC| > 1)
#   - Figure 4A: upregulated ligands (CellTalkDB mouse ligand-receptor pairs, data/rnaseq/mouse_lr_pair.rds) whose
#     receptor is expressed in naive cells (mean logCPM > 6); limma::coolmap scales each gene (row Z-score)
#   - Figure S4G: negative-feedback genes, same heatmap settings
#   - Figure 4J: EGSEA plotPathway (Pathview) of the KEGG NF-kB pathway, contrast polyI:C (FP) vs CpG (FC)
# Figure 4B-C and S4A-B (NicheNet ligand activity) are produced by R/Fig4BC_S4AB_NicheNet.R.
# Packages: limma, edgeR, dplyr, gplots, org.Mm.eg.db, ggplot2, RColorBrewer, EGSEA (Bioconductor).
# Run with the working directory set to the repository root (the folder containing data/).
# Note: not re-run on the analysis computer used for the revision (R not installed there).
# =====================================================================================================
library(limma)
library(edgeR)
library(dplyr)
library(gplots)
library(org.Mm.eg.db)
library(ggplot2)
library("RColorBrewer")
library(EGSEA)

#Preprocessing from count matrix
DGEobj = read.csv("data/rnaseq/merged_FC.csv", header= TRUE)
group=as.factor( c('CF', 'CF', 'CF', 'CL', 'CL', 'CL', 'CT', 'CT', 'CT', 
                   'FC', 'FC', 'FC', 'FF', 'FF', 'FF', 'FL', 'FL', 'FL',
                   'FP', 'FP', 'FP', 'FT', 'FT', 'FT', 'PF', 'PF', 'PF', 
                   'PL', 'PL', 'PL', 'PT', 'PT', 'PT') )

DGEobj = tibble::column_to_rownames(DGEobj[-1], var="gene")
DGEobj = DGEList(counts = DGEobj, group=group)
keep <- filterByExpr(DGEobj, min.count = 25, group=group)
DGEobj <- DGEobj[keep,,keep.lib.sizes=FALSE]
DGEobj <- calcNormFactors(DGEobj)

logCPM <- cpm(DGEobj, log=TRUE, prior.count = 3)
geneid <- rownames(logCPM)
genes <- AnnotationDbi::select(org.Mm.eg.db, keys = geneid, column = c('SYMBOL'), keytype = "ENTREZID")


rownames(logCPM) <- genes$SYMBOL
colnames(logCPM) <-  c('CF', 'CF', 'CF', 'CL', 'CL', 'CL', 'CT', 'CT', 'CT', 
                       'FC', 'FC', 'FC', 'FF', 'FF', 'FF', 'FL', 'FL', 'FL',
                       'FP', 'FP', 'FP', 'FT', 'FT', 'FT', 'PF', 'PF', 'PF', 
                       'PL', 'PL', 'PL', 'PT', 'PT', 'PT') 
write.csv(logCPM, "logCPM.csv")   # optional output

design <- model.matrix(~0+group)
contr.matrix <- makeContrasts(CF = groupCF - groupFF, 
                              CL = groupCL - groupFF, 
                              CT = groupCT - groupFF, 
                              FC = groupFC - groupFF, 
                              FL = groupFL - groupFF, 
                              FP = groupFP - groupFF, 
                              FT = groupFT - groupFF, 
                              PF = groupPF - groupFF, 
                              PL = groupPL - groupFF,
                              PT = groupPT - groupFF,
                              FPFC = groupFP - groupFC,
                              levels = colnames(design))

fit <- lmFit(logCPM, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)
dt <- decideTests(efit, adjust.method = "BH", p.value = 0.01,lfc = 1)


#Ligand-Receptor stuff
mouse_lr_pair = readRDS('data/rnaseq/mouse_lr_pair.rds')

#PIC or CPG
deltas = rownames( dt[ dt[,4] > 0 | dt[,6] > 0, ] )
deltaligs = unique(mouse_lr_pair$ligand_gene_symbol)
deltaligs = deltaligs[deltaligs %in% deltas]
#deltaligs = deltaligs[rowMeans(logCPM[deltaligs, colnames(logCPM) == "FP"] ) > 3 | rowMeans(logCPM[deltaligs, colnames(logCPM) == "FC"] ) > 3]

deltalig_receptors = mouse_lr_pair[mouse_lr_pair$ligand_gene_symbol %in% deltaligs, 3]
deltalig_receptors = deltalig_receptors[deltalig_receptors %in% rownames(logCPM)]
deltalig_receptors = deltalig_receptors[rowMeans(logCPM[deltalig_receptors, colnames(logCPM) == "FF"] ) > 6 ]
deltaligs = intersect( deltaligs, mouse_lr_pair[mouse_lr_pair$receptor_gene_symbol %in% deltalig_receptors, 2] )

meanexp = rowMeans(logCPM[deltaligs, colnames(logCPM) == "FP"] ) + rowMeans(logCPM[deltaligs, colnames(logCPM) == "FC"] ) - rowMeans(logCPM[deltaligs, colnames(logCPM) == "FF"] )
deltaligs = deltaligs[order(meanexp, decreasing = TRUE)]

# Figure 4A: ligand heatmap; columns FP (polyI:C), FC (CpG), FF (media), 3 replicates each
p <- coolmap(logCPM[deltaligs[!( deltaligs %in% c("Mmp13", "Adam17", "Adam9", "Calr", "Col18a1", "Hspa1b", "Mmp9", "Col4a2", "Fbn1", "Col4a5", "Bgn", "Plat") )], 
                    c( 19:21, 10:12, 13:15) ], linkage.row="none", linkage.col="none", show.dendrogram="row", xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))



# Figure S4G: known negative-feedback genes (same samples and colour breaks as Figure 4A)
p_s4g <- coolmap(logCPM[c("Nfkbia", "Nfkbie", "Nfkbib", "Cyld", "Irak3", "Tnfaip3", "Tax1bp1", "Itch", "Nfkb2", "Trim30a", "Socs1", "Socs3"),
                    c( 19:21, 10:12, 13:15) ], linkage.row="none", linkage.col="none", show.dendrogram="row", xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))

# Figure 4J: EGSEA / Pathview (contrast 11 = FP - FC, polyI:C vs CpG)
#EGSEa
foo3 = rownames( dt[ dt[,"FPFC"] != 0, ] )
foo3 = foo3[!is.na(foo3)]
write.csv(foo3, file="FPvsFC.csv")

gs.annots = buildIdx(entrezIDs = rownames(DGEobj), species = "mouse", msigdb.gsets = "none",
                     kegg.exclude = c("Metabolism"))


test = egsea.cnt(counts = DGEobj, group = group, design = design,
                 contrasts = contr.matrix, gs.annots = gs.annots, symbolsMap = genes,
                 baseGSEAs = egsea.base()[-c(2, 12)], sort.by = "avg.rank",
                 num.threads = 4, report = FALSE)

topSets(test, contrast = 11, gs.label = "kegg", number = 15)
plotPathway(test, contrast = 11, "NF-kappa B signaling pathway", gs.label = "kegg", file.name = "NF-kappa B signaling pathway")
plotPathway(test, contrast = 11, "TNF signaling pathway", gs.label = "kegg", file.name = "TNF signaling pathway")
# (not in the paper) heatmap of all FP vs FC differentially expressed genes
coolmap(logCPM[foo3, c( 19:21, 10:12, 13:15)], linkage.row="ward", linkage.col="none", show.dendrogram="row",
        xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))
