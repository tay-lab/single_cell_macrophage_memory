# =====================================================================================================
# Figure 4B-C and Figure S4A-B related to
# Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
# Wang AG, Mercado-Vásquez G, Son M, et al. Corresponding author: tays@uchicago.edu
#
# NicheNet analysis of bulk RNA-seq of BMDMs (GEO: GSE260996; counts in data/rnaseq/merged_FC.csv).
# Candidate ligands are genes upregulated after CpG or polyI:C (CellTalkDB mouse ligand-receptor pairs,
# data/rnaseq/mouse_lr_pair.rds) with a receptor expressed in naive cells (mean logCPM > 6). NicheNet ranks
# them by how well their predicted targets explain the genes upregulated by each stimulus (corrected AUPR).
#   - Figure 4B / 4C: ligand activity (corrected AUPR) of the top 10 ligands, CpG / polyI:C
#   - Figure S4A / S4B: regulatory potential of these ligands on their top 50 predicted targets among the
#     upregulated genes, CpG / polyI:C (scores below the 25% quantile are set to 0). Figure S4B shows the
#     first 89 of the 108 target genes.
# The figures are saved as PDF files in output_nichenet/. The mouse NicheNet v2 networks are downloaded
# from Zenodo (record 7074291) at each run.
# Packages: limma, edgeR, org.Mm.eg.db (Bioconductor), nichenetr 2.2.1.1, tidyverse. Tested with R 4.5.2.
# Run with the working directory set to the repository root (the folder containing data/).
# =====================================================================================================
library(limma)
library(edgeR)
library(org.Mm.eg.db)
library(nichenetr)
library(tidyverse)

out_dir = "output_nichenet"
dir.create(out_dir, showWarnings = FALSE)

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
summary(dt)

#Ligand-Receptor pairs (CellTalkDB)
mouse_lr_pair = readRDS('data/rnaseq/mouse_lr_pair.rds')

#Ligands upregulated after polyI:C (column 6, FP) or CpG (column 4, FC)
deltas = rownames( dt[ dt[,4] > 0 | dt[,6] > 0, ] )
deltaligs = unique(mouse_lr_pair$ligand_gene_symbol)
deltaligs = deltaligs[deltaligs %in% deltas]

#Their receptors, kept if expressed in naive cells (FF)
deltalig_receptors = mouse_lr_pair[mouse_lr_pair$ligand_gene_symbol %in% deltaligs, 3]
deltalig_receptors = deltalig_receptors[deltalig_receptors %in% rownames(logCPM)]
deltalig_receptors = deltalig_receptors[rowMeans(logCPM[deltalig_receptors, colnames(logCPM) == "FF"] ) > 6 ]

#Mouse NicheNet v2 networks
options(timeout = 600)
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))

expressed_ligands <- deltaligs
expressed_receptors <- deltalig_receptors

potential_ligands <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique()

head(potential_ligands)


#### polyI:C (Figure 4C and Figure S4B)
up_genes_FP <- rownames(efit)[dt[, 6] == 1] %>% .[. %in% rownames(ligand_target_matrix)]
background_genes_FP <- rownames(efit)[dt[, 6] != 1] %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = up_genes_FP,
                                               background_expressed_genes = background_genes_FP,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(-aupr_corrected)))
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = up_genes_FP,
         ligand_target_matrix = ligand_target_matrix,
         n = 50) %>% bind_rows()
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
order_targets <- head(order_targets, 89)   # Figure S4B shows the first 89 target genes

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

# Figure S4B: regulatory potential, polyI:C
p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Secreted ligands", "Upregulated ligand targets",
                                               color = "blue", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high ="blue") +
  guides(fill = guide_colorbar(title = "Regulatory potential", title.position = "bottom", title.hjust = 0.5)) +
  theme(legend.position = "top", legend.justification = "right")

p_ligand_target_network
ggsave(file.path(out_dir, "FigS4B_polyIC_regulatory_potential.pdf"), p_ligand_target_network, width = 10, height = 4)

vis_ligand_aupr <- as_tibble(ligand_activities) %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand")
vis_ligand_aupr = vis_ligand_aupr[, "aupr_corrected", drop = FALSE] %>% arrange(aupr_corrected) %>% as.matrix()

# Figure 4C: ligand activity (corrected AUPR), polyI:C
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "polyI:C induced", "",
                                     color = "blue", legend_title = "Ligand activity") +
  theme(axis.text.x.top = element_blank(), legend.position = "bottom") +
  scale_fill_gradient(low = "whitesmoke", high = "blue", limits = c(0.05, 0.2),
                      breaks = seq(0.05, 0.2, by = 0.05), name = "Ligand activity") +
  guides(fill = guide_colorbar(title.position = "bottom", title.hjust = 0.5))
p_ligand_aupr
ggsave(file.path(out_dir, "Fig4C_polyIC_ligand_activity.pdf"), p_ligand_aupr, width = 2.5, height = 4.5)




#### CpG (Figure 4B and Figure S4A)
up_genes_FC <- rownames(efit)[dt[, 4] == 1] %>% .[. %in% rownames(ligand_target_matrix)]
background_genes_FC <- rownames(efit)[dt[, 4] != 1] %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = up_genes_FC,
                                               background_expressed_genes = background_genes_FC,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(-aupr_corrected)))
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = up_genes_FC,
         ligand_target_matrix = ligand_target_matrix,
         n = 50) %>% bind_rows()
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

# Figure S4A: regulatory potential, CpG
p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Secreted ligands", "Upregulated ligand targets",
                                               color = "red", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high ="red") +
  guides(fill = guide_colorbar(title = "Regulatory potential", title.position = "bottom", title.hjust = 0.5)) +
  theme(legend.position = "top", legend.justification = "right")

p_ligand_target_network
ggsave(file.path(out_dir, "FigS4A_CpG_regulatory_potential.pdf"), p_ligand_target_network, width = 10, height = 4)

vis_ligand_aupr <- as_tibble(ligand_activities) %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand")
vis_ligand_aupr = vis_ligand_aupr[, "aupr_corrected", drop = FALSE] %>% arrange(aupr_corrected) %>% as.matrix()

# Figure 4B: ligand activity (corrected AUPR), CpG
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "CpG induced", "",
                                     color = "red", legend_title = "Ligand activity") +
  theme(axis.text.x.top = element_blank(), legend.position = "bottom") +
  scale_fill_gradient(low = "whitesmoke", high = "red", limits = c(0.05, 0.2),
                      breaks = seq(0.05, 0.2, by = 0.05), name = "Ligand activity") +
  guides(fill = guide_colorbar(title.position = "bottom", title.hjust = 0.5))
p_ligand_aupr
ggsave(file.path(out_dir, "Fig4B_CpG_ligand_activity.pdf"), p_ligand_aupr, width = 2.5, height = 4.5)
