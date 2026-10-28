# =============================================================================
# NicheNet ligand activity analysis  -  Figure 4B-C and Figure S4A-B
#
# Wang, Mercado-Vasquez, Son, ..., Tay. "Encoding innate immune memory in single
# macrophages through cell intrinsic and paracrine feedback" (Cell Systems).
#
# Re-analysis of bulk RNA-seq of BMDMs (GEO: GSE260996). Asks which ligands that
# are upregulated after CpG or polyI:C (and whose receptors are expressed in
# naive BMDMs) best explain the genes upregulated by each stimulus.
#
#   Fig 4B  : ligand activity (corrected AUPR), CpG      -> output_nichenet/CpG_ligand_activity.pdf
#   Fig 4C  : ligand activity (corrected AUPR), polyI:C  -> output_nichenet/polyIC_ligand_activity.pdf
#   Fig S4A : regulatory potential, polyI:C              -> output_nichenet/polyIC_regulatory_potential.pdf
#   Fig S4B : regulatory potential, CpG                  -> output_nichenet/CpG_regulatory_potential.pdf
#
# Fixed and commented version of the unsaved RStudio tab "Untitled1"
# (last edited 2026-02-03), kept unchanged as
# nichenet_Fig4BC_S4AB_ORIGINAL_unsaved_20260203.R. With the default settings
# below, the analysis is the same as the original. Changes are marked [FIX].
# See README_nichenet.md for packages and how to run.
#
# Run with the working directory set to the repository root (the folder containing
# data/), as for R/Fig4_RNAseq.R. Inputs are read from data/rnaseq; outputs go to
# output_nichenet/ and the NicheNet networks are cached in nichenet_networks/.
# Packages (versions used): nichenetr 2.2.1.1, edgeR, limma, dplyr, tidyverse, ggplot2; R 4.5.2.
# Note: not re-run on the analysis computer used for the revision (R not installed there).
# =============================================================================


## ---- 0. Settings -----------------------------------------------------------

DATA_DIR    <- file.path("data", "rnaseq")              # merged_FC.csv, mouse_lr_pair.rds
OUT_DIR     <- "output_nichenet"                        # all results are written here
NETWORK_DIR <- "nichenet_networks"                      # cached NicheNet model files

# Analysis choices. The defaults reproduce the original script, so they should
# give the published panels. Change them only on purpose and report it in Methods.

# Ligand set. FALSE (original): every ligand up after CpG or polyI:C (CellTalkDB),
# then kept only if NicheNet's lr_network links it to a receptor that is
# expressed in naive cells. TRUE: also apply the two extra filters used for the
# Fig 4A heatmap (CellTalkDB receptor intersection + manual exclusion list), so
# that 4B-C use exactly "the set of upregulated ligands in (A)" as the legend says.
USE_FIG4A_LIGAND_FILTER <- FALSE
FIG4A_EXCLUDED <- c("Mmp13", "Adam17", "Adam9", "Calr", "Col18a1", "Hspa1b",
                    "Mmp9", "Col4a2", "Fbn1", "Col4a5", "Bgn", "Plat")

# Number of top ligands (by corrected AUPR), as in the published panels: 10 for
# CpG (Fig 4B) and 11 for polyI:C (Fig 4C). Ligands without a target above the
# 0.25 cutoff are not shown in the regulatory potential heatmaps (Fig S4A-B).
N_TOP <- c(CpG = 10, polyIC = 11)

TARGETS_PER_LIGAND    <- 50    # top predicted targets per ligand (S4A-B)
REG_POTENTIAL_CUTOFF  <- 0.25  # quantile cutoff for regulatory potential heatmap
RECEPTOR_MIN_LOGCPM   <- 6     # receptor must have mean logCPM > 6 in naive (FF)

# NicheNet v2 mouse model (Zenodo record 7074291)
NICHENET_URL   <- "https://zenodo.org/record/7074291/files/"
LR_NETWORK_RDS <- "lr_network_mouse_21122021.rds"
LT_MATRIX_RDS  <- "ligand_target_matrix_nsga2r_final_mouse.rds"


## ---- 1. Packages -------------------------------------------------------------
# [FIX] The original loaded library(dp), which does not exist (probably meant
# dplyr), and used org.Mm.eg.db without loading it. So it only ran after the main
# edgeR script had been run in the same session. cowplot/ggpubr were unused.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(nichenetr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(NETWORK_DIR, showWarnings = FALSE)


## ---- 2. Counts -> logCPM (same as edgeR_limma_script_LO.R) --------------------

counts <- read.csv(file.path(DATA_DIR, "merged_FC.csv"), header = TRUE)

# 11 conditions x 3 replicates. Letter 1 = phase I, letter 2 = phase II stimulus.
# F = media, C = CpG, P = polyI:C, L = LPS, T = TNF.  FF = naive, FC = CpG, FP = polyI:C.
sample_groups <- c("CF", "CF", "CF", "CL", "CL", "CL", "CT", "CT", "CT",
                   "FC", "FC", "FC", "FF", "FF", "FF", "FL", "FL", "FL",
                   "FP", "FP", "FP", "FT", "FT", "FT", "PF", "PF", "PF",
                   "PL", "PL", "PL", "PT", "PT", "PT")
group <- factor(sample_groups)

# First column of merged_FC.csv is a row index; "gene" holds Entrez IDs
counts <- tibble::column_to_rownames(counts[-1], var = "gene")
stopifnot(ncol(counts) == length(sample_groups))  # [FIX] guard against wrong file/order

DGEobj <- DGEList(counts = counts, group = group)
keep   <- filterByExpr(DGEobj, min.count = 25, group = group)
DGEobj <- DGEobj[keep, , keep.lib.sizes = FALSE]
DGEobj <- calcNormFactors(DGEobj)                     # TMM
logCPM <- cpm(DGEobj, log = TRUE, prior.count = 3)

# Entrez ID -> gene symbol.
# [FIX] mapIds(multiVals = "first") always returns one symbol per row, in order.
# The original assigned select()$SYMBOL directly, which silently misaligns rows
# (or errors) if any ID maps to more than one symbol. For 1:1 mappings the
# result is identical.
rownames(logCPM) <- mapIds(org.Mm.eg.db, keys = rownames(logCPM),
                           column = "SYMBOL", keytype = "ENTREZID",
                           multiVals = "first")
colnames(logCPM) <- sample_groups

# [FIX] The original wrote logCPM.csv here, overwriting the shared data file.
# Removed. Uncomment to save a copy in the output folder instead:
# write.csv(logCPM, file.path(OUT_DIR, "logCPM.csv"))


## ---- 3. Differential expression (limma on logCPM) ------------------------------
# Same model and cutoffs as the main script: BH FDR < 0.01 and |log2 FC| > 1.

design <- model.matrix(~0 + group)
contr.matrix <- makeContrasts(CF = groupCF - groupFF,
                              CL = groupCL - groupFF,
                              CT = groupCT - groupFF,
                              FC = groupFC - groupFF,   # CpG vs naive
                              FL = groupFL - groupFF,
                              FP = groupFP - groupFF,   # polyI:C vs naive
                              FT = groupFT - groupFF,
                              PF = groupPF - groupFF,
                              PL = groupPL - groupFF,
                              PT = groupPT - groupFF,
                              FPFC = groupFP - groupFC,
                              levels = colnames(design))

fit  <- lmFit(logCPM, design)
fit  <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit)
dt   <- decideTests(efit, adjust.method = "BH", p.value = 0.01, lfc = 1)
print(summary(dt))

# [FIX] Columns are selected by contrast name ("FC", "FP") instead of position
# (4, 6). Same columns, but safe if the contrast list is edited.


## ---- 4. Candidate ligands and receptors (CellTalkDB) ---------------------------

mouse_lr_pair <- readRDS(file.path(DATA_DIR, "mouse_lr_pair.rds"))
# [FIX] check the columns used by position in the original (2 = ligand, 3 = receptor)
stopifnot(all(c("ligand_gene_symbol", "receptor_gene_symbol") %in% names(mouse_lr_pair)))

# Ligands upregulated after CpG or polyI:C
up_any    <- rownames(dt)[dt[, "FC"] > 0 | dt[, "FP"] > 0]
deltaligs <- unique(mouse_lr_pair$ligand_gene_symbol)
deltaligs <- deltaligs[deltaligs %in% up_any]

# Their receptors, kept only if expressed in naive BMDMs (mean logCPM in FF > 6)
deltalig_receptors <- mouse_lr_pair$receptor_gene_symbol[mouse_lr_pair$ligand_gene_symbol %in% deltaligs]
deltalig_receptors <- unique(deltalig_receptors[deltalig_receptors %in% rownames(logCPM)])
naive_mean <- rowMeans(logCPM[deltalig_receptors, colnames(logCPM) == "FF", drop = FALSE])
deltalig_receptors <- deltalig_receptors[naive_mean > RECEPTOR_MIN_LOGCPM]

if (USE_FIG4A_LIGAND_FILTER) {
  # Same filters as the Fig 4A heatmap in edgeR_limma_script_LO.R (lines ~830-836)
  deltaligs <- intersect(deltaligs,
                         mouse_lr_pair$ligand_gene_symbol[mouse_lr_pair$receptor_gene_symbol %in% deltalig_receptors])
  deltaligs <- setdiff(deltaligs, FIG4A_EXCLUDED)
}


## ---- 5. NicheNet v2 model --------------------------------------------------------
# [FIX] Files are downloaded once and cached locally (the original re-downloaded
# them on every run). The timeout is set before the first download (the original
# set it after the first file). weighted_networks was downloaded but never
# used, so it is no longer loaded.

options(timeout = 600)
get_nichenet_file <- function(file) {
  local_file <- file.path(NETWORK_DIR, file)
  if (!file.exists(local_file)) {
    download.file(paste0(NICHENET_URL, file), local_file, mode = "wb")
  }
  readRDS(local_file)
}
lr_network           <- get_nichenet_file(LR_NETWORK_RDS)   # columns: from (ligand), to (receptor)
ligand_target_matrix <- get_nichenet_file(LT_MATRIX_RDS)    # rows = targets, cols = ligands

# Potential ligands = candidate ligands with a NicheNet ligand-receptor link to a
# receptor expressed in naive cells
potential_ligands <- lr_network %>%
  filter(from %in% deltaligs & to %in% deltalig_receptors) %>%
  pull(from) %>% unique()
cat("Potential ligands:", length(potential_ligands), "\n")
writeLines(potential_ligands, file.path(OUT_DIR, "potential_ligands.txt"))


## ---- 6. Ligand activity and regulatory potential, per stimulus ----------------
# The two blocks in the original were copy-pasted. They are now one function.
#
# Note on background: the original passed "all genes NOT in the gene set" as
# background, and NicheNet recommends "all expressed genes". Both give the same
# result, because nichenetr itself removes the gene set from the background
# (convert_gene_list_settings_evaluation, nichenetr 2.2.1.1). No change needed.

run_nichenet <- function(contrast, label, n_top, color) {

  in_model   <- rownames(efit) %in% rownames(ligand_target_matrix)
  is_up      <- dt[, contrast] == 1
  # [FIX] unique(): if two Entrez IDs share one symbol, the original counted that
  # gene twice. The number of such cases is printed. If it is 0, results are identical.
  n_dup <- sum(duplicated(rownames(efit)[in_model]))
  if (n_dup > 0) cat(label, ": duplicated gene symbols collapsed:", n_dup, "\n")
  geneset    <- unique(rownames(efit)[is_up & in_model])    # genes up vs naive
  background <- unique(rownames(efit)[!is_up & in_model])   # as in original (see note above)

  ligand_activities <- predict_ligand_activities(
    geneset                    = geneset,
    background_expressed_genes = background,
    ligand_target_matrix       = ligand_target_matrix,
    potential_ligands          = potential_ligands) %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(-aupr_corrected))

  # [FIX] results are saved (the original only printed them to the console)
  write.csv(ligand_activities, file.path(OUT_DIR, paste0(label, "_ligand_activities.csv")),
            row.names = FALSE)

  # top_n() keeps ties, same as the original
  best_upstream_ligands <- ligand_activities %>%
    top_n(n_top, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand)

  # Regulatory potential: top predicted targets of each ligand within the gene set
  active_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset,
           ligand_target_matrix = ligand_target_matrix,
           n = TARGETS_PER_LIGAND) %>%
    bind_rows()
  write.csv(active_links_df, file.path(OUT_DIR, paste0(label, "_ligand_target_links.csv")),
            row.names = FALSE)

  active_links <- prepare_ligand_target_visualization(
    ligand_target_df     = active_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff               = REG_POTENTIAL_CUTOFF)

  order_ligands <- intersect(best_upstream_ligands, colnames(active_links)) %>% rev()
  order_targets <- active_links_df$target %>% unique() %>% intersect(rownames(active_links))

  # [FIX] report ligands that are in the top list but dropped from the heatmap
  dropped <- setdiff(best_upstream_ligands, colnames(active_links))
  cat(label, ": top ligands =", length(best_upstream_ligands),
      "| shown in regulatory potential heatmap =", length(order_ligands),
      if (length(dropped)) paste("| dropped:", paste(dropped, collapse = ", ")) else "", "\n")

  # Fig S4A/B: regulatory potential heatmap
  vis_ligand_target <- t(active_links[order_targets, order_ligands, drop = FALSE])
  p_target <- make_heatmap_ggplot(vis_ligand_target,
                                  "Secreted ligands", "Upregulated ligand targets",
                                  color = color, legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke", high = color)

  # Fig 4B/C: ligand activity (corrected AUPR) heatmap
  vis_ligand_aupr <- ligand_activities %>%
    filter(test_ligand %in% best_upstream_ligands) %>%
    as.data.frame() %>%
    column_to_rownames("test_ligand")
  vis_ligand_aupr <- vis_ligand_aupr[, "aupr_corrected", drop = FALSE] %>%
    arrange(aupr_corrected) %>% as.matrix()

  list(label = label, color = color, ligand_activities = ligand_activities,
       best_upstream_ligands = best_upstream_ligands,
       vis_ligand_aupr = vis_ligand_aupr, p_target = p_target)
}

res <- list(
  CpG    = run_nichenet("FC", "CpG",    N_TOP[["CpG"]],    "red"),   # Fig 4B, S4B
  polyIC = run_nichenet("FP", "polyIC", N_TOP[["polyIC"]], "blue")   # Fig 4C, S4A
)
# Colors are as in the original (CpG red, polyI:C blue). Note that Fig 4J uses
# the opposite convention (blue = CpG enriched, red = polyI:C enriched). Check
# the colors against the final figure.


## ---- 7. Figures ------------------------------------------------------------------
# [FIX] The original labelled the y axis "Prioritized CAF-ligands" (left over from
# the NicheNet tutorial) and fixed the AUPR color scale to 0.05-0.2, which turns
# values outside that range grey. Here the scale covers the range of both
# stimuli, so 4B and 4C are directly comparable. The range is printed below.

aupr_range <- range(c(res$CpG$vis_ligand_aupr, res$polyIC$vis_ligand_aupr))
cat("Corrected AUPR range in the plotted ligands:", round(aupr_range, 4), "\n")

for (r in res) {
  p_aupr <- make_heatmap_ggplot(r$vis_ligand_aupr,
                                "Prioritized ligands", "Ligand activity",
                                color = r$color, legend_title = "AUPR") +
    theme(axis.text.x.top = element_blank()) +
    scale_fill_gradient(low = "whitesmoke", high = r$color,
                        limits = aupr_range, name = "Corrected\nAUPR")

  ggsave(file.path(OUT_DIR, paste0(r$label, "_ligand_activity.pdf")), p_aupr,
         width = 2.5, height = 4)
  ggsave(file.path(OUT_DIR, paste0(r$label, "_regulatory_potential.pdf")), r$p_target,
         width = 10, height = 4)
  print(p_aupr)
  print(r$p_target)
}


## ---- 8. Record versions --------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("Done. Results in", normalizePath(OUT_DIR), "\n")
