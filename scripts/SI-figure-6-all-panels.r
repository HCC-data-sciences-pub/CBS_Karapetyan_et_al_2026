# Methods for SI Figure 6A-F

library(ggplot2)
library(ggrepel)
library(Seurat)
library(SingleCellExperiment)
library(edgeR)
library(tidyr)
library(tidyverse)
library(dplyr)
library(stringr)

# Starting with Seurat object that was integrated using Harmony, and cell type contaminants were removed

# Gene filtering:
filter_uncharacterized_genes <- function(seurat_obj) {
  
  remove_patterns <- c(
    "^ENSG",
    "^ENS",
    "^LOC",
    "^LINC",
    "-AS[0-9]*$",
    "-OT[0-9]*",
    "-DT$",
    "-IT[0-9]*",
    "^AC[0-9]+\\.[0-9]+",
    "^AL[0-9]+\\.[0-9]+",
    "^AP[0-9]+\\.[0-9]+",
    "^AF[0-9]+\\.[0-9]+",
    "^AJ[0-9]+\\.[0-9]+",
    "^BX[0-9]+\\.[0-9]+",
    "^CR[0-9]+\\.[0-9]+",
    "^CT[0-9]+\\.[0-9]+",
    "^CU[0-9]+\\.[0-9]+",
    "^FP[0-9]+\\.[0-9]+",
    "^Z[0-9]+\\.[0-9]+",
    "^[A-Z]{2}[0-9]+\\.[0-9]+",
    "^KIAA[0-9]+",
    "^HIST[0-9]",
    "^H[1-4][A-Z]",
    "^H[1-4]-",
    "^MT-",
    "^MTRNR",
    "^MTRN",
    "^FAM[0-9]+",
    "^C[0-9]+orf[0-9]+"
  )
  gene_names <- rownames(seurat_obj)
  remove_mask <- rep(FALSE, length(gene_names))
  
  for (pattern in remove_patterns) {
    remove_mask <- remove_mask | grepl(pattern, gene_names)
  }
  
  # Apply IG gene removal pattern (but keep IGF)
  ig_remove_pattern <- "^(IGH|IGL|IGK|IGSF|IGBP|IGHMBP)"
  ig_keep_pattern <- "^IGF"
  
  ig_matches <- grepl(ig_remove_pattern, gene_names)
  ig_keep <- grepl(ig_keep_pattern, gene_names)
  remove_mask <- remove_mask | (ig_matches & !ig_keep)
  
  genes_to_keep <- gene_names[!remove_mask]
  
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)
  return(seurat_obj)
}
seurat_filtered <- filter_uncharacterized_genes(seurat_obj)

# Converted object metadata to dataframe

# Calculated scores: APS, CYT, IFNg, and Tex
ifng_genelist = c("CCL5", "CD27", "CD274","CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
cyt_genelist = c("GZMA", "GZMH", "GZMB", "PRF1", "NKG7")
aps_genelist <- c(
    "HSPA5","PIK3C3","SEC23A","TAPBP","SEC24A","BECN1","ATG14","CANX","HLA-F",
    "SEC31A","SEC24B","SEC24D","SAR1B","SEC13","ERAP1","ERAP2","B2M","PDIA3",
     "TAP1","SEC24C","CALR","BCAP31","PIK3R4","TAP2","HLA-C","HLA-E","HLA-G",
     "HLA-F","TAPBP","TAP2","TAPBP","TAP1","TAP2","HLA-E","HLA-A","HLA-G",
     "HLA-F","TAP2","TAP1","HLA-E","TAP2","TAP1","TAP1","TAP2","HLA-E","HLA-F",
     "HLA-E","HLA-G","TAP1","TAPBP","TAP2","TAP1","HLA-G","HLA-E","HLA-B",
     "HLA-F","HLA-G","HLA-G","TAPBP","HLA-E","HLA-G","HLA-F","TAP2","B2M",
     "HLA-G","CANX")  
t_ex_genelist = c("HAVCR2","KLRC1","SRGAP3","FABP5","ENTPD1","CSF1","PDCD1","TNFSF4","GZMA","MKI67","LGALS3","TOX","TIGIT","FUT8","TNFRSF9","IL2RB","LAG3","NAB1","SLC2A8","CTLA4","ITM2A","GCNT1","PSMB8","PLSCR1","KCNK5","TNIP3","BUB1","PSMA1","TANK","SUB1","ACOT7","NEDD9","HMGB2","SDF2L1","ETFB","EZH2","MCM3","MIS18BP1","RAD21","PSMD9","SLC43A3","PMF1")

# Score calculation, example for Tex score:
pseudo_objT <- AggregateExpression(seu_obj,group.by = c('Patient.Number','time3'),return.seurat = TRUE,slot = "counts")
counts_matrixt <- GetAssayData(object = pseudo_objT , assay = "RNA", layer = "counts")
min_cells <- 20
keep_samples <- names(cell_counts)[cell_counts >= min_cells]
cell_counts = table(seu_obj@meta.data$Patient.Number, seu_obj@meta.data$time3)
keep_samples <- rownames(cell_counts)[cell_counts >= min_cells]
pseudo_objTt <- pseudo_objT[, rownames(pseudo_objT) %in% keep_samples]
dget <- edgeR::DGEList(counts = counts_matrixt)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 5
dget <- dget[keep, , keep.lib.sizes = FALSE]
dget <- edgeR::calcNormFactors(dget, method = "TMM")
logCPMt <- edgeR::cpm(dget, log = TRUE, prior.count = 1)

# Use the specific score genelists, example with Tex:
existing_genes_tex <- intersect(t_ex_genelist, rownames(logCPMt))
logCPM_scaledt <- t(scale(t(logCPMt[existing_genes_tex, ])))
sig_scoret <- colMeans(logCPM_scaledt)
sgscoret = as.data.frame(sig_scoret)
pso8 = as.data.frame(pseudo_objT@meta.data)
pso888 = merge(pso8,sgscoret, by="cell_counts", all.x=TRUE)
# "psl" and "tsl" are the metadata dataframes with the previous 3 scores (IFNg, CYT, APS)
tsl_tex = merge(tsl, pso888, by="cell_counts", all.x=TRUE)
psl_tex2 = merge(psl, pso888, by="cell_counts", all.x=TRUE)
psll_tex4 = psl_tex2 %>%
    pivot_longer(
    cols = c("IFNg.score", "CYT.score", "APS.score", "Tex.score"),
    names_to = "Score.Type",
    values_to = "Score.Value")
tsll_tex2 = tsl_tex %>%
    pivot_longer(
    cols = c("IFNg.score", "CYT.score", "APS.score", "Tex.score"),
    names_to = "Score.Type",
    values_to = "Score.Value")

# Patients that had samples for both Baseline and Week04 timepoints:

til_pat_pairs = c("UOP-002","UOP-003","UOP-004","UOP-005","UOP-006","UOP-007","UOP-008","UOP-009","UOP-014","UOP-015","UOP-018","UOP-023","UOP-024","UOP-026","UOP-028","UOP-029","UOP-030","UOP-031","UOP-034","UOP-036","UOP-038","UOP-042","UOP-043","UOP-044")

pbmc_pat_pairs = c("UOP-002","UOP-003","UOP-004","UOP-005","UOP-006","UOP-007","UOP-008","UOP-009","UOP-014","UOP-015","UOP-017","UOP-018","UOP-020","UOP-026","UOP-028","UOP-029","UOP-030","UOP-033","UOP-034","UOP-035","UOP-036","UOP-038","UOP-040","UOP-041","UOP-042")

til_pp_tsll2 <- tsl_wtexl[tsl_wtexl$Patient.Number %in% til_pat_pairs, ]
pbmc_pp_psll2 <- psl_texl[psl_texl$Patient.Number %in% pbmc_pat_pairs, ]

# 6A baseline only TIL all scores:
# subset data table to have only Baseline timepoint present
base_tsll2_leadin_scm = ggplot(base_tsll2, aes(x = Lead_in, y = Score.Value, fill=Lead_in)) +
geom_boxplot() +
labs(title = "Baseline-only: All Scores by Lead-in (TIL)",
y = "Signature.score",
x = "Lead-in") +facet_wrap(~Score.Type, scales = "free_x", ncol = 4) +theme_minimal()+geom_jitter(width = 0.2, alpha = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_manual(values = c("Combo" = "#DB2EE6", "Nivo" = "#5DCACC", "Rela" = "#737172"))+stat_compare_means(paired = FALSE, method='anova', aes(label = paste0(" ", after_stat(p.format))),size=3, label.y = 1.8, label.x = 1.7)

# 6B baseline only PBMC all scores:
# subsets data table to have only Baseline timepoint present
base_psll2_leadin_scm = ggplot(base_psll2, aes(x = Lead_in, y = Score.Value,    fill=Lead_in)) +
    geom_boxplot() +
    labs(title = "Baseline-only: All Scores by Lead-in (PBMC)",
    y = "Signature.score",
    x = "Lead-in") +facet_wrap(~Score.Type, scales = "free_x", ncol = 4) +theme_minimal()+geom_jitter(width = 0.2, alpha = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_manual(values = c("Combo" = "#DB2EE6", "Nivo" = "#5DCACC", "Rela" = "#737172"))+stat_compare_means(paired = FALSE, method='anova', aes(label = paste0(" ", after_stat(p.format))),size=3, label.y = 1.8, label.x = 1.7)

# 6C PBMC paired all scores:
pbmc_pp_psll_1 = ggplot(data = pbmc_pp_psll2, mapping = aes(x = time3, y = Score.Value, group = time3, color = Lead_in, split_by= time3, scales="free_x")) + geom_line(aes(group = Patient.Number), colour = "gray80")+geom_point(size=3) + facet_wrap(~Lead_in+Score.Type)+theme_minimal()+scale_color_manual(values = point_colors, name = "Lead-in") + labs(title = "Paired Samples: Signatures Scores by Timepoint (PBMC)",y = "Signature.score",x = "Timepoint")

anno_pair_psl = compare_means(Score.Value ~ time3, group.by = c('Lead_in', 'Score.Type'),method = 'wilcox', data = pbmc_pp_psll, p.adjust.method = "holm") %>% mutate(y.position = 1, p.adj = format.pval(p.adj, digits = 2))

pbmc_pp_psll_1_wstats = pbmc_pp_psll_1+stat_pvalue_manual(anno_pair_psl, label='p.format', tip.length = 0.01,y.position=c(2,2,2), size=3)+scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# 6D CYT plot:
pfs_cyt2 = ggplot(tsll_tex2, aes(x = PFS.6month.status, y = CYT.score, fill=Lead_in)) +
geom_boxplot() +
labs(title = "Comparison of CYT Scores by Time and Lead_in (TIL)",
                y = "CYT.score",
                x = "Treatment") +facet_wrap(~time3, scales = "free_x", ncol=4) +theme_minimal()+geom_jitter(width = 0.2, alpha = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_manual(values = c("Combo" = "#DB2EE6", "Nivo" = "#5DCACC", "Rela" = "#737172"))

anno_pfs_cyt6 = compare_means(CYT.score ~ PFS.6month.status, group.by = 'time3',method = 'wilcox', data = tsll_tex2, p.adjust.method = "holm") %>% mutate(y.position = 1, p.adj = format.pval(p.adj, digits = 2))

anno_pfs_cyt6$group2 = c("PFS-free at 6 month","PFS-free at 6 month", "PFS-free at 6 month")
anno_pfs_cyt6$group1 = c("PD/Dead before 6 month","PD/Dead before 6 month", "PD/Dead before 6 month")

pfs_cyt2+stat_pvalue_manual(anno_pfs_cyt6, label='p.format', tip.length = 0.01,y.position=c(2,2,2), size=3)+scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# 6E IFNg plot:
pfs_ifng7 = ggplot(tsll_tex2, aes(x = PFS.6month.status, y = IFNg.score, fill=Lead_in)) +
geom_boxplot() +
labs(title = "Comparison of IFNg Scores by Time and Lead_in (TIL)",
     y = "IFNg.score",
     x = "Treatment") +facet_wrap(~Lead_in, scales = "free_x", ncol=4) +theme_minimal()+geom_jitter(width = 0.2, alpha = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_manual(values = c("Combo" = "#DB2EE6", "Nivo" = "#5DCACC", "Rela" = "#737172"))

anno_pfs_ifng8 = compare_means(IFNg.score ~ PFS.6month.status, group.by = 'Lead_in',method = 'wilcox', data = tsll_tex2, p.adjust.method = "holm") %>% mutate(y.position = 1, p.adj = format.pval(p.adj, digits = 2))

anno_ifng8$group2 = c("PFS-free at 6 month","PFS-free at 6 month", "PFS-free at 6 month")
anno_ifng8$group1 = c("PD/Dead before 6 month","PD/Dead before 6 month", "PD/Dead before 6 month")

pfs_ifng7+stat_pvalue_manual(anno_ifng8, label='p.format', tip.length = 0.01,y.position=c(1.3,1.3,1.3), size=3)+scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# 6F Tex plot:
pfs_tex2 = ggplot(tsll_tex2, aes(x = PFS.6month.status, y = Tex.score, fill=Lead_in)) +
geom_boxplot() +
labs(title = "Comparison of Tex Scores by Time and Lead_in (TIL)",
                y = "Tex.score",
                x = "Treatment") +facet_wrap(~time3, scales = "free_x", ncol=4) +theme_minimal()+geom_jitter(width = 0.2, alpha = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_manual(values = c("Combo" = "#DB2EE6", "Nivo" = "#5DCACC", "Rela" = "#737172"))

anno_pfs_tex6 = compare_means(Tex.score ~ PFS.6month.status, group.by = 'time3',method = 'wilcox', data = tsll_tex2, p.adjust.method = "holm") %>% mutate(y.position = 1, p.adj = format.pval(p.adj, digits = 2))

anno_pfs_tex6$group2 = c("PFS-free at 6 month","PFS-free at 6 month", "PFS-free at 6 month")
anno_pfs_tex6$group1 = c("PD/Dead before 6 month","PD/Dead before 6 month", "PD/Dead before 6 month")

pfs_tex2+stat_pvalue_manual(anno_pfs_tex6, label='p.format', tip.length = 0.01,y.position=c(2,2,2), size=3)+scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# Jennifer Walker, PhD 2026