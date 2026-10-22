Methods for SI Figure 6A-F

R packages:
    R, 4.5.0
    ggplot2, 4.0.3
    ggrepel, 0.9.7
    Seurat, 5.5.1 
    SingleCellExperiment, (BioC 3.21) 1.32.0 
    edgeR, (BioC 32.1) 4.6.3
    tidyr, 1.3.2 
    tidyverse, 2.0.0
    dplyr, 1.2.0
    stringr, 1.6.0 

- Received Seurat object that contained the counts matrix and simple metadata
- Added more metadata based on clustering (Leiden, resolution = 0.25)
- Integration using Seurat with Harmony on Patient.name column
- Cell type annotation with Seurat FindAllMarkers and manual assignment 
    Removed B cell contaminants
    General categories: Macrophages, CD3+ T cells, CD4+ T cells 

- Gene filtering after integration and cell type cleaning, to remove uncharacterized and uninformative genes from the analysis

- Converted object metadata to dataframe, including these columns:
    Patient.name
    Lead-in
    PFS.6month.status
    treatment
    time3 (Baseline, Week04)
    Combined columns:
    time.treatment = time3+treatment
    cell_counts = patient.name+time3

- Calculated scores: APS, CYT, IFNg, and Tex
- Gene lists:
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

- Score calculation, see example in R script on Line 96

- Patients that had samples for both Baseline and Week04 timepoints:
(TIL) til_pat_pairs = c("UOP-002","UOP-003","UOP-004","UOP-005","UOP-006","UOP-007","UOP-008","UOP-009","UOP-014","UOP-015","UOP-018","UOP-023","UOP-024","UOP-026","UOP-028","UOP-029","UOP-030","UOP-031","UOP-034","UOP-036","UOP-038","UOP-042","UOP-043","UOP-044")
(PBMC) pbmc_pat_pairs = c("UOP-002","UOP-003","UOP-004","UOP-005","UOP-006","UOP-007","UOP-008","UOP-009","UOP-014","UOP-015","UOP-017","UOP-018","UOP-020","UOP-026","UOP-028","UOP-029","UOP-030","UOP-033","UOP-034","UOP-035","UOP-036","UOP-038","UOP-040","UOP-041","UOP-042")

- Subset for TIL and PBMC datasets that contained the paired samples above

- 6A baseline only TIL all scores:
# subset data table to have only Baseline timepoint present
    Use ggplot2 to geom_boxplot where the x-axis is the Lead-in (treatment arm), and y-axis is the Score.Value for the TIL data table
    Use facet_wrap on the "~Score.Type"
    Use stat_compare_means for statistics, method="anova"

- 6B baseline only PBMC all scores:
# subsets data table to have only Baseline timepoint present
    Use ggplot2 to geom_boxplot where the x-axis is the Lead-in (treatment arm), and y-axis is the Score.Value for the PBMC data table
    Use facet_wrap on the "~Score.Type"
    Use stat_compare_means for statistics, method="anova"

- 6C PBMC paired all scores:
    Use ggplot2 geom_point, and geom_line
    Map using x as the Time column, y as the Score.Value column
    Group and Split by Time
    Color by Lead-in
    Facet_wrap on "~Lead_in+Score.Type"
    Use compare_means to run Wilcox for statistical analysis
    Add the statistics to the plot using stat_pvalue_manual


- 6D CYT plot:
    Use ggplot2 geom_boxplot
    Plot the tsl/TIL-specific dataset where x is PFS status, and y is CYT score
    Use compare_means to determine statistical significance using Wilcox tests
    Use stat_pvalue_manual to add the stats to the plot

- 6E IFNg plot:
Use ggplot2 geom_boxplot
    Plot the tsl/TIL-specific dataset where x is PFS status, and y is IFNg score
    Use compare_means to determine statistical significance using Wilcox tests
    Use stat_pvalue_manual to add the stats to the plot

- 6F Tex plot:
Use ggplot2 geom_boxplot
    Plot the tsl/TIL-specific dataset where x is PFS status, and y is Tex score
    Use compare_means to determine statistical significance using Wilcox tests
    Use stat_pvalue_manual to add the stats to the plot




# Jennifer Walker, PhD 2026