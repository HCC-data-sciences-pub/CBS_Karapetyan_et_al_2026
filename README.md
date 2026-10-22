# CBS_Karapetyan_et_al_2026

Analysis code from the Hillman Cancer Center Cancer Bioinformatics Services (CBS) core, supporting:

> Karapetyan L, Cillo AR, Wang S, *et al.* **Lead-in PD-1 and LAG-3 blockade alone or in combination reveal distinct immune and clinical effects in metastatic melanoma.** (2026). [LINK TBD]

## Scripts

| Script | Output / Figure |
| --- | --- |
| `scripts/spatial_analysis_minimum.r` | Bivariate Ripley's K spatial analysis (spatialTIME) per sample; produces the tables used by `manuscript_figures.py` |
| `scripts/manuscript_figures.py` | Figure 4: CD8 vs FOXP3 Ripley's K boxplots by lead-in arm (Nivo / Rela / Combo) |
| `scripts/SI-figure-6-all-panels.r` | SI Figure 6A–F (scRNA-seq gene signature scores: APS, CYT, IFNg, Tex) |
| `scripts/SI-figure-6-all-panels.md` | Methods notes and package versions for SI Figure 6 |
