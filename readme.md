# m1A Score Dynamic Change and Its Function in Spinal Cord Injury

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-%3E%3D4.0-orange)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

Code repository for the manuscript investigating the dynamic changes of **m¹A (N1-methyladenosine) modification** and its functional role in **spinal cord injury (SCI)**. This study integrates bulk RNA-seq, single-cell RNA-seq (scRNA-seq), and single-nucleus RNA-seq (snRNA-seq) to characterize m¹A-related gene expression, m¹A scores across cell types, and their association with transcription factors, cellular trajectories, and biological pathways after SCI.

------

## Repository Structure

```
m1A_score_manuscript/
├── Bulk RNA-seq_process_scripts/     # Bulk RNA-seq analysis (human blood, mouse spinal cord)
│   ├── data_processing.R             # Data import, cleaning, normalization
│   ├── GEO_limma.R                   # GEO data download & limma differential expression
│   ├── GO_BP.R                       # GO Biological Process enrichment (compareCluster)
│   ├── ssGSEA.R                      # ssGSEA immune cell infiltration & correlation
│   ├── gsea_groups.R                 # GSEA pathway enrichment analysis
│   ├── for_function_wrap_plot.R      # Multi-gene boxplot wrapper (grouped comparisons)
│   ├── pca.R                         # Principal component analysis
│   ├── ggboxplot.R / ggline.R / ggviolin.R  # Publication-ready visualization helpers
│   ├── transform.R                   # Data transformation utilities
│   └── gsea和cor                      # Additional GSEA & correlation scripts
│
├── scRNA-seq_process_scripts/        # Single-cell RNA-seq (mouse SCI, GSE162610)
│   ├── 1_m1A score and expression.R  # Fig.1: UMAP, m1A regulator expression, UCell scoring
│   ├── 3_SCI_affected_cell_types.R   # Fig.2: SCENIC regulon activity, TF analysis in myeloid
│   ├── 3_m1A_score and tf correlation.R  # m1A score × transcription factor correlation
│   └── 4_microglia_monocle.R         # Fig.4: Monocle2 pseudotime trajectory of microglia
│
├── snRNA-seq_process_scripts/        # Single-nucleus RNA-seq (mouse SCI, GSE172167)
│   ├── 1_m1A score and expression.R  # Fig.1: snRNA-seq m1A scoring & visualization
│   ├── 2_neuro_cellratio_single.R    # Neuron cell ratio analysis
│   ├── 2_neuron_subcluster.R         # Neuron subclustering (SCTransform + integration)
│   ├── 2_neuron_subcluster_using.R   # Neuron subcluster annotation & GO/KEGG enrichment
│   ├── 4_metacell_enhance_cor.R      # SuperCell metacell construction & gene correlation
│   ├── 5_diff_function_m1a_cor.R     # Shared differential pathways (Upset/Venn) × m1A correlation
│   ├── 5_2_diff_metacell_m1a_cor.R   # Metacell-level m1A-pathway correlation
│   └── mysupercell_function.R        # Custom supercell network visualization
│
└── readme                            # Original brief description
```
------

## Analytical Workflow

### 1. Bulk RNA-seq Analysis

- **Data processing**: Import FPKM expression matrices, deduplication, log transformation, group assignment (control, acute, 3/6/12 months post-injury).
- **Differential expression**: `limma` for microarray; between-group comparisons at multiple time points.
- **Functional enrichment**: GO Biological Process (via `clusterProfiler` + MSigDB), GSEA (hallmark & GO gene sets).
- **Immune infiltration**: ssGSEA (`GSVA`) to estimate immune cell proportions and correlate with m¹A gene expression.

### 2. scRNA-seq Analysis (GSE162610)

- **Quality control & clustering**: Seurat pipeline (SCTransform, PCA, UMAP, clustering).
- **m¹A score calculation**:  
  - UCell module scoring with a custom **m¹A gene set** (`m1A_genesets.csv`: writers, erasers, readers).
  - Score comparison across time points (`Uninjured`, `1dpi`, `3dpi`, `7dpi`) and cell types.
- **m¹A regulator expression**: Heatmaps, dot plots, feature plots, and cell-proportion bar plots.
- **Transcription factor analysis**:  
  - pySCENIC → AUCell regulon activity scores in myeloid cells.
  - Spearman correlation between m¹A score and TF regulon activities.
- **Microglia trajectory**: Monocle2 pseudotime ordering with DDRTree, colored by pseudotime, cell subtype, condition, and m¹A score.

### 3. snRNA-seq Analysis (GSE172167)

- **Quality control & clustering**: Seurat SCTransform + integration across samples.
- **m¹A score calculation**: UCell scoring with the same m¹A gene set, compared across conditions (`Uninjured`, `1dpi`, `1wpi`, `3wpi`, `6wpi`).
- **Neuron subclustering**:  
  - Subset neurons → re-integration → subclustering at multiple resolutions.  
  - Marker-based annotation and GO/KEGG enrichment per subcluster.
- **Metacell analysis (SuperCell)**:  
  - Aggregate single nuclei into metacells (gamma=20).  
  - UMAP visualization, purity assessment, and gene-gene correlation at metacell level.
- **Differential pathway × m¹A correlation**:  
  - Identify shared up/down-regulated pathways across time points (Upset + Venn).  
  - Compute pathway scores (UCell) and correlate them with m¹A score within specific neuron subtypes (e.g., Slit2_IN).

------

## Key R Packages Used

| Category                | Packages                                                     |
| ----------------------- | ------------------------------------------------------------ |
| Single-cell / Nucleus   | [Seurat](https://satijalab.org/seurat/) (≥4.0), [UCell](https://github.com/carmonalab/UCell), [AUCell](http://bioconductor.org/packages/AUCell/), [SuperCell](https://github.com/GfellerLab/SuperCell) |
| Trajectory inference    | [monocle](http://cole-trapnell-lab.github.io/monocle-release/) (Monocle2) |
| Gene set / Pathway      | [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [GSVA](https://bioconductor.org/packages/GSVA/) (ssGSEA), [msigdbr](https://github.com/igordot/msigdbr), [enrichplot](https://bioconductor.org/packages/enrichplot/), [GseaVis](https://github.com/junjunlab/GseaVis) |
| Differential expression | [limma](https://bioconductor.org/packages/limma/), [MAST](https://bioconductor.org/packages/MAST/) |
| Transcription factor    | [pySCENIC](https://github.com/aertslab/pySCENIC) (Python) + [AUCell](http://bioconductor.org/packages/AUCell/) (R) |
| Visualization           | [ggplot2](https://ggplot2.tidyverse.org/), [ggpubr](https://rpkgs.datanovia.com/ggpubr/), [pheatmap](https://cran.r-project.org/package=pheatmap), [corrplot](https://cran.r-project.org/package=corrplot), [ggpointdensity](https://cran.r-project.org/package=ggpointdensity), [UpSetR](https://github.com/hms-dbmi/UpSetR), [ggvenn](https://cran.r-project.org/package=ggvenn), [ggrastr](https://cran.r-project.org/package=ggrastr) |
| Utilities               | [tidyverse](https://www.tidyverse.org/) (dplyr, tidyr, ggplot2, readr), [qs](https://github.com/traversc/qs), [cowplot](https://cran.r-project.org/package=cowplot), [export](https://cran.r-project.org/package=export) |

------

## Getting Started

### Prerequisites

- **R** ≥ 4.0 with Bioconductor
- **Python** (for pySCENIC; optional for the regulon analysis in scRNA-seq)
- Sufficient memory (≥ 64 GB recommended for sc/snRNA-seq datasets)

### Installation

```r
# Install core packages
install.packages(c("Seurat", "tidyverse", "ggpubr", "pheatmap", "cowplot",
                   "corrplot", "qs", "ggpointdensity", "UpSetR", "ggvenn",
                   "ggrastr", "export"))

# Bioconductor packages
BiocManager::install(c("limma", "MAST", "clusterProfiler", "org.Mm.eg.db",
                       "GSVA", "msigdbr", "enrichplot", "GEOquery",
                       "UCell", "AUCell", "SuperCell", "monocle"))
```

### Data Preparation

1. Download GEO datasets (GSE226238, GSE162610, GSE172167) to `data/` directories.
2. Place the **m¹A gene set** file (`m1A_genesets.csv`) at the expected paths (project root or script subdirectory).
3. Adjust working directories (`setwd()`) and file paths in each script to match your environment.

### Running the Analysis

Scripts are numbered and named to correspond to manuscript figures. Run them in order within each subdirectory:

```bash
# Example: scRNA-seq analysis
cd scRNA-seq_process_scripts/
Rscript "1_m1A score and expression.R"    # Figure 1
Rscript "3_SCI_affected_cell_types.R"     # Figure 2
Rscript "3_m1A_score and tf correlation.R" # Figure 2 (continued)
Rscript "4_microglia_monocle.R"           # Figure 4
```

------

## Citation

If you use this code in your research, please cite the corresponding manuscript and the original data sources.

------

## License

This project is made available for academic use. Please refer to the repository license file for details.

------

## Contact

For questions or collaborations, please open an issue on this repository or contact the corresponding author.
