# PCB Exposure and Gene Expression in Mouse Brain: snRNA-seq Analysis

## Overview

Single-nuclei RNA-sequencing (snRNA-seq 5') analysis investigating polychlorinated biphenyl (PCB) effects on gene expression in the RTT mouse and human prefrontal cortex.

### Associated Publication

**Female cortical cellular mosaicism underlies shared MeCP2 and PCB impacted gene pathways**

PubMed: https://pubmed.ncbi.nlm.nih.gov/40501678/

## Quick Start

### Installation
```bash
# Clone repository
git clone https://github.com/osmansharifi/PCB_mouse_snRNAseq.git
cd PCB_mouse_snRNAseq

# Create conda environment
conda env create -f environment.yml
conda activate pcb_snrnaseq
```

### Basic Usage
```bash
# Run preprocessing
Rscript scripts/01_preprocessing/qc_and_normalization.R

# Run core analysis
Rscript scripts/04_core_analysis/deg_analysis.R
```

## Repository Structure
```
scripts/
├── 01_preprocessing/          # QC, normalization, ambient RNA correction
├── 02_celltype_annotation/    # Clustering, cell type identification
├── 03_deconvolution/          # Deconvolution methods
├── 04_core_analysis/          # Main DEG & co-expression analysis
│   ├── limma_voom_analysis/         # Limma-Voom DE analysis
│   ├── genotype_comparisons/        # Genotype-specific DE comparisons
│   ├── mosaicism_analysis/          # Cellular mosaicism DE analysis
│   └── summary_analysis/            # Summary and proportions analysis
├── 05_comparative_analyses/   # DE method comparisons
├── 06_human_translation/      # Human data analysis
└── utils/                     # Reusable utility functions

results/
├── figures/                   # Publication-ready figures
├── tables/
│   └── deg_results/          # Differential expression analysis results
│       ├── by_celltype/      # DE results organized by cell type
│       ├── by_genotype/      # DE results organized by genotype
│       ├── mosaicism/        # Cellular mosaicism analysis results
│       ├── method_comparisons/ # DE method comparison results
│       └── summary/          # Summary tables and figures
└── objects/                   # Saved R/Python objects

data/
├── raw/                       # Raw sequencing data
└── processed/                 # Processed count matrices

docs/
├── methods.md                 # Detailed methods & parameters
└── analysis_workflow.md       # Pipeline explanation
```

## Key Results

- **Cell types identified**: Multiple neuronal and glial populations
- **DEGs by treatment**: Identified PCB-responsive genes
- **Co-expression modules**: Functional gene networks

## Dependencies

See `environment.yml` for complete list. Key packages:

- Seurat (v4.3+)
- tidyverse
- ggplot2

## Code Standards

All code follows R standards:
- snake_case naming conventions
- Comprehensive script headers
- Reusable functions in `scripts/utils/`
- Error handling and logging
- Relative paths for portability

## Citation

If you use this code, please cite:
```
Sharifi, O., et al. (2024). "PCB Exposure and Gene Expression in Mouse Brain".
GitHub: https://github.com/osmansharifi/PCB_mouse_snRNAseq
```

## License

MIT License - See LICENSE file for details

## Author

**Osman Sharifi**

---

**Last updated**: 2024-10-31
