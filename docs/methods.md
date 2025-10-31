# Methods and Analysis Parameters

## Quality Control

### Filtering Thresholds
- **Minimum UMI count**: 500
  - Rationale: Removes very low-quality captures
  
- **Minimum genes detected**: 200
  - Rationale: Ensures adequate transcriptome diversity
  
- **Maximum mitochondrial RNA**: 5%
  - Rationale: Indicates cell viability (high MT% = dying cells)
  
- **Maximum UMI count**: 50,000
  - Rationale: Removes potential doublets

## Normalization

### SCTransform (Seurat)
- **Method**: Variance-stabilizing transformation
- **Scale factor**: 10,000
- **Number of variable features**: 2,000

## Clustering & Annotation

### Parameters
- **k-nearest neighbors**: 20
- **Resolution**: 0.8
- **Algorithm**: Leiden (default in Seurat v4)

## Differential Expression Analysis

### Method
- **Test**: Wilcoxon rank-sum test
- **Minimum log2 fold change**: 0.25
- **Minimum percent expressing**: 10%
- **P-value threshold**: 0.05 (adjusted)
- **Adjustment method**: Bonferroni

## Reproducibility Notes

- **Random seed**: 42 (set in all scripts)
- **Software versions**: See environment.yml
- **Total runtime**: Varies by dataset size

---

**Methods version**: 1.0.0  
**Last updated**: 2024-10-31
