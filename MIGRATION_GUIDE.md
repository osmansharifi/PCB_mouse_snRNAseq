# Script Migration Guide

## Quick Reference: Where Scripts Go

### Preprocessing Scripts → scripts/01_preprocessing/
- soupx/PEBBLES_SoupX.R ✅ **Already done**
- soupx/PEBBLES_filter_02.R
- Other QC/filtering scripts

### Core Analysis → scripts/04_core_analysis/
- 02_DEG_analysis/DEG_Limma_analysis_postnatal.R
- 02_DEG_analysis/limmavoom_pcbvsveh_het_DEG_analysis.R
- 06_hdWGCNA/hdWGCNA.R

### Comparative/Specialized → analysis/
- 03_Kari_DEG_analyses/ (entire folder)
- 04_DEG_methods_comparison/ (entire folder)
- 07_mosiacism/ (entire folder)

### Human Translation → scripts/06_human_translation/
- 08_human_hdWGCNA/human_hdWGCNA.R
- 09_human_DEG_analysis/DEanalysis_filtering_01.R
- 09_human_DEG_analysis/limmaVoomCC/DEG_analysis.R

---

## How to Migrate Each Script

### Step 1: Copy Script to New Location
```bash
cp [OLD_PATH]/script.R scripts/[CATEGORY]/script_name.R
```

### Step 2: Update the Header
Use `scripts/SCRIPT_TEMPLATE.R` as a guide. Add:
- Clear script name
- Purpose statement
- Input/output files
- Dependencies list
- Author and date

### Step 3: Update Paths
Replace hardcoded paths with relative paths:
```r
# OLD (hardcoded):
setwd('/share/lasallelab/Osman/2021_PEBBLES_Cortex/...')
load('/share/lasallelab/Osman/...')

# NEW (relative):
seurat_obj <- load_seurat("data/processed/seurat.rds")
saveRDS(seurat_obj, "data/processed/output.rds")
```

### Step 4: Add Logging
Replace `print()` and `cat()` with logging:
```r
# OLD:
print("Processing started")

# NEW:
log_message("Processing started", level = "INFO")
```

### Step 5: Remove `setwd()`
Delete any `setwd()` lines. Use relative paths instead.

### Step 6: Add Error Handling
Wrap risky operations:
```r
tryCatch({
  # Your code here
}, error = function(e) {
  log_message(sprintf("ERROR: %s", e$message), level = "ERROR")
})
```

### Step 7: Commit
```bash
git add scripts/[CATEGORY]/
git commit -m "refactor: standardize [script_name]

- Add professional header
- Use relative paths
- Add logging statements
- Use utility functions
- Add error handling"
```

---

## Priority Order

**High Priority** (do these first):
1. DEG analysis scripts (02_DEG_analysis/)
2. hdWGCNA script (06_hdWGCNA/)

**Medium Priority** (do these next):
3. Human analysis scripts (08_human_hdWGCNA/, 09_human_DEG_analysis/)
4. SoupX preprocessing scripts (soupx/)

**Lower Priority** (do these last):
5. Specialized analyses (03_Kari, 04_methods_comparison, 07_mosaicism)

---

## Quick Template for Migration
```bash
# 1. Copy script
cp [OLD_PATH]/script.R scripts/[CATEGORY]/script_name.R

# 2. Edit file with your editor to:
#    - Add header from SCRIPT_TEMPLATE.R
#    - Replace hardcoded paths with relative paths
#    - Replace print/cat with log_message()
#    - Remove setwd() lines

# 3. Commit
git add scripts/[CATEGORY]/script_name.R
git commit -m "refactor: standardize [script_name] script"

# 4. Move to next script
```

---

## Bulk Operations

### Move all preprocessing scripts:
```bash
cp soupx/*.R scripts/01_preprocessing/
```

### Move all comparative analyses:
```bash
cp -r 03_Kari_DEG_analyses analysis/
cp -r 04_DEG_methods_comparison analysis/
cp -r 07_mosiacism analysis/
```

### Move all human scripts:
```bash
cp 08_human_hdWGCNA/*.R scripts/06_human_translation/
cp 09_human_DEG_analysis/*.R scripts/06_human_translation/
```

---

## After Migration

1. Update headers in each script
2. Replace paths with relative ones
3. Add logging statements
4. Test scripts run without errors
5. Commit all changes

---

**Estimated Time:**
- Quick copy & commit: 2-3 hours
- Full standardization: 6-8 hours
- Your choice!
