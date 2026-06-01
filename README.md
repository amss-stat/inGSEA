https://amss-stat.github.io/inGSEA/

# inGSEA — Improved Gene Set Enrichment Analysis Using a Weighted Integral Statistic

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.3.0-green.svg)]()

A web-based implementation of the inGSEA framework, which extends classical
Gene Set Enrichment Analysis with an Anderson-Darling enrichment score and
Cauchy combination test for improved statistical power and robustness.

---

## Overview

Standard GSEA has limited power when pathways exhibit heterogeneous or
non-concordant expression patterns. inGSEA addresses this by:

- Adopting an **Anderson-Darling (AD) enrichment score** that
  enhances detection of complex signals, paticularly sparse and bidirectional ones
- Combining AD and KS tests via a **Cauchy combination test** for
  robustness across diverse expression patterns
- Approximating permutation null distributions with a **generalised gamma
  distribution**, substantially reducing computational cost

All computation runs entirely in the browser — no data is uploaded to any
server.

---

## Live Demo

[**Launch inGSEA →**](https://amss-stat.github.io/inGSEA)

A built-in demo dataset (200 genes × 20 samples, 3 synthetic pathways) is
available via the **⚡ Load demo dataset** button.

---

## Usage

### Input files

| File | Format | Description |
|------|--------|-------------|
| Expression matrix | CSV or TSV | Rows = genes, columns = samples |
| Gene sets | GMT | MSigDB format (name · URL · genes) |

GMT files can be downloaded from
[MSigDB Collections](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
(free registration required).

### Steps

1. Drop your expression matrix into the **Expression Matrix** zone
2. Drop a GMT file into the **Gene Sets** zone
3. Set the number of **case samples** (first *N* columns)
4. Adjust permutations, weight exponent, and pathway selection as needed
5. Choose a **statistical engine**:
   - *Parametric Approximation* — generalised gamma null fit (faster)
   - *Permutation only* — empirical p-values only
6. Click **▶ Run inGSEA**

### Output

- Sortable results table with NES, NES-AD, p-values, and FDR values
- Interactive enrichment plot (zoom, pan, hover)
- NES bar chart overview (when > 10 pathways)
- CSV export and PNG export of enrichment curves

---

## Methods Summary

| Component | Description |
|-----------|-------------|
| Ranking metric | Welch signal-to-noise ratio |
| KS enrichment score | Weighted running sum (weight exponent *p*, default 1) |
| AD enrichment score | Integral of squared running sum over rank fractions |
| Null approximation | Generalised gamma MLE fit to permutation null |
| p-value combination | Cauchy combination test |
| Multiple testing | Benjamini-Hochberg FDR |

---

## Browser Compatibility

Chrome is strongly recommended for best performance.

---

## References

1. Zhang, Q., & Li, Q. (2026). A rotated multivariate linear mixed model for dual large-scale genome-wide association study. *bioariv*, doi: 
2. Mootha VK, *et al.* (2003). PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. *Nature Genetics*, 34(3), 267–273.
3. Subramanian A, *et al.* (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545–15550.
4. Liberzon A, Subramanian A, Pinchback R, Thorvaldsdóttir H, Tamayo P & Mesirov JP (2011). Molecular Signatures Database (MSigDB) 3.0. *Bioinformatics*, 27(12), 1739–1740.
5. Liberzon A, *et al.* (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems*, 1(6), 417–425.

---

## License

The inGSEA source code is released under the
[MIT License](LICENSE).

GMT files from MSigDB are subject to their own terms of use (CC BY 4.0 for
academic use; commercial use requires a separate agreement with the Broad
Institute of MIT and Harvard).

---

## Contact

For questions or bug reports, please open an issue on
[GitHub](https://github.com/amss-stat/inGSEA/issues).
