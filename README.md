https://amss-stat.github.io/iGSEA-DEMO/

# iGSEA — Improved Gene Set Enrichment Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.2.4-green.svg)]()

A web-based implementation of the iGSEA framework, which extends classical
Gene Set Enrichment Analysis with an Anderson-Darling enrichment score and
Cauchy combination test for improved statistical power and robustness.

> **Note:** The accompanying manuscript is currently under review.
> Citation information will be added upon publication.

---

## Overview

Standard GSEA has limited power when pathways exhibit heterogeneous or
non-concordant expression patterns. iGSEA addresses this by:

- Integrating an **Anderson-Darling (AD) based enrichment score** that
  enhances detection of sparse and bidirectional signals
- Combining KS and AD statistics via a **Cauchy combination test** for
  robustness across diverse expression patterns
- Approximating permutation null distributions with a **generalised gamma
  distribution**, substantially reducing computational cost
- Providing strict **false positive rate control** alongside superior
  statistical power (validated by extensive simulation)

All computation runs entirely in the browser — no data is uploaded to any
server.

---

## Live Demo

[**Launch iGSEA →**](https://amss-stat.github.io/iGSEA-DEMO)

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
6. Click **▶ Run iGSEA**

### Output

- Sortable results table with NES, NES-AD, p-values, and BH-adjusted FDR
- Interactive enrichment score walk plot (zoom, pan, hover)
- NES bar chart overview (when > 10 pathways)
- CSV export and PNG export of enrichment curves

---

## Methods Summary

| Component | Description |
|-----------|-------------|
| Ranking metric | Welch signal-to-noise ratio |
| KS enrichment score | Weighted running sum (weight exponent *p*, default 1) |
| AD statistic | Integral of squared running sum over rank fractions |
| Null approximation | Generalised gamma MLE fit to permutation null |
| p-value combination | Cauchy combination test (Liu & Xie 2020) |
| Multiple testing | Benjamini-Hochberg FDR |

---

## Browser Compatibility

Chrome is strongly recommended for best performance.
Firefox and Edge are supported. Safari may be slower for large datasets.

---

## References

1. *[Our paper — citation to be added upon publication]*
2. Mootha VK, *et al.* (2003). PGC-1α-responsive genes involved in oxidative
   phosphorylation are coordinately downregulated in human diabetes.
   *Nature Genetics*, 34(3), 267–273.
   [doi:10.1038/ng1180](https://doi.org/10.1038/ng1180)
3. Subramanian A, *et al.* (2005). Gene set enrichment analysis: a
   knowledge-based approach for interpreting genome-wide expression profiles.
   *PNAS*, 102(43), 15545–15550.
   [doi:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
4. Liberzon A, Subramanian A, Pinchback R, Thorvaldsdóttir H, Tamayo P &
   Mesirov JP (2011). Molecular Signatures Database (MSigDB) 3.0.
   *Bioinformatics*, 27(12), 1739–1740.
   [doi:10.1093/bioinformatics/btr260](https://doi.org/10.1093/bioinformatics/btr260)
5. Liberzon A, *et al.* (2015). The Molecular Signatures Database (MSigDB)
   hallmark gene set collection. *Cell Systems*, 1(6), 417–425.
   [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)

---

## License

The iGSEA source code is released under the
[MIT License](LICENSE).

GMT files from MSigDB are subject to their own terms of use (CC BY 4.0 for
academic use; commercial use requires a separate agreement with the Broad
Institute of MIT and Harvard).

---

## Contact

For questions or bug reports, please open an issue on
[GitHub](https://github.com/amss-stat/iGSEA-DEMO/issues).
