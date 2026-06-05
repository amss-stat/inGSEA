# inGSEA: Improved Gene Set Enrichment Analysis Using a Weighted Integral Statistic

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-0.3.0-blue.svg)]()
[![Web Application](https://img.shields.io/badge/Web_App-Live-success.svg)](https://amss-stat.github.io/inGSEA/)

Gene Set Enrichment Analysis (GSEA) is one of the most popular methods for transcriptomic analysis, yet its statistical power is limited when biological pathways exhibit heterogeneous or non-concordant expression patterns. 

**inGSEA (integral-based GSEA)** is a powerful and robust extension of the classical GSEA framework. It introduces a novel enrichment score based on the Anderson-Darling (AD) weighted integral statistic and aggregates it with the classic Kolmogorov-Smirnov (KS) statistic via a Cauchy combination test. To overcome the computational bottlenecks of permutation testing, inGSEA utilizes a generalized gamma distribution to approximate the empirical null. 

The tool is accessible as user-friendly, and privacy-preserving web-based software:  
🌐 **[Launch inGSEA Web Application](https://amss-stat.github.io/inGSEA/)**

## Key Features

- **Anderson-Darling (AD) Enrichment Score:** Integrates squared deviations over the entire distribution with heavy tail-weighting. This substantially enhances detection power for complex signals, particularly sparse and bidirectional expression patterns missed by classic GSEA.
- **Robust Cauchy Combination:** Aggregates KS and AD statistics to provide robust sensitivity across diverse enrichment scenarios, balancing their complementary strengths.
- **Fast Generalized Gamma Approximation:** Reduces the computational burden and improves $p$-value resolution by fitting a generalized gamma distribution to the permutation null.
- **100% Client-Side Privacy:** No data is uploaded to any server. All statistical computations are performed locally within your web browser.
- **Interactive Visualization:** Offers fully interactive enrichment plots. Hovering over or clicking specific segments reveals individual genes and their rank metrics. Click the arrow to quickly jump to the pathway database.

## Quick Start

### 1. Launch the Tool
Visit **[https://amss-stat.github.io/inGSEA/](https://amss-stat.github.io/inGSEA/)** (Google Chrome is strongly recommended for optimal JavaScript engine performance).

*Tip: You can click the **"⚡ Load Demo Dataset"** button on the homepage to instantly try inGSEA with synthetic data.*

### 2. Input Data Formats
You only need two files to run the analysis. Simply drag and drop them into the designated zones.

**Expression Matrix (`.csv`, `.tsv`, or `.txt`)**: 
Rows represent genes, and columns represent samples. The first column must contain gene symbols.
| Gene | Case_1 | Case_2 | ... | Ctrl_1 | Ctrl_2 |
|:---:|:---:|:---:|:---:|:---:|:---:|
| TP53 | 2.5 | 1.2 | ... | 4.4 | 2.1 |
| EGFR | 4.2 | 5.1 | ... | 4.0 | 4.3 |
| MYC | 4.8 | 5.3 | ... | 3.5 | 4.9 |

**Gene Sets (`.gmt` or `.txt`)**: 
Standard MSigDB format. Each row represents a pathway: `Name`, `URL/Description`, followed by `Gene symbols`.
| Pathway Name | URL | Gene 1 | Gene 2 | Gene 3 | ... |
|:---:|:---:|:---:|:---:|:---:|:---:|
| HALLMARK_P53_PATHWAY | http://... | TP53 | BAX | CDKN1A | ... |
| HALLMARK_HYPOXIA | http://... | HIF1A | VEGFA | PGK1 | ... |

*Note: GMT files can be downloaded directly from the [MSigDB Collections](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp).*

### 3. Run the Analysis
1. Define the number of **Case samples** (assuming the first $N$ columns are cases, and the rest are controls).
2. Configure **Permutations** (e.g., 1000) and **Weight exponent $p$** (default: 1).
3. Select the **Statistical engine**:
   - *Parametric Approximation* (Recommended): Uses generalized gamma null fit for speed and high-resolution $p$-values.
   - *Permutation only*: Calculates empirical $p$-values strictly based on permutations.
4. Click **▶ Run inGSEA**.

## Output & Interpretation

The software provides dynamic, real-time results directly in the browser:

### Results Table
The output table is sortable and includes the following key metrics:
- **NES / NES-AD:** Normalized Enrichment Scores based on KS and AD statistics.
- **p<sub>KS</sub> / p<sub>AD</sub>:** Significance levels for the KS and AD statistics.
- **p<sub>Cauchy</sub>:** The aggregated omnibus $p$-value via the Cauchy combination test.
- **FDR<sub>KS</sub> / FDR<sub>AD</sub>:** False Discovery Rates based on the normalized enrichment scores.

### Visualizations
- **Interactive Enrichment Walk Plot:** Visualizes the running sum, peak/integral positions, and individual gene hits. Supports zooming, panning, and high-resolution PNG export.
- **NES Bar Chart:** Automatically generated when analyzing > 10 pathways to provide an overview of pathway activation/suppression directions.

## Citation

If you use inGSEA in your research, please cite:

> Zhang, Q., & Li, Q. (2026). *inGSEA: An Improved Method for Gene Set Enrichment Analysis Using a Weighted Integral Statistic.* bioRxiv. doi: [https://doi.org/10.64898/2026.06.02.729106](https://doi.org/10.64898/2026.06.02.729106)

**References for original methodologies:**
1. Subramanian A, *et al.* (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545–15550.
2. Liberzon A, *et al.* (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems*, 1(6), 417–425.

## License & Contact

The inGSEA source code is released under the [MIT License](LICENSE). 

*GMT files from MSigDB are subject to their own terms of use (CC BY 4.0 for academic use; commercial use requires a separate agreement with the Broad Institute of MIT and Harvard).*

For questions, feature requests, or bug reports, please [open an issue](https://github.com/amss-stat/inGSEA/issues) or contact [liqz@amss.ac.cn](mailto:liqz@amss.ac.cn).
