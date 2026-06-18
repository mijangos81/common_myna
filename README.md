# The genomic footprint of myna invasion in Oman showcases desert isolation and urban connectivity

R code and data to reproduce the population-genomic analyses in Al Rawahi et al. (2026), *Scientific Reports* ([10.1038/s41598-026-57674-0](https://doi.org/10.1038/s41598-026-57674-0)).

**Authors:** Qais Al Rawahi¹, Abdullahi Aliyu¹,², Mazen M. Al-Obaidi³, Karim Khalil⁴,⁵, Masooma Al-Lawati¹, Tahani A. Al-Jidaili¹, Basel S. Al-Maskari¹, Ahmed M. Al-Shakili⁶, Jihad A. Al-Toubi¹, Maisaa S. Al-Saadi¹, Adamu Abdul Abubakar¹ & Jose L. Mijangos⁷,⁸

1. Department of Veterinary Medicine, College of Applied and Health Sciences, A'Sharqiyah University, P.O Box 42, 400 Ibra, Sultanate of Oman.
2. Department of Veterinary Pathology, Faculty of Veterinary Medicine, City Campus Complex, Usmanu Danfodiyo University, Sokoto 840212, Sokoto State, Nigeria.
3. Science Department, University of Technology & Applied Sciences, Rustaq, Sultanate of Oman.
4. Department of Animal and Veterinary Sciences, College of Agricultural and Marine Sciences, Sultan Qaboos University, P.O. Box 34, 123 Muscat, Sultanate of Oman.
5. Anatomy and Embryology Department, Faculty of Veterinary Medicine, Cairo University, Giza, Egypt.
6. DG of Natural Conservation, Environment Authority, P.O Box 42, 400 Ibra, Sultanate of Oman.
7. Faculty of Science, Sydney School of Veterinary Science, The University of Sydney, Sydney, NSW 2006, Australia.
8. Diversity Arrays Technology Pty Ltd, University of Canberra, Bruce, Australia.

Corresponding author: Jose L. Mijangos (jose.mijangosaraujo@sydney.edu.au).

---

## Contents

1. [Overview](#overview)
2. [Repository structure](#repository-structure)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Analysis steps](#analysis-steps)
7. [Citation](#citation)
8. [Questions](#questions)
9. [License](#license)

---

## Overview

This repository holds the R code and DArT SNP data to reproduce the analyses in the paper. We genotyped 311 common myna (*Acridotheres tristis*) from four governorates of Oman at 17,460 SNPs. The data resolve two genetic clusters: Dhofar, isolated by desert, and a connected northern group (Muscat, North Batinah and South Sharqiyah) with internal substructure. North Batinah carries a larger effective population size than the other governorates. These patterns point to differentiated management, with localised control feasible in Dhofar and an integrated approach needed across the connected northern governorates.

## Repository structure

```
├── myna_genomic_analysis.R                   # Main analysis script
├── Report_DImy24-9735_SNP_mapping_1.csv.zip  # DArT SNP report (myna; script input)
├── Report_DCrw24-9702_SNP_mapping_1.csv.zip  # DArT SNP report (additional dataset)
├── metadata_myna.csv                         # Sample metadata (IDs, coordinates, populations, morphological traits)
├── sequence_report.csv                       # Assembly chromosome lookup (GCA_037013685.1; used by the GWAS)
├── res_fast_myna.rds                         # Saved fastStructure result
├── common_myna.Rproj                         # RStudio project file
├── LICENSE.txt                               # MIT License
└── README.md
```

The reference genome is not included in the repository (see Requirements).

## Requirements

- R (≥ 4.0)
- The **dartRverse** suite (see Installation)
- For the GWAS: the **GAPIT**, **CMplot** and **colorspace** R packages
- A local NCBI BLAST+ installation, used by `gl.blast` to map loci to the reference genome
- The common myna reference genome FASTA, `GCA_037013685.1_AcTris_vAus2.0_genomic.fna`, downloaded from NCBI (assembly GCA_037013685.1) and placed in the project folder

## Installation

Install the dartRverse meta-package from CRAN, then the sub-packages used here:

```r
install.packages("dartRverse")
library(dartRverse)
dartRverse_install("dartR.base",    rep = "CRAN")  # data import, filtering, BLAST, diversity, PCoA
dartRverse_install("dartR.popgen",  rep = "CRAN")  # fastStructure, LD-Ne
dartRverse_install("dartR.spatial", rep = "CRAN")  # spatial autocorrelation

# GWAS packages: CMplot and colorspace from CRAN, GAPIT from GitHub
install.packages(c("CMplot", "colorspace", "devtools"))
devtools::install_github("jiabowang/GAPIT")
```

Run `dartRverse_install("all")` to list every sub-package, or see the [dartRverse repository](https://github.com/green-striped-gecko/dartRverse) for development versions.

Clone this repository and unzip the SNP report before running the script:

```bash
git clone https://github.com/mijangos81/common_myna.git
cd common_myna
unzip Report_DImy24-9735_SNP_mapping_1.csv.zip
```

## Usage

Open `common_myna.Rproj` in RStudio (or set the working directory to the project folder), then run:

```r
source("myna_genomic_analysis.R")
```

The script reads `Report_DImy24-9735_SNP_mapping_1.csv`, `metadata_myna.csv` and `sequence_report.csv`. It writes summary tables (`res_het.csv`, `res_fst.csv`), figures (`het.pdf`, `grm.pdf`, `faststruc.pdf`, `pca_1_2.pdf`, `autocorr.pdf`, `ibd.pdf`), the fastStructure result (`res_fast_myna.rds`), and the per-trait GAPIT and Manhattan-plot output to the working directory.

## Analysis steps

`myna_genomic_analysis.R` runs the following:

1. Read the DArT SNP report and sample metadata (`gl.read.dart`).
2. Map loci to the reference genome by BLAST and record chromosome and position (`gl.blast`).
3. Remove sex-linked (Z/W) loci (`gl.drop.sexlinked`).
4. Filter individuals by call rate (remove those with < 60% of loci scored).
5. Filter loci by call rate (< 90%), read depth and reproducibility.
6. Summarise diversity and F-statistics (`gl.report.heterozygosity`, `gl.report.fstat`).
7. Compute the genomic relatedness matrix (`gl.grm`).
8. Estimate population structure with fastSTRUCTURE, K = 1–6, 10 replicates (`gl.run.faststructure`).
9. Run Principal Coordinates Analysis (`gl.pcoa`).
10. Test spatial autocorrelation of allele frequencies (`gl.spatial.autoCorr`).
11. Estimate isolation by distance (`gl.ibd`).
12. Estimate effective population size from linkage disequilibrium (`gl.LDNe`).
13. Run genome-wide association for the eight morphological traits (GAPIT, BLINK model, two principal components) and draw Manhattan plots (`CMplot`).

## Citation

If you use this code or data, please cite:

> Al Rawahi, Q., Aliyu, A., Al-Obaidi, M. M., Khalil, K., Al-Lawati, M., Al-Jidaili, T. A., Al-Maskari, B. S., Al-Shakili, A. M., Al-Toubi, J. A., Al-Saadi, M. S., Abubakar, A. A., & Mijangos, J. L. (2026). The genomic footprint of myna invasion in Oman showcases desert isolation and urban connectivity. *Scientific Reports*. https://doi.org/10.1038/s41598-026-57674-0

## Questions

Open an issue for questions about the code, or contact the corresponding author for questions about the study.

## License

Released under the MIT License (see [LICENSE.txt](LICENSE.txt)).
