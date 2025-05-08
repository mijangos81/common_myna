# Desert Isolation and Urban Connectivity: Genomic Footprint of Myna Invasion in Oman

> **Authors:** Qais Al Rawahi¹, Abdullahi Aliyu¹², Mazin Al-Abidi³, Kareem Khalil⁵, Masooma Al-Lawati¹, Tehani A Al-Jidailiy¹, Basel S Al-Maskari, Ahmed M Al-Shakili³, Jihad A Al-Toubi¹, Maisaa S Al-Saadi¹, Adamu Abdul Abubakar¹, Jose L Mijangos⁴⁵

¹ Department of Veterinary Medicine, College of Applied and Health Sciences, A’Sharqiyah University, Ibra, Sultanate of Oman
² Department of Veterinary Pathology, Usmanu Danfodiyo University, Sokoto State, Nigeria
³ Directorate General of Natural Conservation, Environment Authority, Ibra, Sultanate of Oman
⁴ Institute for Applied Ecology, University of Canberra, Australia
⁵ Diversity Arrays Technology Pty Ltd, University of Canberra, Australian Capital Territory, Australia

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Outputs](#outputs)
7. [Citation](#citation)
8. [Contributing](#contributing)
9. [License](#license)

---

## Project Overview

This repository contains R code and supporting files for analyzing SNP data from invasive common myna (*Acridotheres tristis*) populations in Oman. The study explores genomic signatures of isolation in desert landscapes and connectivity within urban environments.

## Repository Structure

```
├── myna_genomic_analysis.R    # Main analysis script
├── Report_DImy24-9735_SNP_mapping_1.csv    # Raw DArT SNP data
├── metadata_myna.csv          # Sample metadata (IDs, locations, populations)
├── GCA_037013685.1_AcTris_vAus2.0_genomic.fna  # Reference genome for BLAST
├── res_fast_myna.rds          # fastStructure results (saved RDS)
├── README.md                  # Project overview and instructions
├── license.txt                # MIT License text
└── LICENSE                    # MIT License file (renamed copy)
```

## Requirements

* **R** version ≥ 4.0
* **dartRverse** R package (and dependencies)

## Installation

1. Install dependencies in R:

   ```r
   install.packages("devtools")
   devtools::install_github("green-striped-gecko/dartRverse")  # or via CRAN
   ```
2. Clone this repository:

   ```bash
   git clone https://github.com/username/myna-genomics-oman.git
   cd myna-genomics-oman
   ```

## Usage

1. Launch R or RStudio.
2. Source the analysis script:

   ```r
   source("myna_genomic_analysis.R")
   ```
3. The script performs the following steps:

   * Data import and annotation
   * Sequential quality-control filtering
   * Diversity and F-statistics summaries
   * Genetic relatedness matrix (GRM)
   * Population structure analysis (fastStructure)
   * Principal Coordinates Analysis (PCoA)
   * Spatial autocorrelation
   * Effective population size estimation (LD-Ne)

## Outputs

* Figures and tables will be generated in your R session.
* `res_fast_myna.rds` contains the fastStructure result object.
* Customize plotting parameters or analysis ranges directly in `myna_genomic_analysis.R` for publication-ready figures.

## Citation

If you use this code or data, please cite:

> Al Rawahi, Q., Aliyu, A., Al-Abidi, M., Khalil, K., Al-Lawatiy, M., Al-Jidailiy, T. A., ... & Mijangos, J. L. (2025). Desert Isolation and Urban Connectivity: The Genomic Footprint of Myna Invasion in Oman.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for bug fixes, feature requests, or enhancements.

## License

This project is licensed under the MIT License. See [LICENSE.txt](LICENSE.txt) for details.
