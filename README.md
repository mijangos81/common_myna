# Desert Isolation and Urban Connectivity: Genomic Footprint of Myna Invasion in Oman

This repository contains code and documentation for analyzing the genomic data of common myna (*Acridotheres tristis*) populations in Oman, with a focus on desert isolation and urban connectivity.

## Project Structure

* **myna\_genomic\_analysis.R**: Main R script performing data import, quality control, diversity summaries, population structure (fastStructure), PCoA, spatial autocorrelation, and Ne estimation.
* **Report\_DImy24-9735\_SNP\_mapping\_1.csv**: Raw SNP genotype data in DArT format.
* **metadata\_myna.csv**: Sample metadata file containing individual IDs, sampling locations, and population labels.
* **GCA\_037013685.1\_AcTris\_vAus2.0\_genomic.fna**: Reference genome file used for BLAST mapping.
* **res\_fast\_myna.rds**: Saved results object from fastStructure analysis.

## Requirements

* R (≥ 4.0)
* [dartRverse](https://cran.r-project.org/package=dartRverse) package and its dependencies

## Installation

```bash
# In R:
install.packages("devtools")
devtools::install_github("green-striped-gecko/dartRverse")  # or CRAN
```

## Usage

1. Clone this repository:

   ```bash
   ```

git clone [https://github.com/username/myna-genomics-oman.git](https://github.com/username/myna-genomics-oman.git) cd myna-genomics-oman

```
2. Open and run the `myna_genomic_analysis.R` script in R or RStudio. It will:
   - Load and annotate the SNP dataset
   - Perform sequential quality-control filtering
   - Generate diversity and F-statistics reports
   - Compute genetic relatedness matrix (GRM)
   - Run population structure analyses (fastStructure)
   - Conduct Principal Coordinates Analysis (PCoA)
   - Evaluate spatial autocorrelation of allele frequencies
   - Estimate effective population size (Ne)

## Output

- Figures and text outputs will be generated in your R session or saved objects (e.g., `res_fast_myna.rds`).
- Modify plotting parameters or K-range in the main script as needed for publication figures.

## Citation

Please cite the following article when using these data or code:

> **Desert Isolation and Urban Connectivity: The Genomic Footprint of Myna Invasion in Oman**  
> Qais Al Rawahi et al. (2025)

## Contributing

Pull requests and issues are welcome. For major changes, please open an issue first to discuss proposed modifications.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

```
