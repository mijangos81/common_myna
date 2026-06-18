# The genomic footprint of myna invasion in Oman showcases desert isolation
# and urban connectivity
#
# Al Rawahi et al. (2026), Scientific Reports. https://doi.org/10.1038/s41598-026-57674-0
#
# Qais Al Rawahi(1), Abdullahi Aliyu(1,2), Mazen M. Al-Obaidi(3), Karim Khalil(4,5),
# Masooma Al-Lawati(1), Tahani A. Al-Jidaili(1), Basel S. Al-Maskari(1),
# Ahmed M. Al-Shakili(6), Jihad A. Al-Toubi(1), Maisaa S. Al-Saadi(1),
# Adamu Abdul Abubakar(1) and Jose L. Mijangos(7,8).
#
# 1 Department of Veterinary Medicine, College of Applied and Health Sciences,
#   A'Sharqiyah University, P.O Box 42, 400 Ibra, Sultanate of Oman.
# 2 Department of Veterinary Pathology, Faculty of Veterinary Medicine,
#   City Campus Complex, Usmanu Danfodiyo University, Sokoto 840212,
#   Sokoto State, Nigeria.
# 3 Science Department, University of Technology & Applied Sciences,
#   Rustaq, Sultanate of Oman.
# 4 Department of Animal and Veterinary Sciences, College of Agricultural and
#   Marine Sciences, Sultan Qaboos University, P.O. Box 34, 123 Muscat,
#   Sultanate of Oman.
# 5 Anatomy and Embryology Department, Faculty of Veterinary Medicine,
#   Cairo University, Giza, Egypt.
# 6 DG of Natural Conservation, Environment Authority, P.O Box 42, 400 Ibra,
#   Sultanate of Oman.
# 7 Faculty of Science, Sydney School of Veterinary Science, The University of
#   Sydney, Sydney, NSW 2006, Australia.
# 8 Diversity Arrays Technology Pty Ltd, University of Canberra, Bruce, Australia.
#
# Corresponding author: Jose L. Mijangos (jose.mijangosaraujo@sydney.edu.au).

# ──────────────────────────────────────────────────────────────────────────────
# 1. Load required packages
# ──────────────────────────────────────────────────────────────────────────────
library(dartRverse)
library(GAPIT)
library(CMplot)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Read SNP dataset and metadata
# ──────────────────────────────────────────────────────────────────────────────
gl_data <- gl.read.dart(
  filename     = "Report_DImy24-9735_SNP_mapping_1.csv",
  ind.metafile = "metadata_myna.csv"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Annotate loci with trimmed sequence information
# ──────────────────────────────────────────────────────────────────────────────
gl_data$other$loc.metrics$TrimmedSequence <-
  gl_data$other$loc.metrics$TrimmedSequenceSnp

# ──────────────────────────────────────────────────────────────────────────────
# 4. Map loci to reference genome and extract chromosome/position
# ──────────────────────────────────────────────────────────────────────────────
gl_data <- gl.blast(
  gl_data,
  ref_genome = "GCA_037013685.1_AcTris_vAus2.0_genomic.fna"
)
# Store BLAST-derived scaffold and start positions
gl_data$chromosome <- as.factor(gl_data$other$loc.metrics$sacc)
gl_data$position   <- gl_data$other$loc.metrics$sstart

# ──────────────────────────────────────────────────────────────────────────────
# 5. Remove sex-linked (Z/W) loci
#    gl.report.sexlinked expects the sex column coded as "F"/"M"; the metadata
#    stores "Female"/"Male", so recode to the first initial (others to NA).
# ──────────────────────────────────────────────────────────────────────────────
sx <- toupper(substr(as.character(gl_data@other$ind.metrics$sex), 1, 1))
sx[!sx %in% c("F", "M")] <- NA
gl_data@other$ind.metrics$sex <- sx

r1      <- gl.report.sexlinked(gl_data, system = "zw")
gl_data <- gl.drop.sexlinked(gl_data, system = "zw")

# ──────────────────────────────────────────────────────────────────────────────
# 6. Quality control – call-rate filtering
#    Remove individuals with < 60% loci scored
# ──────────────────────────────────────────────────────────────────────────────
gl_qc1 <- gl.filter.callrate(gl_data, threshold = 0.6, method = "ind")
gl.report.callrate(gl_qc1, method = "ind")

# ──────────────────────────────────────────────────────────────────────────────
# 7. Quality control – locus-level filtering
#    Remove loci with < 90% individuals genotyped, low read depth and low
#    reproducibility
# ──────────────────────────────────────────────────────────────────────────────
gl_qc2 <- gl.filter.callrate(gl_qc1, threshold = 0.9, method = "loc")
gl_qc2 <- gl.filter.rdepth(gl_qc2)
gl_qc2 <- gl.filter.reproducibility(gl_qc2)
gl.report.callrate(gl_qc2, method = "loc")

# ──────────────────────────────────────────────────────────────────────────────
# 8. Genetic diversity and F-statistics
# ──────────────────────────────────────────────────────────────────────────────
res_het <- gl.report.heterozygosity(
  x               = gl_qc2,
  plot.colors.pop = c("#3283FE", "#FEAF16", "#B00068", "#1CFFCE")
)
write.csv(res_het, "res_het.csv")
ggsave("het.pdf", width = 6, height = 6, units = "in")

res_fst <- gl.report.fstat(gl_qc2)
write.csv(res_fst$Stat_matrices$Fstp, "res_fst.csv")

# ──────────────────────────────────────────────────────────────────────────────
# 9. Genetic relatedness matrix (GRM)
# ──────────────────────────────────────────────────────────────────────────────
gl.grm(
  x                  = gl_qc2,
  palette_convergent = gl.colors("con"),
  palette_discrete   = c("#3283FE", "#FEAF16", "#B00068", "#1CFFCE"),
  legend.title       = "",
  label.size         = 1,
  legendx            = -0.01
)
ggsave("grm.pdf", width = 8, height = 8, units = "in")

# ──────────────────────────────────────────────────────────────────────────────
# 10. Population-structure analysis via fastSTRUCTURE
#     full dataset, K = 1–6, 10 replicates
# ──────────────────────────────────────────────────────────────────────────────
res_fast <- gl.run.faststructure(
  x         = gl_qc2,
  k.range   = 1:6,
  num.k.rep = 10
)
saveRDS(res_fast, file = "res_fast_myna.rds")

gl.plot.faststructure(
  sr          = res_fast,
  k.range     = 2:4,
  x           = gl_qc2,
  den         = FALSE,
  label.size  = 10,
  met_clumpp  = "greedy",
  iter_clumpp = 1000
)
ggsave("faststruc.pdf", width = 10, height = 5, units = "in")

# ──────────────────────────────────────────────────────────────────────────────
# 11. Principal Coordinates Analysis (PCoA)
# ──────────────────────────────────────────────────────────────────────────────
pcoa_res <- gl.pcoa(gl_qc2)
gl.pcoa.plot(
  glPca      = pcoa_res,
  x          = gl_qc2,
  xaxis      = 1,
  yaxis      = 2,
  ellipse    = TRUE,
  pt.size    = 2,
  label.size = 1,
  pt.colors  = c("#3283FE", "#FEAF16", "#B00068", "#1CFFCE")
)
ggsave("pca_1_2.pdf", width = 6, height = 6, units = "in")

# ──────────────────────────────────────────────────────────────────────────────
# 12. Spatial autocorrelation of allele frequencies
# ──────────────────────────────────────────────────────────────────────────────
gl_qc3 <- gl_qc2
pop(gl_qc3) <- rep("Population", nInd(gl_qc3))
res_autocorr <- gl.spatial.autoCorr(
  x    = gl_qc3,
  bins = seq(0, 1e6, by = 2e5)
)
ggsave("autocorr.pdf", width = 6, height = 4, units = "in")

# ──────────────────────────────────────────────────────────────────────────────
# 13. Isolation by distance (IBD)
# ──────────────────────────────────────────────────────────────────────────────
res_ibd <- gl.ibd(gl_qc2, distance = "euclidean", paircols = "pop")
ggsave("ibd.pdf", width = 6, height = 4, units = "in")

# ──────────────────────────────────────────────────────────────────────────────
# 14. Effective population size (Ne) with LD-Ne
# ──────────────────────────────────────────────────────────────────────────────
ne_est <- gl.LDNe(gl_qc2)

# ──────────────────────────────────────────────────────────────────────────────
# 15. Genome-wide association (GWAS) for eight morphological traits
#     GAPIT, BLINK model, two principal components; Z-linked markers excluded
# ──────────────────────────────────────────────────────────────────────────────

# Working copy for the GWAS
data_gapit <- gl_qc2

# Remove loci with missing or invalid chromosome names
chr_raw    <- as.character(data_gapit$chromosome)
keep_chr   <- !is.na(chr_raw) & chr_raw != "<NA>"
data_gapit <- data_gapit[, keep_chr]

# Standardise chromosome labels against the assembly sequence report
chrom_tbl <- read.csv("sequence_report.csv", stringsAsFactors = FALSE)
chrom_tbl <- chrom_tbl[, c("GenBank.seq.accession", "chrom")]
colnames(chrom_tbl) <- c("chr", "chrom")

chroms <- data.frame(
  chr   = as.character(data_gapit$chromosome),
  pos   = data_gapit$position,
  order = seq_len(nLoc(data_gapit)),
  stringsAsFactors = FALSE
)
chroms$chrom_new      <- chrom_tbl$chrom[match(chroms$chr, chrom_tbl$chr)]
data_gapit$chromosome <- as.factor(chroms$chrom_new)

# Exclude unplaced ("Un") and Z-chromosome loci
loc_remove <- which(as.character(data_gapit$chromosome) %in% c("Un", "Z"))
if (length(loc_remove) > 0) {
  data_gapit <- gl.drop.loc(
    data_gapit,
    loc.list = locNames(data_gapit)[loc_remove]
  )
}

# Keep chromosomes with more than 100 loci
chrom_keep_tbl <- table(data_gapit$chromosome)
chrom_keep     <- names(chrom_keep_tbl[chrom_keep_tbl > 100])
loc_keep       <- which(as.character(data_gapit$chromosome) %in% chrom_keep)
data_gapit <- gl.keep.loc(
  data_gapit,
  loc.list = locNames(data_gapit)[loc_keep]
)

gl.plot.snp.density(data_gapit)

# Genotypes in GAPIT format (computed once; shared across traits)
geno_hmp <- as.data.frame(gl2gapit(x = data_gapit))

# Chromosome colour palette for the Manhattan plots
polychrome <- function(n) {
  rep(
    colorspace::darken(
      c("#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF",
        "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356",
        "#85660D", "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF",
        "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C",
        "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF",
        "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3B00FB"),
      0.2
    ),
    length.out = n
  )
}

# Eight traits in the order of Supplementary Figure S2 (a–h). Names are the
# metadata columns; values are the Manhattan-plot titles.
traits <- c(
  Weiht_rams         = "Weight",
  body_width         = "Body width",
  body_circumference = "Body circumference",
  Wing               = "Wing length",
  Tarsus             = "Tarsus length",
  wingspan           = "Wingspan",
  Beak_length        = "Beak length",
  length_body        = "Body length"
)

for (trait in names(traits)) {

  # Phenotype table for this trait; drop individuals with missing values
  pheno_df <- data_gapit$other$ind.metrics[, c("id", trait)]
  colnames(pheno_df) <- c("Taxa", trait)
  pheno_df[[trait]]  <- as.numeric(pheno_df[[trait]])
  pheno_df <- pheno_df[!is.na(pheno_df[[trait]]), ]

  # Run GWAS with the BLINK model and two principal components
  GAPIT(
    Y         = pheno_df,
    G         = geno_hmp,
    PCA.total = 2,
    model     = "Blink"
  )

  # Locate this trait's BLINK result table
  res_file <- list.files(
    pattern = paste0("GWAS_Results.*BLINK.*", trait, ".*\\.csv$")
  )
  if (length(res_file) == 0) {
    warning("No BLINK result file found for trait: ", trait)
    next
  }
  gwas_tab <- read.csv(res_file[1], stringsAsFactors = FALSE)

  # Manhattan plot: genome-wide (7.3) and suggestive (5.0) thresholds
  plot_df <- gwas_tab[, c("SNP", "Chr", "Pos", "P.value")]
  plot_df$Chr     <- as.numeric(plot_df$Chr)
  plot_df$Pos     <- as.numeric(plot_df$Pos)
  plot_df$P.value <- as.numeric(plot_df$P.value)
  plot_df <- plot_df[complete.cases(plot_df), ]

  CMplot(
    plot_df,
    plot.type     = "m",
    LOG10         = TRUE,
    col           = polychrome(length(unique(plot_df$Chr))),
    threshold     = c(5e-8, 1e-5),
    threshold.col = c("red", "blue"),
    amplify       = TRUE,
    main          = traits[[trait]],
    signal.cex    = 1,
    file.name     = trait
  )
}

# End of analysis script
