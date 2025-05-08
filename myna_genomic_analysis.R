# DESERT ISOLATION AND URBAN CONNECTIVITY: THE GENOMIC FOOTPRINT OF MYNA INVASION IN OMAN
# Qais Al Rawahi1, Abdullahi Aliyu1,2, Mazin Al‑Abidi3, Kareem Khalil5,
# Masooma Al‑Lawatiy1, Tehani A Al‑Jidailiy1, Basel S Al‑Maskari, Ahmed M Al‑Shakili3,
# Jihad A Al‑Toubi1, Maisaa S Al‑Saadi1, Adamu Abdul Abubakar1, Jose L Mijangos4,5.
#
# 1 Department of Veterinary Medicine, College of Applied and Health Sciences,
#   A’Sharqiyah University, P.O. Box 42, Postal code 400, Ibra, Sultanate of Oman.
# 2 Department of Veterinary Pathology, Faculty of Veterinary Medicine,
#   City Campus Complex, Usmanu Danfodiyo University, 840212 Sokoto,
#   Sokoto State, Nigeria.
# 3 Directorate General of Natural Conservation, Environment Authority,
#   P.O. Box 42, Postal code 400, Ibra, Sultanate of Oman.
# 4 Institute for Applied Ecology, University of Canberra, Bruce, Australia.
# 5 Diversity Arrays Technology Pty Ltd, University of Canberra,
#   Bruce, Australian Capital Territory, Australia.

# ──────────────────────────────────────────────────────────────────────────────
# 1. Load required package
# ──────────────────────────────────────────────────────────────────────────────
library(dartRverse)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Read SNP dataset and metadata
# ──────────────────────────────────────────────────────────────────────────────
gl_data <- gl.read.dart(
  infile            = "Report_DImy24-9735_SNP_mapping_1.csv",
  ind.metafile      = "metadata_myna.csv"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Annotate loci with trimmed sequence information
# ──────────────────────────────────────────────────────────────────────────────
gl_data$other$loc.metrics$TrimmedSequence <- 
  gl_data$other$loc.metrics$TrimmedSequenceSnp

# ──────────────────────────────────────────────────────────────────────────────
# 4. Map loci to reference genome & extract chromosome/position
# ──────────────────────────────────────────────────────────────────────────────
gl_data <- gl.blast(
  obj         = gl_data,
  ref_genome  = "GCA_037013685.1_AcTris_vAus2.0_genomic.fna"
)
# Store blast‑derived scaffold and start positions
gl_data$chromosome <- as.factor(gl_data$other$loc.metrics$sacc)
gl_data$position   <- gl_data$other$loc.metrics$sstart

# ──────────────────────────────────────────────────────────────────────────────
# 5. Quality control – call‑rate filtering
#    5.1 Remove individuals with ≤ 60% loci scored
# ──────────────────────────────────────────────────────────────────────────────
gl_qc1 <- gl.filter.callrate(
  obj       = gl_data,
  threshold = 0.6,
  method    = "ind"
)
# report remaining individual call‑rates
gl.report.callrate(gl_qc1, method = "ind")

# ──────────────────────────────────────────────────────────────────────────────
# 6. Quality control – locus‑level filtering
#    6.1 Remove loci with ≤ 90% individuals genotyped
#    6.2 Remove loci with low read depth
#    6.3 Remove loci with low reproducibility
# ──────────────────────────────────────────────────────────────────────────────
gl_qc2 <- gl.filter.callrate(
  obj       = gl_qc1,
  threshold = 0.9,
  method    = "loc"
)
gl_qc2 <- gl.filter.rdepth(gl_qc2)
gl_qc2 <- gl.filter.reproducibility(gl_qc2)

# report post‑QC locus call‑rates
gl.report.callrate(gl_qc2, method = "loc")

# ──────────────────────────────────────────────────────────────────────────────
# 7. Summaries of genetic diversity and F‑statistics
# ──────────────────────────────────────────────────────────────────────────────
gl.report.heterozygosity(
  obj            = gl_qc2,
  plot.colors.pop = c("#3283FE","#FEAF16","#B00068","#1CFFCE")
)
gl.report.fstat(gl_qc2)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Genetic Relatedness Matrix (GRM)
# ──────────────────────────────────────────────────────────────────────────────

gl.grm(
  obj                  = gl_qc2,
  palette_convergent   = gl.colors("con"),
  palette_discrete     = c("#3283FE","#FEAF16","#B00068","#1CFFCE"),
  legend.title         = "",
  label.size           = 1
)

# ──────────────────────────────────────────────────────────────────────────────
# 9. Population‑structure analysis via fastStructure
#    • full dataset, K = 1–6, 20 replicates
# ──────────────────────────────────────────────────────────────────────────────
res_fast <- gl.run.faststructure(
  obj       = gl_qc2,
  k.range   = 1:6,
  num.k.rep = 20
)
saveRDS(res_fast, file = "res_fast_myna.rds")

gl.plot.faststructure(
  res          = res_fast,
  k.range      = 2:6,
  x            = gl_qc2,
  den          = TRUE,
  label.size   = 20,
  met_clumpp   = "greedy",
  iter_clumpp  = 1000
)

# ──────────────────────────────────────────────────────────────────────────────
# 10. Principal Coordinates Analysis (PCoA)
# ──────────────────────────────────────────────────────────────────────────────
pcoa_res <- gl.pcoa(gl_qc2)
gl.pcoa.plot(
  pcoa.res    = pcoa_res,
  gl.obj      = gl_qc2,
  xaxis       = 1,
  yaxis       = 2,
  ellipse     = TRUE,
  pt.size     = 3,
  label.size  = 2,
  pt.colors   = c("#3283FE","#FEAF16","#B00068","#1CFFCE")
)

# ──────────────────────────────────────────────────────────────────────────────
# 11. Spatial autocorrelation of allele frequencies
# ──────────────────────────────────────────────────────────────────────────────
pop(gl_qc2) <- rep("Population", nInd(gl_qc2))
gl.spatial.autoCorr(
  obj  = gl_qc2,
  bins = seq(0, 1e6, by = 2e4)
)

# ──────────────────────────────────────────────────────────────────────────────
# 12. Estimate effective population size (Ne) with LD‑Ne
# ──────────────────────────────────────────────────────────────────────────────
ne_est <- gl.LDNe(gl_qc2)

# End of analysis script
