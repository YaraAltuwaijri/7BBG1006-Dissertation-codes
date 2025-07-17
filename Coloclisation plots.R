# Load libraries
library(data.table)
library(tidyr)
library(ggplot2)
library(coloc)

# Define input and output directories
input_dir <- "C:/Users/yaraa/Downloads/plots/coloclization/"
output_dir <- "C:/Users/yaraa/Downloads/plots/coloclization/coloc.plots/"

# Define VCF input files
files <- list(
  "Psoriasis A (GCST90014456)" = file.path(input_dir, "GCST90014456_buildGRCh38.vcf.gz"),
  "Psoriasis B (GCST90014457)" = file.path(input_dir, "GCST90014457_buildGRCh38.vcf.gz"),
  "Lichen sclerosus (LS_meta)" = file.path(input_dir, "LS_meta_v3_ukb_fg_hunt1_buildGRCh38.vcf.gz")
)

# Function to extract data
extract_data <- function(vcf_file) {
  vcf <- fread(vcf_file, skip = "#CHROM")
  colnames(vcf)[1] <- "CHROM"
  info <- vcf[[10]]  # Assuming INFO is in column 10
  
  parts <- tstrsplit(info, ":", fixed = TRUE)
  if (length(parts) < 4) stop("Expecting INFO format: ES:SE:LP:AF")
  
  dt <- data.table(
    POS = vcf$POS,
    BETA = as.numeric(parts[[1]]),
    SE = as.numeric(parts[[2]]),
    LP = as.numeric(parts[[3]]),
    EAF = as.numeric(parts[[4]])
  )
  dt <- na.omit(dt)
  return(dt)
}

# Load data
trait_data <- lapply(files, extract_data)

# Pairwise combinations
combinations <- combn(names(trait_data), 2, simplify = FALSE)

# Run coloc + plot
for (combo in combinations) {
  trait1 <- combo[1]
  trait2 <- combo[2]
  
  dt1 <- trait_data[[trait1]]
  dt2 <- trait_data[[trait2]]
  
  merged <- merge(dt1, dt2, by = "POS", suffixes = c(".1", ".2"))
  if (nrow(merged) < 100) next
  
  # Remove duplicate positions
  merged <- merged[!duplicated(merged$POS), ]
  
  # Run colocalisation
  coloc_res <- coloc.abf(
    dataset1 = list(beta = merged$BETA.1, varbeta = merged$SE.1^2, MAF = merged$EAF.1, snp = merged$POS, N = 10000, type = "quant"),
    dataset2 = list(beta = merged$BETA.2, varbeta = merged$SE.2^2, MAF = merged$EAF.2, snp = merged$POS, N = 10000, type = "quant")
  )
  
  pp4 <- round(coloc_res$summary["PP.H4.abf"], 3)
  
  # Plot
  p <- ggplot(merged, aes(x = LP.1, y = LP.2)) +
    geom_point(alpha = 0.75, size = 0.5) +
    labs(
      title = paste0("Colocalisation: ", trait1, " × ", trait2, "\nPP4 = ", pp4),
      x = paste0("-log10(p): ", trait1),
      y = paste0("-log10(p): ", trait2)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(size = 11)
    )
  
  # Save
  output_file <- paste0(
    gsub("[^a-zA-Z0-9]", "_", trait1),
    "_vs_",
    gsub("[^a-zA-Z0-9]", "_", trait2),
    "_coloc_plot.png"
  )
  
  ggsave(file.path(output_dir, output_file), plot = p, width = 8, height = 6)
}
