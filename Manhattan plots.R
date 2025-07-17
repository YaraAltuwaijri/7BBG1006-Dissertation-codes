# Load required libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)

# Define the input gzipped VCF file and output directory
vcf_file <- "C:/Users/yaraa/Downloads/plots/ALDH2_merged.vcf.gz"
output_dir <- "C:/Users/yaraa/Downloads/plots/"


# Read the gzipped VCF file
vcf_data <- fread(vcf_file, skip = "#CHROM")

# Extract the column headers
column_headers <- colnames(vcf_data)

# Extract the study names from the VCF header (starting from the 10th column)
study_names <- column_headers[10:length(column_headers)]

# Initialize variables to store limits for x and y axes
x_limits <- range(vcf_data$V2, na.rm = TRUE)
y_limits <- c(0, 0)  # Initialize y limits

# Loop through each study to calculate y limits
for (study in study_names) {
  study_data <- vcf_data[, .(CHROM = get("#CHROM"), POS = POS, INFO = get(study))]
  study_data <- study_data %>%
    separate(INFO, into = c("Field1", "Field2", "LogP"), sep = ":") %>%
    mutate(LogP = as.numeric(LogP)) %>%
    na.omit()
  
  # Update y_limits based on the current study
  if (nrow(study_data) > 0) {
    max_logp <- max(study_data$LogP, na.rm = TRUE)
    y_limits[2] <- max(y_limits[2], max_logp)
  }
}

# Loop through each study to generate a Manhattan plot
for (study in study_names) {
  # Extract the data for the current study
  study_data <- vcf_data[, .(CHROM = get("#CHROM"), POS = POS, INFO = get(study))]
  study_data <- study_data %>%
    separate(INFO, into = c("Field1", "Field2", "LogP"), sep = ":") %>%
    mutate(LogP = as.numeric(LogP)) %>%
    na.omit()
  
  # Generate the Manhattan plot using ggplot2
  p <- ggplot(study_data, aes(x = POS, y = LogP)) +
    geom_point(alpha = 0.75, size = 0.5) +
    geom_hline(yintercept = 7.3, color = "red", linetype = "dashed", linewidth = 0.8) +  # genome-wide line
    coord_cartesian(ylim = y_limits) +  # consistent y-axis across plots
    labs(
      title = paste("Regional Manhattan Plot for Study:", study),
      x = "Genomic Position",
      y = "-log10(p-value)"
    ) +
    theme_minimal(base_family = "sans", base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Save the plot to a file
  output_file <- file.path(output_dir, paste0(study, "_manhattan_plot.png"))
  ggsave(output_file, plot = p, width = 10, height = 6)
}
