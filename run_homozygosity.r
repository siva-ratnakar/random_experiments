### Run it as:
#Rscript run_homozygosity.r <path_to_csv> <path_to_output> 2>&1 | tee <path_to_output>/analysis.log
## If all files are in the same directory:
## Rscript run_homozygosity.r ./MyHeritage_raw_dna_data.csv . 2>&1 | tee ./analysis.log
## This script processes MyHeritage DNA data to identify runs of homozygosity (ROH) and estimate inbreeding coefficients and consanguinity relationships of parents.
## The input file should be in CSV format with columns: RSID, CHROMOSOME, POSITION, RESULT. (The same as what you download from MyHeritage DNA result)
## The output directory should be a path where the results and visualizations will be saved.

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ggplot2)
})

# Function to read and process MyHeritage DNA data
process_myheritage_data <- function(file_path) {
  # Read the file, skipping header lines that start with "#"
  raw_data <- fread(file_path, skip = "RSID", header = TRUE)
  
  # Rename columns if needed (adjust based on your file format)
  if (!all(c("RSID", "CHROMOSOME", "POSITION", "RESULT") %in% names(raw_data))) {
    stop("File format doesn't match expected MyHeritage format")
  }
  
  # Convert chromosome column to factor for proper sorting
  raw_data$CHROMOSOME <- factor(raw_data$CHROMOSOME, 
                               levels = c(1:22, "X", "Y", "MT"))
  
  # Extract only SNPs with valid calls (exclude no-calls or indels)
  valid_data <- valid_data[valid_data$RESULT != "--"]
  
  return(valid_data)
}

# Function to determine sex from DNA data
determine_sex <- function(data) {
  # Count X and Y chromosome markers
  x_count <- sum(data$CHROMOSOME == "X")
  y_count <- sum(data$CHROMOSOME == "Y")
  
  cat("X chromosome markers:", x_count, "\n")
  cat("Y chromosome markers:", y_count, "\n")
  
  # If Y chromosome markers present, likely male
  if (y_count > 10) {
    return("male")
  } else {
    return("female")
  }
}

# Function to add homozygosity column, accounting for sex chromosomes
add_homozygosity <- function(data, sex) {
  # For autosomal chromosomes, check if both alleles are the same
  data$homozygous <- substr(data$RESULT, 1, 1) == substr(data$RESULT, 2, 2)
  
  if (sex == "male") {
    # For males, X and Y chromosomes should not be considered homozygous
    # as they only have one copy
    data$homozygous[data$CHROMOSOME %in% c("X", "Y")] <- FALSE
  }
  
  return(data)
}

# Function to identify runs of homozygosity
find_roh <- function(data, sex, min_snps = 25, min_kb = 1000, max_gap = 1000000) {
  roh_list <- list()
  
  # Exclude X chromosome for males
  if (sex == "male") {
    data <- data[!data$CHROMOSOME %in% c("X", "Y", "MT")]
  } else {
    data <- data[!data$CHROMOSOME %in% c("Y", "MT")]
  }
  
  # Ensure homozygous column exists and has no NA values
  if (!"homozygous" %in% names(data)) {
    stop("homozygous column missing from data")
  }
  data[is.na(data$homozygous), "homozygous"] <- FALSE
  
  for (chrom in levels(data$CHROMOSOME)) {
    chrom_data <- data[data$CHROMOSOME == chrom]
    if (nrow(chrom_data) == 0) next
    
    setorder(chrom_data, POSITION)
    
    current_start <- NULL
    current_end <- NULL
    snp_count <- 0
    
    for (i in 1:nrow(chrom_data)) {
      current_row <- chrom_data[i,]
      
      # Safe logical check
      is_homozygous <- isTRUE(current_row$homozygous)
      
      if (is_homozygous) {
        if (is.null(current_start)) {
          current_start <- current_row$POSITION
          current_end <- current_row$POSITION
          snp_count <- 1
        } else if (current_row$POSITION - current_end <= max_gap) {
          current_end <- current_row$POSITION
          snp_count <- snp_count + 1
        } else {
          if (snp_count >= min_snps && (current_end - current_start) >= min_kb) {
            roh_list[[length(roh_list) + 1]] <- list(
              chromosome = chrom,
              start = current_start,
              end = current_end,
              length_kb = (current_end - current_start) / 1000,
              snp_count = snp_count
            )
          }
          current_start <- current_row$POSITION
          current_end <- current_row$POSITION
          snp_count <- 1
        }
      } else {
        if (!is.null(current_start) && snp_count >= min_snps && (current_end - current_start) >= min_kb) {
          roh_list[[length(roh_list) + 1]] <- list(
            chromosome = chrom,
            start = current_start,
            end = current_end,
            length_kb = (current_end - current_start) / 1000,
            snp_count = snp_count
          )
        }
        current_start <- NULL
        current_end <- NULL
        snp_count <- 0
      }
    }
    
    # Check final segment
    if (!is.null(current_start) && snp_count >= min_snps && (current_end - current_start) >= min_kb) {
      roh_list[[length(roh_list) + 1]] <- list(
        chromosome = chrom,
        start = current_start,
        end = current_end,
        length_kb = (current_end - current_start) / 1000,
        snp_count = snp_count
      )
    }
  }
  
  if (length(roh_list) > 0) {
    roh_dt <- rbindlist(roh_list)
    return(roh_dt)
  } else {
    return(data.table(
      chromosome = character(),
      start = numeric(),
      end = numeric(),
      length_kb = numeric(),
      snp_count = numeric()
    ))
  }
}

# Function to estimate inbreeding coefficient based on ROH
# Using only autosomal chromosomes
estimate_inbreeding <- function(roh_data, genome_length = 2881033286) {
  # Calculate total ROH length in bp
  total_roh_length <- sum(roh_data$length_kb) * 1000
  
  # Estimate F_ROH (inbreeding coefficient based on ROH)
  f_roh <- total_roh_length / genome_length
  
  return(f_roh)
}

# Function to estimate consanguinity relationship based on F_ROH
estimate_relationship <- function(f_roh) {
  # Theoretical F values for different relationships
  relationships <- data.table(
    relationship = c("Unrelated", "Third cousins", "Second cousins", 
                     "First cousins once removed", "First cousins", 
                     "Double first cousins/Avuncular", "Half siblings", 
                     "Parent-child/Full siblings"),
    expected_f = c(0, 0.0039, 0.0156, 0.0313, 0.0625, 0.125, 0.125, 0.25)
  )
  
  # Ensure diff is numeric
  relationships[, diff := abs(expected_f - f_roh)]
  
  # Find the closest relationship
  closest <- relationships[which.min(relationships$diff)]
  
  # Get the second closest relationship
  relationships[which.min(relationships$diff), diff := Inf]  # Temporarily set min to Inf
  second_closest <- relationships[which.min(relationships$diff)]
  
  return(list(
    likely_relationship = closest$relationship,
    expected_f = closest$expected_f,
    observed_f = f_roh,
    second_relationship = second_closest$relationship,
    second_expected_f = second_closest$expected_f
  ))
}

# Function to calculate ROH statistics by chromosome
roh_by_chromosome <- function(roh_data) {
  # Calculate statistics by chromosome
  chrom_stats <- roh_data[, .(
    total_length_kb = sum(length_kb),
    average_length_kb = mean(length_kb),
    count = .N
  ), by = chromosome]
  
  return(chrom_stats)
}

# Main analysis function
analyze_consanguinity <- function(file_path, output_dir) {
  genome_length <- 2881033286  # Total genome length in base pairs
  
  cat("Reading and processing data...\n")
  dna_data <- process_myheritage_data(file_path)
  
  cat("\nDetermining biological sex from markers...\n")
  sex <- determine_sex(dna_data)
  cat("Determined biological sex:", sex, "\n")
  
  dna_data <- add_homozygosity(dna_data, sex)
  
  autosomal_data <- dna_data[!dna_data$CHROMOSOME %in% c("X", "Y", "MT")]
  cat("\nData summary (autosomal chromosomes only):\n")
  cat("Total autosomal SNPs analyzed:", nrow(autosomal_data), "\n")
  cat("Homozygous autosomal SNPs:", sum(autosomal_data$homozygous), 
      sprintf("(%.2f%%)\n", 100 * sum(autosomal_data$homozygous)/nrow(autosomal_data)))
  
  cat("\nFinding runs of homozygosity (excluding X if male)...\n")
  roh <- find_roh(dna_data, sex)
  
  if (nrow(roh) == 0) {
    cat("No significant runs of homozygosity found.\n")
    return(NULL)
  }
  
  total_roh_mb <- sum(roh$length_kb) / 1000
  avg_roh_kb <- mean(roh$length_kb)
  total_roh_percentage <- (total_roh_mb * 1000000 / genome_length) * 100
  
  cat("\nRuns of Homozygosity Results:\n")
  cat("Number of ROH segments:", nrow(roh), "\n")
  cat("Total ROH length:", round(total_roh_mb, 2), "Mb\n")
  cat("Total ROH length as percentage of genome:", round(total_roh_percentage, 2), "%\n")
  cat("Average ROH length:", round(avg_roh_kb, 2), "kb\n")
  
  chrom_stats <- roh_by_chromosome(roh)
  cat("\nROH by chromosome (top 5 by total length):\n")
  top_chroms <- chrom_stats[order(-total_length_kb)][1:5]
  print(top_chroms)
  
  f_roh <- estimate_inbreeding(roh)
  cat("\nEstimated inbreeding coefficient (F_ROH):", round(f_roh, 4), "\n")
  
  relationship <- estimate_relationship(f_roh)
  cat("\nConsanguinity estimate:\n")
  cat("Likely parental relationship:", relationship$likely_relationship, "\n")
  cat("Expected F value for this relationship:", relationship$expected_f, "\n")
  cat("Observed F value from your genome:", round(relationship$observed_f, 4), "\n")
  cat("Alternative relationship:", relationship$second_relationship, "\n")
  cat("Expected F value for alternative:", relationship$second_expected_f, "\n")
  
  # Save results to output directory
  write.csv(roh, file.path(output_dir, "roh_segments.csv"), row.names = FALSE)
  
  # Create and save visualizations
  p1 <- plot_roh(roh)
  p2 <- plot_roh_dist(roh)
  
  ggsave(file.path(output_dir, "roh_segments.png"), plot = p1)
  ggsave(file.path(output_dir, "roh_segments_distribution.png"), plot = p2)
  
  return(list(
    sex = sex,
    roh_segments = roh,
    total_roh_mb = total_roh_mb,
    total_roh_percentage = total_roh_percentage,
    inbreeding_coefficient = f_roh,
    likely_relationship = relationship,
    chromosome_stats = chrom_stats
  ))
}

# Plotting functions
plot_roh <- function(roh_data) {
  # Convert chromosome to factor for proper ordering
  roh_data$chromosome <- factor(roh_data$chromosome, 
                               levels = c(1:22))
  
  # Create the plot
  p <- ggplot(roh_data, aes(x = chromosome, y = length_kb)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +
    labs(title = "Runs of Homozygosity by Chromosome (Autosomes Only)",
         x = "Chromosome",
         y = "ROH Length (kb)") +
    theme_minimal()
  
  return(p)
}

plot_roh_dist <- function(roh_data) {
  p <- ggplot(roh_data, aes(x = length_kb)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(title = "Distribution of ROH Segment Lengths",
         x = "ROH Length (kb)",
         y = "Count") +
    theme_minimal()
  
  return(p)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript runs_homozygosity.R <path_to_csv> <path_to_output>")
}

input_file <- args[1]
output_dir <- args[2]

# Run the analysis
analyze_consanguinity(input_file, output_dir)