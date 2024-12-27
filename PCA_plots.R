#!/usr/bin/env Rscript

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: ./pca_plots.R <input_directory> <output_file>")
}

# Assign arguments to variables
input_directory <- args[1]
output_file <- args[2]

# Load required libraries
suppressMessages(lapply(c("DESeq2", "ggplot2", "plotly", "RColorBrewer", "htmlwidgets"), library, character.only = TRUE))

# List all featureCounts output files from the input directory
file_list <- list.files(path = input_directory, 
                        pattern = "*.txt", 
                        full.names = TRUE)

# Function to read featureCounts output
read_featureCounts <- function(file) {
  counts <- read.table(file, header = TRUE, skip = 1)
  counts <- counts[, c(1, (ncol(counts)-2):ncol(counts))]
  colnames(counts) <- c("Geneid", 
                        paste0(gsub(".*(WT|SNF2)(_[0-9]+).*", "\\1\\2", basename(file)), 
                               c("_1", "_2", "_3")))
  return(counts)
}

# Apply the function to all files
count_list <- lapply(file_list, read_featureCounts)

# Merge all count data
merged_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)

# Set row names and remove the Geneid column
rownames(merged_counts) <- merged_counts$Geneid
merged_counts <- merged_counts[, -1]

# Replace NA values with 0
merged_counts[is.na(merged_counts)] <- 0

# Convert to matrix of integers
merged_counts <- as.matrix(merged_counts)
mode(merged_counts) <- "integer"

# Create sample information
sample_info <- data.frame(
  row.names = colnames(merged_counts),
  condition = factor(gsub("_[0-9]+$", "", colnames(merged_counts))),
  replicate = factor(rep(1:3, length.out = ncol(merged_counts)))
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = merged_counts,
                              colData = sample_info,
                              design = ~ condition + replicate)

# Perform variance stabilizing transformation
vst_data <- vst(dds, blind=TRUE)

# Perform PCA
pca_result <- prcomp(t(assay(vst_data)))

# Calculate percentage of variance explained
percent_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Extract PCA coordinates
pca_data <- as.data.frame(pca_result$x)

# Add sample information to PCA data
pca_data$sample <- rownames(pca_data)
pca_data$condition <- dds$condition
pca_data$replicate <- dds$replicate

# Get unique conditions
conditions <- unique(pca_data$condition)
n_conditions <- length(conditions)

# Generate color palette for conditions
if (n_conditions <= 9) {
  color_palette <- brewer.pal(n_conditions, "Set1")
} else if (n_conditions <= 12) {
  color_palette <- brewer.pal(n_conditions, "Paired")
} else {
  color_palette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Paired"))
  color_palette <- color_palette[1:n_conditions]
}

color_vector <- setNames(color_palette, conditions)

# Create interactive PCA plot with Plotly
p <- plot_ly(pca_data, x = ~PC1, y = ~PC2, 
             color = ~condition, 
             colors = color_vector,
             text = ~paste("Sample:", sample, "<br>Condition:", condition, "<br>Replicate:", replicate),
             hoverinfo = "text") %>%
  add_markers(size = 10) %>%
  layout(
    title = "PCA Plot of VST Transformed Data",
    xaxis = list(
      title = paste0("PC1: ", percent_var[1], "% variance"),
      zeroline = FALSE
    ),
    yaxis = list(
      title = paste0("PC2: ", percent_var[2], "% variance"),
      zeroline = FALSE
    ),
    legend = list(title = list(text = "Condition"))
  )

# Save the interactive plot as an HTML file to the specified output file path
htmlwidgets::saveWidget(
  widget = p,
  file = output_file,
  selfcontained = TRUE,
  background = "white",
  title = "PCA Plot"
)
