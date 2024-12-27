# Principal-Component-Analysis
Performs Principal Component Analysis (PCA) on gene expression count data generated by `featureCounts`. It processes multiple `.txt` files from a specified directory, merges the data, and generates an interactive PCA plot using the `DESeq2` and `plotly` libraries. The resulting plot is saved as an HTML file for easy sharing and visualization. The script uses the `vst()` function from the `DESeq2` package to transform the raw count data. VST applies a transformation that approximates a log2 scale for large counts but avoids issues with zero or low counts. This transformation ensures that variance is stabilized across genes, making it easier to compare samples.
# Requirements
Before running the script, ensure you have the following:
	1.	R Installed: The script requires R to be installed on your system.
	•	Download R: https://cran.r-project.org/
	2.	R Libraries: The script uses several R packages. Install them by running the following command in R:

> install.packages(c("ggplot2", "plotly", "RColorBrewer", "htmlwidgets"))
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("DESeq2")

# Usage
1. Clone or Copy the Script
Save the script as `pca_analysis.R`.
2. Make It Executable
Run the following command to make the script executable:

> chmod +x pca_analysis.R

# Run the Script
Execute the script from the terminal with the following syntax:

> ./pca_analysis.R <input_directory> <output_file>

•	`<input_directory>`: Path to the directory containing `.txt` files generated by `featureCounts`.
•	`<output_file>`: Path where the output HTML file (PCA plot) will be saved.

# Example 
> ./pca_analysis.R path/to/your/featurecounts.txt samplePCAplot.html
