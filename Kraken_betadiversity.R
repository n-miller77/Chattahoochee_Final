beta_diversity.py -i /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/kraken_output/bracken_output2/empty_genus_files/*bracken.G.out --type bracken > /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/beta_diversity_S_kreport.bracken.G.tsv

beta_diversity.py \
  -i /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/kraken_output/bracken_output2/empty_genus_files/*bracken.G.out \
  --type bracken \
  --level G \
  > /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/beta_diversity_G_bracken_no_wwtp.tsv



# Load required libraries
library(gplots)
library("pheatmap")
library(ape)

#######################################################################
#################### Beta Diversity Analysis ##########################
#######################################################################

basic_sample_metadata <- read.csv("/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/basic_sample_metadata_no_wwtp.csv")

########################################
########### Create Heatmap #############
########################################

input_data <- read.table("/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/beta_diversity_no_wwtp_kreport.bracken.G.fixed.tsv", quote="\"", comment.char="", check.names=FALSE)

# Replace missing values with a specified value (e.g., 0)
input_data[is.na(input_data)] <- 1

# Create a data frame with sample names and corresponding metadata
metadata_df <- basic_sample_metadata[, c("site", "subsite", "batch")]

# Set row names to be the sample names
rownames(metadata_df) <- basic_sample_metadata$sample

# Specify colors
ann_colors = list(
  subsite = c(REG = "gray", INF = "yellow4", EFF = "yellow1", UP = "coral1", DOWN = "coral4"),
  site = c( DAM = "aliceblue", RL = "lightblue1", SC = "lightblue3", CC = "lightblue4", S1 = "mediumpurple1", S2 = "mediumpurple3", S3 = "mediumpurple4"),
  batch = c("period_1" = "darkolivegreen1", "period_2" = "darkolivegreen3", "period_3" = "darkolivegreen")
)

# Plot heatmap
pheatmap(input_data,
         annotation_col = metadata_df,
         annotation_colors = ann_colors,
         main = "Distance Matrix Heatmap",  # Set the plot title
         xlab = "Samples",  # Set the x-axis label
         ylab = "Samples"   # Set the y-axis label
)


