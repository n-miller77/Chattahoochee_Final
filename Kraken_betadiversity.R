beta_diversity.py -i /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/kraken_output/bracken_output2/empty_genus_files/*bracken.G.out --type bracken > /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/beta_diversity_S_kreport.bracken.G.tsv

beta_diversity.py \
  -i /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/kraken_output/bracken_output2/empty_genus_files/*bracken.G.out \
  --type bracken \
  --level G \
  > /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/beta_diversity_G_bracken_no_wwtp.tsv




# Load required libraries
library(gplots)
library(pheatmap)
library(ape)

#######################################################################
#################### Beta Diversity Analysis ##########################
#######################################################################

basic_sample_metadata <- read.csv("/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/basic_sample_metadata_1124.csv")

########################################
########### Create Heatmap #############
########################################

input_data <- read.table("/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/beta_diversity_no_wwtp_kreport.bracken.G.fixed.tsv",
                         quote="\"", comment.char="", check.names=FALSE)

# Replace missing values with a specified value
input_data[is.na(input_data)] <- 1

# Create annotation dataframe
metadata_df <- basic_sample_metadata[, c("site", "subsite", "batch")]
rownames(metadata_df) <- basic_sample_metadata$sample

# ---- NEW: create mapping for alternate labels ----
custom_labels <- c(
  "KTK04_2_CC_Down2" = "1_C_DOWN",
  "KTK04_2_CC_UP" = "1_C_UP",
  "KTK04_2_DAM" = "1_DAM",
  "KTK04_2_RL_DOWN" = "1_A_DOWN",
  "KTK04_2_RL_UP" = "1_A_UP",
  "KTK04_2_S1" = "1_S1",
  "KTK04_2_S2" = "1_S2",
  "KTK04_2_S3" = "1_S3",
  "KTK04_2_SC_UP" = "1_B_UP",
  "KTK04_3_CC_DOWN" = "2_C_DOWN",
  "KTK04_3_CC_UP" = "2_C_UP",
  "KTK04_3_DAM" = "2_DAM",
  "KTK04_3_RL_DOWN" = "2_A_DOWN",
  "KTK04_3_RL_UP" = "2_A_UP",
  "KTK04_3_S1" = "2_S1",
  "KTK04_3_S2" = "2_S2",
  "KTK04_3_S3" = "2_S3",
  "KTK04_3_SC_UP" = "2_B_UP",
  "KTK04_4_CC_DOWN" = "3_C_DOWN",
  "KTK04_4_CC_UP" = "3_C_UP",
  "KTK04_4_DAM" = "3_DAM",
  "KTK04_4_RL_DOWN" = "3_A_DOWN",
  "KTK04_4_RL_UP" = "3_A_UP",
  "KTK04_4_S1" = "3_S1",
  "KTK04_4_S2" = "3_S2",
  "KTK04_4_S3" = "3_S3",
  "KTK04_4_SC_UP" = "3_B_UP"
  # add all mappings here
)

# Create labels for rows and columns
labels_col <- custom_labels[colnames(input_data)]
labels_row <- custom_labels[rownames(input_data)]

# If any samples lack a custom label → keep original
labels_col[is.na(labels_col)] <- colnames(input_data)[is.na(labels_col)]
labels_row[is.na(labels_row)] <- rownames(input_data)[is.na(labels_row)]

# Colors
ann_colors = list(
  subsite = c(REG = "gray", INF = "yellow4", EFF = "yellow1",
              UP = "coral1", DOWN = "coral4"),
  site = c(DAM="aliceblue", A="lightblue1", B="lightblue3",
           C="lightblue4", S1="mediumpurple1",
           S2="mediumpurple3", S3="mediumpurple4"),
  batch = c("period_1"="darkolivegreen1",
            "period_2"="darkolivegreen3",
            "period_3"="darkolivegreen")
)

# ---- Plot heatmap with custom labels ----
pheatmap(input_data,
         annotation_col = metadata_df,
         annotation_colors = ann_colors,
         labels_col = labels_col,    # x-axis labels
         labels_row = labels_row,    # y-axis labels
         main = "Beta-Diversity Distance Matrix Heatmap",
         xlab = "Samples",
         ylab = "Samples"
)












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


