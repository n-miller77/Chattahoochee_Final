# === Load Required Libraries ===
library(ggplot2)
library(dplyr)
library(readr)
library(Polychrome)

# === STEP 1: Define File List ===
files <- list(
  River = c(
    "Period1_River.csv",
    "Period2_River.csv",
    "Period3_River.csv"
  ),
  WWTP = c(
    "Period1_WWTP.csv",
    "Period2_WWTP.csv",
    "Period3_WWTP.csv"
  )
)

# === STEP 2: Read and Tag Data ===
read_and_tag <- function(file_path, period_label, source_label) {
  df <- read_csv(file_path, show_col_types = FALSE)
  colnames(df)[1:3] <- c("Sample", "Abundance", "Class")
  
  df <- df %>%
    filter(!is.na(Class), Abundance > 0) %>%
    mutate(
      Period = period_label,
      Source = source_label,
      Sample = factor(Sample, levels = Sample)  # Preserve row order
    )
  return(df)
}

# === STEP 3: Load and Combine All Data ===
all_data <- bind_rows(
  read_and_tag(files$River[1], "Period 1", "River"),
  read_and_tag(files$River[2], "Period 2", "River"),
  read_and_tag(files$River[3], "Period 3", "River"),
  read_and_tag(files$WWTP[1], "Period 1", "WWTP"),
  read_and_tag(files$WWTP[2], "Period 2", "WWTP"),
  read_and_tag(files$WWTP[3], "Period 3", "WWTP")
)

# === STEP 4: Create Pastel Color Palette ===
all_classes <- unique(all_data$Class)
set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 5: Factor Levels ===
all_data$Class <- factor(all_data$Class, levels = names(final_color_map))
all_data$Period <- factor(all_data$Period, levels = c("Period 1", "Period 2", "Period 3"))
all_data$Source <- factor(all_data$Source, levels = c("River", "WWTP"))

# === STEP 6: Plot Function ===
plot_faceted_horizontal <- function(df_subset, source_label) {
  df_subset$Class <- factor(df_subset$Class, levels = names(final_color_map))
  
  # 🟡 Preserve sample order within each period
  df_subset <- df_subset %>%
    group_by(Period) %>%
    mutate(Sample = factor(Sample, levels = Sample)) %>%
    ungroup()
  
  ggplot(df_subset, aes(y = Sample, x = Abundance, fill = Class)) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_wrap(~ Period, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = final_color_map, name = "Class") +
    labs(
      title = paste("MAG Distribution by Sample -", source_label),
      x = "Percent Abundance",
      y = "Sample"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# === STEP 7: Generate and Save Plots ===
for (source_type in c("River", "WWTP")) {
  subset_df <- all_data %>% filter(Source == source_type)
  plot <- plot_faceted_horizontal(subset_df, source_type)
  output_file <- paste0("MAG_distribution_", tolower(source_type), ".pdf")
  ggsave(output_file, plot = plot, width = 11, height = 13)
  cat("💾 Saved:", output_file, "\n")
}

cat("✅ Done! Two figures created with correct sample order and pastel colors.\n")






# === Load Required Libraries ===
library(ggplot2)
library(dplyr)
library(readr)
library(Polychrome)

# === STEP 1: Define File List ===
files <- list(
  River = c(
    "Period1_River.csv",
    "Period2_River.csv",
    "Period3_River.csv"
  ),
  WWTP = c(
    "Period1_WWTP.csv",
    "Period2_WWTP.csv",
    "Period3_WWTP.csv"
  )
)

# === STEP 2: Read and Tag All Data ===
read_and_tag <- function(file_path, period_label, source_label) {
  df <- read_csv(file_path, show_col_types = FALSE)
  colnames(df)[1:3] <- c("Sample", "Abundance", "Class")
  
  df <- df %>%
    filter(!is.na(Class), Abundance > 0) %>%
    mutate(
      Period = period_label,
      Source = source_label
    )
  return(df)
}

# Combine all into one data frame
all_data <- bind_rows(
  read_and_tag(files$River[1], "Period 1", "River"),
  read_and_tag(files$River[2], "Period 2", "River"),
  read_and_tag(files$River[3], "Period 3", "River"),
  read_and_tag(files$WWTP[1], "Period 1", "WWTP"),
  read_and_tag(files$WWTP[2], "Period 2", "WWTP"),
  read_and_tag(files$WWTP[3], "Period 3", "WWTP")
)

# === STEP 3: Generate Consistent Pastel Color Palette ===
all_classes <- unique(all_data$Class)
set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Soften color by blending toward white
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 4: Set factor levels ===
all_data$Period <- factor(all_data$Period, levels = c("Period 1", "Period 2", "Period 3"))
all_data$Class <- factor(all_data$Class, levels = names(final_color_map))

# === STEP 5: Plotting Function (for River and WWTP) ===
plot_faceted_horizontal <- function(df_subset, source_label) {
  df_subset$Sample <- factor(df_subset$Sample, levels = unique(df_subset$Sample))
  
  ggplot(df_subset, aes(y = Sample, x = Abundance, fill = Class)) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_wrap(~ Period, ncol = 1, scales = "free_y") +  # stacked vertically by period
    scale_fill_manual(values = final_color_map, name = "Class") +
    labs(
      title = paste("MAG Distribution by Sample -", source_label),
      x = "Percent Abundance",
      y = "Sample"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# === STEP 6: Generate and Save Final Figures ===
for (source_type in c("River", "WWTP")) {
  subset_df <- all_data %>% filter(Source == source_type)
  plot <- plot_faceted_horizontal(subset_df, source_type)
  output_file <- paste0("MAG_distribution_", tolower(source_type), ".pdf")
  ggsave(output_file, plot = plot, width = 11, height = 13)
  cat("💾 Saved:", output_file, "\n")
}

cat("✅ Done! Two figures created with pastel-colored stacked bar charts per sample and period.\n")





































library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)

# === STEP 1: Define File List ===
files <- c(
  "Period1_River.csv",
  "Period1_WWTP.csv",
  "Period2_River.csv",
  "Period2_WWTP.csv",
  "Period3_River.csv",
  "Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 River",
  "Period 1 WWTP",
  "Period 2 River",
  "Period 2 WWTP",
  "Period 3 River",
  "Period 3 WWTP"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and Label Each Dataset ===
all_data <- lapply(seq_along(files), function(i) {
  df <- read.csv(files[i])
  colnames(df)[1:3] <- c("Location", "Percent_Abundance", "Class")

  df <- df %>%
    filter(!is.na(Class), Percent_Abundance > 0) %>%
    mutate(
      Period = gsub("River|WWTP", "", names(files)[i]) |> trimws(),
      Source = ifelse(grepl("River", names(files)[i]), "River", "WWTP")
    )

  return(df)
})

# === STEP 3: Combine All Data ===
combined_data <- bind_rows(all_data)

# === STEP 4: Create Pastel Color Palette ===
all_classes <- unique(combined_data$Class)

set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Soften colors by blending with white
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 5: Factor Levels for Ordering ===
combined_data$Class <- factor(combined_data$Class, levels = names(final_color_map))
combined_data$Period <- factor(combined_data$Period, levels = c("Period 1", "Period 2", "Period 3"))
combined_data$Source <- factor(combined_data$Source, levels = c("River", "WWTP"))
combined_data$Location <- factor(combined_data$Location, levels = unique(combined_data$Location))

# === STEP 6: Plot Function for Each Source ===
plot_by_source <- function(data_subset, source_label) {
  ggplot(data_subset, aes(x = Location, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity", width = 0.5) +  # thinner bars
    facet_wrap(~ Period, ncol = 1, scales = "free_y") +  # vertical layout
    scale_fill_manual(values = final_color_map, name = "MAG Class") +
    labs(
      title = paste("MAG Distribution by Sample -", source_label),
      x = "Sample Location",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# === STEP 7: Generate and Save Both Figures ===
for (src in levels(combined_data$Source)) {
  df_subset <- combined_data %>% filter(Source == src)
  p <- plot_by_source(df_subset, src)
  out_file <- paste0("MAG_distribution_by_sample_", tolower(src), ".pdf")
  cat("💾 Saving:", out_file, "\n")
  ggsave(out_file, plot = p, width = 10, height = 12)
}

cat("✅ Done! Two vertically faceted sample-level figures saved: one for River, one for WWTP.\n")


























# === Load Libraries ===
library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)

# === STEP 1: Define File List ===
files <- c(
  "Period1_River.csv",
  "Period1_WWTP.csv",
  "Period2_River.csv",
  "Period2_WWTP.csv",
  "Period3_River.csv",
  "Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 River",
  "Period 1 WWTP",
  "Period 2 River",
  "Period 2 WWTP",
  "Period 3 River",
  "Period 3 WWTP"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and Label Each Dataset ===
all_data <- lapply(seq_along(files), function(i) {
  df <- read.csv(files[i])
  colnames(df)[1:3] <- c("Location", "Percent_Abundance", "Class")

  df <- df %>%
    filter(!is.na(Class), Percent_Abundance > 0) %>%
    mutate(
      Period = gsub("River|WWTP", "", names(files)[i]) |> trimws(),
      Source = ifelse(grepl("River", names(files)[i]), "River", "WWTP")
    )

  return(df)
})

# === STEP 3: Combine All Data ===
combined_data <- bind_rows(all_data)

# === STEP 4: Create Consistent Soft Pastel Color Palette ===
all_classes <- unique(combined_data$Class)

set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Soften colors by blending with white
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 5: Factor Ordering ===
combined_data$Class <- factor(combined_data$Class, levels = names(final_color_map))
combined_data$Period <- factor(combined_data$Period, levels = c("Period 1", "Period 2", "Period 3"))
combined_data$Source <- factor(combined_data$Source, levels = c("River", "WWTP"))

# === STEP 6: Summarize for Plotting ===
summary_data <- combined_data %>%
  group_by(Source, Period, Class) %>%
  summarise(Total_Abundance = sum(Percent_Abundance), .groups = "drop")

# === STEP 7: Define Plot Function ===
plot_stacked_bar <- function(data_subset, source_label) {
  ggplot(data_subset, aes(x = Period, y = Total_Abundance, fill = Class)) +
    geom_bar(stat = "identity", width = 0.4) +  # thin bars
    scale_fill_manual(values = final_color_map, name = "MAG Class") +
    labs(
      title = paste("Total MAG Abundance -", source_label),
      x = "Period",
      y = "Total Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# === STEP 8: Generate and Save Each Plot ===
for (src in levels(combined_data$Source)) {
  df_plot <- summary_data %>% filter(Source == src)
  p <- plot_stacked_bar(df_plot, src)
  out_file <- paste0("MAG_distribution_", tolower(src), ".pdf")
  cat("💾 Saving:", out_file, "\n")
  ggsave(out_file, plot = p, width = 6, height = 7)
}

cat("✅ Done! Two stacked bar plots saved: one for River, one for WWTP.\n")

















# === Load Libraries ===
library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)

# === STEP 1: Define File List ===
files <- c(
  "Period1_River.csv",
  "Period1_WWTP.csv",
  "Period2_River.csv",
  "Period2_WWTP.csv",
  "Period3_River.csv",
  "Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 River",
  "Period 1 WWTP",
  "Period 2 River",
  "Period 2 WWTP",
  "Period 3 River",
  "Period 3 WWTP"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and Label Each Dataset ===
all_data <- lapply(seq_along(files), function(i) {
  df <- read.csv(files[i])
  colnames(df)[1:3] <- c("Location", "Percent_Abundance", "Class")

  df <- df %>%
    filter(!is.na(Class), Percent_Abundance > 0) %>%
    mutate(
      Period = gsub("River|WWTP", "", names(files)[i]) |> trimws(),
      Source = ifelse(grepl("River", names(files)[i]), "River", "WWTP")
    )

  return(df)
})

# === STEP 3: Combine All Data ===
combined_data <- bind_rows(all_data)

# === STEP 4: Create Consistent Soft Pastel Color Palette ===
all_classes <- unique(combined_data$Class)

set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Soften toward pastel by blending with white
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 5: Summarize Abundance Per Class Per Period ===
summary_data <- combined_data %>%
  group_by(Period, Class) %>%
  summarise(Total_Abundance = sum(Percent_Abundance), .groups = "drop")

# Order factors
summary_data$Class <- factor(summary_data$Class, levels = names(final_color_map))
summary_data$Period <- factor(summary_data$Period, levels = c("Period 1", "Period 2", "Period 3"))

# === STEP 6: Create Final Stacked Bar Plot ===
p <- ggplot(summary_data, aes(x = Period, y = Total_Abundance, fill = Class)) +
  geom_bar(stat = "identity", width = 0.4) +  # thin bars
  scale_fill_manual(values = final_color_map, name = "MAG Class") +
  labs(
    title = "Total MAG Abundance per Period",
    x = "Period",
    y = "Total Percent Abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# === STEP 7: Save Output ===
output_file <- "MAG_distribution_stacked_per_period.pdf"
cat("💾 Saving:", output_file, "\n")
ggsave(output_file, plot = p, width = 6, height = 7)

cat("✅ Done! Stacked bar chart saved with pastel colors, one thin bar per period.\n")




library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)

# === STEP 1: Define File List ===
files <- c(
  "Period1_River.csv",
  "Period1_WWTP.csv",
  "Period2_River.csv",
  "Period2_WWTP.csv",
  "Period3_River.csv",
  "Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 River",
  "Period 1 WWTP",
  "Period 2 River",
  "Period 2 WWTP",
  "Period 3 River",
  "Period 3 WWTP"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and Label Each Dataset ===
all_data <- lapply(seq_along(files), function(i) {
  df <- read.csv(files[i])
  colnames(df)[1:3] <- c("Location", "Percent_Abundance", "Class")
  
  df <- df %>%
    filter(!is.na(Class), Percent_Abundance > 0) %>%
    mutate(
      Period = gsub("River|WWTP", "", names(files)[i]) |> trimws(),
      Source = ifelse(grepl("River", names(files)[i]), "River", "WWTP")
    )
  
  return(df)
})

combined_data <- bind_rows(all_data)

# === STEP 3: Generate Pastel Color Palette ===
all_classes <- unique(combined_data$Class)
set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Soften colors by blending toward white
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 4: Set Factors for Ordering ===
combined_data$Location <- factor(combined_data$Location, levels = unique(combined_data$Location))
combined_data$Class <- factor(combined_data$Class, levels = names(final_color_map))
combined_data$Period <- factor(combined_data$Period, levels = c("Period 1", "Period 2", "Period 3"))
combined_data$Source <- factor(combined_data$Source, levels = c("River", "WWTP"))

# === STEP 5: Function to Create and Save Faceted Plot ===
create_faceted_plot <- function(data_subset, source_label) {
  p <- ggplot(data_subset, aes(x = Location, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Period, nrow = 1, scales = "free_x") +
    coord_flip() +
    scale_fill_manual(values = final_color_map, name = "MAG Class") +
    labs(
      title = paste("Identified MAGs -", source_label),
      x = "Sample Location",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  output_file <- paste0("MAG_distribution_", tolower(source_label), "_periods.pdf")
  cat("💾 Saving:", output_file, "\n")
  ggsave(output_file, plot = p, width = 14, height = 7)
}

# === STEP 6: Generate the Two Figures ===
create_faceted_plot(filter(combined_data, Source == "River"), "River")
create_faceted_plot(filter(combined_data, Source == "WWTP"), "WWTP")

cat("✅ Two faceted figures saved: one for River, one for WWTP.\n")




















# Load required libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# List your 6 file names
file_list <- c("BLA_Period1_River.csv", "BLA_Period1_WWTP.csv", "BLA_Period2_River.csv", 
               "BLA_Period2_WWTP.csv", "BLA_Period3_River.csv", "BLA_Period3_WWTP.csv")

# Define a consistent pastel palette for the 5 categories
category_colors <- c(
  "BlaA" = "#AEC6CF",  # pastel blue
  "BlaC" = "#FFB347",  # pastel orange
  "MBL"  = "#B0E0E6",  # pastel turquoise
  "OXA"  = "#77DD77",  # pastel green
  "TEM"  = "#FF6961"   # pastel red
)

# Loop through each file to generate and save the plot
for (i in seq_along(file_list)) {
  
  # Read the data (change `sep` to "," if it's a CSV)
  df <- read.table(file_list[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Melt the data to long format
  df_long <- melt(df, id.vars = "Sample", 
                  variable.name = "Category", 
                  value.name = "Abundance")
  
  # Generate plot
  p <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = category_colors) +
    labs(title = paste("Stacked Bar Chart:", file_list[i]),
         x = "Sample", y = "Percent Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot as PNG
  ggsave(filename = paste0("stacked_bar_", i, ".png"), plot = p, width = 10, height = 6)
}






library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)

# === STEP 1: Define File List ===
files <- c(
  "BLA_Period1_River.csv",
  "BLA_Period1_WWTP.csv",
  "BLA_Period2_River.csv",
  "BLA_Period2_WWTP.csv",
  "BLA_Period3_River.csv",
  "BLA_Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 River",
  "Period 1 WWTP",
  "Period 2 River",
  "Period 2 WWTP",
  "Period 3 River",
  "Period 3 WWTP"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and Label Each Dataset ===
all_data <- lapply(seq_along(files), function(i) {
  df <- read.csv(files[i])
  colnames(df)[1:3] <- c("Location", "Percent_Abundance", "Class")
  
  df <- df %>%
    filter(!is.na(Class), Percent_Abundance > 0) %>%
    mutate(
      Period = gsub("River|WWTP", "", names(files)[i]) |> trimws(),
      Source = ifelse(grepl("River", names(files)[i]), "River", "WWTP")
    )
  
  return(df)
})

combined_data <- bind_rows(all_data)

# === STEP 3: Generate Pastel Color Palette ===
all_classes <- unique(combined_data$Class)
set.seed(123)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Soften colors by blending toward white
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 4: Set Factors for Ordering ===
combined_data$Location <- factor(combined_data$Location, levels = unique(combined_data$Location))
combined_data$Class <- factor(combined_data$Class, levels = names(final_color_map))
combined_data$Period <- factor(combined_data$Period, levels = c("Period 1", "Period 2", "Period 3"))
combined_data$Source <- factor(combined_data$Source, levels = c("River", "WWTP"))

# === STEP 5: Function to Create and Save Faceted Plot ===
create_faceted_plot <- function(data_subset, source_label) {
  p <- ggplot(data_subset, aes(x = Location, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Period, nrow = 1, scales = "free_x") +
    coord_flip() +
    scale_fill_manual(values = final_color_map, name = "MAG Class") +
    labs(
      title = paste("Identified MAGs -", source_label),
      x = "Sample Location",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  output_file <- paste0("MAG_distribution_", tolower(source_label), "_periods.pdf")
  cat("💾 Saving:", output_file, "\n")
  ggsave(output_file, plot = p, width = 14, height = 7)
}

# === STEP 6: Generate the Two Figures ===
create_faceted_plot(filter(combined_data, Source == "River"), "River")
create_faceted_plot(filter(combined_data, Source == "WWTP"), "WWTP")

cat("✅ Two faceted figures saved: one for River, one for WWTP.\n")
















# === Load Required Libraries ===
library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)

# === STEP 1: Define File List ===
files <- c(
  "Period1_River.csv",
  "Period1_WWTP.csv",
  "Period2_River.csv",
  "Period2_WWTP.csv",
  "Period3_River.csv",
  "Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 River",
  "Period 1 WWTP",
  "Period 2 River",
  "Period 2 WWTP",
  "Period 3 River",
  "Period 3 WWTP"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and Clean Each Dataset ===
all_data <- lapply(files, function(file) {
  df <- read.csv(file)

  # Standardize expected column names
  colnames(df)[1:3] <- c("Location", "Percent_Abundance", "Class")

  # Filter out NA and 0 values
  df <- df %>% 
    filter(!is.na(Class), Percent_Abundance > 0)

  return(df)
})

# === STEP 3: Identify Unique MAG Classes Across All Files ===
all_classes <- unique(unlist(lapply(all_data, function(df) unique(as.character(df$Class)))))
cat("🎨 Total unique classes across datasets:", length(all_classes), "\n")

# === STEP 4: Generate Soft Pastel-Like Color Palette ===
set.seed(123)  # For reproducibility
pastel_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# Optional: Sort palette by perceived brightness for better visual separation
pastel_palette <- pastel_palette[order(colSums(col2rgb(pastel_palette)))]
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 5: Generate and Save Plots ===
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  period_name <- names(files)[i]

  # Preserve sample order
  df$Location <- factor(df$Location, levels = unique(df$Location))

  # Order MAG classes by frequency of occurrence
  class_order <- df %>%
    group_by(Class) %>%
    summarise(freq = n_distinct(Location)) %>%
    arrange(desc(freq)) %>%
    pull(Class)
  df$Class <- factor(df$Class, levels = class_order)

  # === Create Bar Plot ===
  p <- ggplot(df, aes(x = Location, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("Identified MAGs -", period_name),
      x = "Sample Location",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    scale_fill_manual(values = final_color_map, name = "MAG Class")

  # === Save Plot to PDF ===
  output_file <- paste0("MAG_distribution_", gsub(" ", "_", tolower(period_name)), ".pdf")
  cat("💾 Saving plot:", output_file, "\n")
  ggsave(output_file, plot = p, width = 11, height = 6)
}

cat("✅ All plots saved with consistent pastel color mapping across 6 periods.\n")


library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(Polychrome)

cat("📂 Reading and processing wide-format datasets...\n")

# Step 1: CSV files in wide format
files <- list(
  "Period 1 River" = "BLA_Period1_River.csv",
  "Period 1 WWTP" = "BLA_Period1_WWTP.csv",
  "Period 2 River" = "BLA_Period2_River.csv",
  "Period 2 WWTP" = "BLA_Period2_WWTP.csv",
  "Period 3 River" = "BLA_Period3_River.csv",
  "Period 3 WWTP" = "BLA_Period3_WWTP.csv"
)

# Step 2: Read and reshape to long format
all_data <- lapply(files, function(file) {
  df <- read.csv(file, check.names = FALSE)
  colnames(df)[1] <- "Sample"

  df_long <- pivot_longer(df, cols = -Sample, names_to = "Class", values_to = "Percent_Abundance") %>%
    filter(!is.na(Percent_Abundance) & Percent_Abundance > 0)

  return(df_long)
})

# Step 3: Unique classes and pastel palette with Polychrome
all_classes <- unique(unlist(lapply(all_data, \(df) unique(df$Class))))

set.seed(42)
base_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")

# ✅ Softening function: blend toward white to make pastel
soften_color <- function(hex_color, amount = 0.5) {
  rgb_vals <- col2rgb(hex_color) / 255
  pastel_vals <- rgb_vals + (1 - rgb_vals) * amount  # blend toward white
  rgb(pastel_vals[1], pastel_vals[2], pastel_vals[3])
}

# ✅ Apply softening
pastel_palette <- sapply(base_palette, soften_color, amount = 0.5)

# Sort (optional)
pastel_palette <- pastel_palette[order(colSums(col2rgb(pastel_palette)))]
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# Step 4: Generate and save stacked bar plots
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  period_name <- names(files)[i]

  df$Sample <- factor(df$Sample, levels = unique(df$Sample))
  df$Class <- factor(df$Class, levels = names(final_color_map))

  p <- ggplot(df, aes(x = Sample, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = final_color_map, name = "Class") +
    labs(
      title = paste("BLA Distribution -", period_name),
      x = "Sample",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  output_file <- paste0("abundance3_", gsub(" ", "_", tolower(period_name)), ".pdf")
  cat("💾 Saving:", output_file, "\n")
  ggsave(output_file, plot = p, width = 11, height = 6)
}

cat("✅ Done! All plots saved with softened Polychrome pastel colors (no colorspace).\n")









library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(Polychrome)

cat("📂 Reading and processing wide-format datasets...\n")

# Step 1: CSV files in wide format
files <- list(
  "Period 1 River" = "BLA_Period1_River.csv",
  "Period 1 WWTP" = "BLA_Period1_WWTP.csv",
  "Period 2 River" = "BLA_Period2_River.csv",
  "Period 2 WWTP" = "BLA_Period2_WWTP.csv",
  "Period 3 River" = "BLA_Period3_River.csv",
  "Period 3 WWTP" = "BLA_Period3_WWTP.csv"
)

# Step 2: Read and reshape to long format
all_data <- lapply(files, function(file) {
  df <- read.csv(file, check.names = FALSE)
  colnames(df)[1] <- "Sample"

  df_long <- pivot_longer(df, cols = -Sample, names_to = "Class", values_to = "Percent_Abundance") %>%
    filter(!is.na(Percent_Abundance) & Percent_Abundance > 0)

  return(df_long)
})

# Step 3: Unique classes and pastel palette with Polychrome
all_classes <- unique(unlist(lapply(all_data, \(df) unique(df$Class))))

set.seed(42)
pastel_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")
pastel_palette <- pastel_palette[order(colSums(col2rgb(pastel_palette)))]
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# Step 4: Generate and save stacked bar plots
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  period_name <- names(files)[i]

  df$Sample <- factor(df$Sample, levels = unique(df$Sample))
  df$Class <- factor(df$Class, levels = names(final_color_map))

  p <- ggplot(df, aes(x = Sample, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = final_color_map, name = "Class") +
    labs(
      title = paste("BLA Distribution -", period_name),
      x = "Sample",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  output_file <- paste0("abundance2_", gsub(" ", "_", tolower(period_name)), ".pdf")
  cat("💾 Saving:", output_file, "\n")
  ggsave(output_file, plot = p, width = 11, height = 6)
}

cat("✅ Done! All plots saved with consistent Polychrome pastel colors.\n")












# Load required packages
library(ggplot2)
library(reshape2)
library(colorspace)

cat("📂 Reading and processing wide-format datasets...\n")

# === Step 1: List your files ===
files <- list(
  "Period 1" = "BLA_Period1_River.csv",
  "Period 2" = "BLA_Period1_WWTP.csv",
  "Period 3" = "BLA_Period2_River.csv",
  "Period 4" = "BLA_Period2_WWTP.csv",
  "Period 5" = "BLA_Period3_River.csv",
  "Period 6" = "BLA_Period3_WWTP.csv"
)

# === Step 2: Read and reshape each dataset ===
all_data <- lapply(files, function(file) {
  df <- read.csv(file, check.names = FALSE)
  colnames(df)[1] <- "Sample"  # Ensure first column is named 'Sample'

  # Convert from wide to long format
  df_long <- melt(df, id.vars = "Sample", variable.name = "Class", value.name = "Percent_Abundance")

  # Remove zero or NA values
  df_long <- df_long[df_long$Percent_Abundance > 0 & !is.na(df_long$Percent_Abundance), ]

  return(df_long)
})

# === Step 3: Collect all unique Class values and create a pastel color palette ===
all_classes <- unique(unlist(lapply(all_data, function(df) unique(df$Class))))
set.seed(42)
pastel_palette <- createPalette(length(all_classes), seedcolors = "#66C2A5")
pastel_palette <- pastel_palette[order(colSums(col2rgb(pastel_palette)))]
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === Step 4: Generate and save plots ===
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  period_name <- names(files)[i]

  # Preserve sample order
  df$Sample <- factor(df$Sample, levels = unique(df$Sample))
  df$Class <- factor(df$Class, levels = names(final_color_map))

  p <- ggplot(df, aes(x = Sample, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = final_color_map, name = "Class") +
    labs(
      title = paste("Antibiotic Class Distribution -", period_name),
      x = "Sample",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  # Save to PDF
  output_file <- paste0("abundance_", gsub(" ", "_", tolower(period_name)), ".pdf")
  cat("💾 Saving plot:", output_file, "\n")
  ggsave(output_file, plot = p, width = 11, height = 6)
}

cat("✅ All stacked bar charts generated from wide-format CSVs.\n")













# === Load libraries ===
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(Polychrome)

# === STEP 1: Define file list ===
files <- c(
  "Period1_River.csv",
  "Period1_WWTP.csv",
  "Period2_River.csv",
  "Period2_WWTP.csv",
  "Period3_River.csv",
  "Period3_WWTP.csv"
)
names(files) <- c(
  "Period 1 WWTP",
  "Period 1 River",
  "Period 2 WWTP",
  "Period 2 River",
  "Period 3 WWTP",
  "Period 3 River"
)

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and pivot all wide-format datasets ===
all_data <- lapply(files, function(file) {
  df <- read.csv(file)
  colnames(df)[1] <- "Location"  # Rename first column to 'Location'
  df_long <- pivot_longer(df, cols = -Location, names_to = "Class", values_to = "Percent_Abundance")
  df_long <- df_long[df_long$Percent_Abundance > 0 & !is.na(df_long$Class), ]
  df_long
})

# === STEP 3: Build master list of unique classes ===
all_classes <- unique(unlist(lapply(all_data, function(df) unique(as.character(df$Class)))))
cat("🎨 Total unique classes across datasets:", length(all_classes), "\n")

# === STEP 4: Create pastel-style colors with good contrast ===
pastel_palette <- createPalette(length(all_classes), seedcolors = "#88CCEE")
names(pastel_palette) <- all_classes
final_color_map <- pastel_palette

# === STEP 5: Plot and save charts for each period ===
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  period_name <- names(files)[i]

  # Preserve order of sample locations
  df$Location <- factor(df$Location, levels = unique(df$Location))

  # Order stacking by frequency of appearance
  class_order <- df %>%
    group_by(Class) %>%
    summarise(freq = n_distinct(Location)) %>%
    arrange(desc(freq)) %>%
    pull(Class)
  df$Class <- factor(df$Class, levels = class_order)

  # === Create the plot ===
  p <- ggplot(df, aes(x = Location, y = Percent_Abundance, fill = Class)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("Identified MAGs -", period_name),
      x = "Sample Location",
      y = "Percent Abundance"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = final_color_map)

  # === Save to PDF ===
  output_file <- paste0("MAG_distribution_", gsub(" ", "_", tolower(period_name)), ".pdf")
  cat("💾 Saving plot:", output_file, "\n")
  ggsave(output_file, plot = p, width = 10, height = 6)
}

cat("✅ All plots saved with consistent pastel color mapping across 6 periods.\n")

