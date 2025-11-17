# === Load libraries ===
library(ggplot2)
library(dplyr)
library(scales)
library(Polychrome)
library(readr)

# === STEP 1: Define file list (3 periods) ===
files <- c(
  "Period1_classes.csv",
  "Period2_classes.csv",
  "Period3_classes.csv"
)
names(files) <- c("Period 1", "Period 2", "Period 3")

cat("📂 Reading and processing datasets...\n")

# === STEP 2: Read and preprocess ===
all_data <- lapply(names(files), function(name) {
  file <- files[name]
  df <- read_csv(file, show_col_types = FALSE)
  
  colnames(df) <- trimws(colnames(df))
  
  df <- df %>%
    mutate(
      Location = as.character(Location),
      Percent_Abundance = as.numeric(Percent_Abundance),
      Class = as.character(Class),
      Period = name
    ) %>%
    filter(!is.na(Class), Percent_Abundance > 0)
  
  # Combine same class values within each location
  df <- df %>%
    group_by(Period, Location, Class) %>%
    summarise(Percent_Abundance = sum(Percent_Abundance, na.rm = TRUE), .groups = "drop")
  
  df
})

combined_df <- bind_rows(all_data)

# === STEP 3: Specify your preferred location order here! ===
# 👇👇👇  <<< YOU CAN EDIT THIS VECTOR TO CONTROL ORDER OF LOCATIONS >>>
location_order <- c(
  "S3",
  "S2",
  "S1",
  "CC_DOWN",
  "CC_UP",
  "SC_UP",
  "RL_DOWN",
  "RL_UP",
  "DAM"
)
# ^ Add or remove as needed — they must match the "Location" values in your CSVs

# Set factor levels
combined_df$Location <- factor(combined_df$Location, levels = location_order)

# === STEP 4: Create pastel color palette ===
all_classes <- unique(combined_df$Class)
cat("🎨 Total unique classes across all periods:", length(all_classes), "\n")

base_palette <- createPalette(length(all_classes), seedcolors = "#88CCEE")
lighten_color <- function(color, factor = 0.5) {
  rgb_col <- col2rgb(color) / 255
  lightened <- rgb_col + (1 - rgb_col) * factor
  rgb(lightened[1], lightened[2], lightened[3])
}
pastel_palette <- sapply(base_palette, lighten_color)
names(pastel_palette) <- all_classes

# === STEP 5: Plot ===
p <- ggplot(combined_df, aes(x = Location, y = Percent_Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Period, ncol = 1, scales = "fixed") +
  scale_fill_manual(values = pastel_palette) +
  coord_flip() +  # horizontal orientation
  scale_y_continuous(position = "right") +
  labs(
    title = "rMAG Abundance by Class",
    x = "Sample Location",
    y = "Abundance (TAD80/GEQ)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x.right = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# === STEP 6: Save ===
output_file <- "MAG_distribution_three_periods_horizontal_ordered2.pdf"
cat("💾 Saving combined plot:", output_file, "\n")
ggsave(output_file, plot = p, width = 10, height = 10)

cat("✅ Plot saved with custom location order and pastel colors.\n")




































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
















