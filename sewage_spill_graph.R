# Load required packages
library(dplyr)

# Load in data
sewage_spill_data <- read.csv("Sewage_Spills_Report_Chattahoochee2_Basin_07.12.2018_to_08.22.2018.csv")
sample_site_coords <- read.csv("chattahoochee_river_sites.csv")
chattahoochee_coords <- read.table("chattahoochee.geo", header = TRUE)

sample_site_coords_only <- sample_site_coords %>%
  select(lat, long)


plot(chattahoochee_coords, asp=1, type="n", main="Sampling Sites with Sewage Spills", xlab="Longitude", ylab="Latitude")

# Add a blue line to connect the sites and represent the river
lines(chattahoochee_coords, col="light blue", lwd=3)

# Orient the flow of the river with some extra labels
text(-83.0, 34.5, "Upstream", cex=1.2, col="red")
text(-85.25, 29.8, "Downstream", cex=1.2, col="red")

# Plot sample sites with blue dots
points(sample_site_coords$long, sample_site_coords$lat, col="blue", pch=16, cex=1.5)

# Add names of sample sites to the left of the points
text(sample_site_coords$long - 0.02, sample_site_coords$lat, labels = sample_site_coords$site, pos = 2, cex = 0.8, col = "blue")

# Plot sewage spill data with red dots
points(sewage_spill_data$Long, sewage_spill_data$Lat, col="red", pch=16, cex=1.5)

# Add legend
legend("bottomright", legend = c("Sampling Sites", "Sewage Spill Sites"),
       col = c("blue", "red"), pch = c(16, 16), cex = 1.0)

