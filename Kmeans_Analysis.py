import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# === USER INPUTS ===
file = "Period_1_River_1109.csv"
k = 10
spatial_order = ['P1-DAM', 'P1-RL-UP', 'P1-RL-DOWN', 'P1-SC-UP', 'P1-CC-UP', 'P1-CC-DOWN', 'P1-S1', 'P1-S2', 'P1-S3']

# === LOAD DATA ===
df = pd.read_csv(file, index_col=0)

# Enforce spatial order if provided
if spatial_order:
    df = df[spatial_order]

# === SCALE DATA ===
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df)

# === RUN K-MEANS ===
kmeans = KMeans(n_clusters=k, random_state=42)
clusters = kmeans.fit_predict(X_scaled)

# === SAVE RESULTS ===
df['cluster'] = clusters
cluster_members = (
    df[['cluster']]
    .reset_index()
    .rename(columns={'index': 'rMAG'})
    .sort_values('cluster')
)
cluster_members.to_csv("p1_river_clusters10.csv", index=False)

# Compute mean abundance per cluster
cluster_means = df.groupby('cluster').mean()
cluster_means.to_csv("p1_river_cluster_mean_patterns10.csv")

print("✅ K-means clustering complete!")
print(f"Saved: rMAG_clusters.csv (assignments) and cluster_mean_patterns.csv (average abundance per cluster)")







import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.cm import get_cmap

# === USER INPUTS ===
cluster_mean_file = "p1_river_cluster_mean_patterns10.csv"
assignments_file = "p1_river_clusters10.csv"
original_data_file = "Period_1_River_1109.csv"
output_plot = "p1_river_cluster_mean_patterns10.png"

# === LOAD DATA ===
cluster_means = pd.read_csv(cluster_mean_file, index_col=0)
assignments = pd.read_csv(assignments_file)
original = pd.read_csv(original_data_file, index_col=0)

# === PLOT 1: MEAN CLUSTER PATTERNS ===
plt.figure(figsize=(10, 6))

# Use a continuous colormap with enough unique colors
num_clusters = len(cluster_means)
cmap = get_cmap('tab20b') if num_clusters <= 40 else get_cmap('turbo')
colors = cmap(np.linspace(0, 1, num_clusters))

# Plot each cluster with its own color
for i, cluster in enumerate(cluster_means.index):
    plt.plot(cluster_means.columns,
             cluster_means.loc[cluster],
             label=f'Cluster {cluster}',
             color=colors[i],
             alpha=0.9)

# === LABELS & STYLE ===
plt.xlabel("Sampling Sites")
plt.ylabel("Abundance (TAD80/GEQ) log-scale")
plt.title("K-means Clusters of rMAG Abundance Patterns - Period 1")
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.tight_layout()

# === SAVE & SHOW ===
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.show()

print(f"✅ Plot saved as '{output_plot}'")
