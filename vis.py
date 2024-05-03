import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.patches import Circle

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--num_sides", type=int, default=3)
args = parser.parse_args()

num_sides = args.num_sides
df = pd.read_csv(f"./plots/{num_sides}-gon.csv")
df[['X', 'Y']] = df['Coordinates'].str.strip('()').str.split(', ', expand=True).astype(float)

# Define the discrete color map for integer values
num_colors = df['#sides'].max() + 1
cmap = plt.cm.get_cmap('RdBu_r')
new_colors = [cmap(i / (num_colors - 1)) for i in range(num_colors)]
new_cmap = LinearSegmentedColormap.from_list("new_cmap", new_colors, N=num_colors)

# Define the boundaries for discrete colors
boundaries = range(num_colors+1)
norm = BoundaryNorm(boundaries, num_colors)

# Calculate the middle positions for ticks
tick_positions = [(i + 0.5) for i in range(num_colors)]

fig, ax = plt.subplots(figsize=(10, 8))
scatter = ax.scatter(df['X'], df['Y'], c=df['#sides'], cmap=new_cmap, s=10, norm=norm)
cbar = plt.colorbar(scatter, label='Number of Sides', ticks=tick_positions)
cbar.set_ticklabels(range(num_colors))
plt.title('2D Scatter Plot of Polygon Sides')

# Remove x and y sides
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Add x and y axes in the middle
ax.axhline(0, color='black',linewidth=0.5)
ax.axvline(0, color='black',linewidth=0.5)

# Move left y-axis and bottom x-axis to centre, passing through (0,0)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')

circle = Circle((0, 0), 1, edgecolor='black', facecolor='none')
ax.add_patch(circle)

ax.axis('equal')  # Set aspect ratio to be equal
plt.show()
plt.savefig(f"./plots/{num_sides}-gon.pdf", dpi=300)
