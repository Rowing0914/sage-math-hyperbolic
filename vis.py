from sage.all import *
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--num_sides", type=int, default=3)
args = parser.parse_args()

num_sides = args.num_sides
df = pd.read_csv(f"{num_sides}-gon.csv")
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
c = circle((0, 0), 1)
c_matplotlib = c.matplotlib(figure=fig, sub=ax)
ax.axis('equal')  # Set aspect ratio to be equal
plt.show()
plt.savefig(f"{num_sides}-gon.png")
