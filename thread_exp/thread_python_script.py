import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("thread_speedup.csv")
threads = df["threads"]
measured_speedup = df["speedup"]
efficiency = df["efficiency"]
ideal_speedup = threads
gap = ideal_speedup - measured_speedup

# Create figure
fig, ax1 = plt.subplots(figsize=(10, 6))

# Speedup (Left Y-axis)
ax1.plot(threads, measured_speedup, marker='o', color='tab:blue', 
         label="Measured Speedup", linewidth=2, markersize=8)
ax1.plot(threads, ideal_speedup, linestyle='--', color='gray', 
         label="Ideal Speedup", linewidth=2)
ax1.fill_between(threads, measured_speedup, ideal_speedup, 
                color='red', alpha=0.1, label="Parallelization Gap")

# Add speedup labels with adjusted positions
for x, y in zip(threads, measured_speedup):
    ax1.text(x, y-0.2, f'{y:.2f}', ha='center', va='top', 
             color='tab:blue', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.5))

# Efficiency (Right Y-axis)
ax2 = ax1.twinx()
ax2.plot(threads, efficiency, marker='s', color='tab:green', 
         label="Efficiency", linewidth=2, markersize=6)

# Add efficiency labels with adjusted positions
for x, y in zip(threads, efficiency):
    ax2.text(x, y+0.03, f'{y:.2f}', ha='center', va='bottom', 
             color='tab:green', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.5))

# Axis styling
ax1.set_xlabel("Number of Threads", fontsize=12)
ax1.set_ylabel("Speedup", color='tab:blue', fontsize=12)
ax2.set_ylabel("Efficiency", color='tab:green', fontsize=12)
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:green')
ax1.grid(True, linestyle=':', alpha=0.7)
ax1.set_xticks(threads)
ax2.set_ylim(0, 1.1)

# Legend - centered above plot
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
plt.legend(lines1 + lines2, labels1 + labels2, 
           loc='upper center', bbox_to_anchor=(0.5, 1.15),
           ncol=3, frameon=False)

plt.title("Parallel Performance Analysis", fontsize=14, pad=20)
plt.tight_layout()
plt.savefig("speedup_efficiency_plot_v2.png", dpi=150, bbox_inches='tight')
plt.show()
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np

# # Load speedup data from CSV
# df = pd.read_csv("thread_speedup.csv")

# # Extract data
# threads = df["threads"]
# measured_speedup = df["speedup"]
# efficiency = df["efficiency"]
# ideal_speedup = threads  # S(p) = p for ideal linear scaling
# gap = ideal_speedup - measured_speedup  # Calculate gap between ideal and measured

# # Create figure and primary axis for speedup
# fig, ax1 = plt.subplots(figsize=(10, 6))

# # Plot speedup data (left Y-axis)
# ax1.plot(threads, measured_speedup, marker='o', color='tab:blue', 
#          label="Measured Speedup", linewidth=2, markersize=8)
# ax1.plot(threads, ideal_speedup, linestyle='--', color='gray', 
#          label="Ideal Speedup", linewidth=2)

# # Plot gap as shaded area
# ax1.fill_between(threads, measured_speedup, ideal_speedup, 
#                 color='red', alpha=0.1, label="Parallelization Gap")

# # Add data point labels for speedup
# for x, y in zip(threads, measured_speedup):
#     ax1.text(x, y, f'{y:.2f}', ha='center', va='bottom', color='tab:blue')

# # Style primary axis
# ax1.set_xlabel("Number of Threads", fontsize=12)
# ax1.set_ylabel("Speedup", color='tab:blue', fontsize=12)
# ax1.tick_params(axis='y', labelcolor='tab:blue')
# ax1.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.7)
# ax1.set_xticks(threads)

# # Create secondary axis for efficiency
# ax2 = ax1.twinx()
# ax2.plot(threads, efficiency, marker='s', color='tab:green', 
#         label="Efficiency", linewidth=2, markersize=6)

# # Add data point labels for efficiency
# for x, y in zip(threads, efficiency):
#     ax2.text(x, y, f'{y:.2f}', ha='center', va='top', color='tab:green')

# # Style secondary axis
# ax2.set_ylabel("Efficiency", color='tab:green', fontsize=12)
# ax2.tick_params(axis='y', labelcolor='tab:green')
# ax2.set_ylim(0, 1.1)  # Efficiency typically between 0-1

# # Combine legends from both axes
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# plt.title("Parallel Performance: Speedup and Efficiency", fontsize=14, pad=20)
# plt.tight_layout()

# # Save and show
# plt.savefig("speedup_efficiency_plot.png", dpi=150, bbox_inches='tight')
# plt.show()

