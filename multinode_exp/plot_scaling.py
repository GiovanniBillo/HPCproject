import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set style
plt.style.use('seaborn')
plt.rcParams['figure.facecolor'] = 'white'

# =============================================
# STRONG SCALING PLOT (Speedup & Efficiency)
# =============================================
def plot_strong_scaling(strong_file):
    df_strong = pd.read_csv(strong_file)
    print("MULTINODE RESULTS STRONG", df_strong)
    nodes = df_strong["nodes"]
    time = df_strong["time"]
    
    # Calculate metrics
    baseline_time = time.iloc[0]
    speedup = baseline_time / time
    efficiency = speedup / nodes
    ideal_speedup = nodes
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Speedup plot (left axis)
    ax1.plot(nodes, speedup, marker='o', color='tab:blue', 
             label="Measured Speedup", linewidth=2, markersize=8)
    ax1.plot(nodes, ideal_speedup, linestyle='--', color='gray', 
             label="Ideal Speedup", linewidth=2)
    ax1.fill_between(nodes, speedup, ideal_speedup,
                   color='red', alpha=0.1, label="Parallelization Gap")
    
    # Speedup labels
    for x, y in zip(nodes, speedup):
        ax1.text(x, y-0.2, f'{y:.2f}', ha='center', va='top',
                color='tab:blue', fontsize=8, 
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.5))
    
    # Efficiency plot (right axis)
    ax2 = ax1.twinx()
    ax2.plot(nodes, efficiency, marker='s', color='tab:green',
             label="Efficiency", linewidth=2, markersize=6)
    
    # Efficiency labels
    for x, y in zip(nodes, efficiency):
        ax2.text(x, y+0.03, f'{y:.2f}', ha='center', va='bottom',
                color='tab:green', fontsize=8,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.5))
    
    # Axis styling
    ax1.set_xlabel("Number of Nodes", fontsize=12)
    ax1.set_ylabel("Speedup", color='tab:blue', fontsize=12)
    ax2.set_ylabel("Efficiency", color='tab:green', fontsize=12)
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax2.tick_params(axis='y', labelcolor='tab:green')
    ax1.grid(True, linestyle=':', alpha=0.7)
    ax1.set_xticks(nodes)
    ax2.set_ylim(0, 1.1)
    
    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    plt.legend(lines1 + lines2, labels1 + labels2,
               loc='upper center', bbox_to_anchor=(0.5, 1.15),
               ncol=3, frameon=False)
    
    plt.title("Strong Scaling Analysis", fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig('strong_scaling.png', dpi=300, bbox_inches='tight')
    plt.close()

# =============================================
# WEAK SCALING PLOT (Time per unit work)
# =============================================
def plot_weak_scaling(weak_file):
    df_weak = pd.read_csv(weak_file)
    nodes = df_weak["nodes"]
    problem_size = df_weak["x"]*df_weak["y"]
    time = df_weak["time"]
    
    # Calculate time per unit work (normalized)
    time_per_unit = time / problem_size
    normalized_time = time_per_unit / time_per_unit.iloc[0]  # Normalize to 1-node case
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Bar plot for normalized time per unit work
    bars = ax.bar(nodes, normalized_time, color='tab:orange', alpha=0.7, width=0.6)
    
    # Add reference line at 1.0 (ideal weak scaling)
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=2, label="Ideal Scaling")
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.2f}',
                ha='center', va='bottom', fontsize=10)
    
    # Axis styling
    ax.set_xlabel("Number of Nodes", fontsize=12)
    ax.set_ylabel("Normalized Time per Unit Work", fontsize=12)
    ax.set_xticks(nodes)
    ax.grid(True, axis='y', linestyle=':', alpha=0.7)
    
    # Add problem size annotations
    for i, (n, size) in enumerate(zip(nodes, problem_size)):
        ax.text(n, 0.02, f"{size:,}", ha='center', va='bottom', 
               rotation=90, fontsize=8, color='white',
               bbox=dict(facecolor='black', alpha=0.5, pad=1))
    
    plt.legend(loc='upper right')
    plt.title("Weak Scaling Analysis\n(Constant Time per Unit Work = Good Scaling)", 
              fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig('weak_scaling.png', dpi=300, bbox_inches='tight')
    plt.close()

# =============================================
# MAIN EXECUTION
# =============================================
if __name__ == "__main__":
    # Update these filenames to match your actual data files
    plot_strong_scaling("multinode_results_STRONG.csv")
    plot_weak_scaling("multinode_results_WEAK.csv")
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np

# def plot_scaling(metrics_file, scaling_type):
#     df = pd.read_csv(metrics_file)
#     nodes = df["nodes"]
#     time = df["time"]
#     speedup = df["speedup"]
#     efficiency = df["efficiency"]
    
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
#     fig.suptitle(f"{scaling_type} Scaling Analysis", fontsize=16)
    
#     # Speedup plot
#     ax1.plot(nodes, speedup, 'bo-', label="Measured")
#     ax1.plot(nodes, nodes if scaling_type=="Strong" else np.ones_like(nodes), 
#              'r--', label="Ideal")
#     ax1.set_xlabel("Number of Nodes", fontsize=12)
#     ax1.set_ylabel("Speedup", fontsize=12)
#     ax1.set_title("Speedup")
#     ax1.grid(True, linestyle=':')
#     ax1.legend()
    
#     # Efficiency plot
#     ax2.plot(nodes, efficiency, 'go-')
#     ax2.axhline(y=1.0, color='r', linestyle='--')
#     ax2.set_xlabel("Number of Nodes", fontsize=12)
#     ax2.set_ylabel("Efficiency", fontsize=12)
#     ax2.set_title("Parallel Efficiency")
#     ax2.grid(True, linestyle=':')
    
#     # Formatting
#     for ax in [ax1, ax2]:
#         ax.set_xticks(nodes)
#         ax.set_xticklabels([f"{int(n)}" for n in nodes])
    
#     plt.tight_layout()
#     plt.savefig(f"{scaling_type.lower()}_scaling.png", dpi=150)
#     plt.show()

# # Generate plots
# plot_scaling("strong_metrics.csv", "Strong")
# plot_scaling("weak_metrics.csv", "Weak")
