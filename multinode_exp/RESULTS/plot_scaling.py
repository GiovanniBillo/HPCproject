import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set style
# plt.style.use('seaborn')
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

def plot_weak_metrics(weak_metrics_file):
    df = pd.read_csv(weak_metrics_file)
    nodes = df["nodes"]
    problem_size = df["problem_size"] / 1e6  # Convert to millions for readability
    relative_time = df["relative_time"]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Bar plot for relative time
    bars = ax.bar(nodes, relative_time, 
                 color=['tab:green' if x <= 1.05 else 'tab:orange' for x in relative_time],
                 alpha=0.7, width=0.6)
    
    # Perfect scaling line and ±5% tolerance band
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=2, label="Perfect Scaling")
    ax.axhspan(0.95, 1.05, color='green', alpha=0.1, label="±5% Tolerance")
    
    # Value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.2f}',
                ha='center', va='bottom' if height > 1 else 'top',
                fontsize=10, bbox=dict(facecolor='white', alpha=0.7, pad=1))
    
    # Problem size annotations
    for n, size in zip(nodes, problem_size):
        ax.text(n, 0.02, f"{size:.0f}M", ha='center', va='bottom',
               fontsize=8, color='black',
               bbox=dict(facecolor='white', alpha=0.7, pad=1))
    
    # Axis styling
    ax.set_xlabel("Number of Nodes", fontsize=12)
    ax.set_ylabel("Relative Time (Actual/Expected)", fontsize=12)
    ax.set_xticks(nodes)
    ax.set_ylim(0, max(relative_time)*1.15)
    ax.grid(True, axis='y', linestyle=':', alpha=0.7)
    
    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='tab:green', alpha=0.7, label='Good Scaling (≤1.05)'),
        Patch(facecolor='tab:orange', alpha=0.7, label='Suboptimal Scaling (>1.05)'),
        plt.Line2D([0], [0], color='gray', linestyle='--', label='Perfect Scaling'),
        Patch(facecolor='green', alpha=0.1, label='±5% Tolerance')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.title("Weak Scaling: Relative Time to Solution\n"
             "(Values close to 1.0 indicate good scaling)", 
             fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig('weak_scaling_relative.png', dpi=300, bbox_inches='tight')
    plt.close()

# =============================================
# MAIN EXECUTION
# =============================================
if __name__ == "__main__":
    # Update these filenames to match your actual data files

    # plot_strong_scaling("multinode_results_STRONG.csv")
    # plot_weak_metrics("weak_metrics.csv")
    plot_strong_scaling("multinode_results_STRONG_NEW.csv")
    plot_weak_metrics("weak_metrics_NEW.csv")

