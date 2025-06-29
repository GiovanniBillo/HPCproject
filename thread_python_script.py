import pandas as pd
import matplotlib.pyplot as plt

# Load speedup data from CSV
df = pd.read_csv("speedup.csv")

# Extract thread count and measured speedup
threads = df["threads"]
measured_speedup = df["speedup"]

# Compute ideal speedup line: S(p) = p
ideal_speedup = threads  # same as x = y for ideal linear scaling

# Plot
plt.figure(figsize=(8, 5))
plt.plot(threads, measured_speedup, marker='o', label="Measured Speedup", linewidth=2)
plt.plot(threads, ideal_speedup, linestyle='--', color='gray', label="Ideal Speedup (Linear)", linewidth=2)

# Style and labels
plt.title("Speedup vs Number of Threads", fontsize=14)
plt.xlabel("Threads", fontsize=12)
plt.ylabel("Speedup", fontsize=12)
plt.grid(True, which='both', linestyle=':', linewidth=0.5)
plt.legend(loc="upper left")
plt.xticks(threads)  # ensures all thread counts appear as x-ticks

# Save and show
plt.tight_layout()
plt.savefig("speedup_plot.png", dpi=150)
plt.show()

