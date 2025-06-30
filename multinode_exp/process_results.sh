#!/bin/bash

# Process STRONG scaling results (fixed problem size, measure time vs nodes)
awk -F, '
  BEGIN {
    print "nodes,time,speedup,efficiency"
  }
  NR==1 {next}  # skip header
  NR==2 {
    baseline_time = $4  # time for 1 node (sequential reference)
  }
  NR>=2 {
    speedup = baseline_time / $4
    efficiency = speedup / $1
    printf("%d,%.2f,%.2f,%.2f\n", $1, $4, speedup, efficiency)
  }
' multinode_results_STRONG.csv > strong_metrics.csv

# Alternative processing showing relative time growth
awk -F, '
  BEGIN {
    print "nodes,problem_size,time,expected_time,relative_time"
  }
  NR==1 {next}
  NR==2 {
    baseline_time = $4
    baseline_size = $2 * $3
  }
  {
    current_size = $2 * $3
    size_ratio = current_size / baseline_size
    expected_time = baseline_time * size_ratio / $1  # Account for both size and node increase
    relative_time = $4 / expected_time
    printf("%d,%d,%.2f,%.2f,%.2f\n", $1, current_size, $4, expected_time, relative_time)
  }
' multinode_results_WEAK.csv > weak_metrics.csv
# # Process WEAK scaling results (fixed time per node, measure problem size)
# awk -F, '
#   BEGIN {
#     print "nodes,time,speedup,efficiency"
#   }
#   NR==1 {next}  # skip header
#   NR==2 {
#     baseline_time = $4  # time for 1 node (reference)
#   }
#   NR>=2 {
#     # Weak scaling: ideal time should equal baseline (constant)
#     speedup = baseline_time / $4
#     efficiency = speedup / $1
#     printf("%d,%.2f,%.2f,%.2f\n", $1, $4, speedup, efficiency)
#   }
# ' multinode_results_WEAK.csv > weak_metrics.csv
