#!/bin/bash

# Create new file with header
echo "num_threads,func_id,avg_time,max_time,min_time,imbalance_ratio" > combined_thread_stats.csv

# Combine all thread_stats files, remove headers, and sort numerically
grep -vh "num_threads" thread_stats_* | sort -t, -k1n,1 -k2n,2 >> combined_thread_stats.csv

# Verify output
echo "Created combined_thread_stats.csv with:"
head combined_thread_stats.csv
echo "..."
tail combined_thread_stats.csv
