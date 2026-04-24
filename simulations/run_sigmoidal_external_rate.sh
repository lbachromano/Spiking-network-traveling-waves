#!/usr/bin/env bash

 
n=10
master_seed=235

nu_ext=2000
g=3.625
for ((i=1; i<=n; i++)); do
    # Derive a unique, reproducible seed for each run
    rs=$(( master_seed + i*10 ))
    echo "Running simulation i=$i with seed rs=$rs"
    python Simulations_sigmoidal_external_rate.py "$i" "$rs" "$nu_ext" "$g"
    echo "Waiting 2 minutes before next run..."
    sleep 120
done
 
