#!/usr/bin/env bash

external_input=0.9 #in the paper, value of nu_{ext}/\nu_{theta}
homogeneous=0 # 1 to simulate a network with spatially homogeneous connectivity, 1 to simulate a spatially extended network.
anisotropy_EE=0.7 # rho_EE defines the anisotropy in the E to E connectivity
n=1 # number of simulations
trials=1 #number of trials per simulations (each one simulates 1 second of activity)
master_seed=33  # single seed for reproducibility

for ((i=1; i<=n; i++)); do
    # Derive a unique, reproducible seed for each run
    rs=$(( master_seed + i ))

    echo "Running simulation i=$i with seed rs=$rs"
    python Simulations_const_ext_rate.py "$i" "$rs" "$trials" "$homogeneous" "$external_input" "$anisotropy_EE"
    echo "Waiting 2 minutes before next run..."
    sleep 120
done
