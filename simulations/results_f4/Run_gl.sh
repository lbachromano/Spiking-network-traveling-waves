#!/bin/bash


cd Isotropic_g3.3

for seed in 5
do
    for i in 1 2 3 4 5
    do
        python Simulate_dynamics.py "$seed" "$i"
    done
done

 