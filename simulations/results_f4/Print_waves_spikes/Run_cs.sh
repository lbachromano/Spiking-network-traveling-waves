#!/bin/bash


for seed in 2
do
    for i in 1 
    do
        python Mazzoni_2D_exp.py "$seed" "$i"
    done
done
 