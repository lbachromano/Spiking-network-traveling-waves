#!/bin/bash 


for seed in 1 2 3 4 5
do
    for i in 1 2 3 4 5 6 7
    do
        python Mazzoni_2D_exp.py $i $seed
    done
done

 

 