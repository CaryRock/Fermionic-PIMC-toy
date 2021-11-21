#!/bin/sh

# t controls the number of beads
for t in 25 50 75 100; do
    for y in 3000 4000 5000 6000 7000; do
        for z in 50 75 100 125 150; do
            julia 2-julia_sho.jl 0.50 1 $t $y $z 100 &
            julia 2-julia_sho.jl 0.75 1 $t $y $z 100 &
            julia 2-julia_sho.jl 1.00 1 $t $y $z 100 &
            julia 2-julia_sho.jl 1.25 1 $t $y $z 100 &
            julia 2-julia_sho.jl 1.50 1 $t $y $z 100;
        done;
    done;
done;

