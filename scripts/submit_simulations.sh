#!/bin/sh

qsub simulations.pbs

# for D in 1 2 3
# do
#     for S in 1 2 3 4 5 6
#     do
# 	qsub -v "SIM=$S" simulations.pbs
#     done
# done
