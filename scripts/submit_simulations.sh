#!/bin/sh

for S in $(seq 1 $1)
do
    for D in {1..4}
    do
	qsub -v "SIM=$S, DGP=$D" simulations.pbs
    done
done

# for MCMC
# for S in {1..2}
# do
#     qsub -v "SIM=$S" simulations.pbs
# done

# for S in {1..4}
# do
#     qsub -v "SIM=$S" simulations.pbs
# done

# for S in {1..8}
# do
#     qsub -v "SIM=$S" simulations.pbs
# done

# GENERAL SIMS
# for S in {1..6}
# do
#     qsub -v "SIM=$S" simulations.pbs
# done


# for S in {7..11}
# do
#     for D in {1..4}
#     do
# 	qsub -v "SIM=$S, DGP=$D" simulations.pbs
#     done
# done
