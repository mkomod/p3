#!/bin/sh

qsub simulations.pbs

for S in {1..12}
do
    qsub -v "SIM=$S" simulations.pbs
done
