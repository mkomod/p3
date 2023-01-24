#!/bin/sh

for S in {1..6}
do
    qsub -v "SIM=$S" simulations.pbs
done

