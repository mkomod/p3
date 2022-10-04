#!/bin/bash

cd ../r/simulations
export CORES=20

for SIM in {1..12}
do
    export SIM=$SIM
    Rscript 01-simulations.R
done
