#!/bin/bash

export CORES=20

for SIM in {1..12}
do
    export SIM=$SIM
    cd ../r/simulations
    Rscript 01-simulations.R
done
