#!/bin/bash

cd ../r/simulations
export CORES=20

for SIM in {1..15}
do
    export SIM=$SIM
    Rscript 97-bmin.R
done
