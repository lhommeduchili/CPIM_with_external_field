#!/bin/bash

# compile CPIM code
echo "compiling..."
make

# check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "compilation successful"

    # run several instances of the simulation (replicas)
    REPLICAS=10
    
    echo "running $REPLICAS instances of the simulation..."
    
    for i in $(seq 1 $REPLICAS); do
        echo "running instance $i..."
        ./CPIM_external_field "CPIM_extinction_rate_swep_L=256_$i.txt" &
    done
else
    echo "compilation failed"
    exit 1
fi
