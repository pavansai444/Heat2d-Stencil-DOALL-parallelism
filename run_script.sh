#!/bin/bash

# Compile the C programs
GREEN='\033[0;32m'
NC='\033[0m'
echo -e "${GREEN}\nStarting Heat2d\n\n${NC}"
gcc -fopenmp parallel_heat2d.c -o parallel
gcc serial_heat2d.c -o serial

# Define the function to run the program with the given argument
run_parallel() {
    arg=$1
    ./parallel 1e-3 par_${arg}.txt $arg 
}

run_serial(){
    ./serial 1e-3 ser_${1}.txt $1 
}

# Run the program in parallel and for each argument
for arg in 5 50 500 1000 5000; do
    run_serial $arg
    run_parallel $arg  
done

#run gnuplot --persist filename.gp for image visualization
 
echo -e "${GREEN}\nExecution Completed\n${NC}"