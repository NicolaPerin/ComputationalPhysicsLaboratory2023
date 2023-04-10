#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=stevejob
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 24
#SBATCH --time=01:00:00

args="10000 100000"

for procs in {1..24..1}
do
    echo "Running with $procs processes"
    for (( rep=1; rep<=5; rep++ ))
    do
        echo "  Repetition $rep"
        mpirun -np $procs ./parallel_walkers $args
    done
done
