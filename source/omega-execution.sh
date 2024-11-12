#!/bin/bash

# Execution script for cgp3d.x and cgp3d_parallel.x
# This script performs simulations varying omega from 0 to 2 with step size 0.1
# for mesh sizes 100, 150, 200, 250 using the Schwarz preconditioner.

# Default parameter values
MAXIT=20000     # Maximum iterations
RTOL=1e-6     # Relative residual tolerance
PTYPE="as"    # Preconditioner type: as (Schwarz method)
OVERLAP=10    # Schwarz method overlap
VERSION="par" # Use the parallel version
THREADS=1    # Use 20 threads for parallel version

# Mesh sizes to test
MESH_SIZES=(200)

# Omega values from 0 to 2 in steps of 0.1
OMEGA_VALUES=( 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9)

# Output file
OUTPUT_FILE="omega_variation_results.txt"

# Clean and compile the code if necessary
EXECUTABLE="cgp3d_parallel.x"

if [ ! -f "$EXECUTABLE" ]; then
    echo "Executable $EXECUTABLE not found. Compiling..."
    make clean
    make || { echo "Compilation failed"; exit 1; }
fi

# Set the number of threads
export OMP_NUM_THREADS=$THREADS
echo "Setting OpenMP threads to $THREADS"

# Start the simulations
echo "Starting simulations..."
echo "Results will be saved in $OUTPUT_FILE"
echo "" > $OUTPUT_FILE  # Clear the output file

for N in "${MESH_SIZES[@]}"; do
    echo "----------------------------------------" >> $OUTPUT_FILE
    echo "Mesh size: $N" >> $OUTPUT_FILE
    echo "----------------------------------------" >> $OUTPUT_FILE
    for OMEGA in "${OMEGA_VALUES[@]}"; do
        echo "Running with omega = $OMEGA and mesh size = $N"
        
        # Execute the program and capture the output
        OUTPUT=$(./$EXECUTABLE -n $N -M $MAXIT -r $RTOL -p $PTYPE -w $OMEGA -o $OVERLAP)
        
        # Append omega value and output to the results file
        echo "Omega: $OMEGA" >> $OUTPUT_FILE
        echo "$OUTPUT" >> $OUTPUT_FILE
        echo "----------------------------------------" >> $OUTPUT_FILE
    done
done

echo "Simulations completed. Results are saved in $OUTPUT_FILE."
