#!/bin/bash

# Execution script for cgp3d_parallel.x

# Default parameter values
MAXIT=2000     # Maximum iterations
RTOL=1e-6      # Relative residual tolerance
OMEGA=1.9      # SSOR relaxation parameter
OVERLAP=10     # Schwarz method overlap
VERSION="par"  # Set to parallel version
THREADS=20     # Fixed number of threads to 20

# Arrays of preconditioners, RHS types, and mesh sizes
PRECONDITIONERS=("id" "ssor" "as")
RHS_TYPES=("point" "func" "const" "check" "low" "high")
MESH_SIZES=(100 150 200 250)

# Output file to store results
OUTPUT_FILE="simulation_results.txt"
echo "Simulation Results" > $OUTPUT_FILE
echo "==================" >> $OUTPUT_FILE

# Clean and compile if requested
if [ "$CLEAN" == "true" ]; then
    echo "Cleaning and recompiling the code..."
    make clean
    make || { echo "Compilation failed"; exit 1; }
else
    # Set the executable to the parallel version
    EXECUTABLE="cgp3d_parallel.x"

    # Check if the executable exists
    if [ ! -f "$EXECUTABLE" ]; then
        echo "Executable $EXECUTABLE not found. Compiling..."
        make || { echo "Compilation failed"; exit 1; }
    fi
fi

# Set the number of threads for the parallel version
export OMP_NUM_THREADS=$THREADS
echo "Setting OpenMP threads to $THREADS"

# Loop through each preconditioner, RHS type, and mesh size
for PTYPE in "${PRECONDITIONERS[@]}"; do
    echo "----------------------------------------" >> $OUTPUT_FILE
    echo "Preconditioner: $PTYPE" >> $OUTPUT_FILE
    echo "----------------------------------------" >> $OUTPUT_FILE
    for N in "${MESH_SIZES[@]}"; do
        echo "Mesh Size: $N" >> $OUTPUT_FILE
        echo "----------------------------------------" >> $OUTPUT_FILE
        for RHS in "${RHS_TYPES[@]}"; do
            echo "RHS Type: $RHS" >> $OUTPUT_FILE
            echo "----------------------------------------" >> $OUTPUT_FILE
            # Run the executable with the specified parameters
            # Capture all output including important simulation data
            OUTPUT=$(./$EXECUTABLE -n $N -M $MAXIT -r $RTOL -p $PTYPE -w $OMEGA -o $OVERLAP -t $RHS 2>&1)
            # Write parameters and output to the file
            echo "$OUTPUT" >> $OUTPUT_FILE
            echo "" >> $OUTPUT_FILE
        done
        echo "" >> $OUTPUT_FILE
    done
    echo "" >> $OUTPUT_FILE
done

echo "Simulation completed. Results saved in $OUTPUT_FILE."
