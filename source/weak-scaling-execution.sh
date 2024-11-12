#!/bin/bash

# Execution script for cgp3d_parallel.x
# Enhanced for weak scaling tests by looping over threads from 1 to 40
# and adjusting the problem size accordingly

# Default parameter values
BASE_N=100     # Base mesh size for 1 thread
MAXIT=200      # Maximum iterations
RTOL=1e-6      # Relative residual tolerance
PTYPE="id"     # Preconditioner type: id, ssor, as
OMEGA=1.9      # SSOR relaxation parameter
OVERLAP=10     # Schwarz method overlap
VERSION="par"  # Default to parallel version
CLEAN=false    # Default to not clean

# Function to display help message
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h            Print this help message"
    echo "  -M <int>      Maximum iteration count (default: $MAXIT)"
    echo "  -r <float>    Relative residual tolerance (default: $RTOL)"
    echo "  -p <string>   Preconditioner: as, ssor, or id (default: $PTYPE)"
    echo "  -w <float>    SSOR relaxation parameter (default: $OMEGA)"
    echo "  -o <int>      Schwarz method overlap (default: $OVERLAP)"
    echo "  -v <string>   Version to run: 'seq' or 'par' (default: $VERSION)"
    echo "  -c            Clean and recompile the code"
    exit 1
}

# Parse command-line options
while getopts "hM:r:p:w:o:v:c" opt; do
    case $opt in
        h) usage ;;
        M) MAXIT=$OPTARG ;;
        r) RTOL=$OPTARG ;;
        p) PTYPE=$OPTARG ;;
        w) OMEGA=$OPTARG ;;
        o) OVERLAP=$OPTARG ;;
        v) VERSION=$OPTARG ;;
        c) CLEAN=true ;;
        *) usage ;;
    esac
done

# Clean and compile if requested
if [ "$CLEAN" == "true" ]; then
    echo "Cleaning and recompiling the code..."
    make clean
    make || { echo "Compilation failed"; exit 1; }
else
    # Check if the selected executable exists
    EXECUTABLE="cgp3d.x"
    if [ "$VERSION" == "par" ]; then
        EXECUTABLE="cgp3d_parallel.x"
    fi

    if [ ! -f "$EXECUTABLE" ]; then
        echo "Executable $EXECUTABLE not found. Compiling..."
        make || { echo "Compilation failed"; exit 1; }
    fi
fi

# Function to run weak scaling tests
run_weak_scaling_tests() {
    local executable=$1
    local maxit=$2
    local rtol=$3
    local ptype=$4
    local omega=$5
    local overlap=$6

    echo "Starting weak scaling tests for $executable..."

    # Create or clear the results file
    RESULTS_FILE="weak_scaling_results_${executable}.txt"
    echo "" > $RESULTS_FILE

    # Loop over threads from 1 to 40, adjust problem size for weak scaling
    for THREAD in {1..40}; do
        export OMP_NUM_THREADS=$THREAD

        # Adjust problem size n for weak scaling
        N=$(echo "$BASE_N * e(l($THREAD)/3)" | bc -l | awk '{printf "%.0f", $1}')

        echo "Running with $THREAD threads and problem size $N..."

        # Execute the program and capture the output
        OUTPUT=$(./$executable -n $N -M $maxit -r $rtol -p $ptype -w $omega -o $overlap)

        # Append thread count, problem size, and output to results file
        echo "Threads: $THREAD, Problem size: $N" >> $RESULTS_FILE
        echo "$OUTPUT" >> $RESULTS_FILE
        echo "----------------------------------------" >> $RESULTS_FILE

        echo "Completed run with $THREAD threads and problem size $N."
    done

    echo "Weak scaling tests completed. Results saved in $RESULTS_FILE."
}

# Execute weak scaling tests if running the parallel version
if [ "$VERSION" == "par" ]; then
    echo "Initiating weak scaling tests for parallel version..."

    run_weak_scaling_tests "$EXECUTABLE" "$MAXIT" "$RTOL" "$PTYPE" "$OMEGA" "$OVERLAP"

else
    # Set the number of threads to 1 for sequential version
    export OMP_NUM_THREADS=1
    echo "Running sequential version with 1 thread."

    # Execute the selected program with the specified parameters
    echo "Running $EXECUTABLE with the following parameters:"
    echo "  Mesh size (-n): $BASE_N"
    echo "  Maximum iterations (-M): $MAXIT"
    echo "  Relative tolerance (-r): $RTOL"
    echo "  Preconditioner (-p): $PTYPE"
    echo "  SSOR relaxation (-w): $OMEGA"
    echo "  Overlap (-o): $OVERLAP"
    echo

    ./$EXECUTABLE -n $BASE_N -M $MAXIT -r $RTOL -p $PTYPE -w $OMEGA -o $OVERLAP
fi
