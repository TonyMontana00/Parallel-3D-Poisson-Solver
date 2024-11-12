#!/bin/bash

# Execution script for cgp3d_parallel.x
# Enhanced for scaling tests by looping over threads from 1 to 40

# Default parameter values
N=100         # Mesh size
MAXIT=200     # Maximum iterations
RTOL=1e-6     # Relative residual tolerance
PTYPE="id"    # Preconditioner type: id, ssor, as
OMEGA=1.9     # SSOR relaxation parameter
OVERLAP=10    # Schwarz method overlap
VERSION="par" # Default to parallel version
CLEAN=false   # Default to not clean

# Function to display help message
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h            Print this help message"
    echo "  -n <int>      Mesh size (default: $N)"
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
while getopts "hn:M:r:p:w:o:v:c" opt; do
    case $opt in
        h) usage ;;
        n) N=$OPTARG ;;
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

# Function to run scaling tests
run_scaling_tests() {
    local executable=$1
    local n=$2
    local maxit=$3
    local rtol=$4
    local ptype=$5
    local omega=$6
    local overlap=$7

    echo "Starting scaling tests for $executable..."

    # Create or clear the results file
    RESULTS_FILE="scaling_results_${executable}.txt"
    echo "" > $RESULTS_FILE

    # Loop over threads from 1 to 40
    for THREAD in {1..40}; do
        export OMP_NUM_THREADS=$THREAD
        echo "Running with $THREAD threads..."

        # Execute the program and capture the output
        OUTPUT=$(./$executable -n $n -M $maxit -r $rtol -p $ptype -w $omega -o $overlap)

        # Append thread count and output to results file
        echo "Threads: $THREAD" >> $RESULTS_FILE
        echo "$OUTPUT" >> $RESULTS_FILE
        echo "----------------------------------------" >> $RESULTS_FILE

        echo "Completed run with $THREAD threads."
    done

    echo "Scaling tests completed. Results saved in $RESULTS_FILE."
}

# Execute scaling tests if running the parallel version
if [ "$VERSION" == "par" ]; then
    echo "Initiating scaling tests for parallel version..."

    run_scaling_tests "$EXECUTABLE" "$N" "$MAXIT" "$RTOL" "$PTYPE" "$OMEGA" "$OVERLAP"

else
    # Set the number of threads to 1 for sequential version
    export OMP_NUM_THREADS=1
    echo "Running sequential version with 1 thread."

    # Execute the selected program with the specified parameters
    echo "Running $EXECUTABLE with the following parameters:"
    echo "  Mesh size (-n): $N"
    echo "  Maximum iterations (-M): $MAXIT"
    echo "  Relative tolerance (-r): $RTOL"
    echo "  Preconditioner (-p): $PTYPE"
    echo "  SSOR relaxation (-w): $OMEGA"
    echo "  Overlap (-o): $OVERLAP"
    echo

    ./$EXECUTABLE -n $N -M $MAXIT -r $RTOL -p $PTYPE -w $OMEGA -o $OVERLAP
fi