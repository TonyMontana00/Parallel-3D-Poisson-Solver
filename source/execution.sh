#!/bin/bash

# Execution script for cgp3d.x and cgp3d_parallel.x

# Default parameter values
N=100         # Mesh size
MAXIT=2000     # Maximum iterations
RTOL=1e-6     # Relative residual tolerance
PTYPE="id"    # Preconditioner type: id, ssor, as
OMEGA=1.9     # SSOR relaxation parameter
OVERLAP=10    # Schwarz method overlap
VERSION="seq" # Default to sequential version
THREADS=2     # Default to 1 thread for sequential version

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
    echo "  -t <int>      Number of threads for parallel version (default: $THREADS)"
    echo "  -c            Clean and recompile the code"
    exit 1
}

# Parse command-line options
while getopts "hn:M:r:p:w:o:v:t:c" opt; do
    case $opt in
        h) usage ;;
        n) N=$OPTARG ;;
        M) MAXIT=$OPTARG ;;
        r) RTOL=$OPTARG ;;
        p) PTYPE=$OPTARG ;;
        w) OMEGA=$OPTARG ;;
        o) OVERLAP=$OPTARG ;;
        v) VERSION=$OPTARG ;;
        t) THREADS=$OPTARG ;;
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

# Set the number of threads if running the parallel version
if [ "$VERSION" == "par" ]; then
    export OMP_NUM_THREADS=$THREADS
    echo "Setting OpenMP threads to $THREADS"
fi

# Execute the selected program with the specified parameters
# echo "Running $EXECUTABLE with the following parameters:"
# echo "  Mesh size (-n): $N"
# echo "  Maximum iterations (-M): $MAXIT"
# echo "  Relative tolerance (-r): $RTOL"
# echo "  Preconditioner (-p): $PTYPE"
# echo "  SSOR relaxation (-w): $OMEGA"
# echo "  Overlap (-o): $OVERLAP"
# if [ "$VERSION" == "par" ]; then
#     echo "  Number of threads (-t): $THREADS"
# fi
echo

./$EXECUTABLE -n $N -M $MAXIT -r $RTOL -p $PTYPE -w $OMEGA -o $OVERLAP
