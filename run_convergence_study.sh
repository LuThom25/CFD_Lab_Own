#!/bin/bash
# Convergence study: 300x300 Lid-Driven Cavity
# Runs all setups sequentially.

BINARY="./build/fluidchen"
CASES="./example_cases/LidDrivenCavity300"

echo "========================================"
echo " WS3 Convergence Study — 300x300 LDC"
echo "========================================"
echo ""

# 1. Serial (no MPI)
echo "[1/8] Serial (no MPI)"
$BINARY "$CASES/LidDrivenCavity300_serial.dat"
echo ""

# 2. MPI (1x1) — 1 rank
echo "[2/8] MPI (1x1) — 1 rank"
mpirun -np 1 $BINARY "$CASES/LidDrivenCavity300_1_1.dat"
echo ""

# 3. MPI (2x2) — 4 ranks
echo "[3/8] MPI (2x2) — 4 ranks"
mpirun -np 4 $BINARY "$CASES/LidDrivenCavity300_2_2.dat"
echo ""

# 4. MPI (1x4) — 4 ranks
echo "[4/8] MPI (1x4) — 4 ranks"
mpirun -np 4 $BINARY "$CASES/LidDrivenCavity300_1_4.dat"
echo ""

# 5. MPI (4x1) — 4 ranks
echo "[5/8] MPI (4x1) — 4 ranks"
mpirun -np 4 $BINARY "$CASES/LidDrivenCavity300_4_1.dat"
echo ""

# 6. MPI (3x2) — 6 ranks
echo "[6/8] MPI (3x2) — 6 ranks"
mpirun -np 6 $BINARY "$CASES/LidDrivenCavity300_3_2.dat"
echo ""

# 7. MPI (2x3) — 6 ranks
echo "[7/8] MPI (2x3) — 6 ranks"
mpirun -np 6 $BINARY "$CASES/LidDrivenCavity300_2_3.dat"
echo ""

# 8. MPI (3x3) — 9 ranks
echo "[8/8] MPI (3x3) — 9 ranks"
mpirun -np 9 $BINARY "$CASES/LidDrivenCavity300_3_3.dat"
echo ""

echo "========================================"
echo " All runs complete. Generating plots..."
echo "========================================"
echo ""

cd "$(dirname "$0")/example_cases/LidDrivenCavity300"
python3 analyze.py
