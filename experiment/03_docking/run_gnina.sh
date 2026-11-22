#!/bin/bash
# GNINA 도킹 실행 스크립트
# 사용법: wsl -d Ubuntu -- bash /mnt/c/Users/KKY/deeplearning_kgy/experiment/03_docking/run_gnina.sh

WORK_DIR="/mnt/c/Users/KKY/deeplearning_kgy/experiment/03_docking"
cd "$WORK_DIR"

echo "=== GNINA Docking ==="
echo "Receptor: rec.pdb"
echo "Ligands: ligands.sdf (300 compounds)"
echo "Autobox: crystal_ligand.pdb (Erlotinib binding site)"
echo ""

# GNINA 도킹 실행
docker run --rm \
    -v "$WORK_DIR":/scr \
    gnina/gnina \
    gnina \
    -r /scr/rec.pdb \
    -l /scr/ligands.sdf \
    --autobox_ligand /scr/crystal_ligand.pdb \
    --autobox_add 4 \
    -o /scr/docked.sdf \
    --exhaustiveness 8 \
    --num_modes 3 \
    --cpu 6 \
    --seed 42 \
    2>&1 | tee /scr/gnina_log.txt

echo ""
echo "=== Docking Complete ==="
echo "Output: $WORK_DIR/docked.sdf"
echo "Log: $WORK_DIR/gnina_log.txt"
