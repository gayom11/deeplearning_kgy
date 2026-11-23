"""
Step 4-A: GNINA 도킹 준비

1) PDB 1M17 (EGFR + Erlotinib) 다운로드
2) 단백질(receptor)과 리간드(crystal ligand) 분리
3) ADMET 상위 300개 화합물을 SDF로 변환
4) GNINA 실행 스크립트 생성
"""
import os
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
ADMET_FILE = os.path.join(OUTPUT_DIR, "..", "04_admet", "admet_top_hits.csv")
PDB_ID = "3OG7"


def download_pdb(pdb_id, output_path):
    """PDB 구조 다운로드"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading PDB {pdb_id}...")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    with open(output_path, "w") as f:
        f.write(resp.text)
    print(f"  Saved: {output_path}")


def extract_receptor_and_ligand(pdb_path, rec_path, lig_path):
    """PDB에서 단백질(ATOM)과 리간드(HETATM AQ4=Erlotinib) 분리"""
    rec_lines = []
    lig_lines = []

    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                rec_lines.append(line)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname == "032":  # Vemurafenib (PLX4032) in 3OG7
                    lig_lines.append(line)
            elif line.startswith("END"):
                rec_lines.append(line)
                lig_lines.append(line)

    with open(rec_path, "w") as f:
        f.writelines(rec_lines)
    print(f"  Receptor: {rec_path} ({len(rec_lines)} lines)")

    with open(lig_path, "w") as f:
        f.writelines(lig_lines)
    print(f"  Crystal ligand: {lig_path} ({len(lig_lines)} lines)")


def compounds_to_sdf(df, output_path, max_compounds=None):
    """SMILES → 3D SDF 변환 (GNINA 입력용)"""
    if max_compounds:
        df = df.head(max_compounds)

    writer = Chem.SDWriter(output_path)
    success = 0
    fail = 0

    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            fail += 1
            continue
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            fail += 1
            continue

        mol.SetProp("_Name", row.get("molecule_chembl_id", f"mol_{success}"))
        mol.SetProp("predicted_active_prob", f"{row.get('predicted_active_prob', 0):.4f}")
        mol.SetProp("admet_score", f"{row.get('admet_score', 0):.4f}")
        writer.write(mol)
        success += 1

    writer.close()
    print(f"  SDF: {output_path} ({success} molecules, {fail} failed)")
    return success


def create_gnina_script(work_dir, n_compounds):
    """GNINA Docker 실행 스크립트 생성"""
    # WSL 경로로 변환
    script = f"""#!/bin/bash
# GNINA 도킹 실행 스크립트
# 사용법: wsl -d Ubuntu -- bash /mnt/c/Users/KKY/deeplearning_kgy/experiment/03_docking/run_gnina.sh

WORK_DIR="/mnt/c/Users/KKY/deeplearning_kgy/experiment/03_docking"
cd "$WORK_DIR"

echo "=== GNINA Docking ==="
echo "Receptor: rec.pdb"
echo "Ligands: ligands.sdf ({n_compounds} compounds)"
echo "Autobox: crystal_ligand.pdb (Erlotinib binding site)"
echo ""

# GNINA 도킹 실행
docker run --rm \\
    -v "$WORK_DIR":/scr \\
    gnina/gnina \\
    gnina \\
    -r /scr/rec.pdb \\
    -l /scr/ligands.sdf \\
    --autobox_ligand /scr/crystal_ligand.pdb \\
    --autobox_add 4 \\
    -o /scr/docked.sdf \\
    --exhaustiveness 8 \\
    --num_modes 3 \\
    --cpu 6 \\
    --seed 42 \\
    2>&1 | tee /scr/gnina_log.txt

echo ""
echo "=== Docking Complete ==="
echo "Output: $WORK_DIR/docked.sdf"
echo "Log: $WORK_DIR/gnina_log.txt"
"""

    script_path = os.path.join(work_dir, "run_gnina.sh")
    with open(script_path, "w", newline="\n") as f:
        f.write(script)
    print(f"  Script: {script_path}")
    return script_path


def main():
    # 1) PDB 다운로드
    pdb_path = os.path.join(OUTPUT_DIR, f"{PDB_ID}.pdb")
    if not os.path.exists(pdb_path):
        download_pdb(PDB_ID, pdb_path)
    else:
        print(f"PDB already exists: {pdb_path}")

    # 2) 단백질/리간드 분리
    rec_path = os.path.join(OUTPUT_DIR, "rec.pdb")
    lig_path = os.path.join(OUTPUT_DIR, "crystal_ligand.pdb")
    extract_receptor_and_ligand(pdb_path, rec_path, lig_path)

    # 3) ADMET 상위 화합물을 SDF로 변환
    print("\nConverting compounds to SDF...")
    df = pd.read_csv(ADMET_FILE)
    sdf_path = os.path.join(OUTPUT_DIR, "ligands.sdf")
    n = compounds_to_sdf(df, sdf_path)

    # 4) GNINA 실행 스크립트 생성
    print("\nCreating GNINA run script...")
    create_gnina_script(OUTPUT_DIR, n)

    print(f"\n=== Docking Preparation Complete ===")
    print(f"Files in {OUTPUT_DIR}:")
    print(f"  rec.pdb           - EGFR receptor (PDB {PDB_ID})")
    print(f"  crystal_ligand.pdb - Erlotinib (binding site reference)")
    print(f"  ligands.sdf       - {n} compounds for docking")
    print(f"  run_gnina.sh      - Docker execution script")
    print(f"\n실행 방법:")
    print(f'  wsl -d Ubuntu -- bash /mnt/c/Users/KKY/deeplearning_kgy/experiment/03_docking/run_gnina.sh')


if __name__ == "__main__":
    main()
