"""
audit.csv 생성 — 전체 파이프라인 추적

각 화합물이 어느 단계에서 통과/탈락했는지 기록.
교수 확인 사항: "ID와 수가 맞는가", "한 ID로 역추적 가능한가"
"""
import os
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))


def main():
    # 1. 스크리닝 라이브러리 원본
    screen_raw = pd.read_csv(os.path.join(BASE, "00_data_collection", "screening_library.csv"))
    all_ids = set(screen_raw["molecule_chembl_id"])

    # 2. Lipinski+PAINS 통과
    screen_filt = pd.read_csv(os.path.join(BASE, "01_preprocessing", "screening_filtered.csv"))
    lipinski_pass = set(screen_filt["molecule_chembl_id"])

    # 3. Scaffold 선별
    scaffold_file = os.path.join(BASE, "01_preprocessing", "screening_scaffold_selected.csv")
    if os.path.exists(scaffold_file):
        scaffold_sel = pd.read_csv(scaffold_file)
        scaffold_pass = set(scaffold_sel["molecule_chembl_id"])
    else:
        scaffold_pass = lipinski_pass

    # 4. ChemBERTa 상위 1000
    chemberta_file = os.path.join(BASE, "02_model_training", "screening_top_hits.csv")
    if os.path.exists(chemberta_file):
        chemberta = pd.read_csv(chemberta_file)
        chemberta_pass = set(chemberta["molecule_chembl_id"])
    else:
        chemberta_pass = set()

    # 5. ADMET 상위 300
    admet_file = os.path.join(BASE, "04_admet", "admet_top_hits.csv")
    if os.path.exists(admet_file):
        admet = pd.read_csv(admet_file)
        admet_pass = set(admet["molecule_chembl_id"])
    else:
        admet_pass = set()

    # 6. GNINA 도킹 완료
    ranked_file = os.path.join(BASE, "05_final_candidates", "all_ranked.csv")
    if os.path.exists(ranked_file):
        ranked = pd.read_csv(ranked_file)
        docked_pass = set(ranked["molecule_chembl_id"])
    else:
        docked_pass = set()

    # 7. 최종 10개
    final_file = os.path.join(BASE, "05_final_candidates", "final_candidates.csv")
    if os.path.exists(final_file):
        final = pd.read_csv(final_file)
        final_pass = set(final["molecule_chembl_id"])
    else:
        final_pass = set()

    # audit 테이블 생성
    rows = []
    for cid in sorted(all_ids):
        step_lipinski = cid in lipinski_pass
        step_scaffold = cid in scaffold_pass
        step_chemberta = cid in chemberta_pass
        step_admet = cid in admet_pass
        step_docked = cid in docked_pass
        step_final = cid in final_pass

        # 탈락 단계 판정
        if step_final:
            status = "FINAL_CANDIDATE"
            exclude_step = ""
            exclude_reason = ""
        elif step_docked:
            status = "DOCKED_NOT_FINAL"
            exclude_step = "05_final_selection"
            exclude_reason = "lower_combined_score_or_similar_structure"
        elif step_admet:
            status = "ADMET_PASS"
            exclude_step = "04_docking"
            exclude_reason = "docking_failure_or_low_affinity"
        elif step_chemberta:
            status = "CHEMBERTA_HIT"
            exclude_step = "03_admet"
            exclude_reason = "admet_score_below_threshold"
        elif step_scaffold:
            status = "SCAFFOLD_SELECTED"
            exclude_step = "02_chemberta"
            exclude_reason = "predicted_active_prob_below_0.802"
        elif step_lipinski:
            status = "LIPINSKI_PASS"
            exclude_step = "01b_scaffold"
            exclude_reason = "scaffold_group_limit_exceeded"
        else:
            status = "EXCLUDED_EARLY"
            exclude_step = "01_lipinski_pains"
            exclude_reason = "ro5_violation_or_pains"

        rows.append({
            "molecule_chembl_id": cid,
            "step1_lipinski": step_lipinski,
            "step1b_scaffold": step_scaffold,
            "step2_chemberta": step_chemberta,
            "step3_admet": step_admet,
            "step4_docking": step_docked,
            "step5_final": step_final,
            "status": status,
            "exclude_step": exclude_step,
            "exclude_reason": exclude_reason,
        })

    df_audit = pd.DataFrame(rows)
    out = os.path.join(BASE, "audit.csv")
    df_audit.to_csv(out, index=False)

    # 단계별 집계
    print("=== Pipeline Audit Summary ===")
    print(f"Total screening library:   {len(all_ids):>6,}")
    print(f"Step 1  Lipinski+PAINS:    {df_audit['step1_lipinski'].sum():>6,}")
    print(f"Step 1b Scaffold select:   {df_audit['step1b_scaffold'].sum():>6,}")
    print(f"Step 2  ChemBERTa top:     {df_audit['step2_chemberta'].sum():>6,}")
    print(f"Step 3  ADMET pass:        {df_audit['step3_admet'].sum():>6,}")
    print(f"Step 4  Docking complete:  {df_audit['step4_docking'].sum():>6,}")
    print(f"Step 5  FINAL:             {df_audit['step5_final'].sum():>6,}")
    print(f"\nStatus breakdown:")
    print(df_audit["status"].value_counts().to_string())
    print(f"\nOutput: {out}")


if __name__ == "__main__":
    main()
