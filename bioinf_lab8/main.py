"""DNA Pattern Generation from Gene Promoters (Bioinformatics Lab 8).

Uses PromKappa 2.0 (https://github.com/Gagniuc/PromKappa-2.0-makes-Objective-Digital-Stains):
- Test sequence S verification
- Sliding window (L = 30 bp)
- C+G% and Information Content computation
- Pattern plotting (Objective Digital Stain)
- Center of weight (centroid) calculation
- Validation against PromKappa 2.0 benchmark results
"""

import os
import re
from typing import Dict, List, Any
from dna_pattern import (
    sliding_window_pattern,
    calculate_center_of_weight,
    plot_dna_pattern,
    plot_centers_of_weight,
    load_fasta
)

PROMKAPPA_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "PromKappa-2.0-makes-Objective-Digital-Stains"
)


def load_promkappa_benchmark_results(results_path: str) -> Dict[str, Dict[str, float]]:
    """Loads benchmark measurement results from PromKappa's results.htm."""
    if not os.path.exists(results_path):
        return {}

    benchmark = {}
    with open(results_path, "r", encoding="latin1") as f:
        html = f.read()

    rows = re.findall(r"<tr>(.*?)</tr>", html, re.DOTALL)
    for r in rows[1:]:
        cols = [re.sub(r"<.*?>", "", c).strip() for c in re.findall(r"<td>(.*?)</td>", r, re.DOTALL)]
        if len(cols) >= 3:
            gene = cols[0]
            try:
                kic = float(cols[1])
                cg = float(cols[2])
                benchmark[gene] = {"kic": kic, "cg": cg}
            except ValueError:
                continue

    return benchmark


def run_pipeline():
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(output_dir, exist_ok=True)

    print("Assignment: DNA Pattern Generation from Gene Promoters (Lab 8)")
    print(f"PromKappa 2.0 Repository: {PROMKAPPA_DIR}")

    # Step 1: Define test sequence S
    prefix = "CGGACTGATCTATCT"
    middle = "A" * 27
    suffix = "CGTAGCATCTATCGATCTATCTAGCGATCTA"
    s_test = prefix + middle + suffix
    print(f"\nSequence S length: {len(s_test)} bp")

    # Step 2: Sliding window configuration
    L = 30
    step = 1
    num_windows = len(s_test) - L + 1
    print(f"Window size L = {L} bp, step = {step}, windows = {num_windows}")

    # Steps 3 & 4: Process C+G% and Information Content
    pattern_ic = sliding_window_pattern(s_test, window_size=L, step=step, metric="ic")
    pattern_sigma = sliding_window_pattern(s_test, window_size=L, step=step, metric="sigma")

    print("\nSample window computations:")
    print(f"{'Win':<4} | {'Window (30 bp)':<32} | {'C+G (%)':<8} | {'IC (%)':<8} | {'sigma(S)':<8}")
    for i in list(range(5)) + list(range(len(pattern_ic) - 5, len(pattern_ic))):
        w = pattern_ic[i]["window"]
        cg = pattern_ic[i]["cg"]
        ic_val = pattern_ic[i]["ic"]
        sig_val = pattern_sigma[i]["ic"]
        print(f"{i+1:<4} | {w:<32} | {cg:<8.2f} | {ic_val:<8.2f} | {sig_val:<8.2f}")

    # Steps 5 & 6: Center of weight and 2D pattern plotting
    cw_cg_ic, cw_ic_val = calculate_center_of_weight(pattern_ic)
    cw_cg_sig, cw_sig_val = calculate_center_of_weight(pattern_sigma)

    print(f"\nCenter of weight (PromKappa IC): C+G% = {cw_cg_ic:.2f}%, IC = {cw_ic_val:.2f}%")
    print(f"Center of weight (Formula mode): C+G% = {cw_cg_sig:.2f}%, sigma = {cw_sig_val:.2f}%")

    test_pattern_path_ic = os.path.join(output_dir, "01_test_sequence_pattern_promkappa_ic.png")
    test_pattern_path_sig = os.path.join(output_dir, "01_test_sequence_pattern_sigma.png")

    plot_dna_pattern(
        pattern_ic,
        title="DNA Pattern: Test Sequence S (PromKappa IC)",
        center_of_weight=(cw_cg_ic, cw_ic_val),
        metric_label="Information Content (IC %)",
        save_path=test_pattern_path_ic
    )

    plot_dna_pattern(
        pattern_sigma,
        title="DNA Pattern: Test Sequence S (sigma(S) = 100 - IC)",
        center_of_weight=(cw_cg_sig, cw_sig_val),
        metric_label="Information Content sigma(S) = 100 - IC (%)",
        save_path=test_pattern_path_sig
    )

    # Steps 7 & 8: Promoter analysis from PromKappa 2.0 repository
    fasta_path = os.path.join(PROMKAPPA_DIR, "gene_promoters_complete.fa")
    results_path = os.path.join(PROMKAPPA_DIR, "results.htm")

    if not os.path.exists(fasta_path):
        fasta_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gene_promoters_complete.fa")

    all_promoters = load_fasta(fasta_path)
    print(f"\nLoaded {len(all_promoters)} gene promoters from {os.path.basename(fasta_path)}")

    benchmark_data = load_promkappa_benchmark_results(results_path)
    if benchmark_data:
        print(f"Loaded {len(benchmark_data)} benchmark records from PromKappa results.htm")

    selected_genes = [
        "SLC30A8", "FTO_1", "CTLA4_1", "INS_1", "NOC2L_1",
        "PLEKHN1_1", "HES4_1", "ISG15_1", "SEPT6_1", "ABL1_1"
    ]

    promoters = {}
    for g in selected_genes:
        if g in all_promoters:
            promoters[g] = all_promoters[g]["sequence"]
        else:
            base = g.split("_")[0]
            for k, v in all_promoters.items():
                if k == base or k.startswith(base + "_"):
                    promoters[k] = v["sequence"]
                    break

    promoter_centroids = []
    print(f"\nProcessing {len(promoters)} promoter patterns (L = 30 bp):")
    print(f"{'Gene':<12} | {'Length':<7} | {'Windows':<7} | {'Mean CG (%)':<11} | {'Mean IC (%)':<11}")

    for i, (gene_name, seq) in enumerate(promoters.items()):
        pts = sliding_window_pattern(seq, window_size=L, step=1, metric="ic")
        mcg, mic = calculate_center_of_weight(pts)
        promoter_centroids.append({
            "name": gene_name,
            "mean_cg": mcg,
            "mean_ic": mic,
            "length": len(seq),
            "num_windows": len(pts)
        })
        print(f"{gene_name:<12} | {len(seq):<7} | {len(pts):<7} | {mcg:<11.2f} | {mic:<11.2f}")

        if i < 4:
            gene_plot_path = os.path.join(output_dir, f"02_pattern_{gene_name}.png")
            plot_dna_pattern(
                pts,
                title=f"DNA Pattern: {gene_name} (L = {L} bp)",
                center_of_weight=(mcg, mic),
                metric_label="Information Content (IC %)",
                save_path=gene_plot_path
            )

    comp_plot_path = os.path.join(output_dir, "03_promoters_centers_of_weight.png")
    plot_centers_of_weight(
        promoter_centroids,
        title="PromKappa DNA Patterns - Centers of Weight Across Gene Promoters",
        metric_label="Average Information Content (IC %)",
        save_path=comp_plot_path
    )

    # Step 8 Verification: Cross-validation against PromKappa results.htm
    if benchmark_data:
        print("\nCross-Validation against PromKappa 2.0 Benchmark Table (results.htm):")
        print(f"{'Gene':<10} | {'PromKappa CG':<12} | {'Python CG':<10} | {'PromKappa IC':<12} | {'Python IC':<10}")
        tested_count = 0
        for gene_name, b_val in benchmark_data.items():
            matched_seq = None
            for candidate in [gene_name, f"{gene_name}_1", f"{gene_name}_2"]:
                if candidate in all_promoters:
                    matched_seq = all_promoters[candidate]["sequence"]
                    break

            if matched_seq:
                pts_sens = sliding_window_pattern(
                    matched_seq, window_size=L, step=1, metric="ic",
                    cg_mode="sequence_sensitive", include_endpoint=False
                )
                py_cg, py_ic = calculate_center_of_weight(pts_sens)
                print(f"{gene_name:<10} | {b_val['cg']:<12.2f} | {py_cg:<10.2f} | {b_val['kic']:<12.2f} | {py_ic:<10.2f}")
                tested_count += 1
                if tested_count >= 8:
                    break

    print(f"\nOutputs saved to: {output_dir}")


if __name__ == "__main__":
    run_pipeline()
