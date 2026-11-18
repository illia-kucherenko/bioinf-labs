"""
Main Entry Point
Bioinformatics Lab 6 - Transposable Element Simulation, FASTA Analysis & Detection App
Supports both Graphical User Interface (GUI) and Command-Line Interface (CLI).
Allows opening FASTA files from disk (e.g. fasta/ folder) and detecting new candidate
transposons with configurable minimum in-between spacing (default: >= 10 bases).
"""

import sys
import argparse
from typing import Dict, List, Optional

from transposon_simulator import TransposonSimulator, SimulationResult
from transposon_detector import TransposonDetector, DetectionResult
from evaluation import IntersectionAnalyzer, BenchmarkSummary
from fasta_io import read_fasta_file, list_fasta_folder


def print_ascii_genome_map(sim_res: SimulationResult, detections: List[DetectionResult], width: int = 70):
    """Render an ASCII track of the genome and detected elements."""
    total_len = sim_res.final_length
    scale = total_len / width

    print("\n=== ASCII Genomic Map ===")
    print(f"Locus: 0 bp" + " " * (width - 15) + f"{total_len} bp")
    print("-" * width)

    gt_track = ["."] * width
    labels = []
    symbol_chars = ["A", "B", "C", "D", "E", "F"]

    for idx, gt in enumerate(sim_res.ground_truth_intervals):
        x1 = min(width - 1, int(gt.start / scale))
        x2 = min(width - 1, int(gt.end / scale))
        sym = symbol_chars[idx % len(symbol_chars)]
        labels.append((sym, gt.element_name, gt.fragment_type, gt.start, gt.end))
        for x in range(x1, x2 + 1):
            gt_track[x] = sym

    print("GT Track:  " + "".join(gt_track))

    det_track = ["."] * width
    for det in detections:
        x1 = min(width - 1, int(det.start / scale))
        x2 = min(width - 1, int(det.end / scale))
        char = "=" if det.status == "CHAINED_INTERRUPTED" else "#"
        for x in range(x1, x2 + 1):
            det_track[x] = char

    print("Det Track: " + "".join(det_track))
    print("-" * width)
    print("Legend:")
    for sym, el_name, ftype, s, e in labels:
        print(f"  [{sym}] {el_name:<16} ({ftype:<18}) @ [{s}:{e}] ({e - s + 1} bp)")
    print("  [#] Detected element fragment | [=] Chained interrupted locus\n")


def print_ascii_fasta_map(seq_len: int, detections: List[DetectionResult], width: int = 70):
    """Render an ASCII track of a FASTA sequence and newly detected candidate transposons."""
    scale = seq_len / width

    print("\n=== ASCII Genomic Track for FASTA Sequence ===")
    print(f"0 bp" + " " * (width - 12) + f"{seq_len} bp")
    print("-" * width)

    track = ["."] * width
    for det in detections:
        x1 = min(width - 1, int(det.start / scale))
        x2 = min(width - 1, int(det.end / scale))
        for x in range(x1, x2 + 1):
            track[x] = "#"

    print("Candidate TEs: " + "".join(track))
    print("-" * width)


def run_cli_scenario(scenario_key: str, length: int, te_count: int, seed: Optional[int] = None):
    """Execute a single simulation scenario in CLI mode."""
    sim = TransposonSimulator(seed=seed)
    detector = TransposonDetector(min_between=10)

    if scenario_key == "non_intersecting":
        sim_res = sim.simulate_scenario_non_intersecting(length, te_count)
    elif scenario_key == "nested_2":
        sim_res = sim.simulate_scenario_2_intersecting_nested(length, te_count)
    elif scenario_key == "overlapping_2":
        sim_res = sim.simulate_scenario_2_intersecting_overlapping(length, te_count)
    elif scenario_key == "nested_3":
        sim_res = sim.simulate_scenario_3_intersecting_nested(length, te_count)
    else:
        print(f"Unknown scenario: {scenario_key}")
        return

    print("=" * 76)
    print(f" SCENARIO: {sim_res.scenario_name}")
    print("=" * 76)
    print(f"Artificial DNA Length: {sim_res.final_length} bp | Simulated TEs: {len(sim_res.transposon_library)}")
    print(f"Simulation Log:")
    for l in sim_res.simulation_log:
        print(f"  - {l}")

    results = detector.detect_all(sim_res.genome_sequence, sim_res.transposon_library)
    homology_dets = results["homology"]
    structural_dets = results["structural"]

    print("\n--- Detected Transposable Elements (Homology Alignment & Chaining) ---")
    print(f"{'Element Name':<34} {'Status':<22} {'Start':<7} {'End':<7} {'Len':<5} {'Method':<24}")
    print("-" * 105)
    for d in homology_dets:
        print(f"{d.element_name:<34} {d.status:<22} {d.start:<7} {d.end:<7} {d.length:<5} {d.method:<24}")

    print("\n--- Detected Transposable Elements (De Novo Structural TIRs) ---")
    print(f"{'Element Name':<32} {'Start':<7} {'End':<7} {'Len':<5} {'Details'}")
    print("-" * 90)
    for s in structural_dets:
        print(f"{s.element_name:<32} {s.start:<7} {s.end:<7} {s.length:<5} {s.details}")

    print_ascii_genome_map(sim_res, homology_dets)

    bench = IntersectionAnalyzer.evaluate(sim_res, homology_dets)
    for line in bench.intersection_report:
        print(line)

    print("\nPerformance Metrics:")
    print(f"  Precision: {bench.precision:.1f}% | Recall: {bench.recall:.1f}% | F1: {bench.f1_score:.1f}% | Mean IoU: {bench.mean_iou:.1f}%")
    print("=" * 76 + "\n")


def run_cli_fasta(filepath: str, min_between: int = 10, min_tir: int = 8):
    """Execute de novo detection on a FASTA file in CLI mode."""
    records = read_fasta_file(filepath)
    if not records:
        print(f"No valid FASTA records found in: {filepath}")
        return

    print("=" * 76)
    print(f" FASTA ANALYSIS: {filepath}")
    print(f" Found {len(records)} sequence record(s) | Min In-Between Bases: {min_between} bp")
    print("=" * 76)

    detector = TransposonDetector(min_between=min_between, min_tir_length=min_tir)

    for i, rec in enumerate(records, 1):
        print(f"\n--- Record {i}: {rec.id} ({rec.length} bp, GC {rec.gc_content:.1f}%) ---")
        if rec.description:
            print(f"    Description: {rec.description}")

        dets = detector.detect_new_candidates_from_sequence(
            rec.sequence,
            min_between=min_between,
            min_tir_length=min_tir
        )

        print(f"\nDiscovered Candidate Transposons ({len(dets)}):")
        print(f"{'Element Name':<28} {'Start':<7} {'End':<7} {'Len':<5} {'In-Between':<12} {'Left TIR':<14} {'Right TIR':<14} {'TSD':<6}")
        print("-" * 98)
        for d in dets:
            print(f"{d.element_name:<28} {d.start:<7} {d.end:<7} {d.length:<5} {f'{d.in_between_len} bp':<12} {d.left_tir:<14} {d.right_tir:<14} {d.tsd if d.tsd else '-':<6}")

        print_ascii_fasta_map(rec.length, dets)

        report = IntersectionAnalyzer.generate_fasta_report(rec.id, rec.length, dets, min_between=min_between)
        print("\nBehavior & Discovery Report:")
        for line in report:
            print(f"  {line}")


def run_all_benchmarks(length: int = 350, seed: int = 42):
    """Run a comparative benchmark across all 4 scenarios."""
    print("\n" + "=" * 76)
    print("   TRANSPOSABLE ELEMENT INTERSECTION BENCHMARK SUITE")
    print("=" * 76)
    scenarios = [
        ("non_intersecting", "Scenario A: Non-Intersecting Control", 4),
        ("nested_2", "Scenario B: 2-TE Nested Intersection", 4),
        ("overlapping_2", "Scenario C: 2-TE Overlapping Intersection", 4),
        ("nested_3", "Scenario D: 3-TE Nested Chain Intersection", 4),
    ]

    for key, title, count in scenarios:
        run_cli_scenario(key, length, count, seed=seed)


def main():
    parser = argparse.ArgumentParser(
        description="Bioinformatics Lab 6: Transposable Element Simulation, FASTA Analysis & Detection App"
    )
    parser.add_argument("--cli", action="store_true", help="Run in Command-Line Interface mode instead of GUI")
    parser.add_argument("--fasta", type=str, default=None, help="Path to a FASTA file to load and analyze (e.g. fasta/sample1_nested_350bp.fasta)")
    parser.add_argument("--min-between", type=int, default=10, help="Minimum number of bases in between the two TIRs (default: 10 bp)")
    parser.add_argument("--min-tir", type=int, default=8, help="Minimum length of Terminal Inverted Repeats in bp (default: 8 bp)")
    parser.add_argument("--benchmark", action="store_true", help="Run comprehensive automated benchmarks for all intersection scenarios")
    parser.add_argument("--scenario", choices=["non_intersecting", "nested_2", "overlapping_2", "nested_3", "all"],
                        default="nested_2", help="Simulation scenario to execute in CLI mode")
    parser.add_argument("--length", type=int, default=350, help="Genome sequence length in bp (200-500)")
    parser.add_argument("--elements", type=int, default=4, help="Number of transposable elements (3 or 4)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for deterministic output")

    args = parser.parse_args()

    # If --fasta is provided, analyze the FASTA file
    if args.fasta:
        run_cli_fasta(args.fasta, min_between=args.min_between, min_tir=args.min_tir)
        return

    # If --benchmark or --cli is provided, run simulation CLI
    if args.benchmark or args.cli:
        if args.benchmark or args.scenario == "all":
            run_all_benchmarks(length=args.length, seed=args.seed if args.seed is not None else 42)
        else:
            run_cli_scenario(args.scenario, length=args.length, te_count=args.elements, seed=args.seed)
        return

    # Otherwise launch Tkinter GUI
    try:
        from gui import launch_gui
        print("Launching Graphical User Interface (GUI)...")
        launch_gui()
    except Exception as e:
        print(f"Notice: Could not launch GUI ({e}). Falling back to CLI mode.\n")
        run_cli_scenario(args.scenario, length=args.length, te_count=args.elements, seed=args.seed)


if __name__ == "__main__":
    main()
