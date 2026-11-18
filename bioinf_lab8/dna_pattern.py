"""DNA Pattern Generation from Gene Promoters (Bioinformatics Lab 8).

Implements algorithms inspired by PromKappa:
- C+G% (CpG) content calculation
- Information content (Kappa Index of Coincidence / Self-Sequence Alignment)
- Sliding-window pattern (Objective Digital Stain) generation
- Center of weight (centroid) computation
- Comparative multi-promoter visualization
"""

from __future__ import annotations
import os
from typing import Dict, List, Tuple, Any, Optional
import matplotlib.pyplot as plt
import numpy as np


def clean_sequence(seq: str) -> str:
    """Removes whitespace, digits, and line breaks, returning uppercase sequence."""
    return "".join(c for c in seq.upper() if c in {"A", "C", "G", "T", "U", "N"})


def calculate_cg_content(
    sequence: str,
    mode: str = "normal",
    total_cg_pct: Optional[float] = None
) -> float:
    """Computes the C+G (CpG) percentage for a DNA sequence window."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return 0.0

    cg_count = seq.count("C") + seq.count("G")

    if mode == "sequence_sensitive":
        if total_cg_pct is None:
            raise ValueError("total_cg_pct must be provided when mode is 'sequence_sensitive'.")
        return (total_cg_pct / n) * cg_count

    return (cg_count / n) * 100.0


def calculate_information_content(
    sequence: str,
    metric: str = "ic"
) -> float:
    """Computes the Information Content via Self-Sequence Alignment."""
    seq = sequence.upper()
    n = len(seq)
    if n <= 1:
        return 0.0

    total_pct = 0.0
    for u in range(1, n):
        s2 = seq[u:]
        matches = sum(1 for i in range(len(s2)) if seq[i] == s2[i])
        pct_u = (matches / (n - u)) * 100.0
        total_pct += pct_u

    ic_val = total_pct / (n - 1)

    if metric == "sigma":
        return 100.0 - ic_val

    return ic_val


def sliding_window_pattern(
    sequence: str,
    window_size: int = 30,
    step: int = 1,
    metric: str = "ic",
    cg_mode: str = "normal",
    include_endpoint: bool = True
) -> List[Dict[str, Any]]:
    """Scans a DNA sequence with a sliding window, computing (CG%, IC) for each window."""
    clean_seq = clean_sequence(sequence)
    seq_len = len(clean_seq)
    if seq_len < window_size:
        raise ValueError(f"Sequence length ({seq_len}) must be at least window size ({window_size}).")

    total_cg_pct = (clean_seq.count("C") + clean_seq.count("G")) / seq_len * 100.0
    max_idx = (seq_len - window_size + 1) if include_endpoint else (seq_len - window_size)

    pattern_points = []
    for idx in range(0, max_idx, step):
        subseq = clean_seq[idx : idx + window_size]
        cg_val = calculate_cg_content(subseq, mode=cg_mode, total_cg_pct=total_cg_pct)
        ic_val = calculate_information_content(subseq, metric=metric)
        pattern_points.append({
            "index": idx,
            "window": subseq,
            "cg": cg_val,
            "ic": ic_val
        })

    return pattern_points


def calculate_center_of_weight(pattern_points: List[Dict[str, Any]]) -> Tuple[float, float]:
    """Calculates the center of weight (centroid) of the pattern: (mean CG%, mean IC)."""
    if not pattern_points:
        return 0.0, 0.0

    mean_cg = sum(p["cg"] for p in pattern_points) / len(pattern_points)
    mean_ic = sum(p["ic"] for p in pattern_points) / len(pattern_points)
    return mean_cg, mean_ic


def load_fasta(fasta_path: str) -> Dict[str, Dict[str, str]]:
    """Parses a FASTA file and extracts promoter names and sequences."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found at: {fasta_path}")

    promoters: Dict[str, Dict[str, str]] = {}
    current_header = ""
    current_name = ""
    seq_chunks: List[str] = []

    with open(fasta_path, "r", encoding="latin1") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name and seq_chunks:
                    promoters[current_name] = {
                        "header": current_header,
                        "sequence": "".join(seq_chunks)
                    }
                current_header = line
                parts = line.split(";")
                left = parts[0].strip()
                tokens = left.split()
                current_name = tokens[-1] if tokens else f"seq_{len(promoters) + 1}"
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if current_name and seq_chunks:
            promoters[current_name] = {
                "header": current_header,
                "sequence": "".join(seq_chunks)
            }

    return promoters


def plot_dna_pattern(
    pattern_points: List[Dict[str, Any]],
    title: str = "DNA Pattern (Objective Digital Stain)",
    center_of_weight: Optional[Tuple[float, float]] = None,
    metric_label: str = "Information Content (IC %)",
    save_path: Optional[str] = None,
    show: bool = False
) -> None:
    """Plots the 2D DNA pattern trajectory and highlights the Center of Weight."""
    if not pattern_points:
        return

    x_vals = [p["cg"] for p in pattern_points]
    y_vals = [p["ic"] for p in pattern_points]
    indices = np.arange(len(pattern_points))

    if center_of_weight is None:
        center_of_weight = calculate_center_of_weight(pattern_points)

    fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
    fig.patch.set_facecolor("#1e1e24")
    ax.set_facecolor("#2b2b36")

    # Trajectory path and points
    ax.plot(x_vals, y_vals, color="#4ecdc4", alpha=0.35, linewidth=1.0, zorder=1)
    scatter = ax.scatter(
        x_vals, y_vals,
        c=indices,
        cmap="cool",
        s=35,
        alpha=0.85,
        edgecolors="none",
        zorder=2
    )

    # Center of weight marker
    cw_x, cw_y = center_of_weight
    ax.scatter(
        [cw_x], [cw_y],
        color="#ff6b6b",
        s=140,
        edgecolors="white",
        linewidths=1.8,
        marker="o",
        zorder=5,
        label=f"Center of Weight ({cw_x:.2f}%, {cw_y:.2f}%)"
    )

    ax.axvline(cw_x, color="#ff6b6b", linestyle="--", alpha=0.5, linewidth=1.0)
    ax.axhline(cw_y, color="#ff6b6b", linestyle="--", alpha=0.5, linewidth=1.0)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("C + G (%)", color="white", fontsize=12, labelpad=8)
    ax.set_ylabel(metric_label, color="white", fontsize=12, labelpad=8)
    ax.set_title(title, color="white", fontsize=14, pad=12, weight="bold")

    ax.tick_params(colors="white")
    ax.grid(True, linestyle=":", alpha=0.25, color="gray")

    cbar = fig.colorbar(scatter, ax=ax, pad=0.03, fraction=0.046)
    cbar.set_label("Window Index (5' → 3')", color="white", fontsize=10)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

    legend = ax.legend(facecolor="#1e1e24", edgecolor="white", labelcolor="white", loc="upper right")
    legend.get_frame().set_alpha(0.85)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        plt.savefig(save_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        print(f"Saved: {save_path}")

    if show:
        plt.show()
    plt.close(fig)


def plot_centers_of_weight(
    promoter_centroids: List[Dict[str, Any]],
    title: str = "Gene Promoters - Centers of Weight Distribution",
    metric_label: str = "Average Information Content (IC %)",
    save_path: Optional[str] = None,
    show: bool = False
) -> None:
    """Plots comparative centers of weight for multiple gene promoter patterns."""
    if not promoter_centroids:
        return

    x_vals = [p["mean_cg"] for p in promoter_centroids]
    y_vals = [p["mean_ic"] for p in promoter_centroids]

    global_cg = float(np.mean(x_vals))
    global_ic = float(np.mean(y_vals))

    fig, ax = plt.subplots(figsize=(9, 8), dpi=150)
    fig.patch.set_facecolor("#1e1e24")
    ax.set_facecolor("#2b2b36")

    ax.scatter(
        x_vals, y_vals,
        color="#00d2fc",
        s=70,
        alpha=0.85,
        edgecolors="white",
        linewidths=0.8,
        zorder=3
    )

    for p in promoter_centroids:
        ax.annotate(
            p["name"],
            (p["mean_cg"], p["mean_ic"]),
            xytext=(6, 4),
            textcoords="offset points",
            fontsize=8.5,
            color="#e0e0e0",
            weight="medium"
        )

    ax.scatter(
        [global_cg], [global_ic],
        color="#ff4757",
        marker="X",
        s=160,
        edgecolors="white",
        linewidths=1.5,
        zorder=5,
        label=f"Global Mean ({global_cg:.2f}%, {global_ic:.2f}%)"
    )
    ax.axvline(global_cg, color="#ff4757", linestyle=":", alpha=0.6, linewidth=1.2)
    ax.axhline(global_ic, color="#ff4757", linestyle=":", alpha=0.6, linewidth=1.2)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Average C + G (%)", color="white", fontsize=12, labelpad=8)
    ax.set_ylabel(metric_label, color="white", fontsize=12, labelpad=8)
    ax.set_title(title, color="white", fontsize=14, pad=12, weight="bold")

    ax.tick_params(colors="white")
    ax.grid(True, linestyle=":", alpha=0.25, color="gray")

    legend = ax.legend(facecolor="#1e1e24", edgecolor="white", labelcolor="white", loc="upper right")
    legend.get_frame().set_alpha(0.85)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        plt.savefig(save_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        print(f"Saved: {save_path}")

    if show:
        plt.show()
    plt.close(fig)
