"""Dynamic programming sequence alignment using matrices."""

from __future__ import annotations
import os
from typing import List, Tuple, Dict, Any, Optional
import matplotlib.pyplot as plt
import numpy as np


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -1
) -> Dict[str, Any]:
    """
    Performs global sequence alignment using the Needleman-Wunsch algorithm.

    Parameters:
    -----------
    seq1 : str
        First DNA sequence (vertical axis / rows).
    seq2 : str
        Second DNA sequence (horizontal axis / columns).
    match : int
        Score for matching nucleotides.
    mismatch : int
        Penalty/score for mismatching nucleotides.
    gap : int
        Penalty/score for a gap.

    Returns:
    --------
    dict with keys:
        'aligned_seq1': str
        'aligned_seq2': str
        'alignment_marks': str (e.g. '|' for match, '*' for mismatch, ' ' for gap)
        'score': int
        'score_matrix': List[List[int]]
        'traceback_matrix': List[List[str]]
        'path': List[Tuple[int, int]]
        'stats': dict of match/mismatch/gap counts
    """
    s1 = seq1.upper().strip()
    s2 = seq2.upper().strip()
    n = len(s1)
    m = len(s2)

    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[""] * (m + 1) for _ in range(n + 1)]

    # Initialize first column (gaps in seq2)
    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap
        traceback[i][0] = "U"

    # Initialize first row (gaps in seq1)
    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap
        traceback[0][j] = "L"

    # Fill dynamic programming matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub_score = match if s1[i - 1] == s2[j - 1] else mismatch
            diag = score_matrix[i - 1][j - 1] + sub_score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap

            best = max(diag, up, left)
            score_matrix[i][j] = best

            # Priority: Diagonal (match/mismatch), Up (gap in seq2), Left (gap in seq1)
            if best == diag:
                traceback[i][j] = "D"
            elif best == up:
                traceback[i][j] = "U"
            else:
                traceback[i][j] = "L"

    # Traceback path from bottom-right (n, m) to top-left (0, 0)
    i, j = n, m
    al1: List[str] = []
    al2: List[str] = []
    path: List[Tuple[int, int]] = [(i, j)]

    while i > 0 or j > 0:
        step = traceback[i][j]
        if step == "D":
            al1.append(s1[i - 1])
            al2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif step == "U":
            al1.append(s1[i - 1])
            al2.append("-")
            i -= 1
        else:
            al1.append("-")
            al2.append(s2[j - 1])
            j -= 1
        path.append((i, j))

    path.reverse()
    al1_str = "".join(reversed(al1))
    al2_str = "".join(reversed(al2))

    # Alignment visualization string
    marks = []
    matches = 0
    mismatches = 0
    gaps = 0

    for a, b in zip(al1_str, al2_str):
        if a == "-" or b == "-":
            marks.append(" ")
            gaps += 1
        elif a == b:
            marks.append("|")
            matches += 1
        else:
            marks.append("*")
            mismatches += 1

    marks_str = "".join(marks)

    return {
        "seq1": s1,
        "seq2": s2,
        "aligned_seq1": al1_str,
        "aligned_seq2": al2_str,
        "alignment_marks": marks_str,
        "score": score_matrix[n][m],
        "score_matrix": score_matrix,
        "traceback_matrix": traceback,
        "path": path,
        "stats": {
            "matches": matches,
            "mismatches": mismatches,
            "gaps": gaps,
            "length": len(al1_str),
            "identity": (matches / len(al1_str) * 100.0) if al1_str else 0.0
        }
    }


def format_matrix_text(
    score_matrix: List[List[int]],
    seq1: str,
    seq2: str,
    path: Optional[List[Tuple[int, int]]] = None
) -> str:
    """Renders the scoring matrix as a clean formatted table string."""
    header = [" ", "-"] + list(seq2)
    path_set = set(path) if path else set()

    col_width = 6
    lines = []

    # Header line
    lines.append("".join(f"{h:>{col_width}}" for h in header))
    lines.append("-" * (col_width * len(header)))

    row_labels = ["-"] + list(seq1)
    for i, row in enumerate(score_matrix):
        row_str = [f"{row_labels[i]:>{col_width}}"]
        for j, val in enumerate(row):
            marker = "*" if (i, j) in path_set else " "
            cell = f"{val}{marker}"
            row_str.append(f"{cell:>{col_width}}")
        lines.append("".join(row_str))

    return "\n".join(lines)


def plot_alignment_matrix(
    score_matrix: List[List[int]],
    seq1: str,
    seq2: str,
    path: Optional[List[Tuple[int, int]]] = None,
    title: str = "Needleman-Wunsch Alignment Matrix",
    save_path: Optional[str] = None,
    show: bool = False
) -> None:
    """Plots a visual heatmap of the scoring matrix with the traceback path."""
    matrix = np.array(score_matrix)
    n, m = matrix.shape

    fig, ax = plt.subplots(figsize=(max(6, m * 0.8), max(5, n * 0.7)), dpi=150)
    fig.patch.set_facecolor("#1e1e24")
    ax.set_facecolor("#2b2b36")

    cax = ax.imshow(matrix, cmap="Blues", interpolation="nearest", alpha=0.88)

    # Add text numbers into each cell
    for i in range(n):
        for j in range(m):
            val = matrix[i, j]
            color = "white" if val < (matrix.max() + matrix.min()) / 2 else "#111111"
            ax.text(j, i, str(val), ha="center", va="center", color=color, fontsize=10, weight="medium")

    # Traceback path overlay
    if path:
        path_y = [p[0] for p in path]
        path_x = [p[1] for p in path]
        ax.plot(path_x, path_y, color="#ff4757", linewidth=2.5, marker="o", markersize=6, zorder=5, label="Optimal Path")

    col_labels = ["-"] + list(seq2)
    row_labels = ["-"] + list(seq1)

    ax.set_xticks(range(m))
    ax.set_yticks(range(n))
    ax.set_xticklabels(col_labels, color="white", fontsize=11, weight="bold")
    ax.set_yticklabels(row_labels, color="white", fontsize=11, weight="bold")

    ax.set_xlabel(f"Sequence 2 ({seq2})", color="white", fontsize=12, labelpad=8)
    ax.set_ylabel(f"Sequence 1 ({seq1})", color="white", fontsize=12, labelpad=8)
    ax.set_title(title, color="white", fontsize=13, pad=12, weight="bold")

    cbar = fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Alignment Score", color="white", fontsize=10)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

    if path:
        legend = ax.legend(facecolor="#1e1e24", edgecolor="white", labelcolor="white", loc="upper left")
        legend.get_frame().set_alpha(0.85)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        plt.savefig(save_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
        print(f"Saved: {save_path}")

    if show:
        plt.show()
    plt.close(fig)
