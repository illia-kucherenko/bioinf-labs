# Bioinformatics Lab 6: Transposable Element Simulation, FASTA Analysis & Detection System

A software application implemented in **Python 3.12** for:
1. **Simulating artificial DNA sequences (200–500 bp)** with 3 to 4 transposable elements (TEs) to analyze algorithmic behavior under **nested and overlapping intersections**.
2. **Opening and analyzing real or custom FASTA files** (`.fasta`, `.fa`, `.fna`) located in the project's `fasta/` folder (or via browse).
3. **Detecting new putative candidate transposons *de novo*** with configurable minimum internal spacer distance (default: **minimized to 10 bases** between the Terminal Inverted Repeats).

---

## 1. Project Overview & Laboratory Task

As specified in `Task.md`:
> * Make an artificial DNA sequence of 200 to 500 bases in length in which to simulate 3 to 4 transposable elements.
> * Implement a software app which can detect the positions of these transposable elements (start, end) within the created DNA sequence.
> * **Note**: Intersect 2 or 3 transposons to see how the algorithm behaves in these situations.

### Key Biological Concepts
- **Class II Transposable Elements (DNA Transposons)**: Characterized by Terminal Inverted Repeats (TIRs) flanking a core sequence, generating short Target Site Duplications (TSD) upon integration.
- **Nested Intersections**: Transposon B inserts *into* Transposon A, interrupting Transposon A and fragmenting it into 5' and 3' halves.
- **Hierarchical (3-Way) Nesting**: Transposon C inserts into Transposon B, which previously inserted into Transposon A (`[TE1_5'] -> [TE2_5'] -> [TE3] -> [TE2_3'] -> [TE1_3']`).
- **Overlapping Intersections**: Transposon B recombinantly replaces or truncates the terminal repeat of Transposon A.
- **Minimized In-Between Spacing (>= 10 bp)**: Small non-autonomous elements (e.g. MITEs) require as few as 10 internal bases between the left and right TIRs.

---

## 2. Architecture & Modules

| Module | Description |
|---|---|
| [`main.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/main.py) | Unified entry point for Graphical User Interface (GUI) and Command-Line Interface (CLI), supporting simulations and FASTA file scanning. |
| [`gui.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/gui.py) | Interactive desktop application built with Tkinter featuring live sequence track canvas, coordinate tables, FASTA loading, and behavior inspector. |
| [`fasta_io.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/fasta_io.py) | Parser for reading single- and multi-record FASTA files from the `fasta/` directory. |
| [`transposon_simulator.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/transposon_simulator.py) | Generates synthetic DNA (200–500 bp) and simulates 3–4 transposons under independent, nested, and overlapping scenarios. |
| [`transposon_detector.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/transposon_detector.py) | Implements **Homology Alignment & Chaining** and **De Novo Structural TIR/TSD Detection** with configurable `min_between=10` bp spacer constraint. |
| [`evaluation.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/evaluation.py) | Evaluates predictions against ground truth and formats de novo discovery reports for FASTA inputs. |
| [`test_transposons.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/test_transposons.py) | 14 unit tests verifying simulation, detection, coordinates, and intersection edge cases. |
| [`test_fasta.py`](file:///C:/Users/skinc/PycharmProjects/bioinf_lab6/test_fasta.py) | 4 unit tests verifying FASTA parsing and the strict $\ge 10$ bp in-between spacing constraint. |

---

## 3. Working with FASTA Files (`fasta/` Folder)

You can place your FASTA files in the `fasta/` directory located at the project root:
```
bioinf_lab6/
  ├── fasta/
  │    ├── sample1_nested_350bp.fasta
  │    ├── sample2_minimal_10bp_spacer.fasta
  │    └── (add your custom .fasta files here)
```

### Loading FASTA Files in GUI Mode
1. Run `python main.py` to open the GUI.
2. Select the **"Load FASTA File (from fasta/ folder or browse)"** radio button.
3. Choose any file from the **"FASTA File"** dropdown (or click **"Browse..."**).
4. Click **"Refresh fasta/"** anytime you add new files into the `fasta/` folder.
5. Set **"Min Bases In-Between"** (default: **10** bp) and **"Min TIR Length"** (default: **8** bp).
6. Click **"[>] Detect Candidate Transposons"** to scan the sequence and visualize the results!

### Running FASTA Detection in CLI Mode
```powershell
# Analyze sample FASTA with >= 10 bases in-between:
.\.venv\Scripts\python.exe main.py --fasta fasta/sample1_nested_350bp.fasta --min-between 10

# Analyze sequence with minimal 10 bp spacer:
.\.venv\Scripts\python.exe main.py --fasta fasta/sample2_minimal_10bp_spacer.fasta --min-between 10
```

---

## 4. How the Algorithm Behaves When Transposons Intersect

### A. Non-Intersecting (Control Baseline)
- **Homology Detector**: Detects 100% of query sequences with continuous, gap-free alignments (100% recall, 100% precision, 100% IoU).
- **Structural Detector**: Accurately pairs the 5' and 3' TIRs across the standard element length.

### B. Two-Transposon Nested Intersection ($TE_2$ inside $TE_1$)
- **Effect on $TE_1$ (Outer Element)**:
  - $TE_1$ is interrupted and split into two distinct sub-intervals: a 5' prefix fragment and a 3' suffix fragment.
  - Naive global alignment fails because of the large intervening sequence ($TE_2$).
  - Our **Fragment Chainer** detects both fragments, computes their coordinates, and links them into an interrupted compound locus:
    `[start_5p : end_3p]` interrupted by $TE_2$.
- **Effect on $TE_2$ (Inner Element)**:
  - $TE_2$ remains 100% intact and contiguous within the interior of $TE_1$.
  - Detected with full confidence and exact coordinates.
- **Structural (De Novo) Behavior**:
  - Detects the inner TIR pair for $TE_2$.
  - Detects the outer TIR pair for $TE_1$, noticing that the distance between TIRs expanded by the length of $TE_2$, and tags it as `STRUCTURAL_NESTED_PARENT`.

### C. Two-Transposon Overlapping Intersection (Truncation / Recombination)
- **Effect on Upstream Element ($TE_1$)**:
  - Its 3' region is overwritten by $TE_2$.
  - Homology alignment detects only the remaining 5' fragment.
  - Structural detection fails for $TE_1$ because its 3' TIR was destroyed (loss of repeat symmetry).
- **Effect on Downstream Element ($TE_2$)**:
  - Retains both terminal repeats and is detected intact.

### D. Three-Transposon Hierarchical Nesting ($TE_3 \subset TE_2 \subset TE_1$)
- **Architecture**: `[TE1_5']` $\rightarrow$ `[TE2_5']` $\rightarrow$ `[TE3_full]` $\rightarrow$ `[TE2_3']` $\rightarrow$ `[TE1_3']`.
- **Behavior**:
  - $TE_3$ (innermost core): Detected intact with 100% coverage.
  - $TE_2$ (middle layer): Detected as two fragments flanking $TE_3$; chained together into a nested parent.
  - $TE_1$ (outer layer): Detected as two fragments flanking the combined $TE_2+TE_3$ locus; chained together as the outer parent.
  - Structural detector displays a recursive nested hierarchy (`TE1 contains TE2 contains TE3`).

---

## 5. Running the Application & Tests

```powershell
# 1. Launch the Graphical Desktop App (Default)
.\.venv\Scripts\python.exe main.py

# 2. Run CLI Simulation for Nested Intersection
.\.venv\Scripts\python.exe main.py --cli --scenario nested_2 --length 350 --elements 4

# 3. Run CLI FASTA Analysis on any file in fasta/
.\.venv\Scripts\python.exe main.py --fasta fasta/sample2_minimal_10bp_spacer.fasta --min-between 10

# 4. Run All Intersection Benchmarks
.\.venv\Scripts\python.exe main.py --benchmark

# 5. Run All 18 Unit Tests
.\.venv\Scripts\python.exe -m unittest discover -v
```
