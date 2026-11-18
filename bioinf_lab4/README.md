# Bioinformatics Lab 4: De Novo Viral Genome Assembly & Shotgun Sequencing Simulation

This project simulates shotgun sequencing and de novo genome assembly for viral DNA, evaluating the reconstructed sequence against the reference genome and analyzing the fundamental algorithmic bottlenecks of sequence assembly.

---

## 🎯 Task Objectives

1. **NCBI Genome Acquisition**: Download an entire viral genome sequence from the NCBI Entrez API (default: *Escherichia phage PhiX174*, `NC_001422.1`, 5,386 bp).
2. **In-Silico Shotgun Sampling**: Generate 2,000 random short-read samples from the genome of length 100–120 bp and store them in an array/list.
3. **Sequence Reconstruction**: Rebuild the DNA sequence from the random samples using de novo assembly algorithms (Greedy Overlap Assembly / Shortest Superstring heuristic and De Bruijn Graph topological analysis).
4. **Algorithmic Failure Mode Diagnosis**: Identify and empirically prove the primary obstacles and bottlenecks with the algorithmic approach (repeats, coverage dropouts, computational complexity, and circularity).

---

## 📁 Project Structure

```text
bioinf_lab4/
├── main.py                   # Main pipeline: download, sampling, assembly, evaluation & plots
├── ncbi_fetcher.py           # NCBI Entrez API client with local disk caching
├── sampler.py                # In-silico shotgun sequencer & Lander-Waterman coverage profiling
├── assembler.py              # Indexed Greedy Overlap Assembler & De Bruijn Graph analyzer
├── evaluator.py              # Sequence validation (N50, coverage fraction, exact match)
├── experiment_problems.py    # Empirical experiments demonstrating repeat collapse & gaps
├── task.md                   # Original lab task description
├── README.md                 # Project documentation
├── .gitignore                # Git ignore configuration
└── data/                     # Output directory
    ├── NC_001422.1.fasta     # Cached reference genome FASTA
    ├── sampled_reads.fasta   # 2,000 generated shotgun reads
    ├── assembly_report.json  # Machine-readable evaluation report
    ├── coverage_depth.png    # Genome-wide sequencing depth visualization
    └── overlap_distribution.png # Distribution of merged overlap lengths
```

---

## 🚀 Getting Started

### 1. Requirements

- Python 3.10+
- Dependencies: `biopython`, `matplotlib` (both installed in `.venv`)

### 2. Run the Main Pipeline

Executes the complete workflow: downloads genome, draws 2,000 samples, reconstructs the sequence, calculates metrics, and generates visual plots:

```powershell
& .\.venv\Scripts\python.exe main.py
```

### 3. Run the Algorithmic Bottleneck Experiments

Executes targeted empirical demonstrations of:
- **Experiment 1**: The Repeat Resolution Barrier (Repeat Length R < Read Length L vs. R >= L)
- **Experiment 2**: Coverage gaps under varying sequencing depths (Lander-Waterman model)
- **Experiment 3**: Sensitivity to minimum overlap threshold K

```powershell
& .\.venv\Scripts\python.exe experiment_problems.py
```

---

## 📊 Summary of Results

| Metric | Result |
| :--- | :--- |
| **Reference Genome** | *Bacteriophage PhiX174* (`NC_001422.1`, 5,386 bp) |
| **Sampled Reads** | 2,000 reads (100–120 bp, total 220,011 bp) |
| **Mean Sequencing Depth** | **40.85x** |
| **Assembly Algorithm** | Prefix-Indexed Greedy Overlap (Shortest Superstring) |
| **Reconstructed Contig** | **1 single contig** of **5,380 bp** (99.89% genome coverage) |
| **Exact Substring Match** | **True** (100% identity to reference) |
| **Assembly Runtime** | **~0.45 seconds** |

---

## 🧠 The Main Problems with the Algorithmic Approach

As required by the assignment, our experiments identified four primary bottlenecks:

### 1. Repetitive Sequences (The Core Assembly Bottleneck)
- **The Limit**: If a repeat length `R` is greater than or equal to read length `L` (`R >= L`), no read can bridge across the repeat into unique flanking regions.
- **Consequence**: The assembler cannot determine the correct linear arrangement of the flanking segments, resulting in **repeat collapse** (loss of intermediate DNA) or **chimeric joins**. In our test with `R = 140 bp > 120 bp`, the assembly shattered into multiple disconnected contigs.

### 2. Coverage Gaps (Lander-Waterman Island Problem)
- Random sampling does not guarantee 100% coverage (`P(uncovered) = e^(-c)`).
- Any region with 0 reads or overlap shorter than cutoff `K` breaks the sequence into fragmented contigs. In our test at 4.1x coverage, the genome broke into 12 separate contigs.

### 3. NP-Hardness and Local Optima
- Finding the exact Shortest Superstring containing all reads is **NP-hard**.
- Greedy heuristics merge the currently longest overlap first and cannot backtrack. Early erroneous choices on near-repeats permanently trap the algorithm in sub-optimal or chimeric assemblies.

### 4. Circular Viral Genomes & Boundary Overhangs
- Phages like PhiX174 are circular molecules. Linear shotgun assembly cannot bridge the `0 <-> 5,386 bp` origin without circular wrapping logic, leaving small unjoined overhangs at the termini.
