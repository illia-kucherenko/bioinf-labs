"""
FASTA I/O and Sequence Utilities Module
Bioinformatics Lab 6 - Reading, parsing, and managing FASTA genome sequences.
"""

import os
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple


@dataclass
class FastaRecord:
    """Represents a single sequence record parsed from a FASTA file."""
    header: str
    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        g = self.sequence.count('G') + self.sequence.count('g')
        c = self.sequence.count('C') + self.sequence.count('c')
        return (g + c) / len(self.sequence) * 100.0


def parse_fasta_string(content: str) -> List[FastaRecord]:
    """Parse FASTA formatted text into a list of FastaRecord objects."""
    records = []
    current_header = ""
    current_seq_parts = []

    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header or current_seq_parts:
                seq = "".join(current_seq_parts).upper()
                # Clean up any non-alphabetical or whitespace characters
                seq = "".join(b for b in seq if b in "ACGTN")
                header_line = current_header[1:].strip()
                parts = header_line.split(maxsplit=1)
                rec_id = parts[0] if parts else "unnamed_sequence"
                desc = parts[1] if len(parts) > 1 else ""
                records.append(FastaRecord(header=header_line, id=rec_id, description=desc, sequence=seq))
                current_seq_parts = []
            current_header = line
        else:
            current_seq_parts.append(line)

    if current_header or current_seq_parts:
        seq = "".join(current_seq_parts).upper()
        seq = "".join(b for b in seq if b in "ACGTN")
        header_line = current_header[1:].strip() if current_header else "unnamed_sequence"
        parts = header_line.split(maxsplit=1)
        rec_id = parts[0] if parts else "unnamed_sequence"
        desc = parts[1] if len(parts) > 1 else ""
        records.append(FastaRecord(header=header_line, id=rec_id, description=desc, sequence=seq))

    return records


def read_fasta_file(filepath: str) -> List[FastaRecord]:
    """Read and parse a FASTA file from disk."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    return parse_fasta_string(content)


def write_fasta_file(filepath: str, records: List[Tuple[str, str]], line_wrap: int = 60) -> None:
    """
    Write sequence records to a FASTA file.
    records: list of (header, sequence)
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_wrap):
                f.write(seq[i : i + line_wrap] + "\n")


def list_fasta_folder(folder_path: str = "fasta") -> List[str]:
    """Return all FASTA files (.fasta, .fa, .fna, .txt) in the given directory."""
    if not os.path.isdir(folder_path):
        return []

    valid_exts = {".fasta", ".fa", ".fna", ".txt"}
    files = []
    for fname in os.listdir(folder_path):
        ext = os.path.splitext(fname)[1].lower()
        if ext in valid_exts:
            files.append(os.path.join(folder_path, fname))
    return sorted(files)
