"""
ncbi_fetcher.py - Module to download genome sequences from NCBI Entrez API.
"""
import os
import urllib.request
import urllib.error


def fetch_ncbi_genome(accession_id: str = "NC_001422.1", data_dir: str = "data") -> tuple[str, str]:
    """
    Downloads a FASTA sequence from NCBI given an accession ID.
    Caches the file locally in `data_dir`.

    Returns:
        tuple[header: str, sequence: str]
    """
    os.makedirs(data_dir, exist_ok=True)
    file_path = os.path.join(data_dir, f"{accession_id}.fasta")

    if os.path.exists(file_path):
        print(f"[NCBI] Loading cached genome from {file_path}")
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
    else:
        url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=nuccore&id={accession_id}&rettype=fasta&retmode=text"
        )
        print(f"[NCBI] Fetching accession {accession_id} from NCBI Entrez API...")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BioinfLab4/1.0"})
            with urllib.request.urlopen(req) as resp:
                content = resp.read().decode("utf-8")
        except urllib.error.URLError as e:
            raise RuntimeError(f"Failed to download genome {accession_id} from NCBI: {e}")

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"[NCBI] Saved genome to {file_path}")

    # Parse FASTA
    lines = content.strip().splitlines()
    header = lines[0] if lines else ""
    sequence_parts = [line.strip().upper() for line in lines if not line.startswith(">")]
    sequence = "".join(sequence_parts)

    print(f"[NCBI] Genome '{header}' successfully loaded ({len(sequence):,} bp).")
    return header, sequence


if __name__ == "__main__":
    header, seq = fetch_ncbi_genome()
    print("Header:", header)
    print("Length:", len(seq))
    print("First 100 bp:", seq[:100])
