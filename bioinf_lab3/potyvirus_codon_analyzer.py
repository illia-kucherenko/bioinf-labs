import os
import json
import urllib.request
from collections import Counter
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from Bio import SeqIO
from Bio.Seq import Seq

# =====================================================================
# PART 1: GENETIC CODE CONVERSION FUNCTION
# =====================================================================

STANDARD_GENETIC_CODE = {
    # Phenylalanine (Phe / F)
    'UUU': 'Phenylalanine (Phe)', 'UUC': 'Phenylalanine (Phe)',
    # Leucine (Leu / L)
    'UUA': 'Leucine (Leu)', 'UUG': 'Leucine (Leu)',
    'CUU': 'Leucine (Leu)', 'CUC': 'Leucine (Leu)',
    'CUA': 'Leucine (Leu)', 'CUG': 'Leucine (Leu)',
    # Isoleucine (Ile / I)
    'AUU': 'Isoleucine (Ile)', 'AUC': 'Isoleucine (Ile)', 'AUA': 'Isoleucine (Ile)',
    # Methionine (Met / M) - Start
    'AUG': 'Methionine (Met)',
    # Valine (Val / V)
    'GUU': 'Valine (Val)', 'GUC': 'Valine (Val)',
    'GUA': 'Valine (Val)', 'GUG': 'Valine (Val)',
    # Serine (Ser / S)
    'UCU': 'Serine (Ser)', 'UCC': 'Serine (Ser)',
    'UCA': 'Serine (Ser)', 'UCG': 'Serine (Ser)',
    'AGU': 'Serine (Ser)', 'AGC': 'Serine (Ser)',
    # Proline (Pro / P)
    'CCU': 'Proline (Pro)', 'CCC': 'Proline (Pro)',
    'CCA': 'Proline (Pro)', 'CCG': 'Proline (Pro)',
    # Threonine (Thr / T)
    'ACU': 'Threonine (Thr)', 'ACC': 'Threonine (Thr)',
    'ACA': 'Threonine (Thr)', 'ACG': 'Threonine (Thr)',
    # Alanine (Ala / A)
    'GCU': 'Alanine (Ala)', 'GCC': 'Alanine (Ala)',
    'GCA': 'Alanine (Ala)', 'GCG': 'Alanine (Ala)',
    # Tyrosine (Tyr / Y)
    'UAU': 'Tyrosine (Tyr)', 'UAC': 'Tyrosine (Tyr)',
    # Stop Codons
    'UAA': 'Stop (*)', 'UAG': 'Stop (*)', 'UGA': 'Stop (*)',
    # Histidine (His / H)
    'CAU': 'Histidine (His)', 'CAC': 'Histidine (His)',
    # Glutamine (Gln / Q)
    'CAA': 'Glutamine (Gln)', 'CAG': 'Glutamine (Gln)',
    # Asparagine (Asn / N)
    'AAU': 'Asparagine (Asn)', 'AAC': 'Asparagine (Asn)',
    # Lysine (Lys / K)
    'AAA': 'Lysine (Lys)', 'AAG': 'Lysine (Lys)',
    # Aspartic Acid (Asp / D)
    'GAU': 'Aspartic Acid (Asp)', 'GAC': 'Aspartic Acid (Asp)',
    # Glutamic Acid (Glu / E)
    'GAA': 'Glutamic Acid (Glu)', 'GAG': 'Glutamic Acid (Glu)',
    # Cysteine (Cys / C)
    'UGU': 'Cysteine (Cys)', 'UGC': 'Cysteine (Cys)',
    # Tryptophan (Trp / W)
    'UGG': 'Tryptophan (Trp)',
    # Arginine (Arg / R)
    'CGU': 'Arginine (Arg)', 'CGC': 'Arginine (Arg)',
    'CGA': 'Arginine (Arg)', 'CGG': 'Arginine (Arg)',
    'AGA': 'Arginine (Arg)', 'AGG': 'Arginine (Arg)',
    # Glycine (Gly / G)
    'GGU': 'Glycine (Gly)', 'GGC': 'Glycine (Gly)',
    'GGA': 'Glycine (Gly)', 'GGG': 'Glycine (Gly)'
}

def codon_to_amino_acid(codon: str) -> str:
    """
    Converts a 3-nucleotide DNA or RNA codon into its corresponding amino acid
    using the Standard Genetic Code table.
    """
    clean_codon = codon.strip().upper().replace('T', 'U')
    if len(clean_codon) != 3:
        raise ValueError(f"Invalid codon length for '{codon}'. Must be exactly 3 nucleotides.")
    return STANDARD_GENETIC_CODE.get(clean_codon, "Unknown / Ambiguous (X)")


# =====================================================================
# PART 2: GENOME EXTRACTION & ALL-FRAME CODON FREQUENCY ANALYZER
# =====================================================================

def extract_genome_from_fasta(
    fasta_path: str = "potyvirus_thevetiae_OM263475.fasta",
    return_record: bool = False
):
    """
    Extracts the complete viral genome sequence from a local FASTA file
    (e.g., potyvirus_thevetiae_OM263475.fasta).

    Args:
        fasta_path (str): Path to the FASTA file. Defaults to 'potyvirus_thevetiae_OM263475.fasta'.
        return_record (bool): If True, returns the Bio.SeqRecord object.
                              If False (default), returns the uppercase nucleotide sequence as a str.

    Returns:
        str or SeqRecord: The extracted genome sequence string or SeqRecord object.

    Raises:
        FileNotFoundError: If the specified FASTA file cannot be located.
        ValueError: If no valid sequence could be extracted from the file.
    """
    # Check current directory and script directory
    target_path = fasta_path
    if not os.path.exists(target_path):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        alt_path = os.path.join(script_dir, fasta_path)
        if os.path.exists(alt_path):
            target_path = alt_path
        else:
            raise FileNotFoundError(
                f"FASTA file '{fasta_path}' not found at '{os.path.abspath(fasta_path)}' "
                f"or in script directory '{alt_path}'."
            )

    print(f"[*] Extracting genome sequence from '{target_path}'...")
    record = None
    try:
        record = SeqIO.read(target_path, "fasta")
    except Exception:
        # Fallback to SeqIO.parse in case of multiple entries or comments
        try:
            with open(target_path, "r", encoding="utf-8") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
                if records:
                    record = records[0]
        except Exception:
            pass

    if record is not None:
        seq_str = str(record.seq).strip().upper()
        header = record.description
        rec_id = record.id
    else:
        # Fallback raw FASTA parser
        header = ""
        rec_id = ""
        seq_lines = []
        with open(target_path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if not header:
                        header = line[1:].strip()
                        rec_id = header.split()[0] if header else ""
                else:
                    seq_lines.append(line)
        seq_str = "".join(seq_lines).upper()

    if not seq_str:
        raise ValueError(f"No valid sequence data could be extracted from '{target_path}'.")

    # Genomic statistics
    total_len = len(seq_str)
    gc_count = seq_str.count("G") + seq_str.count("C")
    gc_pct = (gc_count / total_len) * 100 if total_len > 0 else 0.0

    print(f"[+] Successfully extracted genome: {header or target_path}")
    print(f"    - Accession / ID: {rec_id}")
    print(f"    - Length: {total_len} bp")
    print(f"    - GC Content: {gc_pct:.2f}%")

    if return_record and record is not None:
        return record
    return seq_str


# Aliases for genome extraction
extract_genome = extract_genome_from_fasta
extract_genome_sequence = extract_genome_from_fasta
fetch_genome_sequence = extract_genome_from_fasta

def extract_all_frames_codons(sequence: str, include_reverse_complement: bool = True) -> dict:
    """
    Extracts codons across all 3 forward reading frames (+1, +2, +3)
    and optionally all 3 reverse-complement frames (-1, -2, -3).
    Returns detailed counts and mathematical metrics.
    """
    seq_obj = Seq(sequence.upper())
    forward_rna = str(seq_obj.transcribe())
    reverse_rna = str(seq_obj.reverse_complement().transcribe())

    frame_codons = {}
    total_counts = Counter()

    # Forward frames (+1, +2, +3)
    for offset in range(3):
        frame_name = f"Frame +{offset + 1}"
        codons = [
            forward_rna[i:i+3]
            for i in range(offset, len(forward_rna) - 2, 3)
            if len(forward_rna[i:i+3]) == 3 and set(forward_rna[i:i+3]).issubset({'A', 'U', 'C', 'G'})
        ]
        frame_codons[frame_name] = len(codons)
        total_counts.update(codons)

    # Reverse complement frames (-1, -2, -3)
    if include_reverse_complement:
        for offset in range(3):
            frame_name = f"Frame -{offset + 1}"
            codons = [
                reverse_rna[i:i+3]
                for i in range(offset, len(reverse_rna) - 2, 3)
                if len(reverse_rna[i:i+3]) == 3 and set(reverse_rna[i:i+3]).issubset({'A', 'U', 'C', 'G'})
            ]
            frame_codons[frame_name] = len(codons)
            total_counts.update(codons)

    total_codons = sum(total_counts.values())

    # Sort frequencies from highest to smallest
    sorted_codons = []
    for codon, count in total_counts.most_common():
        rel_freq = (count / total_codons) * 100 if total_codons > 0 else 0.0
        sorted_codons.append({
            "codon": codon,
            "amino_acid": codon_to_amino_acid(codon),
            "count": count,
            "frequency_pct": round(rel_freq, 4)
        })

    return {
        "sequence_length": len(sequence),
        "frame_breakdown": frame_codons,
        "total_codons_analyzed": total_codons,
        "sorted_frequencies": sorted_codons
    }


# =====================================================================
# PART 3: CHART VISUALIZATION (TOP 10 MOST FREQUENT CODONS)
# =====================================================================

def plot_top_10_codons(analysis_data: dict, accession_id: str = "NC_116972.1", output_path: str = "top10_codons.png"):
    """
    Plots a bar chart of the top 10 most frequent codons labeled with their
    associated amino acids and relative percentages.
    """
    top_10 = analysis_data["sorted_frequencies"][:10]
    labels = [f"{item['codon']}\n({item['amino_acid'].split(' ')[0]})" for item in top_10]
    counts = [item["count"] for item in top_10]
    pcts = [item["frequency_pct"] for item in top_10]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(labels, counts, color="#2b5c8f", edgecolor="black", alpha=0.88)

    plt.title(f"Top 10 Most Frequent Codons Across All 6 Reading Frames\nPotyvirus thevetiae ({accession_id})", fontsize=14, fontweight="bold")
    plt.xlabel("Codon & Associated Amino Acid", fontsize=12)
    plt.ylabel("Total Codon Count (All Frames)", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    # Annotate bars with count and relative frequency percentage
    for bar, count, pct in zip(bars, counts, pcts):
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + (max(counts) * 0.015),
            f"{count}\n({pct:.2f}%)",
            ha="center", va="bottom", fontsize=9, fontweight="bold"
        )

    plt.ylim(0, max(counts) * 1.18)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"[+] Saved Top 10 Codon Chart to '{output_path}'")
    plt.show()


def export_top_10_codons_to_xlsx(analysis_data: dict, output_path: str = "top10_codons.xlsx") -> str:
    """
    Exports the statistical data for only the Top 10 most frequent codons
    into an Excel (.xlsx) spreadsheet.

    Columns included:
      - Rank (1-10)
      - Codon (e.g., AAA, UUU)
      - Amino Acid (e.g., Lysine (Lys))
      - Count
      - Frequency (%)
      - Cumulative Frequency (%)

    Parameters:
        analysis_data (dict): The dictionary returned by extract_all_frames_codons.
        output_path (str): Filepath where the .xlsx file will be saved.

    Returns:
        str: Path to the saved Excel file.
    """
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Top 10 Codons"
    ws.views.sheetView[0].showGridLines = True

    headers = [
        "Rank",
        "Codon",
        "Amino Acid",
        "Count",
        "Frequency (%)",
        "Cumulative Frequency (%)"
    ]
    ws.append(headers)

    header_font = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
    header_fill = PatternFill(start_color="2B5C8F", end_color="2B5C8F", fill_type="solid")
    header_align = Alignment(horizontal="center", vertical="center")
    thin_border = Border(
        left=Side(style='thin', color='D9D9D9'),
        right=Side(style='thin', color='D9D9D9'),
        top=Side(style='thin', color='D9D9D9'),
        bottom=Side(style='thin', color='D9D9D9')
    )

    for col_num in range(1, len(headers) + 1):
        cell = ws.cell(row=1, column=col_num)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = header_align
        cell.border = thin_border

    top_10 = analysis_data["sorted_frequencies"][:10]
    cum_freq = 0.0

    for idx, item in enumerate(top_10, 1):
        freq = float(item["frequency_pct"])
        cum_freq += freq
        row_data = [
            idx,
            item["codon"],
            item["amino_acid"],
            int(item["count"]),
            round(freq, 4),
            round(cum_freq, 4)
        ]
        ws.append(row_data)
        row_num = idx + 1
        row_fill = PatternFill(start_color="F7F9FB", end_color="F7F9FB", fill_type="solid") if idx % 2 == 0 else None

        for col_idx in range(1, len(headers) + 1):
            c = ws.cell(row=row_num, column=col_idx)
            c.border = thin_border
            if row_fill:
                c.fill = row_fill
            if col_idx in (1, 2):
                c.alignment = Alignment(horizontal="center", vertical="center")
            elif col_idx == 3:
                c.alignment = Alignment(horizontal="left", vertical="center")
            else:
                c.alignment = Alignment(horizontal="right", vertical="center")

    # Auto-adjust column widths for optimal display
    for col in ws.columns:
        max_len = max(len(str(cell.value or '')) for cell in col)
        col_letter = get_column_letter(col[0].column)
        ws.column_dimensions[col_letter].width = max(max_len + 5, 12)

    wb.save(output_path)
    print(f"[+] Saved Top 10 Codon Table to '{output_path}'")
    return output_path


# Aliases for Excel export
save_top_10_codons_to_xlsx = export_top_10_codons_to_xlsx
save_top_10_codons_to_excel = export_top_10_codons_to_xlsx


# =====================================================================
# PART 4: LLM PROMPT BUILDER & API / LOCAL MODEL EXECUTION
# =====================================================================

def parse_top_10_codons_from_xlsx(xlsx_path: str = "top10_codons.xlsx") -> list:
    """
    Parses the Top 10 dominant codons from an Excel (.xlsx) file.

    Parameters:
        xlsx_path (str): Path to the Excel file.

    Returns:
        list[dict]: List of parsed codon dictionaries.
    """
    target_path = xlsx_path
    if not os.path.exists(target_path):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        alt_path = os.path.join(script_dir, xlsx_path)
        if os.path.exists(alt_path):
            target_path = alt_path
        else:
            raise FileNotFoundError(f"Excel file '{xlsx_path}' not found.")

    wb = openpyxl.load_workbook(target_path, data_only=True)
    ws = wb.active

    rows = list(ws.iter_rows(values_only=True))
    if not rows or len(rows) < 2:
        raise ValueError(f"No codon data rows found in '{target_path}'.")

    codons_data = []
    for r in rows[1:]:
        if not r or r[1] is None:
            continue
        codons_data.append({
            "rank": r[0],
            "codon": str(r[1]).strip(),
            "amino_acid": str(r[2]).strip(),
            "count": r[3],
            "frequency_pct": r[4],
            "cumulative_frequency_pct": r[5] if len(r) > 5 else None
        })

    print(f"[+] Successfully parsed {len(codons_data)} dominant codons from '{target_path}'")
    return codons_data


LLAMA_SYSTEM_ROLE = (
    "You are an expert computational biologist and clinical nutritionist specializing "
    "in viral genomics, codon translation, and targeted amino acid dietary management. "
    "Provide clear, evidence-based nutritional recommendations with concise justifications."
)


def build_llama_prompt_from_xlsx(
    xlsx_path: str = "top10_codons.xlsx",
    base_instruction: str = (
        "For a person infected with a virus that has the following 10 Codons dominant in its genome, "
        "suggest foods that would restrict or minimize the consumption of the respective amino acids from those 10 dominant Codons."
    )
) -> str:
    """
    Parses the top 10 codons from the Excel spreadsheet and constructs an optimized
    prompt for llama3.2:3b instructing the model to provide good foods to eat, foods to avoid,
    and brief justifications for both.

    Parameters:
        xlsx_path (str): Path to the Excel file containing top 10 codons data.
        base_instruction (str): Primary clinical instruction.

    Returns:
        str: Formatted user prompt.
    """
    top_10 = parse_top_10_codons_from_xlsx(xlsx_path)

    codon_lines = [
        f"  {item['rank']}. Codon: {item['codon']} | Amino Acid: {item['amino_acid']} "
        f"| Count: {item['count']} | Relative Frequency: {item['frequency_pct']}%"
        for item in top_10
    ]
    codon_summary_table = "\n".join(codon_lines)

    # Unique non-stop target amino acids
    target_amino_acids = []
    for item in top_10:
        aa = item['amino_acid']
        if aa and "Stop" not in aa and aa not in target_amino_acids:
            target_amino_acids.append(aa)

    prompt = (
        f"{base_instruction}\n\n"
        f"### Viral Genomic Context (Top 10 Dominant Codons parsed from '{os.path.basename(xlsx_path)}'):\n"
        f"{codon_summary_table}\n\n"
        f"### Target Amino Acids to Restrict/Minimize:\n"
        f"{', '.join(target_amino_acids)}\n\n"
        f"### Instructions & Output Format:\n"
        f"Synthesize your answer concisely into exactly two sections:\n\n"
        f"1. **Good Foods to Eat**:\n"
        f"   - Recommend general categories and specific examples of whole foods naturally lowest in these amino acids (such as fresh fruits, tubers, specific low-protein grains, and watery vegetables).\n"
        f"   - Provide a brief justification explaining why these foods help minimize the pool of amino acids needed for viral replication.\n\n"
        f"2. **Foods to Avoid / Minimize**:\n"
        f"   - List categories and examples of protein-dense foods that are rich in these amino acids (such as animal meats, poultry, dairy/cheeses, seafood, and dense legumes).\n"
        f"   - Provide a brief justification explaining why avoiding these foods restricts viral protein synthesis."
    )
    return prompt


# Aliases for prompt building
build_prompt_from_xlsx = build_llama_prompt_from_xlsx
parse_xlsx_to_prompt = build_llama_prompt_from_xlsx


def query_llama_local(
    prompt: str,
    system_role: str = LLAMA_SYSTEM_ROLE,
    model_name: str = "llama3.2:3b",
    ollama_url: str = "http://localhost:11434",
    timeout: int = 180,
    options: dict = None
) -> str:
    """
    Sends the prompt locally to Ollama using llama3.2:3b.
    100% offline and local (no external cloud APIs or OpenAI dependencies).

    Parameters:
        prompt (str): The user instruction prompt.
        system_role (str): System prompt establishing expertise and behavior.
        model_name (str): Local Ollama model name (default: "llama3.2:3b").
        ollama_url (str): Local Ollama endpoint (default: "http://localhost:11434").
        timeout (int): Request timeout in seconds.
        options (dict): Generation options (e.g., num_predict, temperature).

    Returns:
        str: Generated response from the local model.
    """
    if options is None:
        options = {
            "num_predict": 1500,
            "temperature": 0.3
        }

    chat_url = f"{ollama_url.rstrip('/')}/api/chat"
    chat_payload = json.dumps({
        "model": model_name,
        "messages": [
            {"role": "system", "content": system_role},
            {"role": "user", "content": prompt}
        ],
        "options": options,
        "stream": False
    }).encode("utf-8")

    chat_req = urllib.request.Request(chat_url, data=chat_payload, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(chat_req, timeout=timeout) as response:
            result = json.loads(response.read().decode("utf-8"))
            return result.get("message", {}).get("content", "").strip()
    except Exception as e:
        # Fallback to local /api/generate endpoint if chat encounters an issue
        gen_url = f"{ollama_url.rstrip('/')}/api/generate"
        gen_payload = json.dumps({
            "model": model_name,
            "prompt": f"{system_role}\n\n{prompt}",
            "options": options,
            "stream": False
        }).encode("utf-8")
        gen_req = urllib.request.Request(gen_url, data=gen_payload, headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(gen_req, timeout=timeout) as response:
                result = json.loads(response.read().decode("utf-8"))
                return result.get("response", "").strip()
        except Exception as gen_err:
            return (
                f"[Local Ollama Notice]: Could not reach local Ollama at {chat_url} ({e}) or {gen_url} ({gen_err}).\n"
                f"Ensure Ollama is running locally ('ollama serve') and '{model_name}' is available.\n\n"
                f"Generated Prompt Ready for Local Model:\n{'-'*50}\n{prompt}"
            )


# Aliases for backward compatibility
query_llm = query_llama_local


# =====================================================================
# MAIN EXECUTION PIPELINE
# =====================================================================

if __name__ == "__main__":
    FASTA_FILE = "potyvirus_thevetiae_OM263475.fasta"
    ACCESSION = "OM263475 / NC_116972.1"

    # 1. Extract genome sequence directly from the local FASTA file
    sequence = extract_genome_from_fasta(FASTA_FILE)

    # 2. Compute codon frequencies across all 6 reading frames
    analysis = extract_all_frames_codons(sequence, include_reverse_complement=True)

    # 3. Print Top 10 Summary Table
    print("\n--- TOP 10 MOST FREQUENT CODONS (ALL FRAMES) ---")
    for idx, item in enumerate(analysis["sorted_frequencies"][:10], 1):
        print(f"{idx:2d}. {item['codon']} -> {item['amino_acid']:20s} | Count: {item['count']:5d} ({item['frequency_pct']:.2f}%)")

    # 4. Save Top 10 Codons to Excel Table (.xlsx)
    excel_path = "top10_codons.xlsx"
    export_top_10_codons_to_xlsx(analysis, output_path=excel_path)

    # 5. Plot the Top 10 Chart
    plot_top_10_codons(analysis, ACCESSION)

    # 6. Parse Top 10 Codons from XLSX and Build Prompt for Llama 3.2 (3B)
    llama_prompt = build_llama_prompt_from_xlsx(excel_path)
    print("\n--- PROMPT PARSED FROM XLSX FOR LLAMA 3.2 (3B) ---")
    print(llama_prompt)

    # 7. Query Local Ollama llama3.2:3b model (100% Local)
    print("\n[*] Querying local Ollama with 'llama3.2:3b' model...")
    llm_response = query_llama_local(llama_prompt)
    print("\n--- LLAMA 3.2 (3B) DIETARY RECOMMENDATIONS ---")
    print(llm_response)