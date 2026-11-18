import math


def calculate_both_tms(dna_sequence, na_molar=0.3):
    length = len(dna_sequence)
    frequencies = {}
    for nucleotide in dna_sequence:
        frequencies[nucleotide] = frequencies.get(nucleotide, 0) + 1

    count_g = frequencies.get('G', 0)
    count_c = frequencies.get('C', 0)
    count_a = frequencies.get('A', 0)
    count_t = frequencies.get('T', 0)

    tm_simple = 4 * (count_g + count_c) + 2 * (count_a + count_t)

    percent_gc = ((count_g + count_c) / length) * 100

    tm_salt = 81.5 + (16.6 * math.log10(na_molar)) + (0.41 * percent_gc) - (600 / length)

    print(f"Sequence Length: {length} nucleotides")
    print(f"GC Content: {percent_gc:.1f}%")
    print(f"Salt Concentration ([Na+]): {na_molar} M")
    print("-" * 40)
    print(f"Simple Tm (Wallace): {tm_simple}°C")
    print(f"Salt-Adjusted Tm:    {tm_salt:.2f}°C")

    return tm_simple, tm_salt


sequence_s = input("Enter a DNA sequence: ").upper()
print(f"\nAnalyzing Sequence: {sequence_s}\n")

calculate_both_tms(sequence_s)