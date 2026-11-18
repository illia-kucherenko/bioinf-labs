dna_sequence = "ACGTTGCATTACGTTGCAGT"
frequencies = {}

for nucleotide in dna_sequence:
    frequencies[nucleotide] = frequencies.get(nucleotide, 0) + 1

print("Nucleotide Frequencies:")
for nucleotide, count in frequencies.items():
    print(f"{nucleotide}: {count}")