"""
Transposon Simulator Module
Bioinformatics Lab 6 - Simulation of Transposable Elements in Artificial DNA
"""

import random
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

# Nucleotide complement map
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(COMPLEMENT.get(base.upper(), 'N') for base in reversed(seq))


def random_dna(length: int, seed: Optional[int] = None) -> str:
    """Generate a random DNA sequence of specified length."""
    rng = random.Random(seed)
    return "".join(rng.choice(['A', 'C', 'G', 'T']) for _ in range(length))


@dataclass
class TransposableElement:
    """Represents a transposable element template."""
    name: str
    sequence: str
    tir_length: int = 8  # Length of Terminal Inverted Repeat
    tsd_length: int = 3  # Target site duplication length
    description: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def left_tir(self) -> str:
        """Left Terminal Inverted Repeat."""
        return self.sequence[:self.tir_length]

    @property
    def right_tir(self) -> str:
        """Right Terminal Inverted Repeat (reverse complement of left TIR)."""
        return self.sequence[-self.tir_length:]

    @classmethod
    def create(
        cls,
        name: str,
        total_length: int = 50,
        tir_length: int = 8,
        tsd_length: int = 3,
        seed: Optional[int] = None,
        description: str = ""
    ) -> 'TransposableElement':
        """
        Synthesize a realistic Class II (TIR) transposable element.
        5' - [TIR_L] + [Internal Sequence] + [TIR_R = rev_comp(TIR_L)] - 3'
        """
        rng = random.Random(seed)
        left_tir = "".join(rng.choice(['A', 'C', 'G', 'T']) for _ in range(tir_length))
        right_tir = reverse_complement(left_tir)

        internal_len = max(0, total_length - 2 * tir_length)
        internal_seq = "".join(rng.choice(['A', 'C', 'G', 'T']) for _ in range(internal_len))

        full_seq = left_tir + internal_seq + right_tir
        return cls(
            name=name,
            sequence=full_seq,
            tir_length=tir_length,
            tsd_length=tsd_length,
            description=description
        )


@dataclass
class GenomicInterval:
    """A segment of the simulated genome corresponding to a TE or fragment."""
    element_name: str
    start: int  # 0-based inclusive
    end: int    # 0-based inclusive (length = end - start + 1)
    fragment_type: str = "full"  # "full", "5_prime_fragment", "3_prime_fragment", "internal_fragment"
    parent_element: Optional[str] = None
    original_element_length: int = 0
    nested_in: Optional[str] = None
    interrupted_by: List[str] = field(default_factory=list)
    sequence: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def start_1based(self) -> int:
        return self.start + 1

    @property
    def end_1based(self) -> int:
        return self.end + 1


@dataclass
class SimulationResult:
    """Full outcome of a DNA transposable element simulation."""
    genome_sequence: str
    host_length: int
    final_length: int
    scenario_name: str
    transposon_library: Dict[str, TransposableElement]
    ground_truth_intervals: List[GenomicInterval]
    simulation_log: List[str]


class TransposonSimulator:
    """
    Simulates artificial DNA sequence (200-500 bp) with 3 to 4 transposable elements,
    supporting standard non-overlapping insertions and 2- or 3-transposon intersections.
    """

    def __init__(self, seed: Optional[int] = None):
        self.rng = random.Random(seed)

    def generate_library(self, count: int = 4) -> Dict[str, TransposableElement]:
        """
        Generate a library of 3 or 4 distinct transposable elements.
        Lengths are between 35 and 65 bp, fitting comfortably within 200-500 bp genomes.
        """
        count = max(3, min(4, count))
        configs = [
            ("TE_Mariner_1", 46, 8, 2, "Mariner-like DNA transposon (TA TSD, 8bp TIR)"),
            ("TE_PiggyBac_2", 50, 9, 4, "PiggyBac-like DNA transposon (TTAA TSD, 9bp TIR)"),
            ("TE_IS10_3", 42, 8, 3, "IS-like insertion sequence (3bp TSD, 8bp TIR)"),
            ("TE_Helitron_4", 40, 8, 2, "Helitron/MITE-like element (2bp TSD, 8bp TIR)")
        ]

        library = {}
        for i in range(count):
            name, length, tir_len, tsd_len, desc = configs[i]
            te_seed = self.rng.randint(1000, 999999)
            te = TransposableElement.create(
                name=name,
                total_length=length,
                tir_length=tir_len,
                tsd_length=tsd_len,
                seed=te_seed,
                description=desc
            )
            library[name] = te

        return library

    def simulate_scenario_non_intersecting(
        self,
        genome_length: int = 350,
        te_count: int = 3
    ) -> SimulationResult:
        """
        Scenario A: 3 or 4 transposable elements inserted independently (no intersection).
        Serves as the baseline/control scenario.
        """
        te_count = max(3, min(4, te_count))
        genome_length = max(200, min(500, genome_length))
        library = self.generate_library(te_count)
        elements = list(library.values())

        # Generate background host DNA
        # We need host_len + sum(te_len) = genome_length
        total_te_len = sum(te.length for te in elements)
        if total_te_len >= genome_length - 40:
            # Scale down or expand genome_length to fit
            genome_length = total_te_len + 60

        host_len = genome_length - total_te_len
        host_dna = random_dna(host_len, seed=self.rng.randint(1, 999999))

        # Determine non-overlapping insertion points along the host DNA
        # Split host DNA into te_count + 1 segments
        seg_sizes = [host_len // (te_count + 1)] * (te_count + 1)
        remainder = host_len % (te_count + 1)
        for i in range(remainder):
            seg_sizes[i] += 1

        sim_log = [f"Initial host DNA generated ({host_len} bp)."]
        ground_truth: List[GenomicInterval] = []

        # Assemble sequence with insertions
        current_seq = ""
        current_pos = 0

        for i, te in enumerate(elements):
            chunk = host_dna[current_pos : current_pos + seg_sizes[i]]
            current_seq += chunk
            current_pos += seg_sizes[i]

            start_idx = len(current_seq)
            current_seq += te.sequence
            end_idx = len(current_seq) - 1

            ground_truth.append(GenomicInterval(
                element_name=te.name,
                start=start_idx,
                end=end_idx,
                fragment_type="full",
                original_element_length=te.length,
                sequence=te.sequence
            ))
            sim_log.append(
                f"Inserted {te.name} ({te.length} bp) at position [{start_idx}:{end_idx}] (intact)."
            )

        # Append remaining host chunk
        current_seq += host_dna[current_pos:]
        sim_log.append(f"Final synthetic sequence length: {len(current_seq)} bp.")

        return SimulationResult(
            genome_sequence=current_seq,
            host_length=host_len,
            final_length=len(current_seq),
            scenario_name=f"Non-Intersecting ({te_count} TEs)",
            transposon_library=library,
            ground_truth_intervals=ground_truth,
            simulation_log=sim_log
        )

    def simulate_scenario_2_intersecting_nested(
        self,
        genome_length: int = 350,
        te_count: int = 3
    ) -> SimulationResult:
        """
        Scenario B: 2 Transposons Intersecting via Nested Insertion.
        TE2 inserts into the interior of TE1, splitting TE1 into two fragments:
        - TE1 (5' fragment)
        - TE2 (full, inserted inside TE1)
        - TE1 (3' fragment)
        TE3 (and TE4 if requested) are inserted independently.
        """
        te_count = max(3, min(4, te_count))
        genome_length = max(200, min(500, genome_length))
        library = self.generate_library(te_count)
        elements = list(library.values())

        te1, te2 = elements[0], elements[1]
        other_tes = elements[2:]

        total_te_len = sum(te.length for te in elements)
        if total_te_len >= genome_length - 40:
            genome_length = total_te_len + 60

        host_len = genome_length - total_te_len
        host_dna = random_dna(host_len, seed=self.rng.randint(1, 999999))

        sim_log = [f"Initial host DNA generated ({host_len} bp)."]

        # Step 1: Create nested compound element [TE1_partA] + [TE2] + [TE1_partB]
        # Choose insertion site inside TE1 (e.g., around 40-60% of TE1)
        split_pos = te1.length // 2
        te1_part_a = te1.sequence[:split_pos]
        te1_part_b = te1.sequence[split_pos:]

        nested_compound = te1_part_a + te2.sequence + te1_part_b

        sim_log.append(
            f"Nested intersection event: {te2.name} inserted into {te1.name} at relative position {split_pos}."
        )

        # Place compound element and other_tes
        insertion_units = [("COMPOUND_1_2", nested_compound)] + [(te.name, te.sequence) for te in other_tes]
        num_units = len(insertion_units)

        seg_sizes = [host_len // (num_units + 1)] * (num_units + 1)
        for i in range(host_len % (num_units + 1)):
            seg_sizes[i] += 1

        current_seq = ""
        current_pos = 0
        ground_truth: List[GenomicInterval] = []

        for i, (uname, useq) in enumerate(insertion_units):
            chunk = host_dna[current_pos : current_pos + seg_sizes[i]]
            current_seq += chunk
            current_pos += seg_sizes[i]

            u_start = len(current_seq)

            if uname == "COMPOUND_1_2":
                # TE1 5' fragment
                te1_a_start = u_start
                te1_a_end = te1_a_start + len(te1_part_a) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te1.name,
                    start=te1_a_start,
                    end=te1_a_end,
                    fragment_type="5_prime_fragment",
                    parent_element=te1.name,
                    original_element_length=te1.length,
                    interrupted_by=[te2.name],
                    sequence=te1_part_a
                ))

                # TE2 nested
                te2_start = te1_a_end + 1
                te2_end = te2_start + te2.length - 1
                ground_truth.append(GenomicInterval(
                    element_name=te2.name,
                    start=te2_start,
                    end=te2_end,
                    fragment_type="full",
                    nested_in=te1.name,
                    original_element_length=te2.length,
                    sequence=te2.sequence
                ))

                # TE1 3' fragment
                te1_b_start = te2_end + 1
                te1_b_end = te1_b_start + len(te1_part_b) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te1.name,
                    start=te1_b_start,
                    end=te1_b_end,
                    fragment_type="3_prime_fragment",
                    parent_element=te1.name,
                    original_element_length=te1.length,
                    interrupted_by=[te2.name],
                    sequence=te1_part_b
                ))

                current_seq += useq
                sim_log.append(
                    f"Placed nested locus: {te1.name}[5'] at [{te1_a_start}:{te1_a_end}], "
                    f"{te2.name} at [{te2_start}:{te2_end}], "
                    f"{te1.name}[3'] at [{te1_b_start}:{te1_b_end}]."
                )
            else:
                te_obj = library[uname]
                current_seq += useq
                u_end = len(current_seq) - 1
                ground_truth.append(GenomicInterval(
                    element_name=uname,
                    start=u_start,
                    end=u_end,
                    fragment_type="full",
                    original_element_length=te_obj.length,
                    sequence=useq
                ))
                sim_log.append(f"Placed independent {uname} at [{u_start}:{u_end}].")

        current_seq += host_dna[current_pos:]
        sim_log.append(f"Final synthetic sequence length: {len(current_seq)} bp.")

        return SimulationResult(
            genome_sequence=current_seq,
            host_length=host_len,
            final_length=len(current_seq),
            scenario_name=f"2-Transposon Nested Intersection ({te1.name} intersect {te2.name})",
            transposon_library=library,
            ground_truth_intervals=ground_truth,
            simulation_log=sim_log
        )

    def simulate_scenario_2_intersecting_overlapping(
        self,
        genome_length: int = 350,
        te_count: int = 3
    ) -> SimulationResult:
        """
        Scenario C: 2 Transposons Overlapping (Partial recombination / overwrite).
        TE2 overlaps the 3' end of TE1 by overlap_len base pairs.
        Demonstrates boundary competition and partial truncation.
        """
        te_count = max(3, min(4, te_count))
        genome_length = max(200, min(500, genome_length))
        library = self.generate_library(te_count)
        elements = list(library.values())

        te1, te2 = elements[0], elements[1]
        other_tes = elements[2:]

        overlap_len = 16  # Moderate overlap
        # Construct overlapping region: TE1 prefix + TE2 full (overwriting TE1 suffix)
        te1_prefix = te1.sequence[:-overlap_len]
        overlap_seq = te2.sequence[:overlap_len]
        # In this model, TE2 inserted and displaced the last 16 bp of TE1
        compound = te1_prefix + te2.sequence

        total_te_len = len(compound) + sum(te.length for te in other_tes)
        host_len = max(50, genome_length - total_te_len)
        host_dna = random_dna(host_len, seed=self.rng.randint(1, 999999))

        sim_log = [
            f"Initial host DNA generated ({host_len} bp).",
            f"Overlapping intersection: {te2.name} overlaps last {overlap_len} bp of {te1.name}."
        ]

        insertion_units = [("OVERLAP_1_2", compound)] + [(te.name, te.sequence) for te in other_tes]
        num_units = len(insertion_units)

        seg_sizes = [host_len // (num_units + 1)] * (num_units + 1)
        for i in range(host_len % (num_units + 1)):
            seg_sizes[i] += 1

        current_seq = ""
        current_pos = 0
        ground_truth: List[GenomicInterval] = []

        for i, (uname, useq) in enumerate(insertion_units):
            chunk = host_dna[current_pos : current_pos + seg_sizes[i]]
            current_seq += chunk
            current_pos += seg_sizes[i]

            u_start = len(current_seq)

            if uname == "OVERLAP_1_2":
                te1_start = u_start
                te1_end = u_start + len(te1_prefix) - 1
                te2_start = te1_end + 1
                te2_end = te2_start + te2.length - 1

                ground_truth.append(GenomicInterval(
                    element_name=te1.name,
                    start=te1_start,
                    end=te1_end,
                    fragment_type="5_prime_fragment",
                    parent_element=te1.name,
                    original_element_length=te1.length,
                    interrupted_by=[te2.name],
                    sequence=te1_prefix
                ))

                ground_truth.append(GenomicInterval(
                    element_name=te2.name,
                    start=te2_start,
                    end=te2_end,
                    fragment_type="full",
                    original_element_length=te2.length,
                    sequence=te2.sequence
                ))

                current_seq += useq
                sim_log.append(
                    f"Placed overlapping locus: {te1.name} (truncated prefix) at [{te1_start}:{te1_end}], "
                    f"{te2.name} at [{te2_start}:{te2_end}]."
                )
            else:
                te_obj = library[uname]
                current_seq += useq
                u_end = len(current_seq) - 1
                ground_truth.append(GenomicInterval(
                    element_name=uname,
                    start=u_start,
                    end=u_end,
                    fragment_type="full",
                    original_element_length=te_obj.length,
                    sequence=useq
                ))
                sim_log.append(f"Placed independent {uname} at [{u_start}:{u_end}].")

        current_seq += host_dna[current_pos:]
        sim_log.append(f"Final synthetic sequence length: {len(current_seq)} bp.")

        return SimulationResult(
            genome_sequence=current_seq,
            host_length=host_len,
            final_length=len(current_seq),
            scenario_name=f"2-Transposon Overlapping Intersection ({te1.name} overlap {te2.name})",
            transposon_library=library,
            ground_truth_intervals=ground_truth,
            simulation_log=sim_log
        )

    def simulate_scenario_3_intersecting_nested(
        self,
        genome_length: int = 400,
        te_count: int = 4
    ) -> SimulationResult:
        """
        Scenario D: 3 Transposons Intersecting via Multi-level Nesting.
        TE2 inserts into TE1.
        Then TE3 inserts into TE2!
        Resulting hierarchy:
        [TE1_5'] -> [TE2_5'] -> [TE3_full] -> [TE2_3'] -> [TE1_3']
        TE4 (if present) is inserted independently.
        This provides a rich benchmark for how algorithms behave in complex nested intersections!
        """
        te_count = max(3, min(4, te_count))
        genome_length = max(250, min(500, genome_length))
        library = self.generate_library(te_count)
        elements = list(library.values())

        te1, te2, te3 = elements[0], elements[1], elements[2]
        te4 = elements[3] if len(elements) > 3 else None

        # Split TE2
        split_te2 = te2.length // 2
        te2_part_a = te2.sequence[:split_te2]
        te2_part_b = te2.sequence[split_te2:]

        # Compound of TE2 + TE3
        compound_2_3 = te2_part_a + te3.sequence + te2_part_b

        # Split TE1
        split_te1 = te1.length // 2
        te1_part_a = te1.sequence[:split_te1]
        te1_part_b = te1.sequence[split_te1:]

        # Full triple nested compound
        compound_1_2_3 = te1_part_a + compound_2_3 + te1_part_b

        total_te_len = len(compound_1_2_3) + (te4.length if te4 else 0)
        host_len = max(40, genome_length - total_te_len)
        host_dna = random_dna(host_len, seed=self.rng.randint(1, 999999))

        sim_log = [
            f"Initial host DNA generated ({host_len} bp).",
            f"3-way nested intersection: {te3.name} inserted into {te2.name}, which was inserted into {te1.name}."
        ]

        insertion_units = [("TRIPLE_NESTED", compound_1_2_3)]
        if te4:
            insertion_units.append((te4.name, te4.sequence))

        num_units = len(insertion_units)
        seg_sizes = [host_len // (num_units + 1)] * (num_units + 1)
        for i in range(host_len % (num_units + 1)):
            seg_sizes[i] += 1

        current_seq = ""
        current_pos = 0
        ground_truth: List[GenomicInterval] = []

        for i, (uname, useq) in enumerate(insertion_units):
            chunk = host_dna[current_pos : current_pos + seg_sizes[i]]
            current_seq += chunk
            current_pos += seg_sizes[i]

            u_start = len(current_seq)

            if uname == "TRIPLE_NESTED":
                # TE1 5'
                te1_a_start = u_start
                te1_a_end = te1_a_start + len(te1_part_a) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te1.name,
                    start=te1_a_start,
                    end=te1_a_end,
                    fragment_type="5_prime_fragment",
                    parent_element=te1.name,
                    original_element_length=te1.length,
                    interrupted_by=[te2.name, te3.name],
                    sequence=te1_part_a
                ))

                # TE2 5'
                te2_a_start = te1_a_end + 1
                te2_a_end = te2_a_start + len(te2_part_a) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te2.name,
                    start=te2_a_start,
                    end=te2_a_end,
                    fragment_type="5_prime_fragment",
                    parent_element=te2.name,
                    nested_in=te1.name,
                    original_element_length=te2.length,
                    interrupted_by=[te3.name],
                    sequence=te2_part_a
                ))

                # TE3 Full
                te3_start = te2_a_end + 1
                te3_end = te3_start + te3.length - 1
                ground_truth.append(GenomicInterval(
                    element_name=te3.name,
                    start=te3_start,
                    end=te3_end,
                    fragment_type="full",
                    nested_in=te2.name,
                    original_element_length=te3.length,
                    sequence=te3.sequence
                ))

                # TE2 3'
                te2_b_start = te3_end + 1
                te2_b_end = te2_b_start + len(te2_part_b) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te2.name,
                    start=te2_b_start,
                    end=te2_b_end,
                    fragment_type="3_prime_fragment",
                    parent_element=te2.name,
                    nested_in=te1.name,
                    original_element_length=te2.length,
                    interrupted_by=[te3.name],
                    sequence=te2_part_b
                ))

                # TE1 3'
                te1_b_start = te2_b_end + 1
                te1_b_end = te1_b_start + len(te1_part_b) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te1.name,
                    start=te1_b_start,
                    end=te1_b_end,
                    fragment_type="3_prime_fragment",
                    parent_element=te1.name,
                    original_element_length=te1.length,
                    interrupted_by=[te2.name, te3.name],
                    sequence=te1_part_b
                ))

                current_seq += useq
                sim_log.append(
                    f"Placed triple nested structure: "
                    f"{te1.name}[5'] @ [{te1_a_start}:{te1_a_end}] -> "
                    f"{te2.name}[5'] @ [{te2_a_start}:{te2_a_end}] -> "
                    f"{te3.name}[full] @ [{te3_start}:{te3_end}] -> "
                    f"{te2.name}[3'] @ [{te2_b_start}:{te2_b_end}] -> "
                    f"{te1.name}[3'] @ [{te1_b_start}:{te1_b_end}]."
                )
            else:
                current_seq += useq
                u_end = len(current_seq) - 1
                ground_truth.append(GenomicInterval(
                    element_name=te4.name,
                    start=u_start,
                    end=u_end,
                    fragment_type="full",
                    original_element_length=te4.length,
                    sequence=useq
                ))
                sim_log.append(f"Placed independent {te4.name} at [{u_start}:{u_end}].")

        current_seq += host_dna[current_pos:]
        sim_log.append(f"Final synthetic sequence length: {len(current_seq)} bp.")

        return SimulationResult(
            genome_sequence=current_seq,
            host_length=host_len,
            final_length=len(current_seq),
            scenario_name=f"3-Transposon Nested Intersection ({te1.name} contains {te2.name} contains {te3.name})",
            transposon_library=library,
            ground_truth_intervals=ground_truth,
            simulation_log=sim_log
        )
