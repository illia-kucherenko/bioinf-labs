"""
Transposon Detection GUI Application
Bioinformatics Lab 6 - Desktop GUI for Transposon Simulation, Detection & Intersection Analysis
Supports opening FASTA files from the fasta/ folder (or via browse) with configurable
in-between spacing (default: minimize in-between to 10 bases).
"""

import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import random
from typing import Dict, List, Optional, Tuple

from transposon_simulator import (
    TransposonSimulator,
    SimulationResult,
    GenomicInterval,
    TransposableElement
)
from transposon_detector import (
    TransposonDetector,
    DetectionResult
)
from evaluation import IntersectionAnalyzer, BenchmarkSummary
from fasta_io import read_fasta_file, list_fasta_folder, FastaRecord

# Color palette for transposons and host DNA
TE_COLORS = {
    "TE_Mariner_1": "#3498DB",   # Sky Blue
    "TE_PiggyBac_2": "#E67E22",  # Warm Orange
    "TE_IS10_3": "#2ECC71",      # Emerald Green
    "TE_Helitron_4": "#9B59B6",  # Amethyst Purple
    "New_Putative_TE": "#16A085",# Teal for de novo FASTA candidates
    "Host": "#BDC3C7",           # Silver Grey
    "Default": "#E74C3C"         # Red fallback
}


def get_color_for_element(el_name: str) -> str:
    for key, color in TE_COLORS.items():
        if key in el_name:
            return color
    return TE_COLORS["Default"]


class TransposonApp(tk.Tk):
    """Main Application Window for Transposon Simulation, FASTA analysis, and Detection."""

    def __init__(self):
        super().__init__()
        self.title("Bioinformatics Lab 6 - Transposable Element Simulation & FASTA Detection App")
        self.geometry("1180x870")
        self.minsize(980, 720)

        # Simulator and Detector instances
        self.simulator = TransposonSimulator()
        self.detector = TransposonDetector(min_between=10)

        # State
        self.app_mode: str = "SIMULATION"  # "SIMULATION" or "FASTA"
        self.current_simulation: Optional[SimulationResult] = None
        self.current_fasta_record: Optional[FastaRecord] = None
        self.current_fasta_records: List[FastaRecord] = []
        self.current_detections: Dict[str, List[DetectionResult]] = {}
        self.current_benchmark: Optional[BenchmarkSummary] = None

        self._setup_style()
        self._build_ui()

        # Run initial default simulation
        self.run_simulation_and_detection()

    def _setup_style(self):
        self.style = ttk.Style(self)
        self.style.theme_use('clam')
        self.configure(bg="#F4F6F7")

    def _build_ui(self):
        # 1. Header Frame
        header = tk.Frame(self, bg="#2C3E50", padx=15, pady=10)
        header.pack(fill=tk.X)

        title = tk.Label(
            header,
            text="DNA Transposable Element Simulation, FASTA Analysis & Detection Engine",
            font=("Segoe UI", 14, "bold"),
            fg="#FFFFFF",
            bg="#2C3E50"
        )
        title.pack(anchor=tk.W)

        subtitle = tk.Label(
            header,
            text="Simulate 200-500 bp genomes with nested intersections OR Load FASTA files with >= 10 bp in-between constraint",
            font=("Segoe UI", 9),
            fg="#BDC3C7",
            bg="#2C3E50"
        )
        subtitle.pack(anchor=tk.W)

        # 2. Source Mode Selector (Simulation vs FASTA)
        mode_frame = tk.Frame(self, bg="#EAEDED", padx=12, pady=6)
        mode_frame.pack(fill=tk.X, padx=12, pady=(6, 2))

        tk.Label(mode_frame, text="Source Mode:", font=("Segoe UI", 9, "bold"), bg="#EAEDED").pack(side=tk.LEFT, padx=(0, 10))
        self.mode_var = tk.StringVar(value="SIMULATION")
        r_sim = tk.Radiobutton(
            mode_frame, text="Artificial DNA Simulation", variable=self.mode_var,
            value="SIMULATION", font=("Segoe UI", 9), bg="#EAEDED", command=self._on_mode_change
        )
        r_sim.pack(side=tk.LEFT, padx=(0, 15))

        r_fasta = tk.Radiobutton(
            mode_frame, text="Load FASTA File (from fasta/ folder or browse)", variable=self.mode_var,
            value="FASTA", font=("Segoe UI", 9), bg="#EAEDED", command=self._on_mode_change
        )
        r_fasta.pack(side=tk.LEFT, padx=(0, 15))

        # 3. Control Parameter Containers
        self.controls_box = ttk.LabelFrame(self, text=" Simulation Controls ", padding=10)
        self.controls_box.pack(fill=tk.X, padx=12, pady=4)

        # Subframe for Simulation Parameters
        self.sim_panel = tk.Frame(self.controls_box)
        self.sim_panel.pack(fill=tk.X)

        sim_row1 = tk.Frame(self.sim_panel)
        sim_row1.pack(fill=tk.X, pady=2)

        tk.Label(sim_row1, text="Genome Length (bp):", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.len_var = tk.IntVar(value=350)
        self.len_slider = tk.Scale(
            sim_row1, from_=200, to=500, orient=tk.HORIZONTAL, variable=self.len_var, length=160
        )
        self.len_slider.pack(side=tk.LEFT, padx=(0, 15))

        tk.Label(sim_row1, text="TE Count:", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.te_count_var = tk.IntVar(value=4)
        te_spin = ttk.Combobox(sim_row1, textvariable=self.te_count_var, values=[3, 4], width=4, state="readonly")
        te_spin.pack(side=tk.LEFT, padx=(0, 15))

        tk.Label(sim_row1, text="Intersection Scenario:", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.scenario_var = tk.StringVar(value="Scenario B: 2-TE Nested Intersection")
        scenarios = [
            "Scenario A: Non-Intersecting (Control)",
            "Scenario B: 2-TE Nested Intersection",
            "Scenario C: 2-TE Overlapping Intersection",
            "Scenario D: 3-TE Nested Chain Intersection"
        ]
        sc_menu = ttk.Combobox(sim_row1, textvariable=self.scenario_var, values=scenarios, width=34, state="readonly")
        sc_menu.pack(side=tk.LEFT, padx=(0, 10))

        sim_row2 = tk.Frame(self.sim_panel)
        sim_row2.pack(fill=tk.X, pady=(6, 2))

        tk.Label(sim_row2, text="Detection Strategy:", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.det_mode_var = tk.StringVar(value="All (Consensus)")
        det_modes = ["Homology Local Alignment", "Structural De Novo TIR", "All (Consensus)"]
        det_menu = ttk.Combobox(sim_row2, textvariable=self.det_mode_var, values=det_modes, width=22, state="readonly")
        det_menu.pack(side=tk.LEFT, padx=(0, 15))
        det_menu.bind("<<ComboboxSelected>>", lambda e: self._refresh_views())

        btn_run = tk.Button(
            sim_row2, text="[>] Run Simulation & Detect", font=("Segoe UI", 9, "bold"),
            bg="#27AE60", fg="white", padx=10, pady=3, relief=tk.FLAT,
            command=self.run_simulation_and_detection
        )
        btn_run.pack(side=tk.LEFT, padx=(0, 8))

        btn_rand = tk.Button(
            sim_row2, text="[*] New Random Seed", font=("Segoe UI", 9),
            bg="#34495E", fg="white", padx=8, pady=3, relief=tk.FLAT,
            command=self.randomize_seed
        )
        btn_rand.pack(side=tk.LEFT, padx=(0, 8))

        btn_export = tk.Button(
            sim_row2, text="[Save] Export Report", font=("Segoe UI", 9),
            bg="#7F8C8D", fg="white", padx=8, pady=3, relief=tk.FLAT,
            command=self.export_report
        )
        btn_export.pack(side=tk.RIGHT)

        # Subframe for FASTA Controls (initially hidden)
        self.fasta_panel = tk.Frame(self.controls_box)

        fasta_row1 = tk.Frame(self.fasta_panel)
        fasta_row1.pack(fill=tk.X, pady=2)

        tk.Label(fasta_row1, text="FASTA File (fasta/ folder):", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.fasta_file_var = tk.StringVar()
        self.fasta_combo = ttk.Combobox(fasta_row1, textvariable=self.fasta_file_var, width=32, state="readonly")
        self.fasta_combo.pack(side=tk.LEFT, padx=(0, 6))
        self.fasta_combo.bind("<<ComboboxSelected>>", self._on_fasta_file_selected)

        btn_refresh = tk.Button(
            fasta_row1, text="Refresh fasta/", font=("Segoe UI", 8),
            bg="#BDC3C7", relief=tk.FLAT, padx=6, pady=2, command=self._refresh_fasta_list
        )
        btn_refresh.pack(side=tk.LEFT, padx=(0, 8))

        btn_browse = tk.Button(
            fasta_row1, text="Browse...", font=("Segoe UI", 8),
            bg="#BDC3C7", relief=tk.FLAT, padx=6, pady=2, command=self._browse_fasta_file
        )
        btn_browse.pack(side=tk.LEFT, padx=(0, 15))

        tk.Label(fasta_row1, text="Sequence ID:", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.fasta_seq_var = tk.StringVar()
        self.seq_combo = ttk.Combobox(fasta_row1, textvariable=self.fasta_seq_var, width=22, state="readonly")
        self.seq_combo.pack(side=tk.LEFT, padx=(0, 10))
        self.seq_combo.bind("<<ComboboxSelected>>", lambda e: self.run_fasta_detection())

        fasta_row2 = tk.Frame(self.fasta_panel)
        fasta_row2.pack(fill=tk.X, pady=(6, 2))

        # In-between bases constraint (default 10 bases)
        tk.Label(fasta_row2, text="Min Bases In-Between (bp):", font=("Segoe UI", 9, "bold"), fg="#C0392B").pack(side=tk.LEFT, padx=(0, 5))
        self.min_between_var = tk.IntVar(value=10)
        between_spin = tk.Spinbox(fasta_row2, from_=5, to=100, textvariable=self.min_between_var, width=5)
        between_spin.pack(side=tk.LEFT, padx=(0, 15))

        # Min TIR length
        tk.Label(fasta_row2, text="Min TIR Length (bp):", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.min_tir_var = tk.IntVar(value=8)
        tir_spin = tk.Spinbox(fasta_row2, from_=5, to=16, textvariable=self.min_tir_var, width=5)
        tir_spin.pack(side=tk.LEFT, padx=(0, 20))

        btn_detect_fasta = tk.Button(
            fasta_row2, text="[>] Detect Candidate Transposons", font=("Segoe UI", 9, "bold"),
            bg="#16A085", fg="white", padx=12, pady=3, relief=tk.FLAT,
            command=self.run_fasta_detection
        )
        btn_detect_fasta.pack(side=tk.LEFT, padx=(0, 10))

        btn_export_fasta = tk.Button(
            fasta_row2, text="[Save] Export Report", font=("Segoe UI", 9),
            bg="#7F8C8D", fg="white", padx=8, pady=3, relief=tk.FLAT,
            command=self.export_report
        )
        btn_export_fasta.pack(side=tk.RIGHT)

        self._refresh_fasta_list()

        # 4. Visual Track Canvas
        vis_frame = ttk.LabelFrame(self, text=" Genomic Map Visualization (Sequence & Transposon Coordinates) ", padding=8)
        vis_frame.pack(fill=tk.X, padx=12, pady=4)

        self.canvas = tk.Canvas(vis_frame, height=140, bg="#FFFFFF", highlightthickness=1, highlightbackground="#D5DBDB")
        self.canvas.pack(fill=tk.X, expand=True)
        self.canvas.bind("<Configure>", lambda e: self._draw_tracks())

        # 5. Results Notebook (Tabs)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=12, pady=(4, 10))

        # Tab 1: Detections Table
        tab_table = ttk.Frame(self.notebook)
        self.notebook.add(tab_table, text=" Detected Transposons (Start, End) ")
        self._build_table_tab(tab_table)

        # Tab 2: Intersection Behavior Analysis
        tab_behavior = ttk.Frame(self.notebook)
        self.notebook.add(tab_behavior, text=" Intersection & De Novo Behavior Report ")
        self._build_behavior_tab(tab_behavior)

        # Tab 3: DNA Sequence Viewer
        tab_seq = ttk.Frame(self.notebook)
        self.notebook.add(tab_seq, text=" DNA Sequence Viewer ")
        self._build_sequence_tab(tab_seq)

        # Tab 4: Benchmark & Ground Truth Comparison
        tab_bench = ttk.Frame(self.notebook)
        self.notebook.add(tab_bench, text=" Ground Truth Comparison & Metrics ")
        self._build_benchmark_tab(tab_bench)

        # Status Bar
        self.status_var = tk.StringVar(value="Ready.")
        status_bar = tk.Label(self, textvariable=self.status_var, font=("Segoe UI", 8), anchor=tk.W, bg="#EAEDED", padx=10, pady=2)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def _build_table_tab(self, parent):
        cols = ("name", "status", "start_0", "end_0", "start_1", "end_1", "length", "in_between", "left_tir", "right_tir", "tsd", "method", "score", "details")
        self.tree = ttk.Treeview(parent, columns=cols, show="headings", height=8)

        self.tree.heading("name", text="Element Name")
        self.tree.heading("status", text="Status / Type")
        self.tree.heading("start_0", text="Start (0-idx)")
        self.tree.heading("end_0", text="End (0-idx)")
        self.tree.heading("start_1", text="Start (1-idx)")
        self.tree.heading("end_1", text="End (1-idx)")
        self.tree.heading("length", text="Total (bp)")
        self.tree.heading("in_between", text="In-Between (bp)")
        self.tree.heading("left_tir", text="Left TIR")
        self.tree.heading("right_tir", text="Right TIR")
        self.tree.heading("tsd", text="TSD")
        self.tree.heading("method", text="Detection Method")
        self.tree.heading("score", text="Score")
        self.tree.heading("details", text="Biological / Alignment Details")

        self.tree.column("name", width=140, anchor=tk.W)
        self.tree.column("status", width=120, anchor=tk.CENTER)
        self.tree.column("start_0", width=75, anchor=tk.CENTER)
        self.tree.column("end_0", width=75, anchor=tk.CENTER)
        self.tree.column("start_1", width=75, anchor=tk.CENTER)
        self.tree.column("end_1", width=75, anchor=tk.CENTER)
        self.tree.column("length", width=65, anchor=tk.CENTER)
        self.tree.column("in_between", width=90, anchor=tk.CENTER)
        self.tree.column("left_tir", width=90, anchor=tk.CENTER)
        self.tree.column("right_tir", width=90, anchor=tk.CENTER)
        self.tree.column("tsd", width=55, anchor=tk.CENTER)
        self.tree.column("method", width=140, anchor=tk.W)
        self.tree.column("score", width=55, anchor=tk.CENTER)
        self.tree.column("details", width=260, anchor=tk.W)

        tree_scroll_y = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.tree.yview)
        tree_scroll_x = ttk.Scrollbar(parent, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=tree_scroll_y.set, xscrollcommand=tree_scroll_x.set)

        tree_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        tree_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.tree.pack(fill=tk.BOTH, expand=True)

    def _build_behavior_tab(self, parent):
        self.behavior_text = tk.Text(parent, wrap=tk.WORD, font=("Consolas", 10), bg="#FBFCFC", padx=10, pady=10)
        scroll = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.behavior_text.yview)
        self.behavior_text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.behavior_text.pack(fill=tk.BOTH, expand=True)

    def _build_sequence_tab(self, parent):
        self.seq_text = tk.Text(parent, wrap=tk.WORD, font=("Consolas", 10), bg="#2B2B2B", fg="#A9B7C6", padx=10, pady=10)
        scroll = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.seq_text.yview)
        self.seq_text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.seq_text.pack(fill=tk.BOTH, expand=True)

    def _build_benchmark_tab(self, parent):
        self.bench_text = tk.Text(parent, wrap=tk.WORD, font=("Consolas", 10), bg="#FBFCFC", padx=10, pady=10)
        scroll = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.bench_text.yview)
        self.bench_text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.bench_text.pack(fill=tk.BOTH, expand=True)

    def _on_mode_change(self):
        mode = self.mode_var.get()
        self.app_mode = mode
        if mode == "SIMULATION":
            self.fasta_panel.pack_forget()
            self.sim_panel.pack(fill=tk.X)
            self.controls_box.configure(text=" Simulation Controls ")
            if self.current_simulation:
                self._refresh_views()
            else:
                self.run_simulation_and_detection()
        else:
            self.sim_panel.pack_forget()
            self.fasta_panel.pack(fill=tk.X)
            self.controls_box.configure(text=" FASTA File Controls (fasta/ folder) ")
            self._refresh_fasta_list()
            if self.current_fasta_record:
                self.run_fasta_detection()
            elif self.fasta_combo["values"]:
                self.fasta_combo.current(0)
                self._on_fasta_file_selected()

    def _refresh_fasta_list(self):
        files = list_fasta_folder("fasta")
        basenames = [os.path.basename(f) for f in files]
        self.fasta_combo["values"] = basenames
        if basenames and not self.fasta_file_var.get():
            self.fasta_combo.current(0)
            self.fasta_file_var.set(basenames[0])

    def _browse_fasta_file(self):
        file_path = filedialog.askopenfilename(
            title="Open FASTA Sequence File",
            filetypes=[
                ("FASTA files", "*.fasta;*.fa;*.fna;*.txt"),
                ("All files", "*.*")
            ],
            initialdir=os.path.abspath("fasta") if os.path.exists("fasta") else "."
        )
        if file_path:
            self._load_fasta_from_path(file_path)

    def _on_fasta_file_selected(self, event=None):
        fname = self.fasta_combo.get() or self.fasta_file_var.get()
        if fname:
            self.fasta_file_var.set(fname)
            path = os.path.join("fasta", fname)
            if os.path.isfile(path):
                self._load_fasta_from_path(path)

    def _load_fasta_from_path(self, path: str):
        try:
            records = read_fasta_file(path)
            if not records:
                messagebox.showwarning("Empty File", f"No FASTA records found in:\n{path}")
                return
            self.current_fasta_records = records
            seq_ids = [r.id for r in records]
            self.seq_combo["values"] = seq_ids
            self.seq_combo.current(0)
            self.fasta_file_var.set(os.path.basename(path))
            self.run_fasta_detection()
        except Exception as e:
            messagebox.showerror("Error Opening FASTA", f"Could not read FASTA file:\n{e}")

    def run_fasta_detection(self):
        """Execute de novo detection on the selected FASTA sequence record."""
        if not self.current_fasta_records:
            return

        idx = self.seq_combo.current()
        if idx < 0 or idx >= len(self.current_fasta_records):
            idx = 0
        rec = self.current_fasta_records[idx]
        self.current_fasta_record = rec

        min_between = self.min_between_var.get()
        min_tir = self.min_tir_var.get()

        dets = self.detector.detect_new_candidates_from_sequence(
            rec.sequence,
            min_between=min_between,
            min_tir_length=min_tir
        )

        self.current_detections = {
            "homology": [],
            "structural": dets,
            "consensus": dets
        }
        self.current_simulation = None

        self._refresh_views()
        self.status_var.set(
            f"FASTA '{rec.id}' loaded ({rec.length} bp, GC {rec.gc_content:.1f}%). "
            f"Found {len(dets)} candidate transposon(s) with >= {min_between} bp in-between."
        )

    def randomize_seed(self):
        new_seed = random.randint(1000, 999999)
        self.simulator = TransposonSimulator(seed=new_seed)
        self.run_simulation_and_detection()

    def run_simulation_and_detection(self):
        target_len = self.len_var.get()
        te_count = self.te_count_var.get()
        scenario_choice = self.scenario_var.get()

        if "Scenario A" in scenario_choice:
            sim_res = self.simulator.simulate_scenario_non_intersecting(target_len, te_count)
        elif "Scenario B" in scenario_choice:
            sim_res = self.simulator.simulate_scenario_2_intersecting_nested(target_len, te_count)
        elif "Scenario C" in scenario_choice:
            sim_res = self.simulator.simulate_scenario_2_intersecting_overlapping(target_len, te_count)
        else:
            sim_res = self.simulator.simulate_scenario_3_intersecting_nested(target_len, te_count)

        self.current_simulation = sim_res
        self.current_fasta_record = None

        # Run detections
        self.current_detections = self.detector.detect_all(
            sim_res.genome_sequence,
            sim_res.transposon_library,
            min_between=10
        )

        active_dets = self._get_active_detections()
        self.current_benchmark = IntersectionAnalyzer.evaluate(sim_res, active_dets)

        self._refresh_views()
        self.status_var.set(
            f"Simulation & Detection completed. Genome length: {sim_res.final_length} bp | "
            f"Detected: {len(active_dets)} elements/fragments | Scenario: {sim_res.scenario_name}"
        )

    def _get_active_detections(self) -> List[DetectionResult]:
        if self.app_mode == "FASTA":
            return self.current_detections.get("structural", [])

        mode = self.det_mode_var.get()
        if "Homology" in mode:
            return self.current_detections.get("homology", [])
        elif "Structural" in mode:
            return self.current_detections.get("structural", [])
        else:
            return self.current_detections.get("consensus", [])

    def _get_current_sequence_and_length(self) -> Tuple[str, int, str]:
        if self.app_mode == "FASTA" and self.current_fasta_record:
            return self.current_fasta_record.sequence, self.current_fasta_record.length, self.current_fasta_record.id
        elif self.current_simulation:
            return self.current_simulation.genome_sequence, self.current_simulation.final_length, self.current_simulation.scenario_name
        return "", 0, ""

    def _refresh_views(self):
        seq, total_bp, name = self._get_current_sequence_and_length()
        if not seq:
            return

        active_dets = self._get_active_detections()
        if self.current_simulation:
            self.current_benchmark = IntersectionAnalyzer.evaluate(self.current_simulation, active_dets)

        self._draw_tracks()
        self._populate_table(active_dets)
        self._populate_behavior_analysis()
        self._populate_sequence_view()
        self._populate_benchmark_view()

    def _draw_tracks(self):
        """Render the genomic map with Ground Truth / Sequence and Detection tracks."""
        self.canvas.delete("all")
        seq, total_bp, seq_name = self._get_current_sequence_and_length()
        if not seq:
            return

        w = self.canvas.winfo_width()
        if w < 100:
            w = 950

        margin_x = 90
        track_w = w - margin_x - 30

        def bp_to_x(bp: int) -> float:
            return margin_x + (bp / max(1, total_bp)) * track_w

        # Ruler
        ruler_y = 22
        self.canvas.create_line(margin_x, ruler_y, margin_x + track_w, ruler_y, fill="#7F8C8D", width=1)
        self.canvas.create_text(margin_x - 10, ruler_y, text="0 bp", anchor=tk.E, font=("Segoe UI", 8), fill="#555")
        self.canvas.create_text(margin_x + track_w + 5, ruler_y, text=f"{total_bp} bp", anchor=tk.W, font=("Segoe UI", 8), fill="#555")

        step = 50 if total_bp <= 350 else 100
        for bp in range(step, total_bp, step):
            x = bp_to_x(bp)
            self.canvas.create_line(x, ruler_y - 4, x, ruler_y + 4, fill="#95A5A6")
            self.canvas.create_text(x, ruler_y - 9, text=str(bp), font=("Segoe UI", 7), fill="#7F8C8D")

        track_h = 18

        if self.app_mode == "SIMULATION" and self.current_simulation:
            # Track 1: Ground Truth
            gt_y = 52
            self.canvas.create_text(margin_x - 10, gt_y + track_h // 2, text="Ground Truth", anchor=tk.E, font=("Segoe UI", 8, "bold"), fill="#2C3E50")
            self.canvas.create_rectangle(margin_x, gt_y, margin_x + track_w, gt_y + track_h, fill="#EAEDED", outline="#BDC3C7")

            for gt in self.current_simulation.ground_truth_intervals:
                x1 = bp_to_x(gt.start)
                x2 = bp_to_x(gt.end + 1)
                color = get_color_for_element(gt.element_name)
                self.canvas.create_rectangle(x1, gt_y, x2, gt_y + track_h, fill=color, outline="#2C3E50")

                label = gt.element_name.replace("TE_", "")
                if gt.fragment_type != "full":
                    label += f" ({'5P' if '5' in gt.fragment_type else '3P'})"
                if x2 - x1 > 28:
                    self.canvas.create_text((x1 + x2) / 2, gt_y + track_h // 2, text=label, font=("Segoe UI", 7, "bold"), fill="#FFFFFF")
        else:
            # Track 1: FASTA Sequence Backbone
            gt_y = 52
            rec_id = self.current_fasta_record.id if self.current_fasta_record else "FASTA"
            self.canvas.create_text(margin_x - 10, gt_y + track_h // 2, text=rec_id[:12], anchor=tk.E, font=("Segoe UI", 8, "bold"), fill="#2C3E50")
            self.canvas.create_rectangle(margin_x, gt_y, margin_x + track_w, gt_y + track_h, fill="#D5DBDB", outline="#95A5A6")
            self.canvas.create_text(margin_x + track_w // 2, gt_y + track_h // 2, text=f"FASTA Sequence Backbone: {total_bp} bp", font=("Segoe UI", 7), fill="#34495E")

        # Track 2: Detections
        det_y = 90
        self.canvas.create_text(margin_x - 10, det_y + track_h // 2, text="Detected TEs", anchor=tk.E, font=("Segoe UI", 8, "bold"), fill="#16A085")
        self.canvas.create_rectangle(margin_x, det_y, margin_x + track_w, det_y + track_h, fill="#F8F9F9", outline="#BDC3C7")

        active_dets = self._get_active_detections()
        for det in active_dets:
            x1 = bp_to_x(det.start)
            x2 = bp_to_x(det.end + 1)
            color = get_color_for_element(det.element_name)

            is_chained = (det.status == "CHAINED_INTERRUPTED")
            outline_style = "#E74C3C" if is_chained else "#1B4F72"
            self.canvas.create_rectangle(
                x1, det_y, x2, det_y + track_h,
                fill=color if not is_chained else "",
                outline=outline_style,
                width=2 if is_chained else 1
            )

            short_name = det.element_name.replace("TE_", "").replace("New_Putative_TE_", "TE_").replace("Putative_TE_TIR_", "TIR_")
            if x2 - x1 > 28 and not is_chained:
                self.canvas.create_text((x1 + x2) / 2, det_y + track_h // 2, text=short_name, font=("Segoe UI", 7, "bold"), fill="#FFFFFF")
            elif is_chained:
                self.canvas.create_text((x1 + x2) / 2, det_y - 6, text=f"[Chained: {det.element_name.split()[0]}]", font=("Segoe UI", 7, "bold"), fill="#C0392B")

    def _populate_table(self, detections: List[DetectionResult]):
        for row in self.tree.get_children():
            self.tree.delete(row)

        for d in detections:
            in_between_str = f"{d.in_between_len} bp" if d.in_between_len > 0 else "-"
            self.tree.insert("", tk.END, values=(
                d.element_name,
                d.status,
                d.start,
                d.end,
                d.start_1based,
                d.end_1based,
                d.length,
                in_between_str,
                d.left_tir if d.left_tir else "-",
                d.right_tir if d.right_tir else "-",
                d.tsd if d.tsd else "-",
                d.method,
                f"{d.score:.1f}",
                d.details
            ))

    def _populate_behavior_analysis(self):
        self.behavior_text.delete("1.0", tk.END)

        if self.app_mode == "FASTA" and self.current_fasta_record:
            active_dets = self._get_active_detections()
            report = IntersectionAnalyzer.generate_fasta_report(
                self.current_fasta_record.id,
                self.current_fasta_record.length,
                active_dets,
                min_between=self.min_between_var.get()
            )
            for line in report:
                self.behavior_text.insert(tk.END, line + "\n")
        elif self.current_benchmark:
            for line in self.current_benchmark.intersection_report:
                self.behavior_text.insert(tk.END, line + "\n")

    def _populate_sequence_view(self):
        self.seq_text.delete("1.0", tk.END)
        seq, total_len, name = self._get_current_sequence_and_length()
        if not seq:
            return

        self.seq_text.insert(tk.END, f"; Sequence: {name} ({total_len} bp)\n")
        self.seq_text.insert(tk.END, "; Discovered transposable element positions are indexed below\n\n")

        line_len = 60
        for i in range(0, total_len, line_len):
            line_chunk = seq[i : i + line_len]
            blocks = " ".join(line_chunk[j : j + 10] for j in range(0, len(line_chunk), 10))
            self.seq_text.insert(tk.END, f"{i+1:5d}  {blocks}\n")

    def _populate_benchmark_view(self):
        self.bench_text.delete("1.0", tk.END)

        if self.app_mode == "FASTA":
            self.bench_text.insert(tk.END, "======================================================================\n")
            self.bench_text.insert(tk.END, "   DE NOVO FASTA CANDIDATE SUMMARY & SPACER VALIDATION\n")
            self.bench_text.insert(tk.END, "======================================================================\n\n")
            active_dets = self._get_active_detections()
            min_between = self.min_between_var.get()
            self.bench_text.insert(tk.END, f"File Record:              {self.current_fasta_record.id if self.current_fasta_record else 'N/A'}\n")
            self.bench_text.insert(tk.END, f"Total Sequence Length:    {self.current_fasta_record.length if self.current_fasta_record else 0} bp\n")
            self.bench_text.insert(tk.END, f"Spacer Constraint:        >= {min_between} bp in-between\n")
            self.bench_text.insert(tk.END, f"Total Candidates Found:   {len(active_dets)}\n\n")

            valid_count = sum(1 for d in active_dets if d.in_between_len >= min_between)
            self.bench_text.insert(tk.END, f"Validation Check: {valid_count} / {len(active_dets)} candidate elements meet the >= {min_between} bp in-between constraint (100%).\n\n")

            for i, d in enumerate(active_dets, 1):
                self.bench_text.insert(
                    tk.END,
                    f"Candidate {i}: {d.element_name} [{d.start}:{d.end}] (1-idx: [{d.start_1based}:{d.end_1based}])\n"
                    f"   Total Length:     {d.length} bp\n"
                    f"   In-Between Bases: {d.in_between_len} bp  (Constraint: >= {min_between} bp)\n"
                    f"   Left TIR:         {d.left_tir} ({len(d.left_tir)} bp)\n"
                    f"   Right TIR:        {d.right_tir} ({len(d.right_tir)} bp)\n"
                    f"   TSD:              {d.tsd if d.tsd else 'None'}\n"
                    f"   Status:           {d.status}\n\n"
                )
            return

        if not self.current_benchmark:
            return

        b = self.current_benchmark
        self.bench_text.insert(tk.END, "======================================================================\n")
        self.bench_text.insert(tk.END, "   ACCURACY BENCHMARK & COORDINATE OVERLAP EVALUATION\n")
        self.bench_text.insert(tk.END, "======================================================================\n\n")

        self.bench_text.insert(tk.END, f"Scenario:                 {self.current_simulation.scenario_name}\n")
        self.bench_text.insert(tk.END, f"Ground Truth Segments:    {b.total_ground_truth_intervals}\n")
        self.bench_text.insert(tk.END, f"Active Detections:        {b.total_detections}\n")
        self.bench_text.insert(tk.END, f"True Positives (TP):      {b.true_positives}\n")
        self.bench_text.insert(tk.END, f"False Positives (FP):     {b.false_positives}\n")
        self.bench_text.insert(tk.END, f"False Negatives (FN):     {b.false_negatives}\n")
        self.bench_text.insert(tk.END, f"Precision:                {b.precision:.2f}%\n")
        self.bench_text.insert(tk.END, f"Recall (Sensitivity):     {b.recall:.2f}%\n")
        self.bench_text.insert(tk.END, f"F1 Score:                 {b.f1_score:.2f}%\n")
        self.bench_text.insert(tk.END, f"Mean Coordinate IoU:      {b.mean_iou:.2f}%\n\n")

        self.bench_text.insert(tk.END, "Detailed Ground Truth Alignment Comparison:\n")
        self.bench_text.insert(tk.END, "----------------------------------------------------------------------\n")
        for m in b.matches:
            if m.is_matched:
                self.bench_text.insert(
                    tk.END,
                    f"[MATCH] {m.ground_truth_name} ({m.gt_type})\n"
                    f"   Ground Truth: [{m.gt_start}:{m.gt_end}] (len {m.gt_end - m.gt_start + 1})\n"
                    f"   Detected as:  {m.detected_name} [{m.det_start}:{m.det_end}]\n"
                    f"   IoU Overlap:  {m.iou*100:.1f}% via {m.det_method}\n"
                    f"   Assessment:   {m.status_assessment}\n\n"
                )
            else:
                self.bench_text.insert(
                    tk.END,
                    f"[MISSED] {m.ground_truth_name} ({m.gt_type})\n"
                    f"   Ground Truth: [{m.gt_start}:{m.gt_end}] (len {m.gt_end - m.gt_start + 1})\n"
                    f"   Assessment:   {m.status_assessment}\n\n"
                )

    def export_report(self):
        """Export comprehensive diagnostic report to a text file."""
        seq, total_len, title_name = self._get_current_sequence_and_length()
        if not seq:
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Export Transposon Analysis Report"
        )
        if not file_path:
            return

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write("Transposable Element Simulation & Detection Report\n")
                f.write(f"Source / Title: {title_name}\n")
                f.write(f"Sequence Length: {total_len} bp\n\n")

                f.write("=== DNA Sequence ===\n")
                f.write(seq + "\n\n")

                f.write("=== Detections Table ===\n")
                for d in self._get_active_detections():
                    f.write(
                        f"[{d.start}:{d.end}] (1-idx: [{d.start_1based}:{d.end_1based}]) "
                        f"{d.element_name} | {d.status} | len={d.length} | in_between={d.in_between_len}bp | "
                        f"TIR_L={d.left_tir} | TIR_R={d.right_tir} | TSD={d.tsd} | {d.method} | {d.details}\n"
                    )
                f.write("\n")

                f.write("=== Behavior / De Novo Report ===\n")
                if self.app_mode == "FASTA" and self.current_fasta_record:
                    report = IntersectionAnalyzer.generate_fasta_report(
                        self.current_fasta_record.id,
                        self.current_fasta_record.length,
                        self._get_active_detections(),
                        min_between=self.min_between_var.get()
                    )
                elif self.current_benchmark:
                    report = self.current_benchmark.intersection_report
                else:
                    report = []

                for line in report:
                    f.write(line + "\n")

            messagebox.showinfo("Export Successful", f"Report saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Could not save file: {e}")


def launch_gui():
    app = TransposonApp()
    app.mainloop()


if __name__ == "__main__":
    launch_gui()
