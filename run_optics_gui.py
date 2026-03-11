import tkinter as tk
import tkinter.font as tkFont
import sv_ttk 
from tkinter import ttk, filedialog, scrolledtext, messagebox
import threading
import sys
import os
import random
import queue # Import queue for thread-safe logging
import platform

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter('ignore')

# Attempt to import the run_optics_predictions function
try:
    from optics_predictions import run_optics_predictions
except ImportError:
    # Fallback/Placeholder if file not found during dev
    pass

# Attempt to import the generate_shap_explanation function
try:
    from optics_shap import generate_shap_explanation
except ImportError:
    print("Warning: Could not import 'generate_shap_explanation'. SHAP mode will fail if selected.")

# Attempt to import the run_structural_mapping function
try:
    from optics_structure_map import run_structural_mapping
except ImportError:
    print("Warning: Could not import 'run_structural_mapping'. Structure Mapping mode will fail if selected.")

# Attempt to import the run_structure_annotation function
try:
    from optics_structure_annotations import run_structure_annotation
except ImportError:
    print("Warning: Could not import 'run_structure_annotation'. Annotation mode will fail if selected.")


# --- Eye Animation Class ---
class HumanEyeLoadingScreen(tk.Toplevel):
    """
    A modal loading screen with an animated eye that looks around in multiple
    directions. The loading text has an animated ellipsis.
    """
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.withdraw() # Hide window until it's configured

        # --- Configuration ---
        self.title("Processing...")
        self.geometry("400x300")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()

        # Center the loading screen over the parent window
        # Update: Use parent's geometry logic more safely
        try:
            parent_x = parent.winfo_rootx()
            parent_y = parent.winfo_rooty()
            parent_width = parent.winfo_width()
            parent_height = parent.winfo_height()
            self.geometry(f"400x300+{parent_x + (parent_width // 2) - 200}+{parent_y + (parent_height // 2) - 150}")
        except:
            # Fallback if geometry calc fails
            self.geometry("400x300")

        # --- Widgets ---
        self.canvas = tk.Canvas(self, bg="#1c1c1c", highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.loading_label = ttk.Label(
            self.canvas,
            text="Processing",
            font=("Century Gothic", 14, "bold"),
            background="#1c1c1c",
            foreground="white" # Set a static white color
        )
        self.canvas.create_window(200, 265, window=self.loading_label)

        # --- Animation State ---
        self.animation_running = False
        self.animation_id = None
        self.ellipsis_count = 0

        # --- Eye Animation Specifics ---
        self.gaze_current_x = 0.0
        self.gaze_target_x = 0.0
        self.gaze_current_y = 0.0
        self.gaze_target_y = 0.0
        
        self.gaze_positions = [
            (-40, -25), (45, 20), (0, 0), (30, 25), (-20, -20), (0, 0), 
            (-45, 0), (45, 0), (0, 30), (0, -30), (0, 0)
        ]
        self.gaze_position_index = 0
        self.gaze_timer = 100 

        self.lift()
        self.deiconify()

    def start_animation(self):
        """Starts the eye and text animations."""
        self.animation_running = True
        self.animate_eye()
        self.animate_text()

    def stop_animation(self):
        """Stops the animations and destroys the window."""
        self.animation_running = False
        if self.animation_id:
            self.after_cancel(self.animation_id)
        self.destroy()

    def animate_text(self):
        """Animates an ellipsis on the loading text."""
        if not self.animation_running: return
        self.ellipsis_count = (self.ellipsis_count + 1) % 4
        dots = "." * self.ellipsis_count
        self.loading_label.config(text=f"Processing{dots}")
        self.after(400, self.animate_text) 

    def animate_eye(self):
        """The main animation loop for the eye's movement."""
        if not self.animation_running: return

        self.gaze_timer -= 1
        if self.gaze_timer <= 0:
            self.gaze_position_index = (self.gaze_position_index + 1) % len(self.gaze_positions)
            self.gaze_target_x, self.gaze_target_y = self.gaze_positions[self.gaze_position_index]
            
            if self.gaze_target_x == 0 and self.gaze_target_y == 0:
                self.gaze_timer = random.randint(100, 150) 
            else:
                self.gaze_timer = random.randint(40, 80) 

        self.gaze_current_x += (self.gaze_target_x - self.gaze_current_x) * 0.1
        self.gaze_current_y += (self.gaze_target_y - self.gaze_current_y) * 0.1

        self.canvas.delete("all")
        self.canvas.create_window(200, 265, window=self.loading_label)
        center_x, center_y = 200, 130
        
        self.canvas.create_oval(center_x - 110, center_y - 80, center_x + 110, center_y + 80, fill="#111111", outline="")
        self.canvas.create_oval(center_x - 100, center_y - 70, center_x + 100, center_y + 70, fill="#EAEAEA", outline="")
        
        iris_x = center_x + self.gaze_current_x
        iris_y = center_y + self.gaze_current_y
        iris_radius = 35
        pupil_radius = 15
        self.canvas.create_oval(iris_x - iris_radius, iris_y - iris_radius, iris_x + iris_radius, iris_y + iris_radius, fill="#5DADE2", outline="")
        self.canvas.create_oval(iris_x - pupil_radius, iris_y - pupil_radius, iris_x + pupil_radius, iris_y + pupil_radius, fill="black", outline="")

        highlight_x = iris_x + 10
        highlight_y = iris_y - 10
        self.canvas.create_oval(highlight_x - 5, highlight_y - 5, highlight_x + 3, highlight_y + 3, fill="white", outline="")

        self.animation_id = self.after(30, self.animate_eye)


# --- Text Redirector Class ---
class TextRedirector(object):
    """
    A thread-safe redirector that writes to a Queue instead of directly to a Widget.
    """
    def __init__(self, log_queue, tag="stdout"):
        self.log_queue = log_queue
        self.tag = tag

    def write(self, str_):
        self.log_queue.put((self.tag, str_))

    def flush(self):
        pass


# --- Mode Selector Frame ---
class ModeSelectorFrame(ttk.Frame):
    """
    Selection screen displayed inside the main application window.
    """
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        # Main Layout
        main_frame = ttk.Frame(self, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Use the logo image stored in the controller (App)
        if hasattr(self.controller, 'logo_icon_img') and self.controller.logo_icon_img:
            logo_lbl = ttk.Label(main_frame, image=self.controller.logo_icon_img)
            logo_lbl.pack(pady=(0, 10))

        title_lbl = ttk.Label(main_frame, text="Welcome to OPTICS", font=("Century Gothic", 18, "bold"))
        title_lbl.pack(pady=10)
        
        subtitle_lbl = ttk.Label(main_frame, text="Select your analysis pipeline:", font=("Century Gothic", 12))
        subtitle_lbl.pack(pady=(0, 20))

        # Buttons
        style = ttk.Style()
        style.configure("Big.TButton", font=("Century Gothic", 12))

        pred_btn = ttk.Button(main_frame, text="Standard Predictions\n(λmax & Spectral Tuning)", style="Big.TButton", 
                            command=lambda: self.controller.show_optics_gui('predictions'))
        pred_btn.pack(fill=tk.X, pady=10, ipady=10)

        shap_btn = ttk.Button(main_frame, text="SHAP Interpretation\n(Amino-Acid Importance)", style="Big.TButton", 
                            command=lambda: self.controller.show_optics_gui('shap'))
        shap_btn.pack(fill=tk.X, pady=10, ipady=10)

        struct_btn = ttk.Button(main_frame, text="Structure SHAP Mapping\n(3D Visualization of SHAP)", style="Big.TButton", 
                              command=lambda: self.controller.show_optics_gui('structure'))
        struct_btn.pack(fill=tk.X, pady=10, ipady=10)

        # New Button for Structure Annotations
        annot_btn = ttk.Button(main_frame, text="Structure Annotations\n(Custom 3D Visualization)", style="Big.TButton", 
                              command=lambda: self.controller.show_optics_gui('annotations'))
        annot_btn.pack(fill=tk.X, pady=10, ipady=10)


# --- Main Logic Frame ---
class OpticsGUIFrame(ttk.Frame):
    def __init__(self, parent, controller, mode='predictions'):
        super().__init__(parent)
        self.controller = controller
        self.mode = mode
        
        # Mode-specific settings
        if self.mode == 'predictions':
            self.title_suffix = "Predictions"
        elif self.mode == 'shap':
            self.title_suffix = "SHAP Analysis"
        elif self.mode == 'structure':
            self.title_suffix = "Structure Mapping"
        else: # annotations
            self.title_suffix = "Structure Annotations"
        
        # --- Choices ---
        self.version_choices = ['vpod_1.3']
        self.model_choices = ['whole-dataset', 'wildtype', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'type-one']
        #self.model_choices = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 
        #                      'wildtype-vert', 'type-one', 'whole-dataset-mnm', 
        #                      'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 
        #                      'wildtype-vert-mnm', 'wildtype-mut']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']
        self.viz_ftyp_choices = ['svg', 'png', 'pdf']
        self.software_choices = ['PyMOL', 'ChimeraX'] # New software choices
        
        # SHAP specific choices
        self.shap_mode_choices = ['both', 'comparison', 'single']

        # --- Top Bar ---
        top_bar_frame = ttk.Frame(self)
        top_bar_frame.pack(fill=tk.X, padx=10, pady=(10,0))

        header_text = f"OPTICS: {self.title_suffix}"
        
        # Use logo from controller
        if hasattr(self.controller, 'logo_icon_img') and self.controller.logo_icon_img:
            gui_logo_label = ttk.Label(top_bar_frame, image=self.controller.logo_icon_img, text=header_text, font=("Century Gothic", 16, "bold"), compound="left")
            gui_logo_label.pack(side=tk.LEFT, padx=(0,10))
        else:
            gui_logo_label = ttk.Label(top_bar_frame, text=header_text, font=("Century Gothic", 16, "bold"))
            gui_logo_label.pack(side=tk.LEFT, padx=(0,10))


        self.theme_toggle_button = ttk.Button(top_bar_frame, text="Toggle Theme", command=self.controller.toggle_theme)
        self.theme_toggle_button.pack(side=tk.RIGHT, padx=(0,5))
        
        self.back_button = ttk.Button(top_bar_frame, text="← Back", command=self.go_back)
        self.back_button.pack(side=tk.RIGHT, padx=(0,5))

        # --- Main Scrollable Area ---
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # --- INPUT WIDGETS ---
        current_row = 0
        scrollable_frame.columnconfigure(1, weight=1)

        # 1. Output Directory (Common to all modes)
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Output Directory:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_dir_var = tk.StringVar()
        ttk.Entry(scrollable_frame, textvariable=self.output_dir_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        ttk.Button(scrollable_frame, text="Browse...", command=self.browse_output_dir).grid(row=current_row, column=2, padx=5, pady=5)
        current_row += 1

        if self.mode in ['predictions', 'shap']:
            # -- ML Workflow Inputs --
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Input Sequence/FASTA File:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.input_file_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.input_file_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            ttk.Button(scrollable_frame, text="Browse...", command=self.browse_input_file).grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1

            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Output Filename Prefix:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.prediction_prefix_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.prediction_prefix_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Model Version:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.version_var = tk.StringVar(value=self.version_choices[0])
            ttk.Combobox(scrollable_frame, textvariable=self.version_var, values=self.version_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1

            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Model:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.model_var = tk.StringVar(value=self.model_choices[0])
            ttk.Combobox(scrollable_frame, textvariable=self.model_var, values=self.model_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1

            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Encoding Method:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.encoding_var = tk.StringVar(value=self.encoding_choices[1]) 
            ttk.Combobox(scrollable_frame, textvariable=self.encoding_var, values=self.encoding_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1

        elif self.mode == 'structure':
            # -- Structure Mapping Inputs --
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="SHAP Analysis CSV File (Single or Comparison):").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.shap_csv_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.shap_csv_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            ttk.Button(scrollable_frame, text="Browse...", command=self.browse_shap_csv).grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1

            # New: Pairwise sequence target mapping explicitly added
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Comparison Target Seq:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.comp_target_var = tk.StringVar(value="1")
            comp_target_frame = ttk.Frame(scrollable_frame)
            comp_target_frame.grid(row=current_row, column=1, columnspan=2, sticky=tk.W, padx=5, pady=5)
            ttk.Radiobutton(comp_target_frame, text="Sequence 1", variable=self.comp_target_var, value="1", command=self.toggle_pdb2_state).pack(side=tk.LEFT, padx=(0, 10))
            ttk.Radiobutton(comp_target_frame, text="Sequence 2", variable=self.comp_target_var, value="2", command=self.toggle_pdb2_state).pack(side=tk.LEFT, padx=(0, 10))
            ttk.Radiobutton(comp_target_frame, text="Both", variable=self.comp_target_var, value="both", command=self.toggle_pdb2_state).pack(side=tk.LEFT)
            current_row += 1
            ttk.Label(scrollable_frame, text="(Only applies if mapping a Pairwise Comparison CSV using Query Positions)", font=("Century Gothic", 9)).grid(row=current_row, column=1, columnspan=2, padx=5, pady=0, sticky=tk.W)
            current_row += 1

            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="PDB File/ID 1 (Seq 1 or Default):").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.pdb_input_var = tk.StringVar()
            self.pdb1_entry = ttk.Entry(scrollable_frame, textvariable=self.pdb_input_var, width=60)
            self.pdb1_entry.grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            self.pdb1_btn = ttk.Button(scrollable_frame, text="Browse File...", command=self.browse_pdb_file)
            self.pdb1_btn.grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="PDB File/ID 2 (Optional, for Seq 2):").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.pdb_input_2_var = tk.StringVar()
            self.pdb2_entry = ttk.Entry(scrollable_frame, textvariable=self.pdb_input_2_var, width=60)
            self.pdb2_entry.grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            self.pdb2_btn = ttk.Button(scrollable_frame, text="Browse File...", command=self.browse_pdb_file_2)
            self.pdb2_btn.grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1
            
            self.toggle_pdb2_state() # Initialize the state
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Chain ID:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.chain_var = tk.StringVar(value="A")
            ttk.Entry(scrollable_frame, textvariable=self.chain_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1
            
            self.use_query_pos_var = tk.BooleanVar(value=True)
            self.query_pos_check = ttk.Checkbutton(scrollable_frame, text="Use Query/Target Sequence Numbering", variable=self.use_query_pos_var)
            self.query_pos_check.grid(row=current_row, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)
            current_row += 1
            
            ttk.Label(scrollable_frame, text="(If unchecked, defaults to Reference/Bovine numbering from SHAP file)", font=("Century Gothic", 9)).grid(row=current_row, column=1, padx=5, pady=0, sticky=tk.W)
            current_row += 1

            # New: Map to Bovine Also Checkbox
            self.map_bovine_also_var = tk.BooleanVar(value=False)
            self.map_bovine_also_check = ttk.Checkbutton(scrollable_frame, text="Also map to Bovine Rhodopsin (1U19)", variable=self.map_bovine_also_var)
            self.map_bovine_also_check.grid(row=current_row, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)
            current_row += 1
            
            # New: Top N Labels and Software choice for mapping
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Top 'n' Sites to Label:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.top_n_labels_var = tk.IntVar(value=10)
            ttk.Spinbox(scrollable_frame, from_=0, to=100, textvariable=self.top_n_labels_var, width=58).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Visualization Software:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.struct_software_var = tk.StringVar(value=self.software_choices[1])
            ttk.Combobox(scrollable_frame, textvariable=self.struct_software_var, values=self.software_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1
        
        elif self.mode == 'annotations':
            # -- Structure Annotation Inputs --
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Annotation CSV File:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.annotation_csv_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.annotation_csv_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            ttk.Button(scrollable_frame, text="Browse...", command=self.browse_annotation_csv).grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1
            
            ttk.Label(scrollable_frame, text="(Columns required: 'position'. Optional: 'color', 'style', 'label')", font=("Century Gothic", 9)).grid(row=current_row, column=1, padx=5, pady=0, sticky=tk.W)
            current_row += 1

            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="PDB File or ID (e.g. 1U19):").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.pdb_input_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.pdb_input_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            ttk.Button(scrollable_frame, text="Browse File...", command=self.browse_pdb_file).grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Chain ID:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.chain_var = tk.StringVar(value="A")
            ttk.Entry(scrollable_frame, textvariable=self.chain_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1

            # New Software Selection Dropdown
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Visualization Software:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.software_var = tk.StringVar(value=self.software_choices[1])
            ttk.Combobox(scrollable_frame, textvariable=self.software_var, values=self.software_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1

        
        # --- MODE SPECIFIC OPTIONS ---
        
        if self.mode == 'predictions':
            # -- Prediction specific options --
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Tolerate Non-standard AAs:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.non_standard_aa_var = tk.BooleanVar(value=True) 
            self.non_standard_aa_check = ttk.Checkbutton(scrollable_frame, text="Enable to predict on sequences with non-standard amino-acids", variable=self.non_standard_aa_var)
            self.non_standard_aa_check.grid(row=current_row, column=1, columnspan=2, padx=5, pady=5, sticky=tk.EW)
            current_row += 1
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Tolerate Incomplete Seqs:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.incomp_seqs_var = tk.BooleanVar(value=False) 
            self.incomp_seqs_check = ttk.Checkbutton(scrollable_frame, text="Enable to predict on sequences outside range of 250-650 amino-acids", variable=self.incomp_seqs_var)
            self.incomp_seqs_check.grid(row=current_row, column=1, columnspan=2, padx=5, pady=5, sticky=tk.EW)
            current_row += 1

            # BLASTp Frame
            blastp_frame = ttk.LabelFrame(scrollable_frame, text="BLASTp Options", padding="10")
            blastp_frame.grid(row=current_row, column=0, columnspan=3, padx=5, pady=10, sticky=tk.EW)
            current_row += 1

            self.blastp_enabled_var = tk.BooleanVar(value=True) 
            self.blastp_check = ttk.Checkbutton(blastp_frame, text="Enable BLASTp Analysis", variable=self.blastp_enabled_var, command=self.toggle_blastp_options)
            self.blastp_check.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

            ttk.Label(blastp_frame, font=("Century Gothic", 12), text="BLASTp Report Filename:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
            self.blastp_report_var = tk.StringVar(value="blastp_report.txt")
            self.blastp_report_entry = ttk.Entry(blastp_frame, textvariable=self.blastp_report_var, width=50)
            self.blastp_report_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

            ttk.Label(blastp_frame, font=("Century Gothic", 12), text="Reference Sequence:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
            self.refseq_var = tk.StringVar(value=self.refseq_choices[0]) 
            self.refseq_combo = ttk.Combobox(blastp_frame, textvariable=self.refseq_var, values=self.refseq_choices, state="readonly", width=47)
            self.refseq_combo.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
            self.refseq_combo.bind("<<ComboboxSelected>>", self.toggle_custom_ref_file) 

            ttk.Label(blastp_frame, font=("Century Gothic", 12), text="Custom Reference File:").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
            self.custom_ref_file_var = tk.StringVar()
            self.custom_ref_file_entry = ttk.Entry(blastp_frame, textvariable=self.custom_ref_file_var, width=50)
            self.custom_ref_file_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
            self.custom_ref_file_button = ttk.Button(blastp_frame, text="Browse...", command=self.browse_custom_ref_file)
            self.custom_ref_file_button.grid(row=3, column=2, padx=5, pady=5)

            # Bootstrap Frame
            bootstrap_frame = ttk.LabelFrame(scrollable_frame, text="Bootstrap Options", padding="10")
            bootstrap_frame.grid(row=current_row, column=0, columnspan=3, padx=5, pady=10, sticky=tk.EW)
            current_row += 1

            self.bootstrap_enabled_var = tk.BooleanVar(value=True) 
            self.bootstrap_check = ttk.Checkbutton(bootstrap_frame, text="Enable Bootstrap Predictions", variable=self.bootstrap_enabled_var, command=self.toggle_bootstrap_options)
            self.bootstrap_check.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

            self.visualize_bootstrap_var = tk.BooleanVar(value=True) 
            self.visualize_check = ttk.Checkbutton(bootstrap_frame, text="Visualize Bootstrap Predictions", variable=self.visualize_bootstrap_var, command=self.toggle_bootstrap_options)
            self.visualize_check.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)
            
            ttk.Label(bootstrap_frame, font=("Century Gothic", 12), text="Bootstrap Viz Filename:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
            self.bootstrap_viz_file_var = tk.StringVar(value="bootstrap_viz")
            self.bootstrap_viz_entry = ttk.Entry(bootstrap_frame, textvariable=self.bootstrap_viz_file_var, width=50)
            self.bootstrap_viz_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
            
            ttk.Label(bootstrap_frame, font=("Century Gothic", 12), text="Viz Filetype:").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
            self.bootstrap_fytp_var = tk.StringVar(value=self.viz_ftyp_choices[0]) 
            self.fytp_combo = ttk.Combobox(bootstrap_frame, textvariable=self.bootstrap_fytp_var, values=self.viz_ftyp_choices, state="readonly", width=47)
            self.fytp_combo.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
            
            self.viz_xaxis_scale = tk.BooleanVar(value=False) 
            self.xaxis_check = ttk.Checkbutton(bootstrap_frame, text="Enable Full-Spectrum X-axis (300-650nm)", variable=self.viz_xaxis_scale)
            self.xaxis_check.grid(row=4, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

            # Init options
            self.toggle_blastp_options()
            self.toggle_bootstrap_options()
            self.toggle_custom_ref_file()

        elif self.mode == 'shap':
            # -- SHAP specific options --
            shap_frame = ttk.LabelFrame(scrollable_frame, text="SHAP Interpretation Options", padding="10")
            shap_frame.grid(row=current_row, column=0, columnspan=3, padx=5, pady=10, sticky=tk.EW)
            current_row += 1
            
            # Analysis Mode
            ttk.Label(shap_frame, font=("Century Gothic", 12), text="Analysis Mode:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
            self.shap_mode_var = tk.StringVar(value=self.shap_mode_choices[0])
            ttk.Combobox(shap_frame, textvariable=self.shap_mode_var, values=self.shap_mode_choices, state="readonly", width=47).grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)
            
            # N Positions
            ttk.Label(shap_frame, font=("Century Gothic", 12), text="Top N Features to Show:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
            self.n_positions_var = tk.IntVar(value=10)
            ttk.Spinbox(shap_frame, from_=1, to=100, textvariable=self.n_positions_var, width=48).grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

            # Use Reference Sites
            self.use_ref_sites_var = tk.BooleanVar(value=True)
            ttk.Checkbutton(shap_frame, text="Use Reference Numbering (e.g., Bovine Rhodopsin positions)", variable=self.use_ref_sites_var).grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)
            
            # Filetype
            ttk.Label(shap_frame, font=("Century Gothic", 12), text="Save Visualizations As:").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
            self.shap_save_as_var = tk.StringVar(value=self.viz_ftyp_choices[0])
            ttk.Combobox(shap_frame, textvariable=self.shap_save_as_var, values=self.viz_ftyp_choices, state="readonly", width=47).grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)


        # --- Run Button ---
        if self.mode == 'predictions':
             btn_text = "Run OPTICS Predictions"
        elif self.mode == 'shap':
             btn_text = "Run SHAP Analysis"
        elif self.mode == 'structure':
             btn_text = "Run SHAP Structure Mapping"
        else:
             btn_text = "Run Structure Annotation"

        self.run_button = ttk.Button(scrollable_frame, text=btn_text, command=self.start_run_thread)
        self.run_button.grid(row=current_row, column=0, columnspan=3, padx=5, pady=20)
        current_row += 1

        # --- Status/Output Area ---
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Output Log:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        current_row += 1
        self.output_text = scrolledtext.ScrolledText(scrollable_frame, wrap=tk.WORD, height=15, width=80, font=("Century Gothic", 10))
        self.output_text.grid(row=current_row, column=0, columnspan=3, padx=5, pady=5, sticky=tk.NSEW)
        self.output_text.configure(state='disabled') 

        # Setup Hyperlink Tags for Output Log
        self.output_text.tag_config("hyperlink", foreground="#3498db", underline=1)
        self.output_text.tag_bind("hyperlink", "<Button-1>", self.open_hyperlink)
        self.output_text.tag_bind("hyperlink", "<Enter>", lambda e: self.output_text.config(cursor="hand2"))
        self.output_text.tag_bind("hyperlink", "<Leave>", lambda e: self.output_text.config(cursor=""))
        
        scrollable_frame.columnconfigure(1, weight=1) 
        
        # Loading Screen placeholder
        self.loading_screen = None

    def toggle_pdb2_state(self, *args):
        if hasattr(self, 'comp_target_var') and hasattr(self, 'pdb2_entry'):
            target = self.comp_target_var.get()
            if target in ["2", "both"]:
                self.pdb2_entry.configure(state=tk.NORMAL)
                self.pdb2_btn.configure(state=tk.NORMAL)
            else:
                self.pdb2_entry.configure(state=tk.DISABLED)
                self.pdb2_btn.configure(state=tk.DISABLED)

    def open_hyperlink(self, event):
        """Opens the directory linked in the text widget."""
        try:
            index = self.output_text.index(f"@{event.x},{event.y}")
            # Get the text range of the tag at this index
            tags = self.output_text.tag_names(index)
            if "hyperlink" in tags:
                # Find the bounds of the current hyperlink
                #start = self.output_text.index(f"{index} wordstart")
                # Adjust start if needed to capture full path including slashes/colons
                # Simple approach: grab the whole line and extract the path we inserted
                #line_idx = index.split('.')[0]
                #line_text = self.output_text.get(f"{line_idx}.0", f"{line_idx}.end")
                
                # In write_to_log, we strip the prefix. Let's rely on the text content.
                # Since we insert just the path with the tag, we can try to extract it.
                # Better: The 'current' index points to the char.
                # Let's just grab the whole tagged range.
                ranges = self.output_text.tag_ranges("hyperlink")
                for start_idx, end_idx in zip(ranges[0::2], ranges[1::2]):
                    if self.output_text.compare(start_idx, "<=", index) and self.output_text.compare(index, "<", end_idx):
                        path_to_open = self.output_text.get(start_idx, end_idx).strip()
                        # Handle potential file:/// prefix if we added it, though os.startfile usually handles paths
                        if path_to_open.startswith("file:///"):
                             path_to_open = path_to_open[8:]
                        
                        if platform.system() == "Windows":
                            os.startfile(path_to_open)
                        elif platform.system() == "Darwin":
                            subprocess.call(["open", path_to_open])
                        else:
                            subprocess.call(["xdg-open", path_to_open])
                        return
        except Exception as e:
            messagebox.showerror("Error", f"Could not open path: {e}")

    def write_to_log(self, message, tag):
        """Called by the main controller to update the log widget."""
        try:
            self.output_text.configure(state='normal')
            
            # Check for special Link Token
            link_token = ">>>LINK<<<"
            if link_token in message:
                parts = message.split(link_token)
                # Insert pre-text
                if parts[0]:
                    self.output_text.insert(tk.END, parts[0], (tag,))
                # Insert Link
                link_path = parts[1].strip()
                self.output_text.insert(tk.END, link_path, ("hyperlink",))
                self.output_text.insert(tk.END, "\n") # Ensure newline after link
            else:
                self.output_text.insert(tk.END, message, (tag,))
                
            self.output_text.see(tk.END)
            self.output_text.configure(state='disabled')
        except Exception:
            pass # Widget might be dead, ignore

    def go_back(self):
        """Signals that the user wants to return to the mode selector."""
        self.controller.show_mode_selector()

    def browse_input_file(self):
        filename = filedialog.askopenfilename(title="Select Input FASTA File",
                                               filetypes=(("FASTA files", "*.fasta *.fa *.fna *.faa *.fas"),
                                                          ("Text files", "*.txt"), 
                                                          ("All files", "*.*")))
        if filename:
            self.input_file_var.set(filename)

    def browse_output_dir(self):
        dirname = filedialog.askdirectory(title="Select Output Directory")
        if dirname:
            self.output_dir_var.set(dirname)

    def browse_custom_ref_file(self):
        filename = filedialog.askopenfilename(title="Select Custom Reference File",
                                               filetypes=(("FASTA files", "*.fasta *.fa *.fna *.faa *.fas"),
                                                          ("All files", "*.*")))
        if filename:
            self.custom_ref_file_var.set(filename)

    def browse_shap_csv(self):
        filename = filedialog.askopenfilename(title="Select SHAP Analysis CSV",
                                               filetypes=(("CSV files", "*.csv"),
                                                          ("All files", "*.*")))
        if filename:
            self.shap_csv_var.set(filename)
    
    def browse_annotation_csv(self):
        filename = filedialog.askopenfilename(title="Select Annotation CSV/TSV",
                                               filetypes=(("CSV/TSV files", "*.csv *.tsv *.txt"),
                                                          ("All files", "*.*")))
        if filename:
            self.annotation_csv_var.set(filename)

    def browse_pdb_file(self):
        filename = filedialog.askopenfilename(title="Select PDB File",
                                               filetypes=(("PDB files", "*.pdb"),
                                                          ("All files", "*.*")))
        if filename:
            self.pdb_input_var.set(filename)
            
    def browse_pdb_file_2(self):
        filename = filedialog.askopenfilename(title="Select PDB File 2",
                                               filetypes=(("PDB files", "*.pdb"),
                                                          ("All files", "*.*")))
        if filename:
            self.pdb_input_2_var.set(filename)

    def toggle_blastp_options(self):
        state = tk.NORMAL if self.blastp_enabled_var.get() else tk.DISABLED
        self.blastp_report_entry.configure(state=state)
        self.refseq_combo.configure(state="readonly" if state == tk.NORMAL else tk.DISABLED) 
        self.toggle_custom_ref_file() 

    def toggle_custom_ref_file(self, event=None): 
        if not hasattr(self, 'blastp_enabled_var'): return
        
        blastp_on = self.blastp_enabled_var.get()
        custom_ref_selected = self.refseq_var.get() == "custom"
        
        if blastp_on and custom_ref_selected:
            self.custom_ref_file_entry.configure(state=tk.NORMAL)
            self.custom_ref_file_button.configure(state=tk.NORMAL)
        else:
            self.custom_ref_file_entry.configure(state=tk.DISABLED)
            self.custom_ref_file_button.configure(state=tk.DISABLED)
            if not custom_ref_selected: 
                 self.custom_ref_file_var.set("")

    def toggle_bootstrap_options(self):
        bootstrap_state = tk.NORMAL if self.bootstrap_enabled_var.get() else tk.DISABLED
        self.visualize_check.configure(state=bootstrap_state)
        
        visualize_state = tk.NORMAL if (self.bootstrap_enabled_var.get() and self.visualize_bootstrap_var.get()) else tk.DISABLED
        self.bootstrap_viz_entry.configure(state=visualize_state)
        self.xaxis_check.configure(state=visualize_state)
        self.fytp_combo.configure(state="readonly" if visualize_state == tk.NORMAL else tk.DISABLED)


    def log_message(self, message):
        """Thread-safe logging helper - puts message in queue via print/stdout redirect."""
        print(message) 

    def start_run_thread(self):
        # Validation based on mode
        if self.mode in ['predictions', 'shap']:
            if not self.input_file_var.get():
                messagebox.showerror("Input Error", "Please specify an input sequence/FASTA file.")
                return
        elif self.mode == 'structure':
            if not self.shap_csv_var.get():
                messagebox.showerror("Input Error", "Please specify the SHAP Analysis CSV file.")
                return
            if not self.pdb_input_var.get():
                messagebox.showerror("Input Error", "Please specify at least the primary PDB file path or ID.")
                return
        elif self.mode == 'annotations':
            if not self.annotation_csv_var.get():
                messagebox.showerror("Input Error", "Please specify the Annotation CSV file.")
                return
            if not self.pdb_input_var.get():
                messagebox.showerror("Input Error", "Please specify a PDB file path or ID.")
                return

        # Create default output dir if not specified
        if not self.output_dir_var.get():
            self.output_dir_var.set(os.path.join(os.getcwd(), 'prediction_outputs'))
            
        self.run_button.config(state=tk.DISABLED)
        self.output_text.configure(state='normal')
        self.output_text.delete(1.0, tk.END) 
        self.output_text.configure(state='disabled')
        
        if self.mode == 'predictions':
            log_txt = "Starting OPTICS Predictions..."
        elif self.mode == 'shap':
            log_txt = "Starting SHAP Analysis..."
        elif self.mode == 'structure':
            log_txt = "Starting Structure Mapping..."
        else:
            log_txt = "Starting Structure Annotation..."

        self.log_message(log_txt)

        self.loading_screen = HumanEyeLoadingScreen(self.controller)
        self.loading_screen.start_animation()

        thread = threading.Thread(target=self.run_predictions_logic, daemon=True)
        thread.start()

    def run_predictions_logic(self):
        try:
            # Common arguments
            pred_dir_val = self.output_dir_var.get()
            
            # Ensure path is absolute for linking
            pred_dir_val = os.path.abspath(pred_dir_val)

            if self.mode == 'predictions':
                input_val = self.input_file_var.get()
                output_val = self.prediction_prefix_var.get()
                if not output_val: output_val = "optics_results"
                version_val = self.version_var.get()
                model_val = self.model_var.get()
                encoding_val = self.encoding_var.get()
                tol_non_stan_aa = self.non_standard_aa_var.get()
                tol_incomp_seqs = self.incomp_seqs_var.get()

                blastp_val = self.blastp_enabled_var.get()
                iden_report_val = self.blastp_report_var.get() if blastp_val else None
                refseq_val = self.refseq_var.get() if blastp_val else "bovine" 
                reffile_val = self.custom_ref_file_var.get() if blastp_val and self.refseq_var.get() == "custom" else None
                
                bootstrap_val = self.bootstrap_enabled_var.get()
                visualize_bootstrap_val = self.visualize_bootstrap_var.get() if bootstrap_val else False
                bootstrap_viz_file_val = self.bootstrap_viz_file_var.get() if bootstrap_val and visualize_bootstrap_val else None
                viz_ftyp_val = self.bootstrap_fytp_var.get() if bootstrap_val and visualize_bootstrap_val else None
                xaxis_val = self.viz_xaxis_scale.get() if bootstrap_val and visualize_bootstrap_val else None

                pred_df, output_file_path = run_optics_predictions(
                    input_sequence=input_val,
                    pred_dir=pred_dir_val,
                    output=output_val,
                    model=model_val,
                    encoding_method=encoding_val,
                    blastp=blastp_val,
                    iden_report=iden_report_val,
                    refseq=refseq_val,
                    reffile=reffile_val,
                    bootstrap=bootstrap_val,
                    visualize_bootstrap=visualize_bootstrap_val,
                    bootstrap_viz_file=bootstrap_viz_file_val,
                    save_as=viz_ftyp_val,
                    full_spectrum_xaxis=xaxis_val,
                    model_version=version_val,
                    tolerate_non_standard_aa=tol_non_stan_aa,
                    tolerate_incomplete_seqs=tol_incomp_seqs
                )
                
                self.log_message(f"\n--- Predictions Complete ---")
                if output_file_path:
                    final_dir = os.path.dirname(os.path.abspath(output_file_path))
                    self.log_message(f"Results located at: >>>LINK<<<{final_dir}")
                    self.controller.after(0, lambda: messagebox.showinfo("Success", f"OPTICS predictions completed successfully!\nResults are in: {final_dir}"))
                
            elif self.mode == 'shap':
                input_val = self.input_file_var.get()
                output_val = self.prediction_prefix_var.get()
                if not output_val: output_val = "optics_results"
                version_val = self.version_var.get()
                model_val = self.model_var.get()
                encoding_val = self.encoding_var.get()
                shap_mode = self.shap_mode_var.get()
                n_pos = self.n_positions_var.get()
                use_ref = self.use_ref_sites_var.get()
                save_as = self.shap_save_as_var.get()
                
                generate_shap_explanation(
                    input_file=input_val, 
                    pred_dir=pred_dir_val, 
                    output=output_val, 
                    save_as=save_as, 
                    model=model_val, 
                    encoding_method=encoding_val, 
                    model_version=version_val, 
                    cmd_line="GUI_Execution", 
                    mode=shap_mode, 
                    n_positions=n_pos, 
                    use_reference_sites=use_ref
                )
                self.log_message(f"\n--- SHAP Analysis Complete ---")
                self.log_message(f"Results located at: >>>LINK<<<{pred_dir_val}")
                self.controller.after(0, lambda: messagebox.showinfo("Success", f"SHAP analysis completed successfully!\nResults are in: {pred_dir_val}"))

            elif self.mode == 'structure':
                csv_path = self.shap_csv_var.get()
                
                # Combine PDBs into a comma-separated string if both are provided
                pdb_input_1 = self.pdb_input_var.get()
                pdb_input_2 = self.pdb_input_2_var.get() if hasattr(self, 'pdb_input_2_var') else ""
                
                pdb_inputs = pdb_input_1
                if pdb_input_2:
                    pdb_inputs = f"{pdb_input_1},{pdb_input_2}"
                
                chain_val = self.chain_var.get()
                use_query_pos = self.use_query_pos_var.get()
                map_bovine_also = self.map_bovine_also_var.get()
                comp_target_val = self.comp_target_var.get()
                top_n_val = self.top_n_labels_var.get()
                software_val = self.struct_software_var.get().lower()

                pdb_out = run_structural_mapping(
                    shap_csv=csv_path,
                    pdb_input=pdb_inputs,
                    output_dir=pred_dir_val,
                    use_query_position=use_query_pos,
                    chain=chain_val,
                    map_to_bovine_also=map_bovine_also,
                    comp_target=comp_target_val,
                    top_n_labels=top_n_val,
                    software=software_val
                )
                
                self.log_message(f"\n--- Structure Mapping Complete ---")
                if pdb_out:
                    self.log_message(f"PDB Saved to: {pdb_out}")
                    self.log_message(f"Results located at: >>>LINK<<<{pred_dir_val}")
                    self.controller.after(0, lambda: messagebox.showinfo("Success", f"Mapping completed successfully!\nSaved PDB to: {pdb_out}"))
                else:
                    self.controller.after(0, lambda: messagebox.showerror("Error", "Structure mapping failed. Check log for details."))
            
            elif self.mode == 'annotations':
                csv_path = self.annotation_csv_var.get()
                pdb_input = self.pdb_input_var.get()
                chain_val = self.chain_var.get()
                software_val = self.software_var.get().lower()
                
                run_structure_annotation(
                    annotation_file=csv_path,
                    pdb_input=pdb_input,
                    output_dir=pred_dir_val,
                    chain=chain_val,
                    software=software_val
                )
                
                ext = ".cxc" if software_val == 'chimerax' else ".pml"
                self.log_message(f"\n--- Structure Annotation Complete ---")
                self.log_message(f"Results located at: >>>LINK<<<{pred_dir_val}")
                self.controller.after(0, lambda: messagebox.showinfo("Success", f"Annotation script created successfully!\nCheck output directory for {ext} file."))


        except Exception as e:
            self.log_message(f"\n--- ERROR ---")
            self.log_message(f"An error occurred: {str(e)}")
            import traceback
            self.log_message(f"Traceback:\n{traceback.format_exc()}")
            self.controller.after(0, lambda: messagebox.showerror("Error", f"An error occurred: {e}"))
        finally:
            self.controller.after(0, lambda: self.run_button.config(state=tk.NORMAL))
            if self.loading_screen:
                self.controller.after(0, self.loading_screen.stop_animation)


# --- Main Application Controller ---
class OpticsApp(tk.Tk):
    """
    Main Application Window (Single Root).
    Manages navigation between frames (ModeSelector -> OpticsGUI).
    Keeps the Tcl interpreter alive to prevent Threading/GC errors.
    """
    def __init__(self):
        super().__init__()
        self.title("OPTICS Pipeline")
        self.geometry("1000x950")
        
        # --- Theme ---
        sv_ttk.set_theme("dark")
        self.dark_mode_enabled = True
        
        # --- Define default font ---
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(family="Century Gothic", size=10)
        self.option_add("*Font", default_font)
        
        # --- Logo ---
        self.logo_icon_display_path = "./data/logo/optics_logo_resized.png"
        self.logo_icon_img = None
        try:
            self.logo_icon_img = tk.PhotoImage(file=self.logo_icon_display_path)
            self.iconphoto(True, self.logo_icon_img) 
        except tk.TclError:
            pass

        # --- Logging Setup ---
        self.log_queue = queue.Queue()
        # Redirect stdout/stderr to queue
        sys.stdout = TextRedirector(self.log_queue, "stdout")
        sys.stderr = TextRedirector(self.log_queue, "stderr")

        # --- Frame Container ---
        self.container = ttk.Frame(self)
        self.container.pack(fill="both", expand=True)
        
        self.current_frame = None
        
        # Start with Mode Selector
        self.show_mode_selector()
        
        # Start Log Polling
        self.update_log_from_queue()

    def show_mode_selector(self):
        """Switches to the Mode Selector screen."""
        if self.current_frame:
            self.current_frame.destroy()
        self.current_frame = ModeSelectorFrame(self.container, self)
        self.current_frame.pack(fill="both", expand=True)

    def show_optics_gui(self, mode):
        """Switches to the Main Logic screen for the selected mode."""
        if self.current_frame:
            self.current_frame.destroy()
        self.current_frame = OpticsGUIFrame(self.container, self, mode)
        self.current_frame.pack(fill="both", expand=True)

    def toggle_theme(self):
        if self.dark_mode_enabled:
            sv_ttk.set_theme("light")
            self.dark_mode_enabled = False
        else:
            sv_ttk.set_theme("dark")
            self.dark_mode_enabled = True

    def update_log_from_queue(self):
        """Polls the queue and updates the current frame if it has a log widget."""
        while not self.log_queue.empty():
            try:
                tag, message = self.log_queue.get_nowait()
                # Check if current frame has the logging method
                if hasattr(self.current_frame, 'write_to_log'):
                    self.current_frame.write_to_log(message, tag)
            except queue.Empty:
                break
        
        self.after(100, self.update_log_from_queue)

    def destroy(self):
        """Restore stdout/stderr when app closes."""
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        super().destroy()


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    app = OpticsApp()
    app.mainloop()