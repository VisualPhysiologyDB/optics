import tkinter as tk
import tkinter.font as tkFont
import customtkinter as ctk
from tkinter import filedialog, messagebox
import threading
import sys
import os
import random
import queue
import platform
from PIL import Image

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter('ignore')

# Attempt to import the run_optics_predictions function
try:
    from optics_predictions import run_optics_predictions
except ImportError:
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

# Set CustomTkinter defaults
ctk.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
ctk.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"


# --- Eye Animation Class ---
class HumanEyeLoadingScreen(ctk.CTkToplevel):
    """
    A modal loading screen with an animated eye that looks around in multiple
    directions. The loading text has an animated ellipsis.
    """
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.withdraw()

        # --- Configuration ---
        self.title("Processing...")
        self.geometry("400x300")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()

        # Center the loading screen over the parent window
        try:
            parent_x = parent.winfo_rootx()
            parent_y = parent.winfo_rooty()
            parent_width = parent.winfo_width()
            parent_height = parent.winfo_height()
            self.geometry(f"400x300+{parent_x + (parent_width // 2) - 200}+{parent_y + (parent_height // 2) - 150}")
        except:
            self.geometry("400x300")

        # --- Widgets ---
        bg_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkFrame"]["fg_color"])
        
        self.canvas = tk.Canvas(self, bg=bg_color, highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.loading_label = ctk.CTkLabel(
            self,
            text="Processing",
            font=ctk.CTkFont(family="Century Gothic", size=16, weight="bold"),
            text_color="white"
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
        self.animation_running = True
        self.animate_eye()
        self.animate_text()

    def stop_animation(self):
        self.animation_running = False
        if self.animation_id:
            self.after_cancel(self.animation_id)
        self.destroy()

    def animate_text(self):
        if not self.animation_running: return
        self.ellipsis_count = (self.ellipsis_count + 1) % 4
        dots = "." * self.ellipsis_count
        self.loading_label.configure(text=f"Processing{dots}")
        self.after(400, self.animate_text) 

    def animate_eye(self):
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
    def __init__(self, log_queue, tag="stdout"):
        self.log_queue = log_queue
        self.tag = tag

    def write(self, str_):
        self.log_queue.put((self.tag, str_))

    def flush(self):
        pass


# --- Mode Selector Frame ---
class ModeSelectorFrame(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent, fg_color="transparent")
        self.controller = controller
        
        main_frame = ctk.CTkFrame(self)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=40, pady=40)

        # Logo
        if hasattr(self.controller, 'ctk_logo_img') and self.controller.ctk_logo_img:
            logo_lbl = ctk.CTkLabel(main_frame, image=self.controller.ctk_logo_img, text="")
            logo_lbl.pack(pady=(20, 10))

        title_lbl = ctk.CTkLabel(main_frame, text="Welcome to OPTICS", font=ctk.CTkFont(family="Century Gothic", size=24, weight="bold"))
        title_lbl.pack(pady=10)
        
        subtitle_lbl = ctk.CTkLabel(main_frame, text="Select your analysis pipeline:", font=ctk.CTkFont(family="Century Gothic", size=14))
        subtitle_lbl.pack(pady=(0, 20))

        # Buttons
        btn_font = ctk.CTkFont(family="Century Gothic", size=14)

        pred_btn = ctk.CTkButton(main_frame, text="Standard Predictions\n(λmax & Spectral Tuning)", font=btn_font, height=60,
                            command=lambda: self.controller.show_optics_gui('predictions'))
        pred_btn.pack(fill=tk.X, padx=50, pady=10)

        shap_btn = ctk.CTkButton(main_frame, text="SHAP Interpretation\n(Amino-Acid Importance)", font=btn_font, height=60,
                            command=lambda: self.controller.show_optics_gui('shap'))
        shap_btn.pack(fill=tk.X, padx=50, pady=10)

        struct_btn = ctk.CTkButton(main_frame, text="Structure SHAP Mapping\n(3D Visualization of SHAP)", font=btn_font, height=60,
                              command=lambda: self.controller.show_optics_gui('structure'))
        struct_btn.pack(fill=tk.X, padx=50, pady=10)

        annot_btn = ctk.CTkButton(main_frame, text="Structure Annotations\n(Custom 3D Visualization)", font=btn_font, height=60,
                              command=lambda: self.controller.show_optics_gui('annotations'))
        annot_btn.pack(fill=tk.X, padx=50, pady=10)


# --- Main Logic Frame ---
class OpticsGUIFrame(ctk.CTkFrame):
    def __init__(self, parent, controller, mode='predictions'):
        super().__init__(parent, fg_color="transparent")
        self.controller = controller
        self.mode = mode
        
        self.lbl_font = ctk.CTkFont(family="Century Gothic", size=13)
        self.title_font = ctk.CTkFont(family="Century Gothic", size=18, weight="bold")
        
        if self.mode == 'predictions':
            self.title_suffix = "Predictions"
        elif self.mode == 'shap':
            self.title_suffix = "SHAP Analysis"
        elif self.mode == 'structure':
            self.title_suffix = "Structure Mapping"
        else:
            self.title_suffix = "Structure Annotations"
        
        # --- Choices ---
        self.version_choices = ['vpod_1.3']
        self.model_choices = ['whole-dataset', 'wildtype', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'type-one']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']
        self.viz_ftyp_choices = ['svg', 'png', 'pdf']
        self.software_choices = ['PyMOL', 'ChimeraX'] 
        self.shap_mode_choices = ['both', 'comparison', 'single']

        # --- Top Bar ---
        top_bar_frame = ctk.CTkFrame(self, fg_color="transparent")
        top_bar_frame.pack(fill=tk.X, padx=20, pady=(10, 0))

        header_text = f"OPTICS: {self.title_suffix}"
        
        if hasattr(self.controller, 'ctk_logo_small') and self.controller.ctk_logo_small:
            gui_logo_label = ctk.CTkLabel(top_bar_frame, image=self.controller.ctk_logo_small, text=f"  {header_text}", font=self.title_font, compound="left")
            gui_logo_label.pack(side=tk.LEFT, padx=0)
        else:
            gui_logo_label = ctk.CTkLabel(top_bar_frame, text=header_text, font=self.title_font)
            gui_logo_label.pack(side=tk.LEFT, padx=0)

        self.back_button = ctk.CTkButton(top_bar_frame, text="← Back", width=80, command=self.go_back)
        self.back_button.pack(side=tk.RIGHT, padx=5)

        self.theme_toggle_button = ctk.CTkButton(top_bar_frame, text="Toggle Theme", width=120, command=self.controller.toggle_theme)
        self.theme_toggle_button.pack(side=tk.RIGHT, padx=5)

        # --- Main Scrollable Area ---
        self.scrollable_frame = ctk.CTkScrollableFrame(self)
        self.scrollable_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        self.scrollable_frame.columnconfigure(1, weight=1)

        # --- INPUT WIDGETS ---
        current_row = 0

        # 1. Output Directory (Common)
        ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Output Directory:").grid(row=current_row, column=0, padx=10, pady=(10,5), sticky=tk.W)
        self.output_dir_var = tk.StringVar()
        ctk.CTkEntry(self.scrollable_frame, textvariable=self.output_dir_var).grid(row=current_row, column=1, padx=10, pady=(10,5), sticky=tk.EW)
        ctk.CTkButton(self.scrollable_frame, text="Browse...", width=80, command=self.browse_output_dir).grid(row=current_row, column=2, padx=10, pady=(10,5))
        current_row += 1

        if self.mode in ['predictions', 'shap']:
            # -- ML Workflow Inputs --
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Input Sequence/FASTA File:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.input_file_var = tk.StringVar()
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.input_file_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            ctk.CTkButton(self.scrollable_frame, text="Browse...", width=80, command=self.browse_input_file).grid(row=current_row, column=2, padx=10, pady=5)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Output Filename Prefix:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.prediction_prefix_var = tk.StringVar(value="unnamed")
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.prediction_prefix_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Model Version:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.version_var = ctk.StringVar(value=self.version_choices[0])
            ctk.CTkOptionMenu(self.scrollable_frame, variable=self.version_var, values=self.version_choices).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Model:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.model_var = ctk.StringVar(value=self.model_choices[0])
            ctk.CTkOptionMenu(self.scrollable_frame, variable=self.model_var, values=self.model_choices).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Encoding Method:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.encoding_var = ctk.StringVar(value=self.encoding_choices[1]) 
            ctk.CTkOptionMenu(self.scrollable_frame, variable=self.encoding_var, values=self.encoding_choices).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1

        elif self.mode == 'structure':
            # -- Structure Mapping Inputs --
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="SHAP Analysis CSV File:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.shap_csv_var = tk.StringVar()
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.shap_csv_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            ctk.CTkButton(self.scrollable_frame, text="Browse...", width=80, command=self.browse_shap_csv).grid(row=current_row, column=2, padx=10, pady=5)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Comparison Target Seq:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.comp_target_var = tk.StringVar(value="1")
            
            comp_target_frame = ctk.CTkFrame(self.scrollable_frame, fg_color="transparent")
            comp_target_frame.grid(row=current_row, column=1, columnspan=2, sticky=tk.W, padx=10, pady=5)
            ctk.CTkRadioButton(comp_target_frame, text="Sequence 1", variable=self.comp_target_var, value="1", command=self.toggle_pdb2_state).pack(side=tk.LEFT, padx=(0, 15))
            ctk.CTkRadioButton(comp_target_frame, text="Sequence 2", variable=self.comp_target_var, value="2", command=self.toggle_pdb2_state).pack(side=tk.LEFT, padx=(0, 15))
            ctk.CTkRadioButton(comp_target_frame, text="Both", variable=self.comp_target_var, value="both", command=self.toggle_pdb2_state).pack(side=tk.LEFT)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, text="(Only applies if mapping a Pairwise Comparison CSV using Query Positions)", font=ctk.CTkFont(size=11, slant="italic")).grid(row=current_row, column=1, columnspan=2, padx=10, pady=(0, 10), sticky=tk.W)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="PDB File/ID 1 (Seq 1/Def.):").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.pdb_input_var = tk.StringVar()
            self.pdb1_entry = ctk.CTkEntry(self.scrollable_frame, textvariable=self.pdb_input_var)
            self.pdb1_entry.grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            self.pdb1_btn = ctk.CTkButton(self.scrollable_frame, text="Browse File...", width=80, command=self.browse_pdb_file)
            self.pdb1_btn.grid(row=current_row, column=2, padx=10, pady=5)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="PDB File/ID 2 (Seq 2 Optional):").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.pdb_input_2_var = tk.StringVar()
            self.pdb2_entry = ctk.CTkEntry(self.scrollable_frame, textvariable=self.pdb_input_2_var)
            self.pdb2_entry.grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            self.pdb2_btn = ctk.CTkButton(self.scrollable_frame, text="Browse File...", width=80, command=self.browse_pdb_file_2)
            self.pdb2_btn.grid(row=current_row, column=2, padx=10, pady=5)
            current_row += 1
            
            self.toggle_pdb2_state()
            
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Chain ID:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.chain_var = tk.StringVar(value="A")
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.chain_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1
            
            self.use_query_pos_var = tk.BooleanVar(value=True)
            self.query_pos_check = ctk.CTkCheckBox(self.scrollable_frame, text="Use Query/Target Sequence Numbering", variable=self.use_query_pos_var)
            self.query_pos_check.grid(row=current_row, column=1, columnspan=2, padx=10, pady=5, sticky=tk.W)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, text="(If unchecked, defaults to Reference/Bovine numbering from SHAP file)", font=ctk.CTkFont(size=11, slant="italic")).grid(row=current_row, column=1, padx=10, pady=(0, 10), sticky=tk.W)
            current_row += 1

            self.map_bovine_also_var = tk.BooleanVar(value=False)
            self.map_bovine_also_check = ctk.CTkCheckBox(self.scrollable_frame, text="Also map to Bovine Rhodopsin (1U19)", variable=self.map_bovine_also_var)
            self.map_bovine_also_check.grid(row=current_row, column=1, columnspan=2, padx=10, pady=5, sticky=tk.W)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Top 'n' Sites to Label:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.top_n_labels_var = tk.StringVar(value="10")
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.top_n_labels_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Visualization Software:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.struct_software_var = ctk.StringVar(value=self.software_choices[1])
            ctk.CTkOptionMenu(self.scrollable_frame, variable=self.struct_software_var, values=self.software_choices).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1
        
        elif self.mode == 'annotations':
            # -- Structure Annotation Inputs --
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Annotation CSV File:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.annotation_csv_var = tk.StringVar()
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.annotation_csv_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            ctk.CTkButton(self.scrollable_frame, text="Browse...", width=80, command=self.browse_annotation_csv).grid(row=current_row, column=2, padx=10, pady=5)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, text="(Columns required: 'position'. Optional: 'color', 'style', 'label')", font=ctk.CTkFont(size=11, slant="italic")).grid(row=current_row, column=1, padx=10, pady=(0, 10), sticky=tk.W)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="PDB File or ID:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.pdb_input_var = tk.StringVar()
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.pdb_input_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            ctk.CTkButton(self.scrollable_frame, text="Browse File...", width=80, command=self.browse_pdb_file).grid(row=current_row, column=2, padx=10, pady=5)
            current_row += 1
            
            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Chain ID:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.chain_var = tk.StringVar(value="A")
            ctk.CTkEntry(self.scrollable_frame, textvariable=self.chain_var).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1

            ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Visualization Software:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
            self.software_var = ctk.StringVar(value=self.software_choices[1])
            ctk.CTkOptionMenu(self.scrollable_frame, variable=self.software_var, values=self.software_choices).grid(row=current_row, column=1, padx=10, pady=5, sticky=tk.EW)
            current_row += 1

        
        # --- MODE SPECIFIC OPTIONS ---
        if self.mode == 'predictions':
            # -- Prediction specific options --
            self.non_standard_aa_var = tk.BooleanVar(value=True) 
            ctk.CTkCheckBox(self.scrollable_frame, text="Tolerate Non-standard AAs", variable=self.non_standard_aa_var).grid(row=current_row, column=0, columnspan=2, padx=10, pady=10, sticky=tk.W)
            current_row += 1
            
            self.incomp_seqs_var = tk.BooleanVar(value=False) 
            ctk.CTkCheckBox(self.scrollable_frame, text="Tolerate Incomplete Seqs (outside 250-650 AA range)", variable=self.incomp_seqs_var).grid(row=current_row, column=0, columnspan=2, padx=10, pady=(0, 15), sticky=tk.W)
            current_row += 1

            # BLASTp Frame
            blastp_frame = ctk.CTkFrame(self.scrollable_frame)
            blastp_frame.grid(row=current_row, column=0, columnspan=3, padx=10, pady=10, sticky=tk.EW)
            blastp_frame.columnconfigure(1, weight=1)
            current_row += 1

            self.blastp_enabled_var = tk.BooleanVar(value=True) 
            ctk.CTkCheckBox(blastp_frame, text="Enable BLASTp Analysis", variable=self.blastp_enabled_var, font=self.title_font, command=self.toggle_blastp_options).grid(row=0, column=0, columnspan=2, padx=15, pady=15, sticky=tk.W)

            ctk.CTkLabel(blastp_frame, font=self.lbl_font, text="BLASTp Report Filename:").grid(row=1, column=0, padx=15, pady=5, sticky=tk.W)
            self.blastp_report_var = tk.StringVar(value="blastp_report.txt")
            self.blastp_report_entry = ctk.CTkEntry(blastp_frame, textvariable=self.blastp_report_var)
            self.blastp_report_entry.grid(row=1, column=1, padx=15, pady=5, sticky=tk.EW)

            ctk.CTkLabel(blastp_frame, font=self.lbl_font, text="Reference Sequence:").grid(row=2, column=0, padx=15, pady=5, sticky=tk.W)
            self.refseq_var = ctk.StringVar(value=self.refseq_choices[0]) 
            self.refseq_combo = ctk.CTkOptionMenu(blastp_frame, variable=self.refseq_var, values=self.refseq_choices, command=self.toggle_custom_ref_file)
            self.refseq_combo.grid(row=2, column=1, padx=15, pady=5, sticky=tk.EW)

            ctk.CTkLabel(blastp_frame, font=self.lbl_font, text="Custom Reference File:").grid(row=3, column=0, padx=15, pady=(5, 15), sticky=tk.W)
            self.custom_ref_file_var = tk.StringVar()
            self.custom_ref_file_entry = ctk.CTkEntry(blastp_frame, textvariable=self.custom_ref_file_var)
            self.custom_ref_file_entry.grid(row=3, column=1, padx=15, pady=(5, 15), sticky=tk.EW)
            self.custom_ref_file_button = ctk.CTkButton(blastp_frame, text="Browse...", width=80, command=self.browse_custom_ref_file)
            self.custom_ref_file_button.grid(row=3, column=2, padx=15, pady=(5, 15))

            # Bootstrap Frame
            bootstrap_frame = ctk.CTkFrame(self.scrollable_frame)
            bootstrap_frame.grid(row=current_row, column=0, columnspan=3, padx=10, pady=10, sticky=tk.EW)
            bootstrap_frame.columnconfigure(1, weight=1)
            current_row += 1

            self.bootstrap_enabled_var = tk.BooleanVar(value=True) 
            ctk.CTkCheckBox(bootstrap_frame, text="Enable Bootstrap Predictions", variable=self.bootstrap_enabled_var, font=self.title_font, command=self.toggle_bootstrap_options).grid(row=0, column=0, columnspan=2, padx=15, pady=15, sticky=tk.W)

            self.visualize_bootstrap_var = tk.BooleanVar(value=True) 
            self.visualize_check = ctk.CTkCheckBox(bootstrap_frame, text="Visualize Bootstrap Predictions", variable=self.visualize_bootstrap_var, command=self.toggle_bootstrap_options)
            self.visualize_check.grid(row=1, column=0, columnspan=2, padx=15, pady=5, sticky=tk.W)
            
            ctk.CTkLabel(bootstrap_frame, font=self.lbl_font, text="Bootstrap Viz Filename:").grid(row=2, column=0, padx=15, pady=5, sticky=tk.W)
            self.bootstrap_viz_file_var = tk.StringVar(value="bootstrap_viz")
            self.bootstrap_viz_entry = ctk.CTkEntry(bootstrap_frame, textvariable=self.bootstrap_viz_file_var)
            self.bootstrap_viz_entry.grid(row=2, column=1, padx=15, pady=5, sticky=tk.EW)
            
            ctk.CTkLabel(bootstrap_frame, font=self.lbl_font, text="Viz Filetype:").grid(row=3, column=0, padx=15, pady=5, sticky=tk.W)
            self.bootstrap_fytp_var = ctk.StringVar(value=self.viz_ftyp_choices[0]) 
            self.fytp_combo = ctk.CTkOptionMenu(bootstrap_frame, variable=self.bootstrap_fytp_var, values=self.viz_ftyp_choices)
            self.fytp_combo.grid(row=3, column=1, padx=15, pady=5, sticky=tk.EW)
            
            self.viz_xaxis_scale = tk.BooleanVar(value=False) 
            self.xaxis_check = ctk.CTkCheckBox(bootstrap_frame, text="Enable Full-Spectrum X-axis (300-650nm)", variable=self.viz_xaxis_scale)
            self.xaxis_check.grid(row=4, column=0, columnspan=2, padx=15, pady=(5, 15), sticky=tk.W)

            # Init options
            self.toggle_blastp_options()
            self.toggle_bootstrap_options()
            self.toggle_custom_ref_file()

        elif self.mode == 'shap':
            shap_frame = ctk.CTkFrame(self.scrollable_frame)
            shap_frame.grid(row=current_row, column=0, columnspan=3, padx=10, pady=10, sticky=tk.EW)
            shap_frame.columnconfigure(1, weight=1)
            current_row += 1
            
            ctk.CTkLabel(shap_frame, text="SHAP Options", font=self.title_font).grid(row=0, column=0, columnspan=2, padx=15, pady=15, sticky=tk.W)

            ctk.CTkLabel(shap_frame, font=self.lbl_font, text="Analysis Mode:").grid(row=1, column=0, padx=15, pady=5, sticky=tk.W)
            self.shap_mode_var = ctk.StringVar(value=self.shap_mode_choices[0])
            ctk.CTkOptionMenu(shap_frame, variable=self.shap_mode_var, values=self.shap_mode_choices).grid(row=1, column=1, padx=15, pady=5, sticky=tk.EW)
            
            ctk.CTkLabel(shap_frame, font=self.lbl_font, text="Top N Features to Show:").grid(row=2, column=0, padx=15, pady=5, sticky=tk.W)
            self.n_positions_var = tk.StringVar(value="10")
            ctk.CTkEntry(shap_frame, textvariable=self.n_positions_var).grid(row=2, column=1, padx=15, pady=5, sticky=tk.EW)

            self.use_ref_sites_var = tk.BooleanVar(value=True)
            ctk.CTkCheckBox(shap_frame, text="Use Reference Numbering (e.g., Bovine Rhodopsin)", variable=self.use_ref_sites_var).grid(row=3, column=0, columnspan=2, padx=15, pady=15, sticky=tk.W)
            
            ctk.CTkLabel(shap_frame, font=self.lbl_font, text="Save Visualizations As:").grid(row=4, column=0, padx=15, pady=(5,15), sticky=tk.W)
            self.shap_save_as_var = ctk.StringVar(value=self.viz_ftyp_choices[0])
            ctk.CTkOptionMenu(shap_frame, variable=self.shap_save_as_var, values=self.viz_ftyp_choices).grid(row=4, column=1, padx=15, pady=(5,15), sticky=tk.EW)


        # --- Run Button ---
        btn_text_map = {
            'predictions': "Run OPTICS Predictions",
            'shap': "Run SHAP Analysis",
            'structure': "Run SHAP Structure Mapping",
            'annotations': "Run Structure Annotation"
        }

        self.run_button = ctk.CTkButton(
            self.scrollable_frame, 
            text=btn_text_map[self.mode], 
            font=ctk.CTkFont(size=15, weight="bold"),
            height=45,
            command=self.start_run_thread
        )
        self.run_button.grid(row=current_row, column=0, columnspan=3, padx=10, pady=25)
        current_row += 1

        # --- Status/Output Area ---
        ctk.CTkLabel(self.scrollable_frame, font=self.lbl_font, text="Output Log:").grid(row=current_row, column=0, padx=10, pady=5, sticky=tk.W)
        current_row += 1
        
        self.output_text = ctk.CTkTextbox(self.scrollable_frame, wrap=tk.WORD, height=300, font=ctk.CTkFont(family="Consolas", size=12))
        self.output_text.grid(row=current_row, column=0, columnspan=3, padx=10, pady=(0,20), sticky=tk.NSEW)
        self.output_text.configure(state='disabled') 

        # Setup Hyperlink Tags for Output Log
        self.output_text.tag_config("hyperlink", foreground="#3498db", underline=1)
        self.output_text.tag_bind("hyperlink", "<Button-1>", self.open_hyperlink)
        self.output_text.tag_bind("hyperlink", "<Enter>", lambda e: self.output_text.configure(cursor="hand2"))
        self.output_text.tag_bind("hyperlink", "<Leave>", lambda e: self.output_text.configure(cursor=""))
        
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
        try:
            index = self.output_text.index(f"@{event.x},{event.y}")
            tags = self.output_text.tag_names(index)
            if "hyperlink" in tags:
                ranges = self.output_text.tag_ranges("hyperlink")
                for start_idx, end_idx in zip(ranges[0::2], ranges[1::2]):
                    if self.output_text.compare(start_idx, "<=", index) and self.output_text.compare(index, "<", end_idx):
                        path_to_open = self.output_text.get(start_idx, end_idx).strip()
                        if path_to_open.startswith("file:///"):
                             path_to_open = path_to_open[8:]
                        
                        if platform.system() == "Windows":
                            os.startfile(path_to_open)
                        elif platform.system() == "Darwin":
                            import subprocess
                            subprocess.call(["open", path_to_open])
                        else:
                            import subprocess
                            subprocess.call(["xdg-open", path_to_open])
                        return
        except Exception as e:
            messagebox.showerror("Error", f"Could not open path: {e}")

    def write_to_log(self, message, tag):
        try:
            self.output_text.configure(state='normal')
            
            link_token = ">>>LINK<<<"
            if link_token in message:
                parts = message.split(link_token)
                if parts[0]:
                    self.output_text.insert(tk.END, parts[0], (tag,))
                link_path = parts[1].strip()
                self.output_text.insert(tk.END, link_path, ("hyperlink",))
                self.output_text.insert(tk.END, "\n")
            else:
                self.output_text.insert(tk.END, message, (tag,))
                
            self.output_text.see(tk.END)
            self.output_text.configure(state='disabled')
        except Exception:
            pass

    def go_back(self):
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
        self.refseq_combo.configure(state="normal" if state == tk.NORMAL else tk.DISABLED) 
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
        self.fytp_combo.configure(state="normal" if visualize_state == tk.NORMAL else tk.DISABLED)

    def log_message(self, message):
        print(message) 

    def start_run_thread(self):
        # Validation
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

        if not self.output_dir_var.get():
            self.output_dir_var.set(os.path.join(os.getcwd(), 'prediction_outputs'))
            
        self.run_button.configure(state=tk.DISABLED)
        self.output_text.configure(state='normal')
        self.output_text.delete(1.0, tk.END) 
        self.output_text.configure(state='disabled')
        
        log_txts = {
            'predictions': "Starting OPTICS Predictions...",
            'shap': "Starting SHAP Analysis...",
            'structure': "Starting Structure Mapping...",
            'annotations': "Starting Structure Annotation..."
        }
        self.log_message(log_txts[self.mode])

        self.loading_screen = HumanEyeLoadingScreen(self.controller)
        self.loading_screen.start_animation()

        thread = threading.Thread(target=self.run_predictions_logic, daemon=True)
        thread.start()

    def run_predictions_logic(self):
        try:
            pred_dir_val = os.path.abspath(self.output_dir_var.get())

            if self.mode == 'predictions':
                input_val = self.input_file_var.get()
                output_val = self.prediction_prefix_var.get() or "optics_results"
                
                pred_df, output_file_path = run_optics_predictions(
                    input_sequence=input_val,
                    pred_dir=pred_dir_val,
                    output=output_val,
                    model=self.model_var.get(),
                    encoding_method=self.encoding_var.get(),
                    blastp=self.blastp_enabled_var.get(),
                    iden_report=self.blastp_report_var.get() if self.blastp_enabled_var.get() else None,
                    refseq=self.refseq_var.get() if self.blastp_enabled_var.get() else "bovine",
                    reffile=self.custom_ref_file_var.get() if self.blastp_enabled_var.get() and self.refseq_var.get() == "custom" else None,
                    bootstrap=self.bootstrap_enabled_var.get(),
                    visualize_bootstrap=self.visualize_bootstrap_var.get() if self.bootstrap_enabled_var.get() else False,
                    bootstrap_viz_file=self.bootstrap_viz_file_var.get() if self.bootstrap_enabled_var.get() and self.visualize_bootstrap_var.get() else None,
                    save_as=self.bootstrap_fytp_var.get() if self.bootstrap_enabled_var.get() and self.visualize_bootstrap_var.get() else None,
                    full_spectrum_xaxis=self.viz_xaxis_scale.get() if self.bootstrap_enabled_var.get() and self.visualize_bootstrap_var.get() else None,
                    model_version=self.version_var.get(),
                    tolerate_non_standard_aa=self.non_standard_aa_var.get(),
                    tolerate_incomplete_seqs=self.incomp_seqs_var.get()
                )
                
                self.log_message(f"\n--- Predictions Complete ---")
                if output_file_path:
                    final_dir = os.path.dirname(os.path.abspath(output_file_path))
                    self.log_message(f"Results located at: >>>LINK<<<{final_dir}")
                    self.controller.after(0, lambda: messagebox.showinfo("Success", f"OPTICS predictions completed successfully!\nResults are in: {final_dir}"))
                
            elif self.mode == 'shap':
                output_val = self.prediction_prefix_var.get() or "optics_results"
                
                try:
                    n_pos = int(self.n_positions_var.get())
                except ValueError:
                    n_pos = 10
                    
                generate_shap_explanation(
                    input_file=self.input_file_var.get(), 
                    pred_dir=pred_dir_val, 
                    output=output_val, 
                    save_as=self.shap_save_as_var.get(), 
                    model=self.model_var.get(), 
                    encoding_method=self.encoding_var.get(), 
                    model_version=self.version_var.get(), 
                    cmd_line="GUI_Execution", 
                    mode=self.shap_mode_var.get(), 
                    n_positions=n_pos, 
                    use_reference_sites=self.use_ref_sites_var.get()
                )
                self.log_message(f"\n--- SHAP Analysis Complete ---")
                self.log_message(f"Results located at: >>>LINK<<<{pred_dir_val}")
                self.controller.after(0, lambda: messagebox.showinfo("Success", f"SHAP analysis completed successfully!\nResults are in: {pred_dir_val}"))

            elif self.mode == 'structure':
                pdb_inputs = self.pdb_input_var.get()
                if hasattr(self, 'pdb_input_2_var') and self.pdb_input_2_var.get():
                    pdb_inputs = f"{pdb_inputs},{self.pdb_input_2_var.get()}"
                
                try:
                    top_n_val = int(self.top_n_labels_var.get())
                except ValueError:
                    top_n_val = 10

                pdb_out = run_structural_mapping(
                    shap_csv=self.shap_csv_var.get(),
                    pdb_input=pdb_inputs,
                    output_dir=pred_dir_val,
                    use_query_position=self.use_query_pos_var.get(),
                    chain=self.chain_var.get(),
                    map_to_bovine_also=self.map_bovine_also_var.get(),
                    comp_target=self.comp_target_var.get(),
                    top_n_labels=top_n_val,
                    software=self.struct_software_var.get().lower()
                )
                
                self.log_message(f"\n--- Structure Mapping Complete ---")
                if pdb_out:
                    self.log_message(f"PDB Saved to: {pdb_out}")
                    self.log_message(f"Results located at: >>>LINK<<<{pred_dir_val}")
                    self.controller.after(0, lambda: messagebox.showinfo("Success", f"Mapping completed successfully!\nSaved PDB to: {pdb_out}"))
                else:
                    self.controller.after(0, lambda: messagebox.showerror("Error", "Structure mapping failed. Check log for details."))
            
            elif self.mode == 'annotations':
                software_val = self.software_var.get().lower()
                
                run_structure_annotation(
                    annotation_file=self.annotation_csv_var.get(),
                    pdb_input=self.pdb_input_var.get(),
                    output_dir=pred_dir_val,
                    chain=self.chain_var.get(),
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
            self.controller.after(0, lambda: self.run_button.configure(state=tk.NORMAL))
            if self.loading_screen:
                self.controller.after(0, self.loading_screen.stop_animation)


# --- Main Application Controller ---
class OpticsApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("OPTICS Pipeline")
        self.geometry("1100x950")
        
        # --- Theme ---
        self.dark_mode_enabled = True
        
        # --- Logo Preparation ---
        self.logo_icon_display_path = "./data/logo/optics_logo_resized.png"
        self.ctk_logo_img = None
        self.ctk_logo_small = None
        
        try:
            pil_image = Image.open(self.logo_icon_display_path)
            # Full size for main menu
            self.ctk_logo_img = ctk.CTkImage(light_image=pil_image, dark_image=pil_image, size=(120, 120))
            # Smaller size for the top bar
            self.ctk_logo_small = ctk.CTkImage(light_image=pil_image, dark_image=pil_image, size=(40, 40))
            
            # Application Window Icon (requires standard tkinter PhotoImage in some OSs)
            icon_img = tk.PhotoImage(file=self.logo_icon_display_path)
            self.iconphoto(True, icon_img)
        except Exception as e:
            print(f"Logo not found or could not be loaded: {e}")

        # --- Logging Setup ---
        self.log_queue = queue.Queue()
        sys.stdout = TextRedirector(self.log_queue, "stdout")
        sys.stderr = TextRedirector(self.log_queue, "stderr")

        # --- Frame Container ---
        self.container = ctk.CTkFrame(self, fg_color="transparent")
        self.container.pack(fill="both", expand=True)
        
        self.current_frame = None
        
        self.show_mode_selector()
        self.update_log_from_queue()

    def show_mode_selector(self):
        if self.current_frame:
            self.current_frame.destroy()
        self.current_frame = ModeSelectorFrame(self.container, self)
        self.current_frame.pack(fill="both", expand=True)

    def show_optics_gui(self, mode):
        if self.current_frame:
            self.current_frame.destroy()
        self.current_frame = OpticsGUIFrame(self.container, self, mode)
        self.current_frame.pack(fill="both", expand=True)

    def toggle_theme(self):
        if self.dark_mode_enabled:
            ctk.set_appearance_mode("Light")
            self.dark_mode_enabled = False
        else:
            ctk.set_appearance_mode("Dark")
            self.dark_mode_enabled = True
            
        # Re-apply color for eye animation background if present
        if self.current_frame and hasattr(self.current_frame, 'loading_screen') and self.current_frame.loading_screen:
            bg_color = self.current_frame.loading_screen._apply_appearance_mode(ctk.ThemeManager.theme["CTkFrame"]["fg_color"])
            self.current_frame.loading_screen.canvas.config(bg=bg_color)

    def update_log_from_queue(self):
        while not self.log_queue.empty():
            try:
                tag, message = self.log_queue.get_nowait()
                if hasattr(self.current_frame, 'write_to_log'):
                    self.current_frame.write_to_log(message, tag)
            except queue.Empty:
                break
        
        self.after(100, self.update_log_from_queue)

    def destroy(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        super().destroy()


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    app = OpticsApp()
    app.mainloop()