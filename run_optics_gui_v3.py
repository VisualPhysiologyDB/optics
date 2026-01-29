import tkinter as tk
import tkinter.font as tkFont
import sv_ttk 
from tkinter import ttk, filedialog, scrolledtext, messagebox
import threading
import sys
import os
import math
import random

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter('ignore')

# Attempt to import the run_optics_predictions function
try:
    from optics_predictions import run_optics_predictions
except ImportError:
    messagebox.showerror("Import Error", "Could not import 'run_optics_predictions' from 'optics_predictions.py'. Make sure the script is in the same directory.")
    sys.exit(1)

# Attempt to import the generate_shap_explanation function
try:
    from optics_shap import generate_shap_explanation
except ImportError:
    print("Warning: Could not import 'generate_shap_explanation' from 'optics_shap.py'. SHAP mode will fail if selected.")

# Attempt to import the run_structural_mapping function
try:
    from optics_structure_map import run_structural_mapping
except ImportError:
    print("Warning: Could not import 'run_structural_mapping' from 'optics_structure_map.py'. Structure Mapping mode will fail if selected.")


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
        parent_x = parent.winfo_x()
        parent_y = parent.winfo_y()
        parent_width = parent.winfo_width()
        parent_height = parent.winfo_height()
        self.after(10, lambda: self.geometry(f"400x300+{parent_x + (parent_width // 2) - 200}+{parent_y + (parent_height // 2) - 150}"))

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


class ModeSelector(tk.Tk):
    """
    Startup dialog to select between Predictions, SHAP analysis, and Structure Mapping.
    """
    def __init__(self):
        super().__init__()
        self.title("Select OPTICS Mode")
        self.geometry("500x500")
        self.resizable(False, False)
        self.selected_mode = None
        
        sv_ttk.set_theme("dark")

        # Logo
        self.logo_icon_display_path = "./data/logo/optics_logo_resized.png"
        try:
            self.logo_icon_img = tk.PhotoImage(file=self.logo_icon_display_path)
            self.iconphoto(True, self.logo_icon_img) 
        except tk.TclError:
            self.logo_icon_img = None

        main_frame = ttk.Frame(self, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        if self.logo_icon_img:
            logo_lbl = ttk.Label(main_frame, image=self.logo_icon_img)
            logo_lbl.pack(pady=(0, 10))

        title_lbl = ttk.Label(main_frame, text="Welcome to OPTICS", font=("Century Gothic", 18, "bold"))
        title_lbl.pack(pady=10)
        
        subtitle_lbl = ttk.Label(main_frame, text="Select your analysis pipeline:", font=("Century Gothic", 12))
        subtitle_lbl.pack(pady=(0, 20))

        # Buttons
        style = ttk.Style()
        style.configure("Big.TButton", font=("Century Gothic", 12))

        pred_btn = ttk.Button(main_frame, text="Standard Predictions\n(Spectral Tuning)", style="Big.TButton", command=self.select_predictions)
        pred_btn.pack(fill=tk.X, pady=10, ipady=10)

        shap_btn = ttk.Button(main_frame, text="SHAP Interpretation\n(Feature Importance)", style="Big.TButton", command=self.select_shap)
        shap_btn.pack(fill=tk.X, pady=10, ipady=10)

        struct_btn = ttk.Button(main_frame, text="Structure Mapping\n(3D Visualization)", style="Big.TButton", command=self.select_structure)
        struct_btn.pack(fill=tk.X, pady=10, ipady=10)
        
        # Centering
        self.eval('tk::PlaceWindow . center')

    def select_predictions(self):
        self.selected_mode = 'predictions'
        self.destroy()

    def select_shap(self):
        self.selected_mode = 'shap'
        self.destroy()

    def select_structure(self):
        self.selected_mode = 'structure'
        self.destroy()


class OpticsGUI(tk.Tk):
    def __init__(self, mode='predictions'):
        super().__init__()
        self.mode = mode
        
        if self.mode == 'predictions':
            title_suffix = "Predictions"
        elif self.mode == 'shap':
            title_suffix = "SHAP Analysis"
        else:
            title_suffix = "Structure Mapping"
            
        self.title(f"OPTICS: {title_suffix}")
        self.geometry("1000x950") 
        
        # --- Define a default font ---
        default_font_family = "Century Gothic"
        default_font_size = 10
        style = ttk.Style(self)
        style.configure(".", font=(default_font_family, default_font_size)) 
        
        # --- Logo and Window Icon ---
        self.logo_icon_display_path = "./data/logo/optics_logo_resized.png"
        try:
            self.logo_icon_img = tk.PhotoImage(file=self.logo_icon_display_path)
            self.iconphoto(True, self.logo_icon_img) 
        except tk.TclError:
            self.logo_icon_img = None 

        sv_ttk.set_theme("dark")
        self.dark_mode_enabled = True

        # --- Choices ---
        self.version_choices = ['vpod_1.3']
        self.model_choices = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 
                              'wildtype-vert', 'type-one', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 
                              'wildtype-vert-mnm', 'wildtype-mut']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']
        self.viz_ftyp_choices = ['svg', 'png', 'pdf']
        
        # SHAP specific choices
        self.shap_mode_choices = ['both', 'comparison', 'single']

        # --- Top Frame ---
        top_bar_frame = ttk.Frame(self)
        top_bar_frame.pack(fill=tk.X, padx=10, pady=(10,0))

        header_text = f"OPTICS: {title_suffix}"
        if self.logo_icon_img:
            gui_logo_label = ttk.Label(top_bar_frame, image=self.logo_icon_img, text=header_text, font=("Century Gothic", 16, "bold"), compound="left")
            gui_logo_label.pack(side=tk.LEFT, padx=(0,10))
            gui_logo_label.image = self.logo_icon_img
        else:
            gui_logo_label = ttk.Label(top_bar_frame, text=header_text, font=("Century Gothic", 16, "bold"))
            gui_logo_label.pack(side=tk.LEFT, padx=(0,10))


        self.theme_toggle_button = ttk.Button(top_bar_frame, text="Toggle Light Mode", command=self.toggle_theme)
        self.theme_toggle_button.pack(side=tk.RIGHT, padx=(0,5))
        
        # --- Main Frame ---
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
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="SHAP Analysis CSV File:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.shap_csv_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.shap_csv_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            ttk.Button(scrollable_frame, text="Browse...", command=self.browse_shap_csv).grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1

            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="PDB File or ID (e.g. 1U19):").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.pdb_input_var = tk.StringVar()
            ttk.Entry(scrollable_frame, textvariable=self.pdb_input_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            ttk.Button(scrollable_frame, text="Browse File...", command=self.browse_pdb_file).grid(row=current_row, column=2, padx=5, pady=5)
            current_row += 1
            
            ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Alignment Offset (Integer):").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
            self.offset_var = tk.IntVar(value=0)
            ttk.Spinbox(scrollable_frame, from_=-1000, to=1000, textvariable=self.offset_var, width=58).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
            current_row += 1
            
            ttk.Label(scrollable_frame, text="(Adjusts residue numbering to match PDB. Default 0)", font=("Century Gothic", 9)).grid(row=current_row, column=1, padx=5, pady=0, sticky=tk.W)
            current_row += 1

        
        # --- MODE SPECIFIC OPTIONS ---
        
        if self.mode == 'predictions':
            # -- Prediction specific options --
            
            # Tolerances
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
            self.use_ref_sites_var = tk.BooleanVar(value=False)
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
        else:
             btn_text = "Run Structure Mapping"

        self.run_button = ttk.Button(scrollable_frame, text=btn_text, command=self.start_run_thread)
        self.run_button.grid(row=current_row, column=0, columnspan=3, padx=5, pady=20)
        current_row += 1

        # --- Status/Output Area ---
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Output Log:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        current_row += 1
        self.output_text = scrolledtext.ScrolledText(scrollable_frame, wrap=tk.WORD, height=15, width=80, font=(default_font_family, default_font_size))
        self.output_text.grid(row=current_row, column=0, columnspan=3, padx=5, pady=5, sticky=tk.NSEW)
        self.output_text.configure(state='disabled') 
        
        scrollable_frame.columnconfigure(1, weight=1) 

        # Redirect stdout and stderr
        sys.stdout = TextRedirector(self.output_text, "stdout")
        sys.stderr = TextRedirector(self.output_text, "stderr")
        
        # Loading Screen placeholder
        self.loading_screen = None

    def toggle_theme(self):
        if self.dark_mode_enabled:
            sv_ttk.set_theme("light")
            self.theme_toggle_button.configure(text="Toggle Dark Mode")
            self.dark_mode_enabled = False
        else:
            sv_ttk.set_theme("dark")
            self.theme_toggle_button.configure(text="Toggle Light Mode")
            self.dark_mode_enabled = True

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

    # --- New Browsers for Structure Mode ---
    def browse_shap_csv(self):
        filename = filedialog.askopenfilename(title="Select SHAP Analysis CSV",
                                               filetypes=(("CSV files", "*.csv"),
                                                          ("All files", "*.*")))
        if filename:
            self.shap_csv_var.set(filename)

    def browse_pdb_file(self):
        filename = filedialog.askopenfilename(title="Select PDB File",
                                               filetypes=(("PDB files", "*.pdb"),
                                                          ("All files", "*.*")))
        if filename:
            self.pdb_input_var.set(filename)


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
        self.output_text.configure(state='normal')
        self.output_text.insert(tk.END, message + "\n")
        self.output_text.see(tk.END) 
        self.output_text.configure(state='disabled')
        self.update_idletasks() 

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
        else:
            log_txt = "Starting Structure Mapping..."

        self.log_message(log_txt)

        self.loading_screen = HumanEyeLoadingScreen(self)
        self.loading_screen.start_animation()

        thread = threading.Thread(target=self.run_predictions_logic, daemon=True)
        thread.start()

    def run_predictions_logic(self):
        try:
            # Common arguments
            pred_dir_val = self.output_dir_var.get()
            
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
                    messagebox.showinfo("Success", f"OPTICS predictions completed successfully!\nResults are in: {os.path.dirname(output_file_path)}")
                
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
                messagebox.showinfo("Success", f"SHAP analysis completed successfully!\nResults are in: {pred_dir_val}")

            elif self.mode == 'structure':
                csv_path = self.shap_csv_var.get()
                pdb_input = self.pdb_input_var.get()
                offset = self.offset_var.get()

                pdb_out = run_structural_mapping(
                    shap_csv=csv_path,
                    pdb_input=pdb_input,
                    output_dir=pred_dir_val,
                    offset=offset
                )
                
                self.log_message(f"\n--- Structure Mapping Complete ---")
                if pdb_out:
                    self.log_message(f"PDB Saved to: {pdb_out}")
                    messagebox.showinfo("Success", f"Mapping completed successfully!\nSaved PDB to: {pdb_out}")
                else:
                    messagebox.showerror("Error", "Structure mapping failed. Check log for details.")


        except Exception as e:
            self.log_message(f"\n--- ERROR ---")
            self.log_message(f"An error occurred: {str(e)}")
            import traceback
            self.log_message(f"Traceback:\n{traceback.format_exc()}")
            messagebox.showerror("Error", f"An error occurred: {e}")
        finally:
            self.run_button.config(state=tk.NORMAL)
            if self.loading_screen:
                self.after(0, self.loading_screen.stop_animation)


class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str_):
        self.widget.configure(state='normal')
        self.widget.insert(tk.END, str_, (self.tag,))
        self.widget.see(tk.END) 
        self.widget.configure(state='disabled')
        self.widget.update_idletasks() 

    def flush(self):
        pass


def main():
    # 1. Run Selector
    selector = ModeSelector()
    selector.mainloop()
    
    # 2. Get Selection
    mode = selector.selected_mode
    
    # 3. Launch Main App if selection was made
    if mode:
        app = OpticsGUI(mode=mode)
        app.mainloop()

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    main()