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
# This assumes optics_predictions.py is in the same directory or Python path
try:
    from optics_predictions import run_optics_predictions
except ImportError:
    messagebox.showerror("Import Error", "Could not import 'run_optics_predictions' from 'optics_predictions.py'. Make sure the script is in the same directory.")
    sys.exit(1)
except Exception as e:
    messagebox.showerror("Import Error", f"An unexpected error occurred during import: {e}")
    sys.exit(1)

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
        self.title("Processing Predictions...")
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
        # Gaze control for both X and Y axes
        self.gaze_current_x = 0.0
        self.gaze_target_x = 0.0
        self.gaze_current_y = 0.0
        self.gaze_target_y = 0.0
        
        # A sequence of (x, y) coordinates for the eye to look at.
        # Includes corners, edges, and center to simulate looking around a screen.
        self.gaze_positions = [
            (-40, -25), (45, 20), (0, 0), (30, 25), (-20, -20), (0, 0), 
            (-45, 0), (45, 0), (0, 30), (0, -30), (0, 0)
        ]
        self.gaze_position_index = 0
        
        # Timer for how long to wait before moving to the next gaze position
        self.gaze_timer = 100 # Start with a delay

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
        
        # Cycle through 0, 1, 2, 3 dots
        self.ellipsis_count = (self.ellipsis_count + 1) % 4
        dots = "." * self.ellipsis_count
        self.loading_label.config(text=f"Processing{dots}")
        
        self.after(400, self.animate_text) # Update text every 400ms

    def animate_eye(self):
        """The main animation loop for the eye's movement."""
        if not self.animation_running: return

        # --- Gaze Logic ---
        self.gaze_timer -= 1
        if self.gaze_timer <= 0:
            # Move to the next target position in our sequence
            self.gaze_position_index = (self.gaze_position_index + 1) % len(self.gaze_positions)
            self.gaze_target_x, self.gaze_target_y = self.gaze_positions[self.gaze_position_index]
            
            # Reset the timer. Pause longer when looking at the center.
            if self.gaze_target_x == 0 and self.gaze_target_y == 0:
                self.gaze_timer = random.randint(100, 150) # Longer pause for "staring"
            else:
                self.gaze_timer = random.randint(40, 80) # Shorter pause for glances

        # Smoothly move the current gaze towards the target for a natural effect
        self.gaze_current_x += (self.gaze_target_x - self.gaze_current_x) * 0.1
        self.gaze_current_y += (self.gaze_target_y - self.gaze_current_y) * 0.1

        # --- Drawing Logic ---
        self.canvas.delete("all")
        self.canvas.create_window(200, 265, window=self.loading_label)
        
        center_x, center_y = 200, 130
        
        # 1. Eye shape/socket
        self.canvas.create_oval(center_x - 110, center_y - 80, center_x + 110, center_y + 80, fill="#111111", outline="")

        # 2. Sclera (white of the eye)
        self.canvas.create_oval(center_x - 100, center_y - 70, center_x + 100, center_y + 70, fill="#EAEAEA", outline="")

        # 3. Iris and Pupil (position controlled by both x and y gaze)
        iris_x = center_x + self.gaze_current_x
        iris_y = center_y + self.gaze_current_y
        iris_radius = 35
        pupil_radius = 15
        self.canvas.create_oval(iris_x - iris_radius, iris_y - iris_radius, iris_x + iris_radius, iris_y + iris_radius, fill="#5DADE2", outline="")
        self.canvas.create_oval(iris_x - pupil_radius, iris_y - pupil_radius, iris_x + pupil_radius, iris_y + pupil_radius, fill="black", outline="")

        # 4. Specular Highlight (moves with the iris)
        highlight_x = iris_x + 10
        highlight_y = iris_y - 10
        self.canvas.create_oval(highlight_x - 5, highlight_y - 5, highlight_x + 3, highlight_y + 3, fill="white", outline="")

        # Schedule the next frame
        self.animation_id = self.after(30, self.animate_eye)


class OpticsGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("OPTICS")
        self.geometry("900x825") 
        
        # --- Add a reference for the loading screen ---
        self.loading_screen = None
            
        # --- Define a default font ---
        default_font_family = "Century Gothic"
        default_font_size = 10
        style = ttk.Style(self)
        style.configure(".", font=(default_font_family, default_font_size)) 
        
        # --- Logo and Window Icon ---
        self.logo_icon_display_path = "./data/logo/optics_logo_resized.png"
        try:
            self.logo_icon_img = tk.PhotoImage(file=self.logo_icon_display_path)
            self.iconphoto(True, self.logo_icon_img) # For window icon
        except tk.TclError:
            print(f"Warning: Could not load window icon: '{self.logo_icon_display_path}'. Ensure the file exists and is a valid image type.")
            self.logo_icon_img = None 

        # Apply the Sun Valley theme (default to dark)
        sv_ttk.set_theme("dark")
        self.dark_mode_enabled = True

        # --- Model, Encoding, and Bootstrap Visualization Choices ---
        self.version_choices = ['vpod_1.3']
        self.model_choices = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 
                              'wildtype-vert', 'type-one', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 
                              'wildtype-vert-mnm', 'wildtype-mut']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']
        self.viz_ftyp_choices = ['svg', 'png', 'pdf']

        # --- Top Frame for Logo and Theme Toggle ---
        top_bar_frame = ttk.Frame(self)
        top_bar_frame.pack(fill=tk.X, padx=10, pady=(10,0))

        if self.logo_icon_img: # Display logo in GUI if loaded
            gui_logo_label = ttk.Label(top_bar_frame, image=self.logo_icon_img, text="OPTICS: Opsin Phenotype Tool for Inference of Color Sentivity", font=("Century Gothic", 16, "bold"), compound="left")
            gui_logo_label.pack(side=tk.LEFT, padx=(0,10))
            gui_logo_label.image = self.logo_icon_img

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
        
        # --- Input Widgets within scrollable_frame ---
        current_row = 0
        scrollable_frame.columnconfigure(1, weight=1)

        # Input Sequence/File
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Input Sequence/FASTA File:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_file_var = tk.StringVar()
        ttk.Entry(scrollable_frame, textvariable=self.input_file_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        ttk.Button(scrollable_frame, text="Browse...", command=self.browse_input_file).grid(row=current_row, column=2, padx=5, pady=5)
        current_row += 1

        # Output Directory
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Output Directory:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_dir_var = tk.StringVar()
        ttk.Entry(scrollable_frame, textvariable=self.output_dir_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        ttk.Button(scrollable_frame, text="Browse...", command=self.browse_output_dir).grid(row=current_row, column=2, padx=5, pady=5)
        current_row += 1

        # Prediction Prefix
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Prediction Outputs Prefix:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.prediction_prefix_var = tk.StringVar()
        ttk.Label(scrollable_frame, font=("Century Gothic", 9), text="\n\n\n(No extensions; used to name output folder and files)...").grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.W)
        ttk.Entry(scrollable_frame, textvariable=self.prediction_prefix_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1
        
        # OPTICS Version
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Model Version:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.version_var = tk.StringVar(value=self.version_choices[0])
        ttk.Combobox(scrollable_frame, textvariable=self.version_var, values=self.version_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1

        # Model
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Model:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.model_var = tk.StringVar(value=self.model_choices[0])
        ttk.Combobox(scrollable_frame, textvariable=self.model_var, values=self.model_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1

        # Encoding Method
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Encoding Method:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.encoding_var = tk.StringVar(value=self.encoding_choices[1]) # Default to aa_prop
        ttk.Combobox(scrollable_frame, textvariable=self.encoding_var, values=self.encoding_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1
        
        # Allow predictions on non-standard AA option
        ttk.Label(scrollable_frame, font=("Century Gothic", 12), text="Tolerate Non-standard AAs:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.non_standard_aa_var = tk.BooleanVar(value=True) 
        self.non_standard_aa_check = ttk.Checkbutton(scrollable_frame, text="Enable to predict on sequences with non-standard amino-acids", variable=self.non_standard_aa_var)
        self.non_standard_aa_check.grid(row=current_row, column=1, columnspan=2, padx=5, pady=5, sticky=tk.EW)
        current_row += 1
        
        # --- BLASTp Options Frame ---
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

        ttk.Label(blastp_frame, font=("Century Gothic", 12), text="Reference Sequence (BLASTp):").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.refseq_var = tk.StringVar(value=self.refseq_choices[0]) 
        self.refseq_combo = ttk.Combobox(blastp_frame, textvariable=self.refseq_var, values=self.refseq_choices, state="readonly", width=47)
        self.refseq_combo.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
        self.refseq_combo.bind("<<ComboboxSelected>>", self.toggle_custom_ref_file) 

        ttk.Label(blastp_frame, font=("Century Gothic", 12), text="Custom Reference File (BLASTp):").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        self.custom_ref_file_var = tk.StringVar()
        self.custom_ref_file_entry = ttk.Entry(blastp_frame, textvariable=self.custom_ref_file_var, width=50)
        self.custom_ref_file_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
        self.custom_ref_file_button = ttk.Button(blastp_frame, text="Browse...", command=self.browse_custom_ref_file)
        self.custom_ref_file_button.grid(row=3, column=2, padx=5, pady=5)


        # --- Bootstrap Options Frame ---
        bootstrap_frame = ttk.LabelFrame(scrollable_frame, text="Bootstrap Options", padding="10")
        bootstrap_frame.grid(row=current_row, column=0, columnspan=3, padx=5, pady=10, sticky=tk.EW)
        current_row += 1

        self.bootstrap_enabled_var = tk.BooleanVar(value=True) 
        self.bootstrap_check = ttk.Checkbutton(bootstrap_frame, text="Enable Bootstrap Predictions", variable=self.bootstrap_enabled_var, command=self.toggle_bootstrap_options)
        self.bootstrap_check.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

        self.visualize_bootstrap_var = tk.BooleanVar(value=True) 
        self.visualize_check = ttk.Checkbutton(bootstrap_frame, text="Visualize Bootstrap Predictions", variable=self.visualize_bootstrap_var, command=self.toggle_bootstrap_options)
        self.visualize_check.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)
        
        ttk.Label(bootstrap_frame, font=("Century Gothic", 12), text="Bootstrap Viz Filename Prefix:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.bootstrap_viz_file_var = tk.StringVar(value="bootstrap_viz")
        self.bootstrap_viz_entry = ttk.Entry(bootstrap_frame, textvariable=self.bootstrap_viz_file_var, width=50)
        self.bootstrap_viz_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
        
        
        ttk.Label(bootstrap_frame, font=("Century Gothic", 12), text="Bootstrap Viz Filetype:").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        self.bootstrap_fytp_var = tk.StringVar(value=self.viz_ftyp_choices[0]) # Default to svg
        self.fytp_combo = ttk.Combobox(bootstrap_frame, textvariable=self.bootstrap_fytp_var, values=self.viz_ftyp_choices, state="readonly", width=47)
        self.fytp_combo.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
        
        self.viz_xaxis_scale = tk.BooleanVar(value=False) 
        self.xaxis_check = ttk.Checkbutton(bootstrap_frame, text="Enable Full-Spectrum X-axis (300-650nm). Default = Scaled to Predictions", variable=self.viz_xaxis_scale)
        self.xaxis_check.grid(row=4, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

        # Initial state of options
        self.toggle_blastp_options()
        self.toggle_bootstrap_options()
        self.toggle_custom_ref_file() 

        # --- Run Button ---
        self.run_button = ttk.Button(scrollable_frame, text="Run OPTICS Predictions", command=self.start_run_thread)
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
        filename = filedialog.askopenfilename(title="Select Input FASTA File or Sequence File",
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

    def toggle_blastp_options(self):
        state = tk.NORMAL if self.blastp_enabled_var.get() else tk.DISABLED
        self.blastp_report_entry.configure(state=state)
        self.refseq_combo.configure(state="readonly" if state == tk.NORMAL else tk.DISABLED) 
        self.toggle_custom_ref_file() 

    def toggle_custom_ref_file(self, event=None): 
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
        if not self.input_file_var.get():
            messagebox.showerror("Input Error", "Please specify an input sequence/FASTA file.")
            return
        if not self.output_dir_var.get():
            messagebox.showerror("Input Error", "Please specify an output directory.")
            return
        
        is_file_input = os.path.isfile(self.input_file_var.get())
        is_likely_sequence = not is_file_input and len(self.input_file_var.get()) > 10 # Heuristic for sequence

        if not is_file_input and is_likely_sequence :
             if not messagebox.askyesno("Input Warning", "Input doesn't look like a file path. Is it a direct sequence string?"):
                 return 
        elif not is_file_input and not is_likely_sequence and self.input_file_var.get(): 
             if not messagebox.askyesno("Input Warning", "Input is short and not a file. Is it a direct (short) sequence string?"):
                 return
        elif not is_file_input and not self.input_file_var.get(): 
            messagebox.showerror("Input Error", "Please specify an input sequence/FASTA file.")
            return


        self.run_button.config(state=tk.DISABLED)
        self.output_text.configure(state='normal')
        self.output_text.delete(1.0, tk.END) 
        self.output_text.configure(state='disabled')
        self.log_message("Starting OPTICS predictions...")

        # Create and start the loading screen animation
        self.loading_screen = HumanEyeLoadingScreen(self)
        self.loading_screen.start_animation()

        # Start the prediction logic in a separate thread
        thread = threading.Thread(target=self.run_predictions_logic, daemon=True)
        thread.start()

    def run_predictions_logic(self):
        try:
            input_val = self.input_file_var.get()
            pred_dir_val = self.output_dir_var.get()
            output_val = self.prediction_prefix_var.get()
            version_val = self.version_var.get()
            model_val = self.model_var.get()
            encoding_val = self.encoding_var.get()
            tol_non_stan_aa = self.non_standard_aa_var.get()
            
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
                tolerate_non_standard_aa=tol_non_stan_aa
            )
            
            self.log_message(f"\n--- Predictions Complete ---")
            if output_file_path:
                self.log_message(f"Full results written to files in directory: {os.path.dirname(output_file_path)}")
                self.log_message(f"Primary output file: {output_file_path}")
                messagebox.showinfo("Success", f"OPTICS predictions completed successfully!\nResults are in: {os.path.dirname(output_file_path)}")
            else:
                self.log_message("Output file path not available.")
                messagebox.showinfo("Success", "OPTICS predictions completed, but output path was not returned.")


        except Exception as e:
            self.log_message(f"\n--- ERROR ---")
            self.log_message(f"An error occurred: {str(e)}")
            import traceback
            self.log_message(f"Traceback:\n{traceback.format_exc()}")
            messagebox.showerror("Error", f"An error occurred during prediction: {e}")
        finally:
            # This block runs whether there was an error or not.
            # It's crucial for re-enabling the UI.
            self.run_button.config(state=tk.NORMAL)
            
            # Schedule the loading screen to close.
            # We use self.after() to ensure this GUI operation runs on the main thread.
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


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    app = OpticsGUI()
    app.mainloop()
