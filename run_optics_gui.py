import tkinter as tk
import tkinter.font as tkFont
import sv_ttk 
from tkinter import ttk, filedialog, scrolledtext, messagebox
import threading
import sys
import os

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


class OpticsGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("OPTICS")
        self.geometry("900x825") \
            
        # --- Define a default font ---
        # You can choose any font available on your system
        default_font_family = "Century Gothic"
        default_font_size = 10
        # Set default font for ttk widgets
        style = ttk.Style(self)
        # Set default font for ALL ttk widgets
        style.configure(".", font=(default_font_family, default_font_size)) 
        
        # Optionally, set different fonts for specific ttk widget types
        # style.configure("TButton", font=(default_font_family, default_font_size + 1, "bold")) # Example for Buttons
        # style.configure("TLabel", font=(default_font_family, default_font_size))      # Example for Labels
        # style.configure("TEntry", font=(default_font_family, default_font_size))      # Example for Entry fields
        # style.configure("TCombobox", font=(default_font_family, default_font_size))  # For Combobox list (may need more specific styling for dropdown)
        # style.configure("TCheckbutton", font=(default_font_family, default_font_size))
        # style.configure("TLabelFrame.Label", font=(default_font_family, default_font_size, "bold")) # For LabelFrame titles
        
        # --- Logo and Window Icon ---
        self.logo_icon_display_path = "./data/logo/optics_logo_resized.png"
        try:
            self.logo_icon_img = tk.PhotoImage(file=self.logo_icon_display_path)
            self.iconphoto(True, self.logo_icon_img) # For window icon
        except tk.TclError:
            print(f"Warning: Could not load window icon: '{self.logo_icon_display_path}'. Ensure the file exists and is a valid image type.")
            self.logo_icon_img = None # To prevent errors if used later and not loaded

        # Apply the Sun Valley theme (default to dark)
        sv_ttk.set_theme("dark")
        self.dark_mode_enabled = True

        # --- Model and Encoding Choices ---
        self.model_choices = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 
                              'wildtype-vert', 'type-one', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 
                              'wildtype-vert-mnm', 'wildtype-mut']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']

        # --- Top Frame for Logo and Theme Toggle ---
        top_bar_frame = ttk.Frame(self)
        top_bar_frame.pack(fill=tk.X, padx=10, pady=(10,0))

        if self.logo_icon_img: # Display logo in GUI if loaded
            gui_logo_label = ttk.Label(top_bar_frame, image=self.logo_icon_img, text="OPTICS: Opsin Phenotype Tool for Inference of Color Sentivity", font=("Century Gothic", 16, "bold"), compound="left")
            gui_logo_label.pack(side=tk.LEFT, padx=(0,10))
             # Keep a reference if not already done by being an attribute of self
            gui_logo_label.image = self.logo_icon_img

        self.theme_toggle_button = ttk.Button(top_bar_frame, text="Toggle Light Mode", command=self.toggle_theme)
        self.theme_toggle_button.pack(side=tk.RIGHT, padx=(0,5))
        
        # --- Main Frame ---
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Create a canvas and a vertical scrollbar
        canvas = tk.Canvas(main_frame) # tk.Canvas might not be fully themed by sv_ttk
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
        self.output_text = scrolledtext.ScrolledText(scrollable_frame, wrap=tk.WORD, height=15, width=80, font=(default_font_family, default_font_size))        # Note: scrolledtext.ScrolledText is not a ttk widget, so its theming by sv_ttk might be limited.
        # For better theme consistency, you might need a ttk-based text widget or style this manually.
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
        # You might need to manually update colors for non-ttk widgets if sv_ttk doesn't handle them
        # For example, the tk.Canvas background if it doesn't change.
        # And the ScrolledText background/foreground if they look out of place.
        # Example for ScrolledText (adjust colors as needed for the theme):
        # if self.dark_mode_enabled:
        #     self.output_text.config(bg="gray15", fg="white")
        # else:
        #     self.output_text.config(bg="white", fg="black")


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
                 return # User said no, it's not a sequence
        elif not is_file_input and not is_likely_sequence and self.input_file_var.get(): # If it's short and not a file
             if not messagebox.askyesno("Input Warning", "Input is short and not a file. Is it a direct (short) sequence string?"):
                 return
        elif not is_file_input and not self.input_file_var.get(): # No input at all
            messagebox.showerror("Input Error", "Please specify an input sequence/FASTA file.")
            return


        self.run_button.config(state=tk.DISABLED)
        self.output_text.configure(state='normal')
        self.output_text.delete(1.0, tk.END) 
        self.output_text.configure(state='disabled')
        self.log_message("Starting OPTICS predictions...")

        thread = threading.Thread(target=self.run_predictions_logic, daemon=True)
        thread.start()

    def run_predictions_logic(self):
        try:
            input_val = self.input_file_var.get()
            pred_dir_val = self.output_dir_var.get()
            output_val = self.prediction_prefix_var.get()
            model_val = self.model_var.get()
            encoding_val = self.encoding_var.get()
            
            blastp_val = self.blastp_enabled_var.get()
            iden_report_val = self.blastp_report_var.get() if blastp_val else None
            refseq_val = self.refseq_var.get() if blastp_val else "bovine" 
            reffile_val = self.custom_ref_file_var.get() if blastp_val and self.refseq_var.get() == "custom" else None
            
            bootstrap_val = self.bootstrap_enabled_var.get()
            visualize_bootstrap_val = self.visualize_bootstrap_var.get() if bootstrap_val else False
            bootstrap_viz_file_val = self.bootstrap_viz_file_var.get() if bootstrap_val and visualize_bootstrap_val else None

            self.log_message(f"Parameters:\n"
                             f"  Input: {input_val}\n"
                             f"  Output Dir: {pred_dir_val}\n"
                             f"  Prediction Prefix: {output_val}\n"
                             f"  Model: {model_val}\n"
                             f"  Encoding: {encoding_val}\n"
                             f"  BLASTp: {blastp_val}\n"
                             f"    Report: {iden_report_val}\n"
                             f"    RefSeq: {refseq_val}\n"
                             f"    Custom Ref File: {reffile_val}\n"
                             f"  Bootstrap: {bootstrap_val}\n"
                             f"    Visualize: {visualize_bootstrap_val}\n"
                             f"    Viz File Prefix: {bootstrap_viz_file_val}\n"
                             f"------------------------------------")

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
                bootstrap_viz_file=bootstrap_viz_file_val
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
            self.run_button.config(state=tk.NORMAL)


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