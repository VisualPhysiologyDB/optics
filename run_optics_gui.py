import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import threading
import sys
import os

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
        self.title("OPTICS Predictions GUI")
        self.geometry("750x700") # Adjusted for more widgets

        # --- Model and Encoding Choices ---
        self.model_choices = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 
                              'wildtype-vert', 'type-one', 'whole-dataset-mnm', 
                              'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 
                              'wildtype-vert-mnm', 'wildtype-mut']
        self.encoding_choices = ['one_hot', 'aa_prop']
        self.refseq_choices = ['bovine', 'squid', 'microbe', 'custom']

        # --- Main Frame ---
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Create a canvas and a vertical scrollbar
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

        # Input Sequence/File
        ttk.Label(scrollable_frame, text="Input Sequence/FASTA File:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_file_var = tk.StringVar()
        ttk.Entry(scrollable_frame, textvariable=self.input_file_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        ttk.Button(scrollable_frame, text="Browse...", command=self.browse_input_file).grid(row=current_row, column=2, padx=5, pady=5)
        current_row += 1

        # Output Directory
        ttk.Label(scrollable_frame, text="Output Directory:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_dir_var = tk.StringVar(value=".") # Default to current directory
        ttk.Entry(scrollable_frame, textvariable=self.output_dir_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        ttk.Button(scrollable_frame, text="Browse...", command=self.browse_output_dir).grid(row=current_row, column=2, padx=5, pady=5)
        current_row += 1

        # Prediction Prefix
        ttk.Label(scrollable_frame, text="Prediction File Prefix:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.prediction_prefix_var = tk.StringVar(value="optics_predictions")
        ttk.Entry(scrollable_frame, textvariable=self.prediction_prefix_var, width=60).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1

        # Model
        ttk.Label(scrollable_frame, text="Model:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.model_var = tk.StringVar(value=self.model_choices[0])
        ttk.Combobox(scrollable_frame, textvariable=self.model_var, values=self.model_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1

        # Encoding Method
        ttk.Label(scrollable_frame, text="Encoding Method:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        self.encoding_var = tk.StringVar(value=self.encoding_choices[1]) # Default to aa_prop
        ttk.Combobox(scrollable_frame, textvariable=self.encoding_var, values=self.encoding_choices, state="readonly", width=57).grid(row=current_row, column=1, padx=5, pady=5, sticky=tk.EW)
        current_row += 1
        
        # --- BLASTp Options Frame ---
        blastp_frame = ttk.LabelFrame(scrollable_frame, text="BLASTp Options", padding="10")
        blastp_frame.grid(row=current_row, column=0, columnspan=3, padx=5, pady=10, sticky=tk.EW)
        current_row += 1

        self.blastp_enabled_var = tk.BooleanVar(value=True) # Default based on function signature
        self.blastp_check = ttk.Checkbutton(blastp_frame, text="Enable BLASTp Analysis", variable=self.blastp_enabled_var, command=self.toggle_blastp_options)
        self.blastp_check.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

        ttk.Label(blastp_frame, text="BLASTp Report Filename:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.blastp_report_var = tk.StringVar(value="blastp_report.txt")
        self.blastp_report_entry = ttk.Entry(blastp_frame, textvariable=self.blastp_report_var, width=50)
        self.blastp_report_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Label(blastp_frame, text="Reference Sequence (BLASTp):").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.refseq_var = tk.StringVar(value=self.refseq_choices[0]) # Default bovine
        self.refseq_combo = ttk.Combobox(blastp_frame, textvariable=self.refseq_var, values=self.refseq_choices, state="readonly", width=47)
        self.refseq_combo.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
        self.refseq_combo.bind("<<ComboboxSelected>>", self.toggle_custom_ref_file) # Corrected line

        ttk.Label(blastp_frame, text="Custom Reference File (BLASTp):").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        self.custom_ref_file_var = tk.StringVar()
        self.custom_ref_file_entry = ttk.Entry(blastp_frame, textvariable=self.custom_ref_file_var, width=50)
        self.custom_ref_file_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
        self.custom_ref_file_button = ttk.Button(blastp_frame, text="Browse...", command=self.browse_custom_ref_file)
        self.custom_ref_file_button.grid(row=3, column=2, padx=5, pady=5)


        # --- Bootstrap Options Frame ---
        bootstrap_frame = ttk.LabelFrame(scrollable_frame, text="Bootstrap Options", padding="10")
        bootstrap_frame.grid(row=current_row, column=0, columnspan=3, padx=5, pady=10, sticky=tk.EW)
        current_row += 1

        self.bootstrap_enabled_var = tk.BooleanVar(value=True) # Default based on function signature
        self.bootstrap_check = ttk.Checkbutton(bootstrap_frame, text="Enable Bootstrap Predictions", variable=self.bootstrap_enabled_var, command=self.toggle_bootstrap_options)
        self.bootstrap_check.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)

        self.visualize_bootstrap_var = tk.BooleanVar(value=True) # Default based on function signature
        self.visualize_check = ttk.Checkbutton(bootstrap_frame, text="Visualize Bootstrap Predictions", variable=self.visualize_bootstrap_var, command=self.toggle_bootstrap_options)
        self.visualize_check.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W)
        
        ttk.Label(bootstrap_frame, text="Bootstrap Viz Filename Prefix:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.bootstrap_viz_file_var = tk.StringVar(value="bootstrap_viz")
        self.bootstrap_viz_entry = ttk.Entry(bootstrap_frame, textvariable=self.bootstrap_viz_file_var, width=50)
        self.bootstrap_viz_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)

        # Initial state of options
        self.toggle_blastp_options()
        self.toggle_bootstrap_options()
        self.toggle_custom_ref_file() # Call once to set initial state

        # --- Run Button ---
        self.run_button = ttk.Button(scrollable_frame, text="Run OPTICS Predictions", command=self.start_run_thread)
        self.run_button.grid(row=current_row, column=0, columnspan=3, padx=5, pady=20)
        current_row += 1

        # --- Status/Output Area ---
        ttk.Label(scrollable_frame, text="Output Log:").grid(row=current_row, column=0, padx=5, pady=5, sticky=tk.W)
        current_row += 1
        self.output_text = scrolledtext.ScrolledText(scrollable_frame, wrap=tk.WORD, height=15, width=80)
        self.output_text.grid(row=current_row, column=0, columnspan=3, padx=5, pady=5, sticky=tk.NSEW)
        self.output_text.configure(state='disabled') # Make read-only initially
        
        scrollable_frame.columnconfigure(1, weight=1) # Allow entry fields to expand

        # Redirect stdout and stderr
        sys.stdout = TextRedirector(self.output_text, "stdout")
        sys.stderr = TextRedirector(self.output_text, "stderr")


    def browse_input_file(self):
        filename = filedialog.askopenfilename(title="Select Input FASTA File or Sequence File",
                                               filetypes=(("FASTA files", "*.fasta *.fa *.fna *.faa"),
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
                                               filetypes=(("FASTA files", "*.fasta *.fa *.fna *.faa"),
                                                          ("All files", "*.*")))
        if filename:
            self.custom_ref_file_var.set(filename)

    def toggle_blastp_options(self):
        state = tk.NORMAL if self.blastp_enabled_var.get() else tk.DISABLED
        self.blastp_report_entry.configure(state=state)
        self.refseq_combo.configure(state=state if state == tk.NORMAL else "disabled") # Combobox needs "disabled"
        self.toggle_custom_ref_file() # Update custom ref based on blastp state too

    def toggle_custom_ref_file(self, event=None): # event is passed by combobox command
        blastp_on = self.blastp_enabled_var.get()
        custom_ref_selected = self.refseq_var.get() == "custom"
        
        if blastp_on and custom_ref_selected:
            self.custom_ref_file_entry.configure(state=tk.NORMAL)
            self.custom_ref_file_button.configure(state=tk.NORMAL)
        else:
            self.custom_ref_file_entry.configure(state=tk.DISABLED)
            self.custom_ref_file_button.configure(state=tk.DISABLED)
            if not custom_ref_selected: # Clear if not custom, even if blastp is on
                 self.custom_ref_file_var.set("")


    def toggle_bootstrap_options(self):
        bootstrap_state = tk.NORMAL if self.bootstrap_enabled_var.get() else tk.DISABLED
        self.visualize_check.configure(state=bootstrap_state)
        
        visualize_state = tk.NORMAL if (self.bootstrap_enabled_var.get() and self.visualize_bootstrap_var.get()) else tk.DISABLED
        self.bootstrap_viz_entry.configure(state=visualize_state)


    def log_message(self, message):
        self.output_text.configure(state='normal')
        self.output_text.insert(tk.END, message + "\n")
        self.output_text.see(tk.END) # Scroll to the end
        self.output_text.configure(state='disabled')
        self.update_idletasks() # Ensure GUI updates

    def start_run_thread(self):
        # Basic Validation
        if not self.input_file_var.get():
            messagebox.showerror("Input Error", "Please specify an input sequence/FASTA file.")
            return
        if not self.output_dir_var.get():
            messagebox.showerror("Input Error", "Please specify an output directory.")
            return
        if not os.path.exists(self.input_file_var.get()) and len(self.input_file_var.get()) < 10: # crude check for sequence vs file
             if not messagebox.askyesno("Input Warning", "Input doesn't look like a file path. Is it a direct sequence string?"):
                 return # User said no, it's not a sequence
        elif not os.path.isfile(self.input_file_var.get()) and len(self.input_file_var.get()) >=10: # if it's long and not a file
            messagebox.showerror("Input Error", f"Input file not found: {self.input_file_var.get()}")
            return


        self.run_button.config(state=tk.DISABLED)
        self.output_text.configure(state='normal')
        self.output_text.delete(1.0, tk.END) # Clear previous output
        self.output_text.configure(state='disabled')
        self.log_message("Starting OPTICS predictions...")

        # Run in a separate thread
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
            refseq_val = self.refseq_var.get() if blastp_val else "bovine" # Pass default if blastp off
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

            # Call the main prediction function
            # The run_optics_predictions function is expected to print to stdout/stderr
            # which will be redirected to our text widget.
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
            self.log_message(f"Main results dataframe (first 5 rows):\n{pred_df.head().to_string()}")
            self.log_message(f"Full results written to files in directory: {os.path.dirname(output_file_path)}")
            self.log_message(f"Primary output file: {output_file_path}")
            messagebox.showinfo("Success", f"OPTICS predictions completed successfully!\nResults are in: {os.path.dirname(output_file_path)}")

        except Exception as e:
            self.log_message(f"\n--- ERROR ---")
            self.log_message(f"An error occurred: {str(e)}")
            import traceback
            self.log_message(f"Traceback:\n{traceback.format_exc()}")
            messagebox.showerror("Error", f"An error occurred during prediction: {e}")
        finally:
            self.run_button.config(state=tk.NORMAL)


class TextRedirector(object):
    """A class to redirect stdout/stderr to a Tkinter text widget."""
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str_):
        self.widget.configure(state='normal')
        self.widget.insert(tk.END, str_, (self.tag,))
        self.widget.see(tk.END) # Scroll to the end
        self.widget.configure(state='disabled')
        self.widget.update_idletasks() # Force GUI update

    def flush(self):
        # Required for file-like object interface
        pass


if __name__ == "__main__":
    # It's good practice to ensure the script directory is in sys.path
    # if optics_predictions.py might not be found by default import mechanisms.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    app = OpticsGUI()
    app.mainloop()