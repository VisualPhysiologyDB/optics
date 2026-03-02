import argparse
import os
import sys
import shutil
import tempfile
import glob
import logging
from PIL import Image

# Set up logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# --- IMPORTANT ---
# This wrapper script assumes that 'optics_shap.py' is in the same directory
try:
    import optics_shap
except ImportError:
    log.error("Fatal Error: Could not import 'optics_shap.py'.")
    log.error("Make sure 'optics_shap.py' is in the same directory as this wrapper.")
    sys.exit(1)
# --- /IMPORTANT ---


def main():
    parser = argparse.ArgumentParser(description="Galaxy wrapper for OPTICS SHAP comparison.")

    # --- Input Files ---
    parser.add_argument("--input", required=True, help="Path to input FASTA file (must contain exactly 2 sequences).")

    # --- Output Files ---
    parser.add_argument("--output_csv", required=True, help="Path for the final output CSV with all SHAP values.")
    parser.add_argument("--output_viz", required=True, help="Path for the final SHAP comparison visualization.")
    parser.add_argument("--output_log", required=True, help="Path for the SHAP script log file.")

    # --- Core OPTICS SHAP Args ---
    parser.add_argument("--prediction_prefix", default="optics_shap_run", help="Prefix for output files.")
    parser.add_argument("--model_version", default="vpod_1.3", help="Model version (currently only vpod_1.3).")
    parser.add_argument("--model", default="whole-dataset", help="Prediction model to use.")
    parser.add_argument("--encoding", default="aa_prop", help="Encoding method.")
    parser.add_argument("--save_viz_as", default="svg", help="Format for SHAP visualization (svg, png, pdf).")
    parser.add_argument("--mode", 
                        help="Analysis mode: 'comparison' (pairwise), 'single' (individual), or 'both'.", 
                        type=str, default="both", choices=['comparison', 'single', 'both'])
    parser.add_argument("--use_reference_sites", help="Use reference site numbering, instead of feature names", action="store_true")


    args = parser.parse_args()

    # Create a temporary directory for optics_shap.py to write into.
    temp_run_dir = tempfile.mkdtemp()
    log.info(f"Created temporary run directory: {temp_run_dir}")

    # Capture the full command for logging
    cmd_line = " ".join(sys.argv)

    try:
        # --- Run the main OPTICS SHAP script ---
        log.info("Starting optics_shap.generate_shap_explanation...")
        optics_shap.generate_shap_explanation(
            input_file=args.input,
            pred_dir=temp_run_dir,  # Tell script to use our temp dir
            output=args.prediction_prefix,
            save_as=args.save_viz_as,
            model=args.model,
            encoding_method=args.encoding,
            model_version=args.model_version,
            cmd_line=cmd_line,
            mode=args.mode,
            use_reference_sites=args.use_reference_sites,
            n_jobs=-1
        )
        log.info("optics_shap.generate_shap_explanation finished.")

        # --- Find the output directory ---
        # The script creates a timestamped subdirectory like 'optics_shap_on_...'.
        report_dirs = glob.glob(os.path.join(temp_run_dir, f'optics_shap_on_{args.prediction_prefix}_*'))
        if not report_dirs:
            # Check for the fallback directory if prefix was 'unnamed' and script defaulted
            if args.prediction_prefix == "unnamed":
                report_dirs = glob.glob(os.path.join(temp_run_dir, 'optics_on_unamed_*'))
            
            if not report_dirs:
                log.error(f"Fatal Error: Could not find output directory from SHAP script inside {temp_run_dir}")
                log.error(f"Expected pattern: optics_shap_on_{args.prediction_prefix}_*")
                # Log directory contents for debugging
                log.error(f"Contents of temp_run_dir ({temp_run_dir}): {os.listdir(temp_run_dir)}")
                sys.exit(1)
        
        actual_report_dir = report_dirs[0]
        log.info(f"Found report directory: {actual_report_dir}")

        # --- Define paths to generated files ---
        gen_csv = os.path.join(actual_report_dir, f"{args.prediction_prefix}_all_shap_values.csv")
        gen_viz = os.path.join(actual_report_dir, f"{args.prediction_prefix}_viz_shap_comparison.{args.save_viz_as}")
        gen_log = os.path.join(actual_report_dir, "shap_script_log.txt")

        # --- List of output datasets that are NOT optional ---
        # We expect all three files to be created.
        always_expected_outputs = {
            args.output_csv,
            args.output_viz,
            args.output_log
        }

        def copy_file(src, dst):
            if dst and dst != 'null':
                if os.path.exists(src):
                    shutil.copy(src, dst)
                    log.info(f"Copied {src} to {dst}")
                else:
                    if dst in always_expected_outputs:
                        log.warning(f"Warning: Expected output file '{src}' was not found.")
                        log.warning(f"Creating empty file at {dst} to satisfy Galaxy.")
                        with open(dst, 'w') as f:
                            f.write(f'File "{os.path.basename(src)}" was expected but not generated by optics_shap.py.\n')
                    else:
                        log.info(f"Optional file '{src}' was not generated. Skipping copy to {dst}.")
            else:
                log.info(f"Skipping copy for {src} (destination is 'null').")

        # --- Move files to final Galaxy paths ---
        log.info("Copying generated files to final Galaxy paths...")
        copy_file(gen_csv, args.output_csv)
        copy_file(gen_log, args.output_log)
        
        # --- New logic for compiling PNGs into a single PDF ---
        # We replace the direct copy of gen_viz with this conditional block.
        if args.save_viz_as == 'png':
            # Target files by the '.png' extension as requested, ignoring prefix patterns
            viz_files = sorted(glob.glob(os.path.join(actual_report_dir, "*.png")))
            
            if not viz_files:
                log.warning(f"Visualization was set to PNG, but no .png files were found in {actual_report_dir}")
                # Fallback to creating an empty file if needed to satisfy Galaxy
                if args.output_viz in always_expected_outputs:
                     open(args.output_viz, 'a').close()
            else:
                try:
                    # Open all PNGs found
                    image_list = [Image.open(f).convert('RGB') for f in viz_files]
                    if image_list:
                        temp_pdf_name = os.path.join(actual_report_dir, 'compiled_shap_viz.pdf')
                        image_list[0].save(
                            temp_pdf_name, 
                            save_all=True, 
                            append_images=image_list[1:],
                            title="SHAP Comparison Visualizations"
                        )
                        copy_file(temp_pdf_name, args.output_viz)
                        log.info(f"Successfully compiled {len(image_list)} PNGs into single PDF: {args.output_viz}")
                    else:
                         log.warning("Found PNG files but failed to open them with Pillow.")
                         open(args.output_viz, 'a').close()

                except Exception as e:
                    log.error(f"Error compiling PNGs into PDF: {e}")
                    # Fallback: Just copy the first part found to avoid a crash
                    copy_file(viz_files[0], args.output_viz)
        else:
            # For SVG or PDF, direct copy of the main expected file
            copy_file(gen_viz, args.output_viz)

    except Exception as e:
        log.error(f"An error occurred during the OPTICS SHAP run: {e}")
        # Write the error to the log file if possible
        try:
            with open(args.output_log, 'w') as f:
                f.write("A fatal error occurred during the wrapper's 'try' block.\n")
                f.write(str(e))
                f.write("\n")
                f.write(f"temp_run_dir was: {temp_run_dir}")
        except Exception:
            pass # Failed to write to log
        sys.exit(1)
    finally:
        # Clean up the temporary directory
        if os.path.exists(temp_run_dir):
            shutil.rmtree(temp_run_dir)
            log.info(f"Cleaned up temporary directory: {temp_run_dir}")

if __name__ == "__main__":
    main()