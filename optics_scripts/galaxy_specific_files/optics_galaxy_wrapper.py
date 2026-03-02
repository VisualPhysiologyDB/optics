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
# This wrapper script assumes that 'optics_predictions.py' is in the
# same directory or available in the PYTHONPATH.
try:
    import optics_predictions
except ImportError:
    log.error("Fatal Error: Could not import 'optics_predictions.py'.")
    log.error("Make sure 'optics_predictions.py' is in the same directory as this wrapper.")
    sys.exit(1)
# --- /IMPORTANT ---


def main():
    parser = argparse.ArgumentParser(description="Galaxy wrapper for OPTICS predictions.")

    # --- Input Files ---
    parser.add_argument("--input", required=True, help="Path to input FASTA file.")

    # --- Output Files ---
    # These paths are provided by Galaxy and are where we must move the final files.
    parser.add_argument("--output_tsv", required=True, help="Path for the final output TSV.")
    parser.add_argument("--output_excel", required=True, help="Path for the final output XLSX.")
    parser.add_argument("--output_blastp", required=True, help="Path for the final BLASTp CSV (or 'null').")
    parser.add_argument("--output_bootstrap_viz", required=True, help="Path for the final bootstrap visualization (or 'null').")
    parser.add_argument("--output_removed_seqs", required=True, help="Path for the removed sequences log.")
    parser.add_argument("--output_fig_tree", required=True, help="Path for the FigTree annotation file.")
    parser.add_argument("--output_itol", required=True, help="Path for the iTOL annotation file.")

    # --- Core OPTICS Args ---
    parser.add_argument("--prediction_prefix", default="optics_run", help="Prefix for output files.")
    parser.add_argument("--model_version", default="vpod_1.3", help="Model version.")
    parser.add_argument("--model", default="whole-dataset", help="Prediction model to use.")
    parser.add_argument("--encoding", default="aa_prop", help="Encoding method.")
    parser.add_argument("--tolerate_non_standard_aa", action="store_true", help="Tolerate non-standard AAs.")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of parallel processes.")

    # --- BLASTp Args ---
    parser.add_argument("--blastp", action="store_true", help="Enable BLASTp.")
    parser.add_argument("--refseq", default="bovine", help="Reference sequence for BLASTp.")
    parser.add_argument("--custom_ref_file", default=None, help="Path to custom reference FASTA.")

    # --- Bootstrap Args ---
    parser.add_argument("--bootstrap", action="store_true", help="Enable bootstrap.")
    parser.add_argument("--bootstrap_num", type=int, default=100, help="Number of bootstrap replicates.")
    parser.add_argument("--visualize_bootstrap", action="store_true", help="Enable bootstrap visualization.")
    parser.add_argument("--save_viz_as", default="svg", help="Format for bootstrap visualization.")
    parser.add_argument("--full_spectrum_xaxis", action="store_true", help="Use full spectrum x-axis.")

    args = parser.parse_args()

    # Create a temporary directory for optics_predictions.py to write into.
    # This avoids cluttering the main directory and makes finding outputs easy.
    temp_run_dir = tempfile.mkdtemp()
    log.info(f"Created temporary run directory: {temp_run_dir}")

    # Define the fixed filenames that optics_predictions.py will be instructed to use.
    blastp_report_name = "galaxy_blastp_report"
    bootstrap_viz_name = "galaxy_bootstrap_viz" # The script will add the extension

    try:
        # --- Run the main OPTICS script ---
        # We tell it to write all its output into the temporary directory.
        log.info("Starting optics_predictions.run_optics_predictions...")
        optics_predictions.run_optics_predictions(
            input_sequence=args.input,
            pred_dir=temp_run_dir,  # Tell script to use our temp dir as the *parent*
            output=args.prediction_prefix,
            model=args.model,
            encoding_method=args.encoding,
            blastp=args.blastp,
            iden_report=blastp_report_name,
            refseq=args.refseq,
            reffile=args.custom_ref_file,
            bootstrap=args.bootstrap,
            bootstrap_num=args.bootstrap_num,
            visualize_bootstrap=args.visualize_bootstrap,
            bootstrap_viz_file=bootstrap_viz_name,
            save_as=args.save_viz_as,
            full_spectrum_xaxis=args.full_spectrum_xaxis,
            model_version=args.model_version,
            n_jobs=args.n_jobs,
            tolerate_non_standard_aa=args.tolerate_non_standard_aa
        )
        log.info("optics_predictions.run_optics_predictions finished.")

        # --- Find the output directory ---
        # The script creates a timestamped subdirectory like 'optics_on_...'.
        # We need to find it.
        report_dirs = glob.glob(os.path.join(temp_run_dir, f'optics_on_{args.prediction_prefix}_*'))
        if not report_dirs:
            log.error(f"Fatal Error: Could not find output directory from OPTICS script inside {temp_run_dir}")
            log.error(f"Expected pattern: optics_on_{args.prediction_prefix}_*")
            sys.exit(1)
        
        actual_report_dir = report_dirs[0]
        log.info(f"Found report directory: {actual_report_dir}")

        # --- Define paths to generated files ---
        gen_tsv = os.path.join(actual_report_dir, f"{args.prediction_prefix}_predictions.tsv")
        gen_excel = os.path.join(actual_report_dir, f"{args.prediction_prefix}_predictions_for_excel.xlsx")
        gen_blastp = os.path.join(actual_report_dir, blastp_report_name + '.csv')
        # We handle bootstrap viz separately below, so we just define the base name pattern
        gen_bootstrap_viz_pattern = os.path.join(actual_report_dir, f"{bootstrap_viz_name}_part*") 
        gen_removed = os.path.join(actual_report_dir, "removed_sequences.txt")
        gen_fig_tree = os.path.join(actual_report_dir, "fig_tree_color_annotation.txt")
        gen_itol = os.path.join(actual_report_dir, "itol_color_annotation.txt")

        # --- Move files to final Galaxy paths ---
        log.info("Copying generated files to final Galaxy paths...")
        
        # --- List of output datasets that are NOT optional in the XML ---
        # These MUST exist, even if empty, to prevent a Galaxy error.
        always_expected_outputs = {
            args.output_tsv,
            args.output_excel,
            args.output_removed_seqs,
            args.output_fig_tree,
            args.output_itol
        }

        def copy_file(src, dst):
            if dst and dst != 'null':
                if os.path.exists(src):
                    shutil.copy(src, dst) # <-- Changed from move to copy
                    log.info(f"Copied {src} to {dst}")
                else:
                    # File doesn't exist.
                    if dst in always_expected_outputs:
                        # This file is NOT optional. Galaxy will fail if it's missing.
                        # Create an empty file to satisfy Galaxy.
                        log.warning(f"Warning: Expected output file '{src}' was not found.")
                        log.warning(f"Creating empty file at {dst} to satisfy Galaxy.")
                        with open(dst, 'w') as f:
                            f.write('No sequences were removed, so this file is empty. :)\nGalaxy demands that this file exists for the tool to run, so here we are.')
                    else:
                        # This was an optional file (like blastp) that just wasn't created.
                        # This is fine.
                        log.info(f"Optional file '{src}' was not generated. Skipping copy to {dst}.")
            else:
                log.info(f"Skipping copy for {src} (destination is 'null').")

        copy_file(gen_tsv, args.output_tsv)
        copy_file(gen_excel, args.output_excel)
        copy_file(gen_removed, args.output_removed_seqs)
        copy_file(gen_fig_tree, args.output_fig_tree)
        copy_file(gen_itol, args.output_itol)

        if args.blastp:
            copy_file(gen_blastp, args.output_blastp)
        
        # --- New logic for compiling bootstrap visualizations ---
        if args.bootstrap and args.visualize_bootstrap and args.output_bootstrap_viz != 'null':
            if args.save_viz_as == 'png':
                # This is the ideal case. We find all PNGs and compile them into a single PDF.
                viz_files = sorted(glob.glob(f"{gen_bootstrap_viz_pattern}.png"))
                if not viz_files:
                    log.warning(f"Bootstrap visualization was enabled, but no PNG parts were found with pattern: {gen_bootstrap_viz_pattern}.png")
                else:
                    try:
                        image_list = [Image.open(f).convert('RGB') for f in viz_files]
                        if image_list:
                            image_list[0].save(
                                f'{actual_report_dir}/pdf_placeholder.pdf', 
                                save_all=True, 
                                append_images=image_list[1:],
                                title="OPTICS Bootstrap Visualizations"
                            )
                            copy_file(f'{actual_report_dir}/pdf_placeholder.pdf', args.output_bootstrap_viz)
                            log.info(f"Successfully compiled {len(image_list)} PNGs into single PDF: {args.output_bootstrap_viz}")
                        else:
                             log.warning("Found PNG files but failed to open them with Pillow.")
                             # Fallback to empty file
                             open(args.output_bootstrap_viz, 'a').close()

                    except Exception as e:
                        log.error(f"Error compiling PNGs into PDF: {e}")
                        # Fallback: Just copy the first part to avoid a crash, even though it's not a PDF
                        copy_file(viz_files[0], args.output_bootstrap_viz)
            
            else:
                # For SVG or PDF, merging is not supported by current libraries.
                # We'll log this and just copy the first part as a fallback.
                log.warning(f"Merging visualization parts for format '{args.save_viz_as}' is not supported.")
                log.warning("Only the first part of the visualization will be provided.")
                
                fallback_file = os.path.join(actual_report_dir, f"{bootstrap_viz_name}_part1.{args.save_viz_as}")
                
                if args.save_viz_as == 'pdf':
                    # If format is already PDF, just copy it
                    copy_file(fallback_file, args.output_bootstrap_viz)
                else:
                    # If it's SVG, we can't provide a PDF.
                    # We create an empty file to satisfy Galaxy, since the output MUST be a PDF.
                    log.warning(f"Output format is 'pdf' but input was 'svg'. Creating empty PDF.")
                    with open(args.output_bootstrap_viz, 'w') as f:
                        f.write('Cannot merge SVGs into a PDF with current tool. Please re-run selecting PNG format.')

    except Exception as e:
        log.error(f"An error occurred during the OPTICS run: {e}")
        sys.exit(1)
    finally:
        # Clean up the temporary directory
        if os.path.exists(temp_run_dir):
            shutil.rmtree(temp_run_dir)
            log.info(f"Cleaned up temporary directory: {temp_run_dir}")

if __name__ == "__main__":
    main()