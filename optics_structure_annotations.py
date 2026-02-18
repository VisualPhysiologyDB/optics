import os
import sys
import argparse
import pandas as pd
import urllib.request
import pathlib

# --- Helper: Fetch PDB ---
def fetch_pdb_if_needed(pdb_input, save_dir):
    """Retrieves PDB from RCSB if input is a 4-char code, otherwise verifies local path."""
    if os.path.isfile(pdb_input):
        return pdb_input
    
    if len(pdb_input) == 4 and pdb_input.isalnum():
        pdb_id = pdb_input.upper()
        dest = os.path.join(save_dir, f"{pdb_id}.pdb")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Downloading PDB structure {pdb_id}...")
        try:
            urllib.request.urlretrieve(url, dest)
            return dest
        except Exception as e:
            print(f"Error downloading PDB: {e}")
            return None
    return None

def run_structure_annotation(annotation_file, pdb_input="1U19", output_dir=None, chain="A", software="pymol"):
    """
    Generates a script to visualize arbitrary protein annotations in PyMOL or ChimeraX.
    
    Args:
        annotation_file (str): Path to CSV/TSV with columns: 'position', 'color', 'style'.
        pdb_input (str): Path to PDB or PDB ID (default 1U19).
        output_dir (str): Where to save the output script.
        chain (str): Chain ID to apply annotations to.
        software (str): 'pymol' or 'chimerax'.
    """
    software = software.lower()
    if software not in ['pymol', 'chimerax']:
        print(f"Error: Unsupported software '{software}'. Choose 'pymol' or 'chimerax'.")
        return

    if output_dir is None:
        output_dir = os.path.dirname(annotation_file)
    if not os.path.exists(output_dir) and output_dir != '':
        os.makedirs(output_dir, exist_ok=True)

    # 1. Load Annotation Data
    try:
        if annotation_file.endswith('.tsv') or annotation_file.endswith('.txt'):
            df = pd.read_csv(annotation_file, sep='\t')
        else:
            df = pd.read_csv(annotation_file)
        
        df.columns = [c.lower().strip() for c in df.columns]
        if 'position' not in df.columns:
            raise ValueError("Input file must have a 'position' column.")
    except Exception as e:
        print(f"Error reading annotation file: {e}")
        return

    # 2. Get PDB
    pdb_path = fetch_pdb_if_needed(pdb_input, output_dir)
    if not pdb_path:
        print("Could not locate or download PDB file.")
        return
    
    pdb_name = os.path.basename(pdb_path)
    script_base_name = os.path.splitext(os.path.basename(annotation_file))[0]

    # --- Generate Script based on Software ---
    
    if software == 'pymol':
        # --- PYMOL GENERATION ---
        ext = ".pml"
        out_path = os.path.join(output_dir, f"visualize_{script_base_name}_pymol{ext}")
        
        print(f"Generating PyMOL script for {len(df)} annotations...")
        with open(out_path, 'w') as f:
            f.write(f"# OPTICS Structure Annotation Script (PyMOL)\n")
            f.write(f"load {pdb_name}\n")
            f.write("hide everything\n")
            f.write("bg_color black\n")
            
            # Base Cartoon
            f.write(f"show cartoon, chain {chain}\n")
            f.write(f"color gray50, chain {chain}\n") 
            f.write(f"set cartoon_transparency, 0.2, chain {chain}\n")

            f.write("\n# --- Annotations ---\n")
            for idx, row in df.iterrows():
                pos = str(row['position'])
                sel_name = f"site_{pos}"
                f.write(f"select {sel_name}, chain {chain} and resi {pos}\n")
                
                # Color Handling (PyMOL accepts 0x for hex)
                color = row.get('color', 'red')
                if str(color).startswith('#'):
                    color = "0x" + str(color).lstrip('#')
                f.write(f"color {color}, {sel_name}\n")
                
                # Style Handling
                style = row.get('style', 'sphere').lower()
                if 'sphere' in style or 'dot' in style:
                    f.write(f"show spheres, {sel_name}\n")
                    f.write(f"set sphere_scale, 1.2, {sel_name}\n")
                elif 'stick' in style:
                    f.write(f"show sticks, {sel_name}\n")
                elif 'cartoon' in style:
                    f.write(f"set cartoon_transparency, 0.0, {sel_name}\n")
                
                if 'label' in row and pd.notna(row['label']):
                    f.write(f"label {sel_name} and n. ca, '{str(row['label'])}'\n")

            f.write("\ndeselect\norient\n")

    elif software == 'chimerax':
        # --- CHIMERAX GENERATION ---
        ext = ".cxc"
        out_path = os.path.join(output_dir, f"visualize_{script_base_name}_chimerax{ext}")
        
        print(f"Generating ChimeraX script for {len(df)} annotations...")
        with open(out_path, 'w') as f:
            f.write(f"# OPTICS Structure Annotation Script (ChimeraX)\n")
            f.write(f"open {pdb_name}\n")
            f.write("hide atoms\n")
            f.write("hide cartoons\n")
            f.write("set bgColor black\n")
            
            # Base Cartoon
            # Note: ChimeraX selector syntax is /Chain:Residue
            f.write(f"show /{chain} cartoons\n")
            f.write(f"color /{chain} gray\n")
            f.write(f"transparency /{chain} 20 target c\n") # 20% transparent

            f.write("\n# --- Annotations ---\n")
            for idx, row in df.iterrows():
                pos = str(row['position'])
                sel = f"/{chain}:{pos}"
                
                # Color Handling (ChimeraX accepts # for hex)
                color = str(row.get('color', 'red'))
                # ChimeraX handles #RRGGBB directly, no change needed usually
                f.write(f"color {sel} {color}\n")
                
                # Style Handling
                style = row.get('style', 'sphere').lower()
                if 'sphere' in style or 'dot' in style:
                    f.write(f"show {sel} atoms\n")
                    f.write(f"style {sel} sphere\n")
                elif 'stick' in style:
                    f.write(f"show {sel} atoms\n")
                    f.write(f"style {sel} stick\n")
                elif 'cartoon' in style:
                    # Make opaque
                    f.write(f"transparency {sel} 0 target c\n")
                
                if 'label' in row and pd.notna(row['label']):
                    f.write(f"label {sel} text \"{str(row['label'])}\"\n")
            
            f.write("\nview\n")

    print(f"Done! Script saved to: {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize arbitrary annotations on a protein structure.")
    parser.add_argument("-a", "--annotation_file", required=True, help="CSV/TSV file.")
    parser.add_argument("-p", "--pdb", default="1U19", help="PDB ID or path.")
    parser.add_argument("-o", "--output_dir", default=None, help="Output directory.")
    parser.add_argument("--chain", default="A", help="Chain identifier.")
    parser.add_argument("--software", default="chimerax", choices=['pymol', 'chimerax'], help="Target software.")
    
    args = parser.parse_args()
    
    run_structure_annotation(args.annotation_file, args.pdb, args.output_dir, args.chain, args.software)