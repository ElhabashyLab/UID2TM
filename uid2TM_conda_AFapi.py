import os
import pandas as pd
import requests
import subprocess
import numpy as np
from Bio.PDB import PDBParser
import io

# --- CONFIGURATION ---
# Ensure these are installed: pip install pandas requests biopython numpy
# TMalign must be in your system PATH or provide full path below
TMALIGN_CMD = "TMalign" 
# ---------------------

def fetch_pdb_url_via_api(uniprot_id):
    """
    Queries the AlphaFold API to get the correct PDB URL 
    (Handles v3/v4 versions automatically).
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        r = requests.get(api_url)
        if r.status_code == 200:
            data = r.json()
            if data and isinstance(data, list):
                # Return the URL for the PDB format
                return data[0]['pdbUrl']
    except Exception:
        pass
    return None

def fetch_protein_structure(uniprot_id, output_path):
    # 1. Get the dynamic URL from API
    pdb_url = fetch_pdb_url_via_api(uniprot_id)
    
    if not pdb_url:
        # Fallback: Try constructing the URL manually if API fails (rare)
        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

    # 2. Download the file
    response = requests.get(pdb_url)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded: {uniprot_id}")
    else:
        raise ValueError(f"Error 404: Could not find PDB for {uniprot_id}. Check ID.")

def calculate_radius_of_gyration(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_path)
        atom_coords = np.array([atom.coord for atom in structure.get_atoms() if atom.name == 'CA']) # Use CA for speed, or remove if atom for all
        
        if len(atom_coords) == 0:
            # Fallback to all atoms if CA missing
            atom_coords = np.array([atom.coord for atom in structure.get_atoms()])

        center_of_mass = np.mean(atom_coords, axis=0)
        rg = np.sqrt(np.mean(np.sum((atom_coords - center_of_mass)**2, axis=1)))
        return rg
    except Exception as e:
        print(f"Rg Error: {e}")
        return np.nan

def run_tmalign(pdb1, pdb2, output_prefix):
    """
    Runs TMalign and parses STDOUT directly (more robust than reading the PDB header).
    """
    cmd = f"{TMALIGN_CMD} {pdb1} {pdb2} -o {output_prefix}"
    
    # Run command
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise Exception(f"TM-align failed: {result.stderr}")

    # Parse STDOUT for metrics
    # Looking for lines like:
    # "Aligned length=  186, RMSD=   2.45, Seq_ID=n_identical/n_aligned= 0.151"
    # "TM-score= 0.65432 (if normalized by length of Chain_1)"
    
    aligned_len = np.nan
    rmsd = np.nan
    tm_score = np.nan
    
    lines = result.stdout.split('\n')
    for line in lines:
        if line.startswith("Aligned length="):
            parts = line.split(',')
            for part in parts:
                if "Aligned length=" in part:
                    aligned_len = int(part.split('=')[1].strip())
                elif "RMSD=" in part:
                    rmsd = float(part.split('=')[1].strip())
        
        # We usually prioritize the TM-score normalized by the FIRST chain (Chain_1)
        elif line.startswith("TM-score=") and "Chain_1" in line:
            tm_score = float(line.split('=')[1].split('(')[0].strip())

    return aligned_len, rmsd, tm_score

# --- MAIN ---

csv_path = input("Enter the path to the CSV file: ").strip()
work_dir = input("Enter the path to the working directory: ").strip()

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

# 1. Load CSV and STRIP WHITESPACE (Crucial Fix)
df = pd.read_csv(csv_path)
df['uid1'] = df['uid1'].astype(str).str.strip()
df['uid2'] = df['uid2'].astype(str).str.strip()
df['prefix'] = df['prefix'].astype(str).str.strip()

# Initialize columns
df['Rg1'] = np.nan
df['Rg2'] = np.nan
df['RMSD'] = np.nan
df['TM-Score'] = np.nan
df['Aligned_Length'] = np.nan

for index, row in df.iterrows():
    prefix = row['prefix']
    uid1 = row['uid1']
    uid2 = row['uid2']

    print(f"Processing {prefix}: {uid1} vs {uid2}")

    pair_dir = os.path.join(work_dir, prefix)
    os.makedirs(pair_dir, exist_ok=True)

    pdb1_path = os.path.join(pair_dir, f"{uid1}.pdb")
    pdb2_path = os.path.join(pair_dir, f"{uid2}.pdb")
    aligned_out_prefix = os.path.join(pair_dir, f"{prefix}_aligned")

    # A. Download
    try:
        if not os.path.exists(pdb1_path): fetch_protein_structure(uid1, pdb1_path)
        if not os.path.exists(pdb2_path): fetch_protein_structure(uid2, pdb2_path)
    except ValueError as e:
        print(f"  Skipping {prefix}: {e}")
        continue

    # B. Radius of Gyration
    df.at[index, 'Rg1'] = calculate_radius_of_gyration(pdb1_path)
    df.at[index, 'Rg2'] = calculate_radius_of_gyration(pdb2_path)

    # C. Alignment
    try:
        aln_len, rmsd, tm = run_tmalign(pdb1_path, pdb2_path, aligned_out_prefix)
        df.at[index, 'RMSD'] = rmsd
        df.at[index, 'TM-Score'] = tm
        df.at[index, 'Aligned_Length'] = aln_len
    except Exception as e:
        print(f"  Alignment failed: {e}")

# Save
output_csv_path = os.path.join(work_dir, os.path.basename(csv_path).replace('.csv', '_TMscore.csv'))
df.to_csv(output_csv_path, index=False)
print(f"Done! Saved to {output_csv_path}")