import os
import pandas as pd
import requests
import subprocess
import numpy as np
from Bio.PDB import PDBParser

# Ask user for input CSV file path and working directory
csv_path = input("Enter the path to the CSV file (containing prefix, uid1, uid2): ")
work_dir = input("Enter the path to the working directory where results will be saved: ")

# Load CSV input
df = pd.read_csv(csv_path)
df['Rg1'] = np.nan
df['Rg2'] = np.nan
df['RMSD'] = np.nan
df['TM-Score1'] = np.nan
df['TM-Score2'] = np.nan
df['Aligned Length'] = np.nan
df['Seq_ID'] = np.nan

# Create working directory if it doesn't exist
if not os.path.exists(work_dir):
    os.makedirs(work_dir)

# AlphaFold base URL for downloading PDB files
alphafold_url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v1.pdb"

# Function to fetch protein structure
def fetch_protein_structure(uniprot_id, output_path):
    url = alphafold_url.format(uniprot_id)
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            f.write(response.content)
    else:
        raise ValueError(f"Error fetching PDB for {uniprot_id}, status code: {response.status_code}")

# Function to calculate radius of gyration using Bio.PDB
def calculate_radius_of_gyration(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)

    # Extract all atom coordinates
    atom_coords = np.array([atom.coord for atom in structure.get_atoms()])

    if len(atom_coords) == 0:
        raise ValueError(f"No atoms found in structure {pdb_path}")

    # Calculate center of mass
    center_of_mass = np.mean(atom_coords, axis=0)

    # Calculate radius of gyration
    rg = np.sqrt(np.mean(np.sum((atom_coords - center_of_mass)**2, axis=1)))
    return rg

# Function to align structures using TM-align and save output to TMalign.out
def align_with_tmalign(pdb1_path, pdb2_path, output_dir):
    # Path to TMalign.out file
    tmalign_output_path = os.path.join(output_dir, "TMalign.out")
    # Run TM-align with output redirected to TMalign.out
    cmd = f"/media/hadeer/Elements/seq2TM/TMalign {pdb1_path} {pdb2_path} > {tmalign_output_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise Exception(f"TM-align failed: {result.stderr}")

    return tmalign_output_path

# Function to parse TM-align output file
def parse_tmalign_output(tmalign_output_path):
    aligned_length = rmsd = seq_id = tm_score1 = tm_score2 = None
    with open(tmalign_output_path, 'r') as file:
        for line in file:
            if line.startswith("Aligned length="):
                parts = line.split(',')
                # Extract values after the last '=' for each part
                aligned_length = int(parts[0].split('=')[-1].strip())
                rmsd = float(parts[1].split('=')[-1].strip())
                seq_id = float(parts[2].split('=')[-1].strip())
            elif line.startswith("TM-score="):
                # Extract only the first numeric value after '='
                score_value = line.split('=')[1].split()[0]  # First part after '=' before any space
                if "normalized by length of Chain_1" in line:
                    tm_score1 = float(score_value)
                elif "normalized by length of Chain_2" in line:
                    tm_score2 = float(score_value)

    # Ensure all values are found, else raise an exception
    if None in (aligned_length, rmsd, seq_id, tm_score1, tm_score2):
        raise Exception("Failed to find relevant information in TM-align output.")
    
    return aligned_length, rmsd, seq_id, tm_score1, tm_score2



# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    prefix = row['prefix']
    uid1 = row['uid1']
    uid2 = row['uid2']

    # Create a directory for the current protein pair
    pair_dir = os.path.join(work_dir, prefix)
    os.makedirs(pair_dir, exist_ok=True)

    # Paths for PDB files
    pdb1_path = os.path.join(pair_dir, f"{uid1}.pdb")
    pdb2_path = os.path.join(pair_dir, f"{uid2}.pdb")

    # Fetch the protein structures from AlphaFold
    try:
        fetch_protein_structure(uid1, pdb1_path)
        fetch_protein_structure(uid2, pdb2_path)
    except ValueError as e:
        print(e)
        continue

    # Calculate radius of gyration for both proteins
    try:
        Rg1 = calculate_radius_of_gyration(pdb1_path)
        Rg2 = calculate_radius_of_gyration(pdb2_path)
        df.at[index, 'Rg1'] = Rg1
        df.at[index, 'Rg2'] = Rg2
    except Exception as e:
        print(f"Error calculating radius of gyration for {uid1} or {uid2}: {e}")
        continue

    # Align protein structures and calculate TM-score and RMSD
    try:
        tmalign_output_file = align_with_tmalign(pdb1_path, pdb2_path, pair_dir)
        aligned_length, rmsd, seq_id, tm_score1, tm_score2 = parse_tmalign_output(tmalign_output_file)
        
        # Update DataFrame with parsed values
        df.at[index, 'RMSD'] = rmsd
        df.at[index, 'TM-Score1'] = tm_score1
        df.at[index, 'TM-Score2'] = tm_score2
        df.at[index, 'Aligned Length'] = aligned_length
        df.at[index, 'Seq_ID'] = seq_id
    except Exception as e:
        print(f"Error aligning {uid1} and {uid2}: {e}")
        continue

# Save the updated DataFrame to a new CSV file with "_TMscore" appended to the original filename
output_csv_path = os.path.join(work_dir, os.path.basename(csv_path).replace('.csv', '_TMscore.csv'))
df.to_csv(output_csv_path, index=False)

print(f"Results saved to {output_csv_path}")

