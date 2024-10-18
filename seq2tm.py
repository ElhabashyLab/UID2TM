import pandas as pd
import requests
from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer

# Function to fetch protein sequence from UniProt
def fetch_sequence(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text.strip().split('\n')[1:]  # Skip the header line
        return ''.join(fasta_data)
    else:
        print(f"Error fetching sequence for {uniprot_id}: {response.status_code}")
        return None

# Function to fetch protein structure from AlphaFold
def fetch_structure(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v1.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error fetching structure for {uniprot_id}: {response.status_code}")
        return None

# Function to calculate TM-score
def calculate_tm_score(pdb1, pdb2):
    parser = PDBParser()
    structure1 = parser.get_structure("structure1", pdb1)
    structure2 = parser.get_structure("structure2", pdb2)

    # Assuming the first model is used for comparison
    model1 = structure1[0]
    model2 = structure2[0]

    # Align the structures and calculate TM-score
    superimposer = Superimposer()
    superimposer.set_atoms(list(model1.get_atoms()), list(model2.get_atoms()))
    tm_score = superimposer.get_rms()
    return tm_score

# Main function
def main(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)

    results = []

    for index, row in df.iterrows():
        uid1, uid2 = row['uid1'], row['uid2']
        
        sequence1 = fetch_sequence(uid1)
        sequence2 = fetch_sequence(uid2)

        structure1 = fetch_structure(uid1)
        structure2 = fetch_structure(uid2)

        if structure1 and structure2:
            tm_score = calculate_tm_score(structure1, structure2)
            results.append({'uid1': uid1, 'uid2': uid2, 'TM_score': tm_score})
        else:
            results.append({'uid1': uid1, 'uid2': uid2, 'TM_score': None})

    # Save results to a new CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv('tm_scores.csv', index=False)

if __name__ == "__main__":
    main('your_file.csv')  # Replace with your CSV file path
