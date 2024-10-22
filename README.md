# uid2TM
This Python script automates the process of comparing protein structures by downloading their PDB files from the AlphaFold protein structure database, computing their radius of gyration (Rg), aligning the structures using TM-align, and extracting key metrics like RMSD, TM-score, and aligned length. The results are saved in a CSV file for each protein pair provided in the input CSV.

# How to run the script:
python3 uid2TM.py

# Input:
A CSV file containing columns: prefix, uid1, and uid2.
prefix: An index for each protein pair.
uid1 and uid2: UniProt IDs of the two proteins to be aligned.
A working directory path to save the downloaded PDB files and results.

# Output:
A CSV file containing the following columns inaddition to the columns in the input file:
Rg1: Radius of gyration of protein 1.
Rg2: Radius of gyration of protein 2.
RMSD: Root Mean Square Deviation from the TM-align output.
TM-Score: TM-align score for protein alignment.
Aligned Length: Number of aligned residues between the two proteins.

# Library Requirements:
BioPython:
pip3 install biopython
Pandas:
pip3 install pandas
Requests:
pip3 install requests
