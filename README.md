# uid2TM
This Python script automates the process of comparing protein structures by downloading their PDB files from the AlphaFold protein structure database, computing their radius of gyration (Rg), aligning the structures using TM-align, and extracting key metrics like RMSD, TM-score, and aligned length. The results are saved in a CSV file for each protein pair provided in the input CSV. Based on statistics:
- 0.0 < TM-score < 0.30, random structural similarity              
- 0.5 < TM-score < 1.00, in about the same fold                   

## How to run the script:
> python3 uid2TM.py

## Input:
A CSV file containing columns: prefix, uid1, and uid2.
- prefix: An index for each protein pair.
- uid1 and uid2: UniProt IDs of the two proteins to be aligned.
A working directory path to save the downloaded PDB files and results.

## Output:
A CSV file named <youinputname_TMscore.csv> containing the following columns inaddition to the columns in the input file:
- Rg1: Radius of gyration of protein 1.
- Rg2: Radius of gyration of protein 2.
- RMSD: Root Mean Square Deviation from the TM-align output.
- TM-Score: TM-align score for protein alignment.
- Aligned Length: Number of aligned residues between the two proteins.

## Requirements:
-TMalign
The uid2TM.py script uses TMalign which is a python module that provides wrappers to TMalign, TMscore and MMalign.
The executables can be downloaded from http://zhanglab.ccmb.med.umich.edu/TM-align/ and should be saved to any directory in PATH. Also see https://pymolwiki.org/index.php/TMalign
Click [here](https://zhanggroup.org/TM-align/TMalign.cpp) to downlowd TMalign 
Install TM-align as follows:
> g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp 

Find where TM-align is installed on your system and make sure the executable is available in your PATH in the bashrc file:
> vi ~/.bashrc

Add the following line to the file, replacing /path/to/TMalign with the actual path to the TM-align executable
> export TMalign="/path/to/TMalign"
Save the file and exit the editor

Reload the shell configuration
> source ~/.bashrc

Verify the installtion as follows
> TMalign -h

-BioPython:
> pip3 install biopython

-Pandas:
> pip3 install pandas

-Requests:
> pip3 install requests

-Numpy
> pip3 install numpy

MDAnalysis
>pip3 install MDAnalysis --only-binary :all:
