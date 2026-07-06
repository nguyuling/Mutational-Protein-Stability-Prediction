import pandas as pd
import requests
from Bio import SeqIO
import io
import time

# 1. Load the s669 dataset
df = pd.read_csv("S669.csv")

# 2. Function to fetch the absolute Wild-Type sequence using the PDB ID
def fetch_wildtype_sequence(pdb_id):
    pdb_id = pdb_id.lower()[:4]  # first 4 char is the pdb id
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            fasta_io = io.StringIO(response.text)
            for record in SeqIO.parse(fasta_io, "fasta"):
                return str(record.seq).upper()
        return None
    except Exception:
        return None

# 3. Function to modify the Wild-Type sequence to create the Mutant sequence
def generate_mutant_sequence(row):
    wild_seq = row['Wt_seq']
    mutation = str(row['Mut']) # e.g., "S11A"
    
    if not wild_seq or pd.isna(wild_seq):
        return None
        
    try:
        wt_aa = mutation[0].upper()     # Original Amino Acid (e.g., 'S')
        mut_aa = mutation[-1].upper()   # Target Mutant Amino Acid (e.g., 'A')
        
        # Convert 1-based PDB position index into 0-based Python string index
        pos = int(mutation[1:-1]) - 1 
        
        # Verification Check: Ensure the WT residue actually matches the position
        if pos < len(wild_seq) and wild_seq[pos] == wt_aa:
            seq_list = list(wild_seq)
            seq_list[pos] = mut_aa
            return "".join(seq_list)
        else:
            return None # Skip due to offset or indexing mismatch
    except Exception:
        return None

# --- EXECUTE THE PIPELINE ---

print("Step 1: Fetching original sequences from RCSB PDB...")
# Map unique PDB codes first to minimize API calls
unique_pdbs = df['Protein'].unique()
pdb_seq_map = {}

for pdb in unique_pdbs:
    pdb_seq_map[pdb] = fetch_wildtype_sequence(pdb)
    time.sleep(0.1) # Soft polite delay for the server

# Create and map initial required structural columns
df['Mut'] = df['Mut_seq']
df['Wt_seq'] = df['Protein'].map(pdb_seq_map)

print("Step 2: Deriving mutated sequences...")
df['Mut_seq_derived'] = df.apply(generate_mutant_sequence, axis=1)

# Rename experimental ground-truth column to match target specification
df['DDG'] = df['DDG_checked_dir']

# 4. Final Dataset Structuring & Cleaning
# Drop rows where the sequence could not be fetched or mutation mapping failed
df_cleaned = df.dropna(subset=['Wt_seq', 'Mut_seq_derived'])

# Select and order exactly according to your layout request
final_columns = ['Protein', 'Mut', 'Wt_seq', 'Mut_seq', 'DDG']
# Replace old text tracking with our validated derived sequence column
df_cleaned = df_cleaned.rename(columns={'Mut_seq_derived': 'Mut_seq'})
df_final = df_cleaned[final_columns]

# Save output to a clean CSV
df_final.to_csv("S669_cleaned.csv", index=False)
print("Pipeline complete. Processed output saved to 's669_ordered_cleaned.csv'")

# View the final formatted layout matrix matching your target example
print(df_final.head(1))