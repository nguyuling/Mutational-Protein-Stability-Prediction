import os
import io
import time
import pandas as pd
import requests
import torch
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBList

# Load your s669 dataset
df = pd.read_csv("S669.csv")

# Initialize tools
pdbl = PDBList()
cif_parser = MMCIFParser(QUIET=True)

# Create absolute folder paths
os.makedirs("structures", exist_ok=True)
os.makedirs("processed_pt_coords", exist_ok=True)

def fetch_wildtype_sequence(pdb_id):
    """Fetches full FASTA sequence from RCSB PDB."""
    pdb_id = pdb_id.lower()[:4]
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

def extract_heavy_atom_coordinates(pdb_id, chain_id='A'):
    """Downloads mmCIF and safely maps heavy atom coordinates to numerical matrix arrays."""
    try:
        pdb_id = pdb_id.lower()
        # Retrieve the file safely 
        pdbl.retrieve_pdb_file(pdb_id, pdir="structures", file_format="mmCif")
        
        # Resolve target structure folders
        expected_path_1 = os.path.join("structures", f"{pdb_id}.cif")
        expected_path_2 = os.path.join("structures", pdb_id[1:3], f"{pdb_id}.cif")
        
        if os.path.exists(expected_path_1):
            file_path = expected_path_1
        elif os.path.exists(expected_path_2):
            file_path = expected_path_2
        else:
            print(f"   ⚠️ Could not find downloaded file for {pdb_id} in local directories.")
            return None
            
        structure = cif_parser.get_structure(pdb_id, file_path)
        model = structure[0]
        
        # Ensure the chain ID is standardized
        chain_id = chain_id.upper()
        if chain_id not in model:
            available_chains = [c.get_id() for c in model.get_chains()]
            chain_id = available_chains[0] # Fallback to first chain
            
        chain = model[chain_id]
        all_residue_coords = []
        
        for residue in chain:
            if residue.id[0] == ' ':  # Standard residues only
                residue_atoms = []
                for atom in residue:
                    # FIX: Biopython uses atom.element (property) or checking atom name prefix
                    # Strip away Hydrogen ('H') and Deuterium ('D') atoms
                    atom_name = atom.get_name().strip()
                    if not (atom_name.startswith('H') or atom_name.startswith('D')):
                        residue_atoms.append(atom.get_coord())
                if residue_atoms:
                    all_residue_coords.append(residue_atoms)
                    
        return all_residue_coords
    except Exception as e:
        print(f"   ❌ Error processing {pdb_id}: {e}")
        return None

# --- EXECUTE PIPELINE ---

print("Step 1: Fetching wild-type sequences from RCSB...")
unique_pdbs = df['Protein'].unique()
pdb_seq_map = {pdb: fetch_wildtype_sequence(pdb) for pdb in unique_pdbs}

# Standardize data framing columns explicitly
df['Wt_seq'] = df['Protein'].map(pdb_seq_map)
df['Mut'] = df['Mut_seq']  
df['DDG'] = df['DDG_checked_dir']

# Filter out rows that failed sequence retrieval
df = df.dropna(subset=['Wt_seq']).copy()

print(f"\nStep 2: Processing {len(unique_pdbs)} unique structures...")
coordinate_registry = {}

for pdb in unique_pdbs:
    pdb_code = pdb[:4]
    chain_code = pdb[4] if len(pdb) > 4 else 'A'
    
    print(f" -> Downloading & parsing: {pdb_code} (Chain {chain_code})")
    coords_list = extract_heavy_atom_coordinates(pdb_code, chain_code)
    
    if coords_list and len(coords_list) > 0:
        residue_tensors = [torch.tensor(res, dtype=torch.float32) for res in coords_list]
        pt_filename = f"processed_pt_coords/{pdb}_coords.pt"
        torch.save(residue_tensors, pt_filename)
        coordinate_registry[pdb] = pt_filename
        print(f"   ✅ Saved coordinates to: {pt_filename}")
    else:
        coordinate_registry[pdb] = None

# Apply the coordinate paths back to your spreadsheet
df['Coord_File_Path'] = df['Protein'].map(coordinate_registry)

# Filter out rows where the structural coordinate generation failed
df_final = df.dropna(subset=['Coord_File_Path']).copy()

# Enforce your requested ordered layout
final_columns = ['Protein', 'Mut', 'Wt_seq', 'Mut_seq', 'DDG', 'Coord_File_Path']
df_output = df_final[final_columns]

# Save spreadsheet
df_output.to_csv("s669_pipeline_summary.csv", index=False)
print(f"\nSuccess! Processed {len(df_output)} rows. Check 's669_pipeline_summary.csv'")