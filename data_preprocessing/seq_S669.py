import io
import pandas as pd
import requests
from Bio import SeqIO

#! load dataset
df = pd.read_csv("raw_data/S669.csv")
df["pdb_id"] = df["Protein"].str[:4]
df["chain"] = df["Protein"].str[4:]

# Fetch unique wild-type sequences with global tracking
pdb_seq_map = {}
successful_wt_count = 0

for _, row in df[["pdb_id", "chain"]].drop_duplicates().iterrows():
    pdb_id, chain_id = row["pdb_id"], row["chain"]
    try:
        response = requests.get(
            f"https://www.rcsb.org/fasta/entry/{pdb_id}", timeout=10
        )
        if response.status_code == 200:
            # Parse FASTA records to find exact matching chain
            chains = {
                rec.id.split(":")[-1].strip().upper(): str(rec.seq).upper()
                for rec in SeqIO.parse(io.StringIO(response.text), "fasta")
            }
            if chain_id in chains:
                pdb_seq_map[(pdb_id, chain_id)] = chains[chain_id]
                successful_wt_count += 1
                continue
    except Exception:
        pass

    # If logic reaches here, fetching failed for this protein combo
    print(
        f"Failed to fetch wild type sequence {pdb_id} (supposed chain {chain_id})"
    )

print(f"Successfully fetch {successful_wt_count} wild-type sequences")

# Map sequences back and drop missing wildtypes
df["wildtype_seq"] = df.set_index(["pdb_id", "chain"]).index.map(pdb_seq_map)
df = df.dropna(subset=["wildtype_seq"]).copy()

# Derive mutant sequences
mutant_seqs = []
successful_mut_count = 0

for _, row in df.iterrows():
    wt_seq = row["wildtype_seq"]
    mut = str(row["PDB_Mut"])
    try:
        wt_aa, mut_aa, pos = mut[0].upper(), mut[-1].upper(), int(mut[1:-1]) - 1
        if 0 <= pos < len(wt_seq) and wt_seq[pos] == wt_aa:
            mut_seq = wt_seq[:pos] + mut_aa + wt_seq[pos + 1 :]
            mutant_seqs.append(mut_seq)
            successful_mut_count += 1
            continue
    except Exception:
        pass

    print(f"Failed to derive mutant sequence {row['Protein']} with mutation {mut}")
    mutant_seqs.append(None)

print(f"Successfully derived {successful_mut_count} mutant sequences")

# Assign mutant sequences and drop failures
df["mutant_seq"] = mutant_seqs
df = df.dropna(subset=["mutant_seq"])

# Select, rename and save final columns according to the schema
final_df = pd.DataFrame(
    {
        "pdb_id": df["pdb_id"],
        "chain": df["chain"],
        "wildtype_seq": df["wildtype_seq"],
        "mutant_seq": df["mutant_seq"],
        "temperature": df["TEMP"],
        "pH": df["pH"],
        "ddG_dir": df["DDG_checked_dir"],
        "ddG_inv": df["DDG_checked_inv"],
    }
)

final_df.to_csv("seq_S669.csv", index=False)