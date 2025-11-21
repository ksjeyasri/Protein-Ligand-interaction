from Bio.PDB import PDBParser, PDBList
import numpy as np
from collections import defaultdict
import pandas as pd
import os

# === Step 1: Input PDB ID and download ===
pdb_id = input("Enter PDB ID (e.g., 3VMT): ").strip().upper()
pdb_dir = "pdb_files"
os.makedirs(pdb_dir, exist_ok=True)

pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format="pdb")
print(f"✅ PDB file downloaded: {pdb_file}")

# === Step 2: Load structure ===
parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdb_file)

# === Step 3: Interaction parameters ===
covalent_radii = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07}
bond_tolerance = 0.45
nonbond_cutoff = 5.0  # Å

# === Step 4: Separate protein and ligand atoms ===
protein_atoms = []
ligand_atoms = []

for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] == " ":
                # Standard amino acid → protein
                protein_atoms.extend(list(residue.get_atoms()))
            elif residue.id[0].startswith("H_") or residue.id[0] == "HETATM":
                # Ligand (exclude water)
                if residue.get_resname() != "HOH":
                    ligand_atoms.extend(list(residue.get_atoms()))

print(f"Protein atoms: {len(protein_atoms)}, Ligand atoms: {len(ligand_atoms)}")

# === Step 5: Interaction counters ===
covalent_count = defaultdict(int)
hydrogen_count = defaultdict(int)
vdw_count = defaultdict(int)

# === Step 6: Calculate protein-ligand interactions ===
for atom1 in protein_atoms:
    for atom2 in ligand_atoms:
        elem1, elem2 = atom1.element, atom2.element
        if not elem1 or not elem2:
            continue

        dist = np.linalg.norm(atom1.coord - atom2.coord)
        if dist < 0.4:  # skip unrealistic
            continue

        bond_type = "-".join(sorted([elem1, elem2]))

        # Covalent bonds
        if elem1 in covalent_radii and elem2 in covalent_radii:
            max_dist = covalent_radii[elem1] + covalent_radii[elem2] + bond_tolerance
            if dist <= max_dist:
                covalent_count[bond_type] += 1
                continue

        # Nonbonded interactions
        if dist <= nonbond_cutoff:
            if dist <= 3.5 and (
                (elem1 == "N" and elem2 == "O")
                or (elem1 == "O" and elem2 == "N")
                or (elem1 == "N" and elem2 == "N")
            ):
                hydrogen_count[bond_type] += 1
            else:
                vdw_count[bond_type] += 1
# === Step 7: Print results ===
print("\n=== Bonded Interactions ===")
print("Covalent bonds:", sum(covalent_count.values()))
for k, v in covalent_count.items():
    print(f"{k}: {v} interactions")

print("\n=== Non-Bonded Interactions ===")
print("Hydrogen bonds:", sum(hydrogen_count.values()))
for k, v in hydrogen_count.items():
    print(f"{k}: {v} interactions")

print("Van der Waals / Nonbonded:", sum(vdw_count.values()))
for k, v in vdw_count.items():
    print(f"{k}: {v} interactions")

# === Step 8: Save to Excel ===
summary_data = []

# Covalent
summary_data.append(["Covalent bonds", "TOTAL", sum(covalent_count.values())])
for k, v in covalent_count.items():
    summary_data.append(["Covalent bonds", k, v])

# Hydrogen
summary_data.append(["Hydrogen bonds", "TOTAL", sum(hydrogen_count.values())])
for k, v in hydrogen_count.items():
    summary_data.append(["Hydrogen bonds", k, v])

# Van der Waals
summary_data.append(["Van der Waals / Nonbonded", "TOTAL", sum(vdw_count.values())])
for k, v in vdw_count.items():
    summary_data.append(["Van der Waals / Nonbonded", k, v])

df = pd.DataFrame(summary_data, columns=["Interaction Type", "Bond", "Count"])
output_file = f"{pdb_id}_interaction_summary.xlsx"
df.to_excel(output_file, index=False)
print(f"\n✅ Summary saved to {output_file}")
