import numpy as np
from Bio import PDB
import matplotlib.pyplot as plt

def calculate_distance(atom1, atom2):

  # Oblicza odległość euklidesową między dwoma atomami CA.

    diff_vector = atom1.get_vector() - atom2.get_vector()
    distance = np.linalg.norm(diff_vector)
    return distance

def contact_map(structure, threshold=8.0):

   # Generuje mapę kontaktów dla danej struktury białka.

    ca_atoms = [atom for atom in structure.get_atoms() if atom.get_id() == 'CA']
    num_atoms = len(ca_atoms)
    contact_matrix = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = calculate_distance(ca_atoms[i], ca_atoms[j])
            if distance <= threshold:
                contact_matrix[i, j] = 1
                contact_matrix[j, i] = 1

    return contact_matrix

def plot_contact_map(contact_matrix):

   # Wyświetla scatter plot mapy kontaktów.

    plt.imshow(contact_matrix, cmap='viridis', interpolation='none', origin='upper')
    plt.colorbar(label='Contact')
    plt.title('Contact Map')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()

if __name__ == "__main__":
    # Załaduj strukturę białka
    pdb_id = "102m"
    pdb_filename = f"{pdb_id}.pdb"
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(pdb_id, pdb_filename)

    # Wyznacz mapę kontaktów
    contact_matrix = contact_map(structure)

    # Wizualizuj mapę kontaktów
    plot_contact_map(contact_matrix)
