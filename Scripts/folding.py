import AminoAcidClass
from fasta_parser import fasta_read
import argparse
import math
import random
import matplotlib.pyplot as plt
from tqdm import tqdm

pa = argparse.ArgumentParser(description=("Program to fold an HP protein by minimizing \
its energy"))
pa.usage = "folding.py fasta_input iterations init_method -s"
pa.add_argument("fasta_input", type=str, help="Name of the starting protein")
pa.add_argument("iterations", type=int, help="Number of tries to move an amino acid")
pa.add_argument("init_method", type=str, 
                help="Nombre de segments de la dimension z")
pa.add_argument("-s", "--sample", action="store_true",
                help=("show the first 4 frames to track the moves"))

args = pa.parse_args()

if __name__ == "__main__":
    prot_seq = fasta_read(args.fasta_input)
    amino_acids = [AminoAcidClass.AminoAcid(i + 1, char) for i, char in enumerate(prot_seq)]
    num_iterations = args.iterations

    AminoAcidClass.AminoAcid.initialize(args.init_method)
    
    print(f"\nStarting energy is {AminoAcidClass.AminoAcid.calculate_total_energy()}\n")
    AminoAcidClass.AminoAcid.visualize_molecule()
    
    frame = 1
    for i in tqdm(range(num_iterations), desc="Working on the moves, please wait ..."):
        random_amino_acid = random.choice(AminoAcidClass.AminoAcid._registry)
        move_accepted = random_amino_acid.movement()
        if move_accepted and frame < 5 and args.sample:
            AminoAcidClass.AminoAcid.visualize_molecule()
            frame += 1

    AminoAcidClass.AminoAcid.visualize_molecule()
    print(f"\nFinal energy is {AminoAcidClass.AminoAcid.calculate_total_energy()}\n")
