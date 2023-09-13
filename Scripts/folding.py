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
                help="Initialization method, should be either linear or random")
pa.add_argument("-s", "--sample", action="store_true",
                help=("show the first 4 frames to track the moves"))
pa.add_argument("-e", "--energy", action="store_true",
                help=("Show the energy plot"))
pa.add_argument("--temp", nargs=1, const=None, type=int, help=("value to get out of energy wells"))

args = pa.parse_args()

if __name__ == "__main__":
    prot_seq = fasta_read(args.fasta_input)
    amino_acids = [AminoAcidClass.AminoAcid(i + 1, char) for i, char in enumerate(prot_seq)]
    num_iterations = args.iterations

    AminoAcidClass.AminoAcid.initialize(args.init_method)
    
    print(f"\nStarting energy is {AminoAcidClass.AminoAcid.calculate_total_energy()}\n")
    AminoAcidClass.AminoAcid.visualize_molecule("Start")
    
    energy_values = []
    frame = 1
    for i in tqdm(range(num_iterations), desc="Working on the moves, please wait ..."):
        random_amino_acid = random.choice(AminoAcidClass.AminoAcid._registry)
        move_accepted = random_amino_acid.movement(args.temp)
        
        if move_accepted:
            energy = AminoAcidClass.AminoAcid.calculate_total_energy()
            energy_values.append(energy)
        
            if frame < 5 and args.sample:
                AminoAcidClass.AminoAcid.visualize_molecule(f"Frame_{frame}")
                frame += 1

    AminoAcidClass.AminoAcid.visualize_molecule("End")
    print(f"\nFinal energy is {AminoAcidClass.AminoAcid.calculate_total_energy()}\n")
    
    
    plt.close('all')
    # Plot the energy values
    if args.energy:
        plt.plot(energy_values)
        plt.xlabel("Iteration")
        plt.ylabel("Energy")
        plt.title("Energy per Iteration")
        plt.savefig("./Results/energy.png")
