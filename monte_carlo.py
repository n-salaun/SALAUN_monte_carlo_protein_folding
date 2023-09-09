import numpy as np
import random
import matplotlib.pyplot as plt


class AminoAcidDictionary(type):
    """This is a metaclass that allows you to create a dictionary of all Amino Acid instances."""
    def __iter__(cls):
        return iter(cls._registry)


class AminoAcid(object, metaclass=AminoAcidDictionary):
    # The summary holds a dictionary with the info of all Amino Acid
    summary = {}
    _registry = []

    def __init__(self, number, class_hp):
        self._registry.append(self)
        self.number = number
        self.class_hp = class_hp
        self.xcoord = 0
        self.ycoord = 0

        AminoAcid.summary[f"aa{self.number}"] = (self.class_hp, self.xcoord, self.ycoord)

    def set_coordinates(self, x, y):
        self.xcoord = x
        self.ycoord = y

    def get_coordinates(self):
        return self.xcoord, self.ycoord

    @staticmethod
    def initialize(mode="linear"):
        if mode == "random":
            for attempt in range(5):  # Number of attempts to distribute amino acids randomly
                used_coordinates = []
                neighbors = [(1, 0), (-1, 0), (0, 1), (0, -1)]

                for aa_object in AminoAcid:
                    if aa_object == AminoAcid._registry[0]:
                        aa_object.set_coordinates(0, 0)
                        used_coordinates.append((0, 0))
                    else:
                        loop_breaker = False
                        iteration = 0
                        while not loop_breaker:
                            iteration += 1
                            aa_key = f"aa{aa_object.number}"
                            last_x, last_y = used_coordinates[-1]
                            random_neighbor = random.choice(neighbors)
                            new_x, new_y = last_x + random_neighbor[0], last_y + random_neighbor[1]
                            if (new_x, new_y) not in used_coordinates:
                                aa_object.set_coordinates(new_x, new_y)
                                used_coordinates.append((new_x, new_y))
                                AminoAcid.summary[aa_key] = (aa_object.class_hp, new_x, new_y)
                                loop_breaker = True

                            if iteration > 50:
                                print(f"Unable to find a random coordinate for amino acid {aa_object.number}")
                                break
                        if iteration > 50:
                            break

                if len(used_coordinates) == len(AminoAcid._registry):
                    break  # Successfully distributed all amino acids
                else:
                    print("Resetting initialization and trying again...")
                    used_coordinates = []  # Reset used_coordinates for the next attempt

            if len(used_coordinates) < len(AminoAcid._registry):
                print("Failed to distribute all amino acids after multiple attempts.")
        else:
            if mode != "linear":
                print("Invalid mode format, should be 'random' or 'linear', proceeding with 'linear'")

            for aa_object in AminoAcid:
                aa_key = f"aa{aa_object.number}"
                aa_object.set_coordinates(aa_object.number - 1, 0)
                AminoAcid.summary[aa_key] = (aa_object.class_hp, aa_object.number - 1, 0)

    @staticmethod
    def calculate_total_energy():
        total_energy = 0
        for i, aa_object in enumerate(AminoAcid):
            if aa_object.class_hp == 'H':
                neighbors = [
                    (aa_object.xcoord + 1, aa_object.ycoord),
                    (aa_object.xcoord - 1, aa_object.ycoord),
                    (aa_object.xcoord, aa_object.ycoord + 1),
                    (aa_object.xcoord, aa_object.ycoord - 1)
                ]

                for neighbor_x, neighbor_y in neighbors:
                    neighbor_coords = (neighbor_x, neighbor_y)
                    neighbor_aa = None

                    # Find the neighbor amino acid based on its coordinates
                    for aa_obj in AminoAcid:
                        if aa_obj.get_coordinates() == neighbor_coords:
                            neighbor_aa = aa_obj
                            break

                    if (
                            neighbor_aa is not None and
                            neighbor_aa.class_hp == 'H' and
                            (abs(i - AminoAcid._registry.index(neighbor_aa)) != 1)
                    ):
                        total_energy -= 1
        return total_energy / 2

    def visualize(self):
        if self.class_hp == 'P':
            marker_color = 'red'
            marker_style = 'o'  # Use a circle for P class amino acids
            label = "P Class"
        else:
            marker_color = 'blue'
            marker_style = 's'  # Use a square for H class amino acids
            label = "H Class"

        plt.plot([self.xcoord], [self.ycoord], marker_style, color=marker_color, markersize=5, zorder=2)

    @staticmethod
    def visualize_molecule():
        plt.figure(figsize=(8, 8))
        legend_added = False  # Variable to keep track of whether the legend has been added

        for aa_object in AminoAcid:
            aa_object.visualize()

        # Draw lines connecting amino acids
        for i in range(len(AminoAcid._registry) - 1):
            aa1 = AminoAcid._registry[i]
            aa2 = AminoAcid._registry[i + 1]
            plt.plot([aa1.xcoord, aa2.xcoord], [aa1.ycoord, aa2.ycoord], 'black', zorder=1)

        plt.title('Amino Acid Molecule Visualization')
        plt.plot([], [], 'ro', markersize=5, label='P Class')
        plt.plot([], [], 'bs', markersize=5, label='H Class')
        plt.legend()
        plt.show()

    # Moveset
    def end_move(self):
        return ()

    def corner_move(self):
        return ()

    def crankshaft_move(self):
        return ()

    def pull_move(self):
        return ()


###Functions

def fasta_read(name):
    """The function read a FASTA file, and return the proteic sequence in HP
    model type"""
    SEQ = []
    # read the file and select only the sequence
    with open(name, "r") as filout:
        for lines in filout:
            if not lines.startswith(">"):
                SEQ.extend(lines.strip())

    # Change the sequence from amino acid to HP model
    final_seq = []
    for char in SEQ:
        if char in "VIFLMCW":
            final_seq.append("H")
        else:
            final_seq.append("P")

    # return string of the final sequence in the HP model
    return ''.join(final_seq)

###Main

prot_seq = "PHPHPPPHPHPPHPHHPPHPHPPPHPHPHPHPPH"
amino_acids = [AminoAcid(i + 1, char) for i, char in enumerate(prot_seq)]

AminoAcid.initialize("random")
print(AminoAcid.calculate_total_energy())
AminoAcid.visualize_molecule()
