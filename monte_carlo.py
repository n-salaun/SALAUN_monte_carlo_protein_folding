import numpy as np
import random

class AminoAcidDictionary(type):
    def __iter__(cls):
        return iter(cls._registry)


class AminoAcid(object, metaclass=AminoAcidDictionary):
    # The summary holds a dictionary with the info of all Amino Acid instances
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
        coordinates = {}

        if mode == "random":
            aa_key = "aa1"
            coordinates[aa_key] = (0, 0)
            AminoAcid.summary[aa_key] = (AminoAcid(1).class_hp, 0, 0)

            for number in range(2, len(AminoAcid.summary) + 1):
                aa_key = f"aa{number}"
                neighbors = [(1, 0), (-1, 0), (0, 1), (0, -1)]

                while True:
                    random_neighbor = random.choice(neighbors)
                    x, y = AminoAcid.summary[f"aa{number - 1}"][1], AminoAcid.summary[f"aa{number - 1}"][2]
                    new_x, new_y = x + random_neighbor[0], y + random_neighbor[1]
                    if (new_x, new_y) not in coordinates.values():
                        coordinates[aa_key] = (new_x, new_y)
                        AminoAcid.summary[aa_key] = (AminoAcid(number, 0).class_hp, new_x, new_y)
                        break

        else:
            if mode != "linear":
                print("Format de mode invalide, devrait être 'random' ou 'linear', procédant avec 'linear'")

            for aa_object in AminoAcid:
                aa_key = f"aa{aa_object.number}"
                coordinates[aa_key] = (aa_object.number - 1, 0)
                AminoAcid.summary[aa_key] = (aa_object.class_hp, aa_object.number - 1, 0)

        for aaobject in AminoAcid:
            aa_key = f"aa{aaobject.number}"
            new_x, new_y = coordinates[aa_key]
            aaobject.set_coordinates(new_x, new_y)


    @staticmethod
    def calculate_total_energy():
        total_energy = 0
        for key, (class_hp, x, y) in AminoAcid.summary.items():
            if class_hp == "H":
                neighbors = [
                    (x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)
                ]
                for neighbor_x, neighbor_y in neighbors:
                    neighbor_key = f"aa{neighbor_x + 1}{neighbor_y}"
                    if neighbor_key in AminoAcid.summary and AminoAcid.summary[neighbor_key][0] == 'H':
                        total_energy -= 1
        return total_energy / 2

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
    #read the file and select only the sequence
    with open(name, "r") as filout:
        for lines in filout:
            if not lines.startswith(">"):
                SEQ.extend(lines.strip())
    
    #Change the sequence from amino acid to HP model
    FINAL_SEQ = []
    for char in SEQ:
        if char in "VIFLMCW":
            FINAL_SEQ.append("H")
        else:
            FINAL_SEQ.append("P")
            
    #return string of the final sequence in the HP model
    return ''.join(FINAL_SEQ)

###Main

# Créer des instances d'AminoAcid
aa1 = AminoAcid(1, "H")
aa2 = AminoAcid(2, "P")
aa3 = AminoAcid(3, "H")
aa4 = AminoAcid(4, "H")
aa5 = AminoAcid(5, "P")
aa6 = AminoAcid(6, "H")

# Afficher les informations des instances
print(aa1.number, aa1.class_hp, aa1.get_coordinates())
print(aa2.number, aa2.class_hp, aa2.get_coordinates())
print(aa3.number, aa3.class_hp, aa3.get_coordinates())
print(aa4.number, aa4.class_hp, aa4.get_coordinates())
print(aa5.number, aa5.class_hp, aa5.get_coordinates())
print(aa6.number, aa6.class_hp, aa6.get_coordinates())

# Initialiser les coordonnées de manière aléatoire
AminoAcid.initialize("linear")

# Afficher les nouvelles coordonnées
print(aa1.number, aa1.class_hp, aa1.get_coordinates())
print(aa2.number, aa2.class_hp, aa2.get_coordinates())
print(aa3.number, aa3.class_hp, aa3.get_coordinates())
print(aa4.number, aa4.class_hp, aa4.get_coordinates())
print(aa5.number, aa5.class_hp, aa5.get_coordinates())
print(aa6.number, aa6.class_hp, aa6.get_coordinates())

print(AminoAcid.summary)
print(AminoAcid.calculate_total_energy())
