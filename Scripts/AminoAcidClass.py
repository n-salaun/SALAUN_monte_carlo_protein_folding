import math
import random
import matplotlib.pyplot as plt
import numpy as np

class AminoAcidDictionary(type):
    """This is a metaclass of AminoAcid to iterate through instances, or get infos about all the instances"""
    def __iter__(cls):
        return iter(cls._registry)


class AminoAcid(object, metaclass=AminoAcidDictionary):
    """
        Represents individual amino acids in the molecular system.

        Attributes:
        - number: A unique identifier for the amino acid.
        - class_hp: Class type ('H' for hydrophobic or 'P' for polar).
        - xcoord: X-coordinate of the amino acid's position.
        - ycoord: Y-coordinate of the amino acid's position.

        Methods:
        - set_coordinates(x, y): Set the coordinates of the amino acid.
        - get_coordinates(): Get the current coordinates of the amino acid.
        - initialize(mode="linear"): Initialize the positions of amino acids, either randomly or linearly.
        - movement(): Simulate movement of the amino acid and calculate energy changes.
        - calculate_total_energy(): Calculate the total energy of the system based on the HP model.
        - visualize(): Visualize the amino acid's position.
        - visualize_molecule(): Visualize the entire molecular structure.
        - used_coordinates(): Get a list of coordinates that are currently occupied by amino acids.
        - distance(point1, point2): Calculate the distance between two points.

        """
    summary = {}
    _registry = []

    def __init__(self, number, class_hp):
        """Initialize an AminoAcid instance with a unique number and class type."""
        self._registry.append(self)
        self.number = number
        self.class_hp = class_hp
        self.xcoord = 0
        self.ycoord = 0

        AminoAcid.summary[f"aa{self.number}"] = (self.class_hp, self.xcoord, self.ycoord)

    def set_coordinates(self, x, y):
        """Set the new coordinates of the amino acid."""
        self.xcoord = x
        self.ycoord = y

    def get_coordinates(self):
        """Get the current coordinates of the amino acid."""
        return self.xcoord, self.ycoord

    @staticmethod
    def initialize(mode="linear"):
        """
                Initialize the positions of amino acids.

                Args:
                - mode (str): The initialization mode, either 'random' or 'linear'.

                """
        if mode == "random":
            for attempt in range(10000):  # Number of attempts to distribute amino acids randomly
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
                                break # Unable to find a next random coordinate for amino acid
                        if iteration > 50:
                            break

                if len(used_coordinates) == len(AminoAcid._registry):
                    break  # Successfully distributed all amino acids
                else:
                    #Resetting initialization and trying again...
                    used_coordinates = []  # Reset used_coordinates for the next attempt

            if len(used_coordinates) < len(AminoAcid._registry):
                print("Failed to distribute all amino acids after multiple attempts. You could try linear")
        else:
            if mode != "linear":
                print("Invalid mode format, should be 'random' or 'linear', proceeding with 'linear'")

            for aa_object in AminoAcid:
                aa_key = f"aa{aa_object.number}"
                aa_object.set_coordinates(aa_object.number - 1, 0)
                AminoAcid.summary[aa_key] = (aa_object.class_hp, aa_object.number - 1, 0)

    def movement(self, threshold):
        """Simulate movement of the amino acid and accept or reject the new position"""
        possible_moves = []
        end_move = self.end_move()
        corner_move = self.corner_move()
        crankshaft_move = self.crankshaft_move()
        prev_total_energy = AminoAcid.calculate_total_energy()

        if end_move:
            possible_moves.append(("end", end_move))
        if corner_move:
            possible_moves.append(("corner", corner_move))
        if crankshaft_move:
            possible_moves.append(("crankshaft", crankshaft_move))

        if possible_moves:
            if len(possible_moves) == 1:
                selected_move_type, selected_move = possible_moves[0]
            else:
                selected_move_type, selected_move = random.choice(possible_moves)

            if selected_move_type == "end":
                prev_coords = self.get_coordinates()
                if len(selected_move) > 1:
                    new_coords = random.choice(selected_move)
                else:
                    new_coords = selected_move[0]
                self.set_coordinates(*new_coords)
            elif selected_move_type == "corner":
                prev_coords = self.get_coordinates()
                new_coords = selected_move
                self.set_coordinates(*new_coords)
            elif selected_move_type == "crankshaft":
                prev_coords = {}
                for key, val in selected_move.items():
                    prev_coords[key] = AminoAcid._registry[int(key)-1].get_coordinates()
                    AminoAcid._registry[int(key)-1].set_coordinates(*val)

            # Calculate the new total energy
            new_total_energy = AminoAcid.calculate_total_energy()
            
            if threshold is not None:
                K_b = 0.0019872041
                if (new_total_energy <= prev_total_energy or
                random.random() > math.exp(-(new_total_energy - prev_total_energy) / (threshold[0] * K_b))):
                    return True
                else:
                    return False
            
            # If the new energy is higher, revert the move
            elif new_total_energy > prev_total_energy:
                if selected_move_type == "crankshaft":
                    for key, val in prev_coords.items():
                        AminoAcid._registry[int(key)-1].set_coordinates(val[0], val[1])
                else:
                    self.set_coordinates(*prev_coords)
                return False
            else:
                return True
        else:
            return False

    @staticmethod
    def calculate_total_energy():
        """Calculate the total energy of the system based on the HP model."""
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
        """Visualize the amino acid's position."""
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
    def visualize_molecule(name):
        """Visualize the entire molecular structure."""
        plt.figure(figsize=(15, 15))

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
        plt.savefig("./Results/" + name + ".png")

    @staticmethod
    def used_coordinates():
        """Get a list of coordinates that are currently occupied by amino acids."""
        used_coordinates = []
        for aa_object in AminoAcid:
            used_coordinates.append((aa_object.xcoord, aa_object.ycoord))
        return used_coordinates

    @staticmethod
    def distance(point1, point2):
        """Function to calculate the distance between two points"""
        return math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)

    # Moveset
    def end_move(self):
        """Function to check if an end move is possible on the current amino acid
        returns the new coords possible if true"""
        index = AminoAcid._registry.index(self)
        if index == 0:
            next_aa = AminoAcid._registry[1]
        elif index == len(AminoAcid._registry) - 1:
            next_aa = AminoAcid._registry[-2]
        else:
            return None

        neighbors = [
            (next_aa.xcoord + 1, next_aa.ycoord),
            (next_aa.xcoord - 1, next_aa.ycoord),
            (next_aa.xcoord, next_aa.ycoord + 1),
            (next_aa.xcoord, next_aa.ycoord - 1)
        ]
        valid_neighbors = [(x, y) for x, y in neighbors if (x, y) not in AminoAcid.used_coordinates()]
        return valid_neighbors

    def corner_move(self):
        """Function to check if a corner move is possible, by trying to check if i-1 and i+1 aa can form a square,
        which gives us the coordinates of the new coords to move the aa into"""

        index = AminoAcid._registry.index(self)
        if index == 0 or index == len(AminoAcid._registry) - 1:
            # Corner moves are not possible at the chain ends
            return None

        current_coords = (self.xcoord, self.ycoord)
        prev_aa_coords = AminoAcid._registry[index - 1].get_coordinates()
        next_aa_coords = AminoAcid._registry[index + 1].get_coordinates()
        dist1_2 = AminoAcid.distance(prev_aa_coords, current_coords)
        dist1_3 = AminoAcid.distance(prev_aa_coords, next_aa_coords)
        dist2_3 = AminoAcid.distance(current_coords, next_aa_coords)

        # Check if two of the distances are equal (forming a corner angle)
        if dist1_2 == dist1_3:
            # Calculate the fourth point's coordinates
            new_x = current_coords[0] + next_aa_coords[0] - prev_aa_coords[0]
            new_y = current_coords[1] + next_aa_coords[1] - prev_aa_coords[1]
        elif dist1_2 == dist2_3:
            new_x = prev_aa_coords[0] + next_aa_coords[0] - current_coords[0]
            new_y = prev_aa_coords[1] + next_aa_coords[1] - current_coords[1]
        elif dist1_3 == dist2_3:
            new_x = prev_aa_coords[0] + current_coords[0] - next_aa_coords[0]
            new_y = prev_aa_coords[1] + current_coords[1] - next_aa_coords[1]
        else:
            return None
        if (new_x, new_y) not in AminoAcid.used_coordinates():
            return new_x, new_y

    def crankshaft_move(self):
        """Function to check if a crankshaft move is possible, by trying to check if the aa can move
        in a hinge motion by being in a U shape"""
        def move_test(aa1_coords, aa2_coords, aa3_coords, aa4_coords):
            new_aa2 = None
            new_aa3 = None

            if aa1_coords[0] == aa2_coords[0] and aa3_coords[0] == aa4_coords[0]:
                if aa1_coords[1] - aa2_coords[1] > 0:
                    new_aa2 = aa1_coords[0], aa1_coords[1] + 1
                    new_aa3 = aa4_coords[0], aa4_coords[1] + 1
                elif aa1_coords[1] - aa2_coords[1] < 0:
                    new_aa2 = aa1_coords[0], aa1_coords[1] - 1
                    new_aa3 = aa4_coords[0], aa4_coords[1] - 1

            if aa1_coords[1] == aa2_coords[1] and aa3_coords[1] == aa4_coords[1]:
                if aa1_coords[0] - aa2_coords[0] > 0:
                    new_aa2 = aa1_coords[0] + 1, aa1_coords[1]
                    new_aa3 = aa4_coords[0] + 1, aa4_coords[1]
                elif aa1_coords[0] - aa2_coords[0] < 0:
                    new_aa2 = aa1_coords[0] - 1, aa1_coords[1]
                    new_aa3 = aa4_coords[0] - 1, aa4_coords[1]
            return new_aa2, new_aa3

        index = AminoAcid._registry.index(self)
        if index <= 1 or index >= len(AminoAcid._registry)-2:
            return None

        aa_minus_2 = AminoAcid._registry[index - 2].get_coordinates()
        aa_minus_1 = AminoAcid._registry[index - 1].get_coordinates()
        current_aa = AminoAcid._registry[index].get_coordinates()
        aa_plus_1 = AminoAcid._registry[index + 1].get_coordinates()
        aa_plus_2 = AminoAcid._registry[index + 2].get_coordinates()

        if (AminoAcid.distance(current_aa, aa_minus_1) == AminoAcid.distance(current_aa, aa_plus_1) ==
                AminoAcid.distance(aa_plus_1, aa_plus_2) == AminoAcid.distance(aa_minus_1, aa_plus_2)):
            new_coords = move_test(aa_minus_1, current_aa, aa_plus_1, aa_plus_2)
            new_current_aa = new_coords[0]
            new_aa_plus_1 = new_coords[1]
            if new_current_aa not in AminoAcid.used_coordinates() and new_aa_plus_1 not in AminoAcid.used_coordinates():
                return {f"{index+1}": new_current_aa, f"{index+2}": new_aa_plus_1}

        elif (AminoAcid.distance(aa_plus_1, aa_minus_2) == AminoAcid.distance(current_aa, aa_minus_1) ==
              AminoAcid.distance(current_aa, aa_plus_1) == AminoAcid.distance(aa_minus_2, aa_plus_1)):
            new_coords = move_test(aa_minus_2, aa_minus_1, current_aa, aa_plus_1)
            new_aa_minus_1 = new_coords[0]
            new_current_aa = new_coords[1]
            if new_aa_minus_1 not in AminoAcid.used_coordinates() and new_current_aa not in AminoAcid.used_coordinates():
                return {f"{index}": new_aa_minus_1, f"{index+1}": new_current_aa}

        else:
            return None

    def pull_move(self):
        """Should be the most efficient move, capable of minimizing with efficiency the energy
        THIS FUNCTION IS NOT WORKING AS INTENDED AND CAUSES PROBLEMS SOMETIMES"""
        def check_neighbors(coords, empty=False):
            valid_neighbors = []
            neighbors = [
                (coords[0] + 1, coords[1]),
                (coords[0] - 1, coords[1]),
                (coords[0], coords[1] + 1),
                (coords[0], coords[1] - 1)
            ]
            for neighbor in neighbors:
                if empty and neighbor not in AminoAcid.used_coordinates():
                    valid_neighbors.append(neighbor)
                elif neighbor in AminoAcid.used_coordinates() and empty is False:
                    valid_neighbors.append(neighbor)
            return valid_neighbors

        index = AminoAcid._registry.index(self)-1

        if index <= 1 or index >= len(AminoAcid._registry)-1:
            # Pull moves can only be initiated from residue i >= 2
            return None

        # Get the current coordinates of the amino acid
        current_aa = AminoAcid._registry[index]
        prev_aa = AminoAcid._registry[index - 1]
        next_aa = AminoAcid._registry[index + 1]

        neighbors = check_neighbors(next_aa.get_coordinates(), empty=True)
        if len(neighbors) < 1:
            return None
        elif len(neighbors) == 1:
            pos_L = neighbors[0]
        else:
            pos_L = random.choice(neighbors)

        vector_i = (pos_L[0] - current_aa.get_coordinates()[0], pos_L[1] - current_aa.get_coordinates()[1])
        pos_C = (prev_aa.get_coordinates()[0] + vector_i[0], prev_aa.get_coordinates()[1] + vector_i[1])

        if pos_C in AminoAcid.used_coordinates():
            return None

        if pos_L in check_neighbors(prev_aa.get_coordinates()):
            current_aa.set_coordinates(*pos_L)
            return "Done"

        j = AminoAcid._registry[index - 2]
        print(f"j coord : {j.get_coordinates()}")

        if j.get_coordinates() in check_neighbors(pos_C):
            current_aa.set_coordinates(*pos_L)
            prev_aa.set_coordinates(*pos_C)
            return "Done"

        else:
            j.set_coordinates(*current_aa.get_coordinates())
            current_aa.set_coordinates(*pos_L)
            prev_aa.set_coordinates(*pos_C)
            increment = 2
            flag = True
            while flag:
                increment += 1
                new_j = AminoAcid._registry[index - increment]
                if new_j == AminoAcid._registry[0]:
                    new_j.set_coordinates(*AminoAcid._registry[index - increment + 1].get_coordinates())
                    flag = False
                elif new_j in check_neighbors(AminoAcid._registry[index - increment + 1].get_coordinates()):
                    flag = False
                else:
                    new_j.set_coordinates(*AminoAcid._registry[index - increment + 2].get_coordinates())

