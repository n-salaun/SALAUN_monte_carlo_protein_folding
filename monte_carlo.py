import numpy as np
import random

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

def energy_calculator():
    """Function to calculate the total energy of the protein"""
    return(energy)

def check_available_moves():
    return()

#Movesets
def corner_move():
    return()
    
def crankshaft_move():
    return()

def pull_move():
    return()

###Main

#Select a random amino acid

#check possible moves

#select a move

#apply the move

#check if the molecule is more stable, or accepted

