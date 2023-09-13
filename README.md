# SALAUN_monte_carlo_protein_folding

Student project submitted for the M2 of BioInformatics at University of Paris (Paris 7). The goal is to use a Monte Carlo method, specifically REMC, to find the best folded version of a protein in a 2D lattice. Issued 02/09/2023, due 13/09/2023

## Setup the conda environment

Clone the repository:

```bash
git clone git@github.com:n-salaun/SALAUN_monte_carlo_protein_folding.git
```

Move to the new directory:

```bash
cd SALAUN_monte_carlo_protein_folding
```

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Create the `env-montecarlo` conda environment:

```
conda env create -f env-montecarlo.yml
```

Load the `monte-carlo` conda environment:

```
conda activate env-montecarlo
```

## Run the program

### Main command 'folding.py'

```bash
python ./Scripts/folding.py fasta_input iterations init_method -s -e
```

### Dictionnary

```
fasta_input : Name of the starting sample protein
iterations : Number of tries to move an amino acid
init_method : Initialization method, should be either linear or random
-s : show the first 4 frames to track the moves
-e : Show a plot of the energy 
```

All results will be put in the Results folder
