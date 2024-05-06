# Structure-aided computational design of triazole-based targeted covalent inhibitors fo Cruzipain.

## Repository description:

In this repository, the source datafiles and Python scripts used as part of the scientific research reported in the manuscript are included.

The corresponding data is organized as part of the following folders:

0-`Database`: this folder contains a .db file which includes the SMILES, inchi_key and some other structural and/or drug-like properties of the building blocks, 4FPMK-, Ester- and Aldehyde-based triazole derivatives evaluated against CZP. Additional datasets (e.g.: Emolecules compounds, full lists of amines, aldehydes and amino acids, etc.) were not included due to size limiting aspects; but are available from the authors upon request.

1-`Training_set`: this folder contains source `Training_Set.smi` file including all compounds in included as training set for the molecular docking workflow. Compounds are provided in the SMILES notation within this file. Additionally, this folder contains the corresponding ligands in the .pdbtq format as required by the molecular docking software (AutoDock-GPU).

2- `Receptor-grids`: this folder contains the receptor files as used in the molecular docking assays with AutoDock-GPU.

3- `3-4FPMK-mostrelevant`: this folder contains the .pdbqt structures of the hit compuounds derived from the 4FPMK warhead.

4- `Aldehydes`: this folder contains the .pdbqt structures of the hit compuounds derived from the Ald warhead.

5- `Esters`: this folder contains the .pdbqt structures of the hit compuounds derived from the Es warhead.

6- `Python_scripts`: this folder contains general Python scripts used during the virtual screening workflow, including input files to perform molecular docking, MM-MD and MMPBSA energetics analysis.

7- `Python_scripts`: this folder contains all .prmtop and .inpcrd files required for MM-MD simulations of the training set and the Es- and Ald-based triazole derivatives synthesised.
