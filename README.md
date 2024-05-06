# Structure-aided computational design of triazole-based targeted covalent inhibitors fo Cruzipain.

## Repository description:

In this repository, the source datafiles and Python scripts used as part of the scientific research reported in the manuscript are included.

The corresponding data is organized as part of the following folders:

0-`Database`: this folder contains a database (.db file) that includes the SMILES, inchi keys and properties of the filtered building blocks and 4FPMF-, Es- and Ald-based triazole derivatives virtually synthesised, screened and analysed. Due to size limitations, the complete library of potential reagents extracted from eMolecules has not been included in this repository, but is available upon request to the authors.

1-`Training_set`: this folder contains source `Training_Set.smi` file including all compounds considered as training set for the molecular docking workflow. Compounds are provided in the SMILES notation within this file. Additionally, this folder contains the corresponding ligands in the .pdbtq format as required by the molecular docking software (AutoDock-GPU).

2- `Receptor-grids`: this folder contains the receptor files as used in the molecular docking assays with AutoDock-GPU.

3- `3-4FPMK-mostrelevant`: this folder contains the .pdbqt structures of the hit compounds derived from the 4FPMK warhead.

4- `Aldehydes`: this folder contains the .pdbqt structures of the hit compounds derived from the Ald warhead.

5- `Esters`: this folder contains the .pdbqt structures of the hit compounds derived from the Es warhead.

6- `Python_scripts`: this folder contains general Python scripts used during the virtual screening workflow, including molecular docking, MM-MD and MMPBSA studies.

7- `MM_MD`: this folder contains the parameter/topology (.prmtop) and coordinate (.inpcrd) files of the Training set and of the synthesised Es- and Ald-based triazole derivatives, required to perform MM-MD simulations.
