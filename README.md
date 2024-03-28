# Structure-aided computational design of triazole-based targeted covalent inhibitors fo Cruzipain.

## Repository description:

In this repository, the source datafiles and Python scripts used as part of the scientific research reported in the manuscript are included.

The corresponding data is organized as part of the following folders:

1-`Training_set`: this folder containes source `Training_Set.smi` file including all compounds in included as training set for the molecular docking workflow. Compounds are provided in the SMILES notation within this file. Additionally, this folder contains the corresponding ligands in the .pdbtq format as required by the molecular docking software (AutoDock-GPU).

2- `Receptor-grids`: this folder contains the receptor files as used in the molecular docking assays with AutoDock-GPU.

3- `3-4FPMK-mostrelevant`: this folder contains the .pdbqt structures of the hit compuounds derived from the 4FPMK warhead.

4- `Aldehydes`: this folder contains the .pdbqt structures of the hit compuounds derived from the Ald warhead.

5- `Esters`: this folder contains the .pdbqt structures of the hit compuounds derived from the Es warhead.

6- `Python_scripts`: this folder contains general Python scripts used during the virtual screening workflow.
