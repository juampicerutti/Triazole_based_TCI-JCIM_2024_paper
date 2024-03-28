import os
import shutil
from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem
from sys import exit
from meeko import MoleculePreparation
from glob import glob

def create_variables(system):
    ligands_file = system[0]
    ligands_prefix = ligands_file.split('.')[0]

    return ligands_file, ligands_prefix

def create_ligands_folder(general_dir,ligands_file,ligands_prefix):
    if not os.path.isdir(f'{general_dir}/2-Prepare_ligands/{ligands_prefix}'):
        os.makedirs(f'{general_dir}/2-Prepare_ligands/{ligands_prefix}')
        shutil.copyfile(f'{general_dir}/1-Classify_Ligands/{ligands_file}',f'{general_dir}/2-Prepare_ligands/{ligands_prefix}/{ligands_file}')

        ligands_dir = f'{general_dir}/2-Prepare_ligands/{ligands_prefix}'

    else:
        print("The ligands folder already exists... EXITING.")
        exit()

    return ligands_dir

def create_pdbqt_files(ligands_dir,ligands_file):
    suppl = Chem.SmilesMolSupplier(f'{ligands_dir}/{ligands_file}')

    for mol in suppl:
        name = mol.GetProp("_Name")
        # Create the RDKit object
        m3d = Chem.AddHs(mol)
        # Prepare the conformer generation
        props = AllChem.ETKDG()
        props.pruneRmsThresh = 0.25
        props.useRandomCoords = True
        props.numThreads = 0
        # Create the object containing the conformers
        confs = AllChem.EmbedMultipleConfs(m3d,100,props)
        ps = AllChem.MMFFGetMoleculeProperties(m3d)

        # Prepare the confromers
        AllChem.MMFFOptimizeMoleculeConfs(m3d, mmffVariant='MMFF94', maxIters=10)

        conformers_energies_dict={} # Empty dictionary to store conformers

        for conf in confs:
            filename='conformer-' + str(conf + 1) + '.pdb'
            ff = AllChem.MMFFGetMoleculeForceField(m3d,ps,confId=conf)
            ff.Minimize()
            #AllChem.MMFFOptimizeMoleculeConfs(m3d, confId=conf, maxIters=10)
            energy_value = ff.CalcEnergy()
            conformers_energies_dict[conf] = energy_value
            # Generates a sorted dictionary containing the conformers
            conformers_energies_dict_sorted=sorted(conformers_energies_dict.items(), key=lambda x: x[1])


        lowest_energy_conformer = conformers_energies_dict_sorted[0][0]
        print(lowest_energy_conformer)

        # Activate the line below if writing the lowest energy conformer in .pdb format is desired.
        #AllChem.MolToPDBFile(m3d,f'{ligands_dir}/{name}_LEC.pdb',confId=lowest_energy_conformer)

        # These two lines will create an RDKit mol of the lowest energy conformed, in order to advance MEEKO processing
        selected = Chem.MolToMolBlock(m3d,confId=lowest_energy_conformer)
        selected_mol = Chem.MolFromMolBlock(selected, removeHs=False)
        AllChem.MolToPDBFile(selected_mol,f'{ligands_dir}/temp.pdb')

        new_mol = Chem.MolFromPDBFile(f'{ligands_dir}/temp.pdb',removeHs=False)

        meeko_processing(new_mol,ligands_dir,name)

def meeko_processing(mol,ligands_dir,name):

    # Create a dictionary with atom typing
    dictio = {
    "ATOM_PARAMS": {
    "defaults+4H_triazole": [
    {"smarts": "[#1]", "atype": "H",},
    {"smarts": "[#1][#7,#8,#9,#15,#16]","atype": "HD"},
    {"smarts": "[#5]", "atype": "B"},
    {"smarts": "[C]", "atype": "C"},
    {"smarts": "[c]", "atype": "A"},
    {"smarts": "[#7]", "atype": "NA"},
    {"smarts": "[#8]", "atype": "OA"},
    {"smarts": "[#9]", "atype": "F"},
    {"smarts": "[#12]", "atype": "Mg"},
    {"smarts": "[#14]", "atype": "Si"},
    {"smarts": "[#15]", "atype": "P"},
    {"smarts": "[#16]", "atype": "S"},
    {"smarts": "[#17]", "atype": "Cl"},
    {"smarts": "[#20]", "atype": "Ca"},
    {"smarts": "[#25]", "atype": "Mn"},
    {"smarts": "[#26]", "atype": "Fe"},
    {"smarts": "[#30]", "atype": "Zn"},
    {"smarts": "[#35]", "atype": "Br"},
    {"smarts": "[#53]", "atype": "I"},
    {"smarts": "[#7X3v3][a]", "atype": "N", "comment": "pyrrole, aniline"},
    {"smarts": "[#7X3v3][#6X3v4]", "atype": "N", "comment": "amide"},
    {"smarts": "[#7+1]", "atype": "N", "comment": "ammonium, pyridinium"},
    {"smarts": "[SX2]", "atype": "SA", "comment": "sulfur acceptor"},
    {"smarts": "[#1][#6X3]:[#6X3][#7]:[#7][#7]", "atype": "HD", "comment": "4,5-H in 1,2,3-triazole"},
    ]
    }
    }

    preparator = MoleculePreparation(atom_type_smarts=dictio)
    preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()
    pdbqt_string_ready = pdbqt_string.replace('UNL','MOL')

    with open(f'{ligands_dir}/{name}.pdbqt','w') as pdbqt_file:
        pdbqt_file.write(pdbqt_string_ready)

# Do not change these variables
general_dir = 'PUT HERE THE FULL $PATH TO THE PROJECT'

### Systems to simulate
# Notation: [0]: The ligands file to prepare
system = ['aldehydes.smi']

# Execute the workflow
ligands_file, ligands_prefix = create_variables(system)
ligands_dir = create_ligands_folder(general_dir,ligands_file,ligands_prefix)
create_pdbqt_files(ligands_dir,ligands_file)
