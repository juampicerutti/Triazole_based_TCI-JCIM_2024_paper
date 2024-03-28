import os
import shutil
from sys import exit
import glob
import subprocess

docking_path = "ADD HERE THE CORRESPONDING DOCKING $PATH"
autodock_GPU_path = "ADD HERE THE $PATH TO AUTODOCK_GPU EXECUTABLE"

def create_variables(system):
    ligands_file = system[0]
    ligands_prefix = ligands_file.split('.')[0]
    ligands_folder = f'{general_dir}/2-Prepare_ligands/{ligands_prefix}'
    receptor_folder = f'{general_dir}/3-Receptor-grids/grids/{system[1]}'
    number_of_runs = system[2]

    return ligands_file, ligands_prefix, ligands_folder, receptor_folder, number_of_runs

def create_docking_results_folder(general_dir,ligands_prefix):

    results_dir = f'{general_dir}/4-Docking_results/{ligands_prefix}/{system[1]}'

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    else:
        print("The results folder already exists... EXITING.")
        exit()

    return results_dir

def docking_execution(ligands_folder,receptor_folder, results_dir, number_of_runs):
    ligands_pdbqt_files = glob.glob(ligands_folder + '/*.pdbqt')
    #receptor_fld_file = glob.glob(receptor_folder + '/grilla-nativa/*.fld')[0]
    receptor_fld_file = glob.glob(receptor_folder + '/*.fld')[0]

    # Create the folder for the .dlg files obtained from the docking
    os.makedirs(results_dir+'/docked_dlg_files/')
    dlg_output_dir = results_dir+'/docked_dlg_files/'

    for ligand in ligands_pdbqt_files:
        ligand_prefix = ligand.split('/')[-1].replace('.pdbqt','')
        shutil.copyfile(ligand, dlg_output_dir+f'/{ligand_prefix}.pdbqt')
        ligand_dlg_file = ligand_prefix+'.dlg'
        print(ligand_prefix)
        # # Execute the docking
        bash_command1 = f'{autodock_GPU_path}/autodock_gpu_64wi -lfile {dlg_output_dir}/{ligand_prefix}.pdbqt -ffile {receptor_fld_file} --xmloutput 0 --contact_analysis 1 --nrun {number_of_runs}'
        process = subprocess.Popen(bash_command1.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        #print(output, error)


        bash_command2 = f'gzip {dlg_output_dir}/{ligand_dlg_file}'
        process = subprocess.Popen(bash_command2.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

# Do not change these variables
general_dir = 'PUT HERE THE FULL PATH TO THE PROJECT FILES'

### Systems to simulate
# Notation: [0]: The ligands file to prepare, [1]: receptor folder, [2]: number_of_runs
system = ['aldehydes.smi','PUT HERE HE FULL $PATH TO THE CORRESPONDING GRID FILES',100]


# Execute the workflow
ligands_file, ligands_prefix, ligands_folder, receptor_folder, number_of_runs = create_variables(system)
results_dir = create_docking_results_folder(general_dir,ligands_prefix)
docking_execution(ligands_folder,receptor_folder, results_dir, number_of_runs)