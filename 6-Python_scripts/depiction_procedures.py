
import sqlite3
import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D

docking_path = "DEFINE HERE THE CORRESPONDING DOCKING $PATH"

import sys
sys.path.append(f'{docking_path}}/drug_screening_package')
import general_procedures.objects_processing as obj_proc

def generate_depiction_grid(db_name, table_name, smiles_column_name, label_column_name, max_mols_ppage,output_path):
    """
    This function will create .png documents depicting grids of molecules corresponding to the table and column of smiles read as input
    ------
    Parameters
    ------
    - db_name: the full path of the database cotaining the SMILES to be depicted in the grid.
    - table_name: the table included in the 'db-name' that contains the SMILES strings to depict.
    - column_name: the field name containing the SMILES strings.
    - max_mols_ppage: the maximum number of molecules that will be depipected per page.
    - output_path: the full path were the generated .png files of the depicted images will be saved.    
    """

    conn = sqlite3.connect(db_name)
    sql = f'''SELECT {smiles_column_name},{label_column_name}
            FROM {table_name};'''

    df_query = pd.read_sql_query(sql,conn)

    smiles_list = df_query[smiles_column_name].to_list()
    labels_list = df_query[label_column_name].to_list()
    indexes_list = [f'cpd_index: {x}' for x in range(len(smiles_list))]
        
    smiles_list_chunks = obj_proc.split_list_in_chunks(smiles_list, max_mols_ppage)
    labels_list_chunks = obj_proc.split_list_in_chunks(labels_list, max_mols_ppage)
    indexes_list_chunks = obj_proc.split_list_in_chunks(indexes_list, max_mols_ppage)
    
    counter = 0 # used to number figures with chunks
    for index, list_item in enumerate(smiles_list_chunks):
        mols_item = [Chem.MolFromSmiles(smile) for smile in list_item]
        img = MolsToGridImage(mols=mols_item, legends=indexes_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
        img.save(f"{output_path}/{table_name}_{counter}.png")
        counter+=1



if __name__ == "__main__":

    #db_name = "/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/level_1_projects/sample_level1.db"
    db_name = f"{docking_path}/docked.db"

    generate_depiction_grid(db_name,'REACTANT',"SMILES","inchi_key",25,f"{docking_path}/compounds_depictions")