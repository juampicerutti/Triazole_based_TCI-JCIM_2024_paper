from rdkit.Chem import rdmolops
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def generate_clean_smiles(smiles):
    sanitized_smiles = sanitize_smiles_errors(smiles)
    smiles_free_stereo = remove_smiles_stereo(sanitized_smiles)
    largest_mol_smiles = get_largest_fragment(smiles_free_stereo)
    # Add here potential processing steps for cleaning/sanitizing the SMILES notation 
    #####
    #
    clean_smiles = largest_mol_smiles # The last cleaning step should be assigned to 'clean_smiles'
    inchi_key = compute_inchi_key(largest_mol_smiles)
    return clean_smiles, inchi_key

def remove_smiles_stereo(smiles):
    #smiles_free_stereo = smiles.replace("@","").replace("[","").replace("]","").replace("H","")
    smiles_free_stereo = smiles.replace("@","")
    if "[CH]" in smiles_free_stereo:
        smiles_free_stereo = smiles_free_stereo.replace("[CH]","C")
    return smiles_free_stereo

def get_largest_fragment(smiles):   
    mol = Chem.MolFromSmiles(smiles)
    mol_frags = rdmolops.GetMolFrags(mol, asMols = True)
    largest_mol = max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())
    largest_mol_smiles = Chem.MolToSmiles(largest_mol)    
    return largest_mol_smiles

def compute_inchi_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    inchi_key = Chem.MolToInchiKey(mol)
    return inchi_key

def depict_smiles(smiles):
    try: 
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol)
        return img
    except:
        print(f"Error procesando: {smiles}")

def sanitize_smiles_errors(smiles):
    
    if "[nH]" in smiles:
        sanitized_smiles = smiles.replace("[nH]","[NH]")
        return sanitized_smiles
    else:
        return smiles