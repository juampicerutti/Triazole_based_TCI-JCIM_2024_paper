import pandas as pd
from rdkit import Chem
from dimorphite_dl import DimorphiteDL

dimorphite_dl = DimorphiteDL(

    min_ph=5.4,
    max_ph=5.6,
    max_variants=128,
    label_states=False,
    pka_precision=0.5
)

files_path = "PUT HERE THE FULL $PATH TO THE PROJECT"

df = pd.read_csv(f'{files_path}/aldehydes.smi')
info_file_ion = 'List_of_smiles_triazoles_ion_ALL.smi'
for smi, name in  zip(df['isosmiles'], df['code']):
        tupla = (tuple(dimorphite_dl.protonate(smi)))
        print(tupla[0])
        for count, smi_ion in enumerate(tupla):
            if "[nH+]" not in smi_ion and "[n-]" not in smi_ion:
                print(name+'_'+str(count), smi_ion)
                with open(info_file_ion, 'a') as info_file:
                    info_file.write(f'{smi_ion},{name+"_"+str(count)}\n')


header_list = ["isosmiles", "code"]
df1 = pd.read_csv(f'{files_path}/aldehydes.smi',names=header_list)
print(df1)

files = ['aldehydes.smi']

header_list = ["isosmiles", "code", "NO"]
for file in files:
    df2_pre = pd.read_csv(f'{files_path}/{file}',sep=' ',names=header_list)
    df2 = df2_pre.iloc[:, 0:2]

    test = pd.merge(df1, df2, how="inner", on=["isosmiles"])
    test.to_csv(f'pH5_5-{file}', index=False)
    print(test)
