import os
import pandas as pd
import shutil

#Script to move top 3 PDB files from TFII Predictions
#Inputs: 
#  - Model Directory with AF out Dirs
#  - Input Key file (same file used in process_inputs scripts)
#  - Output Directory 
#Warning - Contains lines that are hardcoded to the 6NMI PDB


dir = "/Users/ryanmckeown/High_Throughput_AF/data/move_rename/models"
end_dir = "/Users/ryanmckeown/High_Throughput_AF/data/move_rename/out"
key_path = "/Users/ryanmckeown/High_Throughput_AF/data/move_rename/key_files/03132022_missed.csv"

# Make list of model folders
models_dir = os.listdir(dir)

#Make a dictionary from the keys 
d = pd.read_csv(key_path, index_col= "pair_group").to_dict('index')

#Hardcodded PDB dictionary #### ONLY FOR 6NMI

pdb_chain_ids = {}
pdb_chain_ids = {
"P19447" : "A" , 
"P18074" : "B", 
"P32780" : "C", 
"Q92759" : "D", 
"Q13888" : "E",
"Q13889" : "F", 
"Q6ZYL4" : "G", 
"P51948" : "H"}

## Create a new name for a model

def make_dir_name(dir, input_dict, pdb_dict):
    dir_name_split = str.split(dir, sep = "_")

    chain_a_id = int(dir_name_split[0])
    chain_b_id = int(dir_name_split[1])

    chain_a_name = input_dict[chain_a_id]["Name"]
    chain_b_name = input_dict[chain_b_id]["Name"]

    #Pull the monomer chain relative to 6NMI
    a_m_id = input_dict[chain_a_id]["Monomer_id"]
    b_m_id = input_dict[chain_b_id]["Monomer_id"]

    a_pdb_id = pdb_dict[a_m_id]
    b_pdb_id = pdb_dict[b_m_id]

    new_name = str(chain_a_name + "_" + chain_b_name + "(" + a_pdb_id + "_" + b_pdb_id + ")")
    print(new_name)
    return(new_name)

## Pull top 3 models ### Make this less hard coded 

for model in models_dir: 
    name = make_dir_name(model, d , pdb_chain_ids)
    pbd_0 = str(dir + "/" + model + "/" + "ranked_0.pdb")
    pdb_1 = str(dir + "/" + model + "/" + "ranked_1.pdb")
    pdb_2 = str(dir + "/" + model + "/" + "ranked_2.pdb")

    shutil.copyfile(pbd_0, str(end_dir + "/" + name + "_" + "0" + ".pdb"))
    shutil.copyfile(pbd_0, str(end_dir + "/" + name + "_" + "1" + ".pdb"))
    shutil.copyfile(pbd_0, str(end_dir + "/" + name + "_" + "3" + ".pdb"))