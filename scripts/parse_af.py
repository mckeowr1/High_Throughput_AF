import numpy as np
from Bio.PDB import *
import Bio.PDB
import pandas as pd
import json
import sys
import glob
import itertools
import string


def open_rankings(path): 
    """Depreciated - contained in pull top model"""
    file_list = glob.glob(path + "/ranking_debug.json")
    
    print(file_list)
    if len(file_list) > 1: 
        print("There is more than 1 rankings file in:" + path)
        sys.exit()
    rankings_file = file_list[0]
    return(rankings_file)


def pull_top_model(file):
    file_list = glob.glob(file + "/ranking_debug.json")
    #print(file_list)
    if len(file_list) > 1: 
        print("There is more than 1 rankings file in:" + file)
        sys.exit()
    rankings_file = file_list[0]
    f = open (rankings_file , "r")
    data = json.load(f)
    top_model = data['order'][0]
    print("The top model is: " + top_model)
    return(top_model)


def open_pkl(path, model):
    name = "/result_" + model + ".pkl"
    file_list = glob.glob(path + name)
    if len(file_list) > 1: 
        print("There is more than 1 models file in: " + path)
        sys.exit()
    pkl_file = file_list[0]
    data = pd.read_pickle(pkl_file)
    return(data)

def open_features(path): 
    name = "/features.pkl"
    file_list = glob.glob(path + name)
    if len(file_list) > 1: 
        print("There is more than 1 features file in: " + path)
        sys.exit()
    features_file = file_list[0]
    data = pd.read_pickle(features_file)
    return(data)

def open_pdb(path, model): #Maybe want to update this function to actually parse the file 
    name = "/relaxed_" + model + ".pdb"
    file_list = glob.glob(path + name)
    if len(file_list) > 1: 
        print("There is more than 1 models file in: " + path)
        sys.exit()
    pdb  = file_list[0]
    
    #Read in as a PDB object
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb)
    model = structure[0] #We want the model out of the structure file
    return(model)

#Pull the scoring metrics we need out of the pkl files 
def pull_array(pkl , data_name): 
    data = pkl[data_name]
    return(data)

def pull_interfaces(chain_a, chain_b, distance_max): 
    # chain_a = pdb['A'] 
    # chain_b = pdb['B'] 
    i_residue_a = []
    i_residue_b = []
    dist = []
    for residue_a in chain_a:
        for residue_b in chain_b:
            try:
                distance = residue_a['CA'] - residue_b['CA']
            except KeyError:
            ## no CA atom, e.g for H_NAG
                continue
            if distance <= distance_max: 
                index_a = residue_a.id[1] #Pull Residue index from the residue obj
                index_b = residue_b.id[1]

                i_residue_a.append(index_a)
                i_residue_b.append(index_b)
                dist.append(distance)
    data = [i_residue_a , i_residue_b , dist]
    return(data)

def convert_chains(chain_a, chain_b , chain_index_map, interfaces_data , seq_ids):
    """
    Function to convert 1 based residue postion to 0 based residue position
    
    Since the PDB Residue ID is specific to a chain ex) Chain_A is residues 1-100 and Chain_B is residues 1-100
    we convert to the Numpy PAE array standard of labeleing Chain_A is residue 0-99 and Chain_B is resiudes 100-199
    
    The function does this in a way using seqIDs to be able to handle protein fragments of different lengths 
    
    Args:
        chains: chain from parsed pdb file
        chain_index_map: a dictionary that assoicates chain id (from model) to entity id 
            (from af)
        interfaces_data: interfaces array 
        seq_ids: af output that shows where chain is in the full sequence
    Returns: 
        interfaces data where PDB residues have converted to proper af index

    """
    chain_a_base = chain_a.full_id[2] # Get the base id
    chain_b_base = chain_b.full_id[2]
    
    a_entity_id = chain_index_map.get(chain_a_base)
    b_entity_id = chain_index_map.get(chain_b_base)

    if a_entity_id == 1:
        def convert_chain_a(chain_a, a_entity_id):
            chain_a_id = chain_a -1 
            return(chain_a_id)
    else: 
        def convert_chain_a(chain_a, a_entity_id): 
            chain_a_id = chain_a # -1 # Convert to 0 based
            sub = np.where(seq_ids == a_entity_id)[0][0]  #Get the value we need to subtract from all chain, B
            new_id = chain_a_id + (sub -1)
            return(new_id)
    if b_entity_id == 1:
        def convert_chain_b(chain_b, b_entity_id):
            chain_b_id = chain_b -1 
            return(chain_b_id)
    else: 
        def convert_chain_b(chain_b, b_entity_id): 
            chain_b_id = chain_b # -1 # Convert to 0 based
            sub = np.where(seq_ids == b_entity_id)[0][0]  #Get the value we need to subtract from all chain, B
            new_id = chain_b_id + (sub -1)
            return(new_id)
    
    #Loop through all elements of interfaces_data
    x = range(0 , len(interfaces_data[1]))

    for i in x: 
        interfaces_data[0][i] = convert_chain_a(interfaces_data[0][i], a_entity_id)
        interfaces_data[1][i] = convert_chain_b(interfaces_data[1][i], b_entity_id)
    return(interfaces_data)

def check_pae(interfaces_data , pae): 
    x = range(0 , len(interfaces_data[1]))
    pae_scores = []
    for i in x: 
        residue_a = interfaces_data[0][i]
        residue_b = interfaces_data[1][i]
        pae_score = pae[residue_a][residue_b]
        pae_scores.append(pae_score)
    interfaces_data.append(pae_scores)
    return(interfaces_data)

def check_plddt(interfaces_data , plddt): 
    x = range(0 , len(interfaces_data[1]))
    avg_plddt = []
    max_plddt = []
    for i in x: 
        residue_a = interfaces_data[0][i]
        residue_b = interfaces_data[1][i]
        avg = (plddt[residue_a] + plddt[residue_b]) / 2
        
        avg_plddt.append(avg)

        mp = max(plddt[residue_a] , plddt[residue_b]) 
        max_plddt.append(mp)
    interfaces_data.append(avg_plddt)
    interfaces_data.append(max_plddt)
    return(interfaces_data)

def pull_domain_start(pfam_json_path, domain_name):
    """Pulls the start of a single domain"""
    domains_json = open(pfam_json_path)
    domains = json.load(domains_json)
    regions = domains['regions']
    for domain in regions: 
        if domain['text'] == domain_name: 
            start = domain["start"]
            #end = domain["end"]
            #print(start)
            #print(end)
        else:
            continue
    return(start)

def pull_chain_start(key, chain):
    """Pulls the start of the chain from the key for each input"""
    chain_base = chain.full_id[2] # Get the base id
    chain_name = chain_base.split(".")[0]
    chain_monomer = chain_base.split(".")[1]
    chain_df = key[(key['Name'] == chain_name) & (key['Monomer_id'] == chain_monomer)]
    chain_start = int(chain_df["Protein_start"])
    return(chain_start)

# def create_chain_key(key_path): 
#     df = pd.read_csv(key_path)
#     bet = list(string.ascii_uppercase) #create an alphabetvector 
#     chains= bet[0:len(df)]
#     df["chain_id"] = chains
#     return(df)

def create_file_dict(key_file, pair_id):
    """Creates a dictionary for an alphafold output folder based on input keys 
    - Only works for two chain output 
    - dictionary key is the chain ID
    """
    df = pd.read_csv(key_file, index_col= "pair_group") #Read in with the pair group id
    chain_a_id = int(pair_id.split(sep = "_")[0])
    chain_b_id = int(pair_id.split(sep = "_")[1])
    pair_df = df.filter(items  = [chain_a_id, chain_b_id], axis= 0) #filter the df to proteins in the file
    bet = list(string.ascii_uppercase) #create an alphabetvector 
    chains= bet[0:len(pair_df)]
    pair_df["chain_id"] = chains #add chain id to the dataframe
    d = pair_df.set_index('chain_id').T.to_dict() #Create a dictionary from df
    return(d)


def convert_full_protein(chain_a, chain_b, chain_index_map, interfaces_data, seq_id, chain_a_start, chain_b_start): 

    """Converts alphafold domain interface data to full protein index
    Args: 
        chains: pbd model chains 
        chain_index_map: a dictionary that assoicates chain id (from model) to entity id 
            (from af)
        interfaces_data: interfaces that have been coverted to AF index
        seq_ids: af output that shows where chain is in the full sequence
        chain_start: where this fragment starts relative to the whole protein
    
    """
    chain_a_base = chain_a.full_id[2] # Get the base id
    chain_b_base = chain_b.full_id[2]
    
    a_entity_id = chain_index_map.get(chain_a_base) #Pull entity ID from the dictionary 
    b_entity_id = chain_index_map.get(chain_b_base)

    if a_entity_id == 1:
        def convert_chain_a(chain_a, start, entity_id): 
            chain_id = chain_a + start 
            return(chain_id)
    else:
        def convert_chain_a(chain_a, start, entity_id):
            sub = np.where(seq_id == entity_id)[0][0] #Get the value we need to subtract from all chain, B
            chain_id = chain_a - sub + start 
            return(chain_id)
    if b_entity_id == 1:
        def convert_chain_b(chain_b , start, entity_id): 
            chain_id = chain_b + start 
            return(chain_id)
    else: 
        def convert_chain_b(chain_b, start, entity_id):
            sub = np.where(seq_id == entity_id)[0][0] #Get the value we need to subtract from all chain, B
            chain_id = chain_b - sub + start 
            return(chain_id)

    x = range(0 , len(interfaces_data[1])) 
    for i in x: 
        interfaces_data[0][i] = convert_chain_a(interfaces_data[0][i], chain_a_start, a_entity_id )
        interfaces_data[1][i] = convert_chain_b(interfaces_data[1][i], chain_b_start, b_entity_id)
    return(interfaces_data)

def rename_chains(pdb, model_dict):
    """Renames chains in model to the user fragment name and the full uniprot ID"""
    for chain in pdb.child_list: 
        chain_base = chain.full_id[2] # Get the base id
        chain_name = model_dict[chain_base]['Name']
        chain_monomer = model_dict[chain_base]['Monomer_id']
        #chain_name = key.loc[key['chain_id'] == chain_base, 'Name'].values[0] #pull sting out of df
        #chain_monomer = key.loc[key['chain_id'] == chain_base, 'Monomer_id'].values[0]
        chain_id = str(chain_name + "." + chain_monomer)
        print(chain_id)
        chain.id = chain_id
    
def get_chain_pairs(model):
    """Pairwise Combination of Chains
        Args: Parsed PDB Model 
        Returns: List of chain pairs [(chain_a, chain_b)]
    """
    chains = model.child_list
    pairs = list(itertools.combinations(chains, 2))
    return(pairs) #list of pairs first pair can be accessed by p[1].. p[n]

#Create a class for chain interface
def map_chain_index(model): 
    """Create a distionary based on chain id and entity id for later conversions"""
    #Could be improved with model.child_dict
    ps = model.child_list #List of chains in model
    base_id = []
    for chain in ps: #Get the base ID which is a string we can recognize later
        #print(chain)
        f_id = chain.full_id
        #print(f_id)
        base_id.append(f_id[2])
    n = list(range(1, len(base_id)+1))

    id_to_ent = dict(zip(base_id, n)) 
    return(id_to_ent) #Can be searched later with a id_to_ent.get("B") call

class chain_interfaces:
    def __init__(self, chain_a, chain_b,  array):
        self.chain_a = chain_a 
        self.chain_b = chain_b 
        self.array = array 
    
    def print_pair(self):
        print(self.chain_a.full_id[2])








