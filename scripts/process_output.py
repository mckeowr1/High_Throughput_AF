from posixpath import basename
import sys , getopt
import argparse
import glob
import os
import json
import re

from sympy import residue
from parse_af import * 

print('Number of arguments:', len(sys.argv), 'arguments.')
#print('Argument List:', str(sys.argv))
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
def file_path(string): 
    if os.path.isfile(string): 
        return string
    else: 
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser(description="Process an Alphafold Output Directory")
parser.add_argument('-p', '--directory', action = 'store', nargs= 1, type=dir_path, help = "Path to the output directory containg pairwise af domain predictions")
parser.add_argument('-d', '--max_distance', action= 'store', nargs = 1 , type=int)
#parser.add_argument('-da', '--chain_A_domains', action = 'store', nargs= 1, type=file_path,
                   # help="must be .json downloaded from PFAM database. Named UNIPROTID.json")
#parser.add_argument('-db', '--chain_B_domains', action = 'store', nargs= 1, type=file_path)
parser.add_argument('-k', '--directory_key_file', action = 'store', nargs= 1, type=file_path, help = "The path to the key file that breaks down what is in this pairwise dir")

args = parser.parse_args()
dir = args.directory[0]
distance_max  = (args.max_distance[0])
#a_domains_file = args.chain_A_domains[0] 
#b_domains_file = args.chain_B_domains[0]
key_file = args.directory_key_file[0]


#Get the names of the files in the directory 

def listdir_nohidden(path): #Remove hidden files from the directory
    return glob.glob(os.path.join(path, '*'))

file_list = listdir_nohidden(dir)

dir_key = pd.read_csv(key_file)

def process_domain_pred(file, max_distance, dir_key):
    #print(file)
    file_name = os.path.basename(file)
    print("Processing AF out dir:", file_name)
    key_file = dir_key.loc[dir_key['dir_name'] == file_name, 'key_file'].values[0] #Get the key file path for the file
    


    top_model = pull_top_model(file) 
    print("The top model is" , top_model)
    
    #Define the two domains being analyzed
    # def split_at(string, sep, n):
    #     words = string.split(sep)
    #     return sep.join(words[:n]), sep.join(words[n:])
    
    # split = split_at(file_name, ".", 1)
    # a_name = split[0]
    # b_name = split[1]

    # domain_a = split_at(a_name, "_", 1)[1]
    # #domain_a = re.sub(r"_\d+$", "", domain_a)
    # print(domain_a)

    # domain_b = split_at(b_name, "_", 1)[1]
    # #domain_b = re.sub(r"_\d+$", "", domain_b)
    # print(domain_b)
    # #Open the top model files
    #Open a the alpha fold outputs
    pkl = open_pkl(file, top_model)
    ft = open_features(file)
    #Pull data out of pkls
    plddt = pull_array(pkl , 'plddt')
    pae = pull_array(pkl , "predicted_aligned_error")
    seq_id = pull_array(ft , "entity_id") 

    chain_key = create_chain_key(key_file) #Give each chain an ID based on its order in FASTA - corresponds to PDB chain in model
    pdb = open_pdb(file, top_model)


    #Rename Chains here - Based on ID in key file
    rename_chains(pdb, chain_key)

    #Create a List of the all the chain pairs
    chain_pairs = get_chain_pairs(pdb)

    #Create a Chain Dictionary
    chain_dict = map_chain_index(pdb)

    all_chains = []
    for pair in chain_pairs: 
        chain_a = pair[0]
        chain_b = pair[1]
        print("Pulling interfaces for:" , chain_a.id, chain_b.id)
        pair_id = str(chain_a.id + chain_b.id)

        i = pull_interfaces(chain_a = chain_a, chain_b = chain_b, distance_max = 6) #interfaces
        c = convert_chains(chain_a, chain_b, chain_dict , i , seq_id)
        i_pae = check_pae(c , pae)
        i_plddt = check_plddt(i_pae , plddt)
        chain_a_start = pull_chain_start(chain_key, chain_a)
        chain_b_start = pull_chain_start(chain_key, chain_b)
        #print(chain_a, "starts at:" ,chain_a_start)
        #print(chain_b, "starts at:" ,chain_b_start)
        final = convert_full_protein(i_plddt, chain_a_start, chain_b_start, seq_id)
        
        #add chain ids to pair names 
        id_list = [pair_id] * len(final[0])
        final.append(id_list)    

        #Convert the list of lists to DF    
        chain_df = pd.DataFrame(final)
        chain_df = chain_df.transpose()
        chain_df.colnames = ["residue_a", "residue_b", "distance", "pae", "avg_plddt", "max_plddt", "domain_pair"]
        all_chains.append(chain_df)
    all_chains_df = pd.concat(all_chains)
    return(all_chains_df)

all_interfaces = []
for file in file_list:
    print(file)
    pair_data = process_domain_pred(file, max_distance = distance_max, dir_key = dir_key)
    all_interfaces.append(pair_data)
all_interfaces_df = pd.concat(all_interfaces)

df_name = dir + "/all_af_domain_interfaces.csv"
all_interfaces_df.colnames = ["residue_a", "residue_b", "distance", "pae", "avg_plddt", "max_plddt", "domain_pair"]


all_interfaces_df.to_csv(path_or_buf = df_name , index = False)









# #Pull the scoring metrics we need out of the pkl files 

# plddt = pull_array(pkl , 'plddt')
# pae = pull_array(pkl , "predicted_aligned_error")
# seq_id = pull_array(ft , "entity_id") 



# # Now We'll start pulling interface data

# #IMPORT THESE FUNCTIONS





# interfaces = pull_interfaces(pdb, distance_max)
# chains = convert_chains(interfaces, seq_id)
# scored_pae = check_pae(chains, pae)
# scored_plddt = check_plddt(chains, plddt)


#Write out the file - inefficent at the moment

# file = open("test.txt", "w")
# writer = csv.writer(file)

# names = ['residue_a', 'residue_b', 'pae', 'avg_plddt', 'max_plddt']
# writer.writerow(names)

# for i in range(len(scored_plddt[1])):
#     writer.writerow([scored_plddt[0][i], 
#                     scored_plddt[1][i], #, #Write resiude 2 
#                     scored_plddt[2][i], #Write pair pae
#                     scored_plddt[3][i], #Write pair avg plddt
#                     scored_plddt[4][i]]) #Write pair max plddt


# file.close()  





## Need to modularize at some point but here we start to parse for interacting residues










