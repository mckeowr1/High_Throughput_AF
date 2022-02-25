import os 
import argparse
from parse_pdb import *
import csv
import glob
import itertools
import pandas as pd


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

parser = argparse.ArgumentParser(description="Process a PDB file for comparison with domain AF prediction")
parser.add_argument('-p', '--pdb_file', action = 'store', nargs= 1, type=file_path)
parser.add_argument('-a', '--chain_A_fasta', action = 'store', nargs= 1, type=file_path, 
                    help = 'must be a fasta downloaded from uniprot with header format:  sp|UNIPROT_ID|ALIAS')
parser.add_argument('-b', '--chain_B_fasta', action= 'store', nargs = 1 , type=file_path,
                    help = 'must be a fasta downloaded from uniprot with header format:  sp|UNIPROT_ID|ALIAS')
parser.add_argument('-d', '--max_distance', action = 'store', nargs= 1, type=int)
parser.add_argument('-da', '--chain_A_domains', action = 'store', nargs= 1, type=file_path,
                    help="must be .json downloaded from PFAM database. Named UNIPROTID.json")
parser.add_argument('-db', '--chain_B_domains', action = 'store', nargs= 1, type=file_path)
parser.add_argument('-o', '--output_dir', action= 'store', nargs = 1 , type=dir_path)
parser.add_argument('-r', '--pad_size', action = 'store', nargs = 1 ,  type=int, help = "The percent of domain length you would like to pad with, ex 10% = 10")

args = parser.parse_args()

pdb_file = args.pdb_file[0]
a_fasta_file = args.chain_A_fasta[0]
b_fasta_file = args.chain_B_fasta[0]
a_domains_file = args.chain_A_domains[0] 
b_domains_file = args.chain_B_domains[0]
distance_max  = args.max_distance[0]
out_dir = args.output_dir[0]
pad_size = args.pad_size[0] 

if pad_size > 0: 
    pad_size = pad_size / 100 #convert pad size to decimal

#print(a_fasta_file)
#print(type(a_fasta_file))

#### Determine ID of Each chain

mappings = define_chain_mapping(pdb_file, a_fasta_file, b_fasta_file)


#### Output Residues not included in PDB structure

chain_a = load_chain(pdb_path = pdb_file, chain_id = 'A')
chain_b = load_chain(pdb_path = pdb_file, chain_id = 'B')

a_crystal = get_crystal_resiudes(chain_a)
b_crystal = get_crystal_resiudes(chain_b)

a_absent = compare_fasta_crystal(a_fasta_file , a_crystal)
b_absent = compare_fasta_crystal(b_fasta_file, b_crystal)

a_absent.insert(0 , "A") #Add some type of ID to the list
b_absent.insert(0 , "B")

absent = [a_absent , b_absent]

file_name = str(out_dir + "/absent_residues.txt")
absent_file = open(file_name, "w+")
writer = csv.writer(absent_file)
writer.writerows(absent)
absent_file.close()

#### Pull pfam domains from the PDB structure

if pad_size == 0:
    a_domains = pull_pfam_domains(a_domains_file, pad_decimal= pad_size, pad = False)
    b_domains = pull_pfam_domains(b_domains_file, pad_decimal= pad_size, pad = False)
else:
    print("calculating pad size")
    print(pad_size)
    a_domains = pull_pfam_domains(a_domains_file, pad_decimal = pad_size, pad = True)
    b_domains = pull_pfam_domains(b_domains_file, pad_decimal = pad_size, pad = True)


extract_pdb_domains(pdb_path = pdb_file, chain_id = 'A', domains = a_domains, out_dir= out_dir)
extract_pdb_domains(pdb_path = pdb_file, chain_id = 'B', domains = b_domains, out_dir= out_dir)

#### Now lets pull interfaces from the PDB files

#Get a list of files we want to work with

#search_name = str()
file_list = glob.glob(out_dir + "/*.pdb")
try: 
    file_list.remove(pdb_file) #Get rid of the original pdb file
except:
    print("PDB_file not in out dir")
#Split the list into lists based on original protein



print(mappings['chain_a'])
chain_a_files = []
chain_b_files = []
problem_files = []
for file in file_list: 
    fname = os.path.basename(file)
    protein_name = fname.split('_')[0]
    #Check to see if the file contains anything 
    parser = PDBParser()
    structure = parser.get_structure(id = "structure", file = file)
    if len(structure) == 0:
        continue
    if protein_name == mappings['chain_a']:
        chain_a_files.append(file)
    else:
        chain_b_files.append(file)
#print(chain_a_files)
#print(chain_b_files)



#Generate all pairwise combinations of files

file_pairs = list(itertools.product(chain_a_files , chain_b_files))

#Pull interfaces from the data
all_pairs_list = [] #empty list to store all pairs DF
for pair in file_pairs: 
    file_a = pair[0]
    file_b = pair[1]
    
    fname_a = os.path.basename(file_a)
    fname_b = os.path.basename(file_b)

    p_d_a = os.path.splitext(fname_a)[0] #get protein and domain name
    p_d_b = os.path.splitext(fname_b)[0]


    interfaces = pull_pdb_domain_interfaces(file_a , file_b , distance_max = distance_max)

    #Transform into pair dataframe 
    domain_a = [fname_a] * len(interfaces[0]) #Produce a name vector the same length as the number of interfaces
    domain_b = [fname_b] * len(interfaces[0])
    residue_a = interfaces[0]
    residue_b = interfaces[1]
    distance = interfaces[2]

    pair_df = pd.DataFrame(list(zip(domain_a, domain_b, residue_a, residue_b, distance)), columns= ["domain_a","domain_b",'residue_a', 'residue_b', 'distance'])
    all_pairs_list.append(pair_df)
all_pairs_df = pd.concat(all_pairs_list)
    
print(all_pairs_df)

df_name = out_dir + "/all_PDB_domain_interfaces.csv"

all_pairs_df.to_csv(path_or_buf = df_name , index = False)

