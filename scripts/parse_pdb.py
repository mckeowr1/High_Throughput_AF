from telnetlib import SE
from Bio import SeqIO
import Bio.PDB
from Bio.PDB import *
import json 

def get_chain_mapping(pdb_file): 
    """Creates dictionary of what chain is actually named in PDB file"""
    parser = PDBParser()
    structure = parser.get_structure(id = "strcutre", file = pdb_file)
    chain_a_id = structure.header['compound']['1']['molecule']
    chain_b_id = structure.header['compound']['2']['molecule']
    mappings = {'chain_a': chain_a_id, 'chain_b' : chain_b_id}
    return(mappings)

def define_chain_mapping(pdb_file, fasta_a_path, fasta_b_path):
    """Assigns chain ID based on FASTA A and FASTA B Names"""
    seq_a = SeqIO.read(fasta_a_path, "fasta")
    seq_b = SeqIO.read(fasta_b_path, "fasta") 
    name_a = seq_a.id
    name_b = seq_b.id

    uni_a = name_a.split('|')[1] #get uniprot name
    uni_b = name_b.split('|')[1]

    mappings = {'chain_a': uni_a, 'chain_b' : uni_b}
    return(mappings)



def get_crystal_resiudes(pdb_chain):
    """Pull list of residue numbers included in crystal stucture"""
    residues = []
    for residue in pdb_chain:
        if residue.id[0] == 'W': #Filter out any water molecules
            continue
        else: 
            id = residue.id[1] - 1 #Subtract to make the index 0 based 
            residues.append(id)
    return(residues)

def compare_fasta_crystal(fasta_path , crystal_residues):
    """Return a list of residues that do not have a predicted
    postion in the PDB file we are comparing to"""
    
    for record in SeqIO.parse(fasta_path, "fasta"): 
        #Should have a check here that the fasta only contains one record
        fasta_seq = record.seq #Create a sequence object
    
    residue_range = range(len(fasta_seq))
    fasta_residues = list(residue_range)
    #return(fasta_residues)

    #Get the differences of the two list s 
    non_matching = list(set(fasta_residues) - set(crystal_residues))

    return(non_matching)

def pull_pfam_domains(pfam_json_path, pad_decimal , pad = False):
    """Reads in Pfam JSON and produces a domain List 
    [UNIPROT_ID, DOMAIN_NAME, START, END]
    """
    #Get Protein Name
    domains_json = open(pfam_json_path)
    domains = json.load(domains_json)

    protein = domains['metadata']['accession']

    regions = domains['regions']
    domain_data = []
    for domain in regions: 
        #print(domain)
        data = []
        name = domain['text']
        if pad == True:
            start = domain["start"]
            end = domain["end"]
            print("Pre-Pad", name)
            print("start:", start, " end:", end)
            length = end - start
            pad_size = round(length * pad_decimal)
            start = domain['start'] - pad_size
            if start <= 0: 
                start = 1
            end = domain['end'] + pad_size
            print("start:", start, " end:", end)
        else:
            start = domain['start']
            end = domain['end']
        #print(name, start, end)
        data.append(protein)
        data.append(name)
        data.append(start)
        data.append(end)
        domain_data.append(data)
    return(domain_data)

def extract_pdb_domains(pdb_path, chain_id, domains, out_dir):
    """
    Takes in a  PDB structure and Domain Data
    Returns PDB file for just the fragment
    Outputs to Specific Directory
    """
    #Load structure
    parser = PDBParser()
    structure = parser.get_structure(id = "strcutre", file = pdb_path)

    for domain in domains:
        protein = domain[0] 
        name = domain[1]
        start = domain[2]
        end = domain[3]
        file_name = out_dir + "/" + protein + '_' + name + ".pdb"
        Bio.PDB.extract(structure = structure, chain_id = chain_id , start = start , end = end , filename = file_name)


def load_chain(pdb_path, chain_id):
    """Load chain from Domain PDB Path - Must note original chain"""
    parser = PDBParser()
    structure = parser.get_structure(id = "structure", file = pdb_path)
    model = structure[0] 
    #print(model)
    chain = model[chain_id]
    #chain = model['A']
    return(chain)



def pull_pdb_domain_interfaces(domain_a_pdb_path , domain_b_pdb_path, distance_max):
    """ Load two PDB Domain Files and pull their interfacing residues """



    a = load_chain(domain_a_pdb_path, chain_id='A') #Must match with what chain they were subset from
    b = load_chain(domain_b_pdb_path, chain_id='B')
    
    #print(len(a))
    #print(len(b))
    

    i_residue_a = []
    i_residue_b = []
    dist = []
    for residue_a in a: 
        #print(residue_a)
        for residue_b in b: 
            try:
                #print(residue_a, residue_b)
                distance = residue_a['CA'] - residue_b['CA'] 
                #print(distance)
            except KeyError:
                print("distance calculation error")
                ## no CA atom, e.g for H_NAG
                continue
            if distance <= distance_max: 
                index_a = residue_a.id[1] #Pull Residue index from the residue obj
                index_b = residue_b.id[1]

                i_residue_a.append(index_a)
                i_residue_b.append(index_b)
                dist.append(distance)
                #print(distance)
    data = [i_residue_a , i_residue_b , dist]
    return(data)

