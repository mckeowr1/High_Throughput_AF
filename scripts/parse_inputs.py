from Bio import SeqIO
import pandas as pd
import numpy as np


def make_fasta_dict(fasta_file): 
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file , "fasta"))
    for key in list(fasta_dict.keys()): #Rename to just uniprot ID
        fasta_id = key 
        up_id = fasta_id.split(sep = "|")[1]
        fasta_dict[up_id] = fasta_dict.pop(fasta_id)
    return(fasta_dict)

# def search_fasta(fasta_iter , id):
#     """"Searches through a Fasta iteration object for a sequence of interest
#         Returns a seq record object """
    
#     soi = (r for r in fasta_iter) if r.id 
#         if id in uniprot_id: #Check if the Fasta ID has a id of interest
#             print("Found" + uniprot_id)
#             soi = seq
#             return(soi)
#         #else: 
#             #print("Not Found")
#             #continue
#         #print(soi)
#         #return(soi)

# records = (r for r in SeqIO.parse(input_file, "sff") if r.id in wanted)


def parse_seq(seq_rec , start , stop, name): 
    stop = stop + 1 #Adj for slicing where last num excluded
    sequence = seq_rec.seq[start : stop]

    #Edit the sequence id to note what regions it contains 
    new_id = seq_rec.id + "_" + name + "_"+ str(start) + "_" + str(stop)
    
    new_rec = SeqIO.SeqRecord(sequence, id =  new_id)
    #print(new_rec)

    return(new_rec)


def parse_pair_group(dict, fasta_dict):
    group_seq_recs = []
    for record in dict:  #Looks like running into problem since length is 1 
        name = record['Name']
        monomer_id = record['Monomer_id']
        start = record['Protein_start']
        end = record['Protein_end']
        print("Processing_" + name + monomer_id)

        monomer_record = fasta_dict[monomer_id]

        subset = parse_seq(monomer_record, start , end, name)
        group_seq_recs.append(subset)

    return(group_seq_recs)
        
        #print("Processing:" + dict['Name'][i] + "_" + dict['Monomer_id'][i])
        #name = dict['Name'][i]
        #print(name)
        #id = dict['Monomer_id'][i]
        #print(id)
        #start = dict['Protein_start'][i]
        #print(start)
        #stop = dict['Protein_end'][i]
        #print(stop)

        #monomer_record = fasta_dict[id]
        #subset = parse_seq(monomer_record, start , stop, name)
        #group_seq_recs.append(subset)
    #return(monomer_record) 
    #return(group_seq_recs)