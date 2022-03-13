import pandas as pd 
from parse_inputs import * 

pairings_file = "/projects/b1059/projects/Ryan/protein_structure/High_Throughput_AF/data/03082021_ben1_tags_multimer/key_files/03082021_ben1_tags_multimer.csv"
pairings_table = pd.read_csv(pairings_file)
print(pairings_table)

fasta_file = "/projects/b1059/projects/Ryan/protein_structure/High_Throughput_AF/data/03082021_ben1_tags_multimer/key_files/ben.fa"
fasta_dict = make_fasta_dict(fasta_file)
fasta_dict.keys()
pair_scheme = ["1_5", "2_5", "3_5", "4_5"]

out_dir = "/projects/b1059/projects/Ryan/protein_structure/High_Throughput_AF/data/03082021_ben1_tags_multimer/fastas"

for pair in pair_scheme:
    p_1 = pair.split(sep= "_")[0]
    print(p_1)
    p_2 = pair.split(sep= "_")[1]
    print(p_2)

    pt1 = pairings_table.loc[pairings_table['pair_group'] == int(p_1)]
    #print(pt1)
    pd1 = pt1.to_dict('records')
    p1_seqs = parse_pair_group(pd1, fasta_dict)
    #print(p1_seqs)
    

    pt2 = pairings_table.loc[pairings_table['pair_group'] == int(p_2)]
    #print(pt2)
    pd2 = pt2.to_dict('records')
    p2_seqs = parse_pair_group(pd2, fasta_dict)
    #print()

    full_pair = p1_seqs + p2_seqs

    print(full_pair)
    out_file = out_dir + "/" + p_1 + "_" + p_2 + ".fa"
    SeqIO.write(full_pair, out_file, "fasta")


#Need to add code to write out a file key 
