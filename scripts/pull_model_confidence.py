

from parse_af import * 
import os 
import argparse

path = "/projects/b1059/projects/Ryan/protein_structure/High_Throughput_AF/data/03052022_TF/initial_models"

def listdir_nohidden(path): #Remove hidden files from the directory
    return glob.glob(os.path.join(path, '*'))

file_list = listdir_nohidden(path)

def split_at(string, sep, n):
    words = string.split(sep)
    return sep.join(words[:n]), sep.join(words[n:])

#Get list of the files
def pull_iptm(af_dir):
    file_name = os.path.basename(af_dir)
    split = split_at(file_name, "_", 1)
    #a_name = split[0]
    #b_name = split[1]
    #domain_a = split_at(a_name, "_", 1)[1]
    #domain_b = split_at(b_name, "_", 1)[1]
    
    pair_name = split[0] + "." + split[1]



    top_model = pull_top_model(af_dir)
    pkl = open_pkl(path = af_dir, model = top_model)
    iptm = float(pkl['iptm'])

    data = [pair_name , iptm]
    return data


all_interfaces = []
for file in file_list:
    if os.path.isdir(file):
        print(file)
        data = pull_iptm(af_dir = file)
        all_interfaces.append(data)
    else: 
        continue
#all_interfaces_df = pd.concat(all_interfaces)
print(all_interfaces)

import csv

with open("/projects/b1059/projects/Ryan/protein_structure/High_Throughput_AF/data/03052022_TF/initial_models/model_iptm.csv", "w+", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(all_interfaces)