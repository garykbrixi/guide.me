import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

with open('prime_configs.csv', mode='r') as infile:
    reader = csv.reader(infile)
    config_dict = {rows[0]:rows[1].split(' ') for rows in reader}

folder_names = config_dict['folder_names']
norm_ext_list = config_dict['ctrl_ext_list']
data_ext_list = config_dict['data_ext_list']
norm_trunc_list = config_dict['ctrl_trunc_list']
data_target = config_dict['data_target'][0]
ctrl_target = config_dict['ctrl_target'][0]
data_trunc_list = config_dict['data_trunc_list']

failures = []
file_pos = []
sample_names = []
for fname in folder_names:
    file_pos.extend([fname + s for s in os.listdir(os.getcwd()+'/data/'+fname)])
    sample_names.extend(os.listdir(os.getcwd()+'/data/'+fname))

names = []
for i in range(len(file_pos)):
    filename = file_pos[i]
    sample_name = sample_names[i]
    if not '.DS_Store' in filename:
        try:
            df = pd.read_csv('data/{}/CRISPResso_quantification_of_editing_frequency.txt'.format(filename), sep='\t')
            df.insert(0, 'Batch', filename)
            names.append(sample_name[14:])
        except:
            print("FAILED: {}",format(sample_name[14:]))
            failures.append(sample_names)
            continue

merged_dict = {}
for name in names:
    merged_dict[name] = {}
    found = False
    for fname in folder_names:
        try:
            data = pd.read_csv('data/'+fname + '/CRISPResso_on_{}/CRISPResso_quantification_of_editing_frequency.txt'.format(name), sep='\t')
            for index, row in data.iterrows():
                merged_dict[name]["Reads_total"] = row["Reads_total"]
                merged_dict[name][row["Amplicon"]] = 100*row["Reads_aligned"]/row["Reads_total"]
                if(row["Amplicon"]==data_target or row["Amplicon"]==ctrl_target):
                    merged_dict[name][row["Amplicon"]+"_unmodified"] = 100*row["Unmodified"]/row["Reads_total"]
            found = True
        except:
            continue
        if(not found):
            print("FAILED: {}",format(name))

merged_Data = pd.DataFrame.from_dict(merged_dict)

df_concat = pd.concat((merged_Data['noPE2nogd1rep1_S114'], merged_Data['noPE2nogd2rep1_S115']))
by_row_index = df_concat.groupby(df_concat.index)
normalizer_full = by_row_index.mean()

large_dict = {}

#match appropraite controls with sample type. Will be transferred to dictionary
names = []
names1 = pd.read_csv('s1_r6CAGi_names.csv')
names2 = pd.read_csv('s2_r6CAGi_names.csv')
names = pd.concat([names1,names2])
names["index1"] = names["skip"].str[13:]
names.sort_values(by="index1", inplace=True)
del names["index1"]


for index, name in names.iterrows():
    name = name[0]
    name = name.replace("-", "")
    if name not in failures:
        large_dict[name] = {}
        if(merged_Data[name]["Reads_total"] < 2000):
            large_dict[name]["filtered"] = "Fail"
        elif(merged_Data[name]["Reads_total"] < 10000):
            large_dict[name]["filtered"] = "Warning"
        else:
            large_dict[name]["filtered"] = "Pass"
        large_dict[name]["reads"]=merged_Data[name]["Reads_total"]
        large_dict[name]["edits_unmodified"] = merged_Data[name][data_target+"_unmodified"]-normalizer_full[ctrl_target+"_unmodified"]
        large_dict[name]["edits_modified"] = (merged_Data[name][data_target] - merged_Data[name][data_target+"_unmodified"])- (normalizer_full[ctrl_target] - normalizer_full[ctrl_target+"_unmodified"])
        large_dict[name]["extensions"] = merged_Data[name][data_ext_list].to_numpy().sum()-normalizer_full[norm_ext_list].to_numpy().sum()
        large_dict[name]["truncation"] = merged_Data[name][data_trunc_list].to_numpy().sum()-normalizer_full[norm_trunc_list].to_numpy().sum()
        large_dict[name]["indels"] = large_dict[name]["extensions"]+large_dict[name]["truncation"]

big_frame = pd.DataFrame.from_dict(large_dict)

big_frame.to_csv('{}.csv'.format(config_dict['save_name'][0]), index = True)
