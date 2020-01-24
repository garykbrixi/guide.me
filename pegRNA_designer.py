import csv
import pandas as pd
import numpy as np
import sys

# highlights = ['\033[95m','\033[94m','\033[92m','\033[93m','\033[91m','\033[0m','\033[1m','\033[4m']

highlights = ['\33[31m','\33[32m','\33[33m','\33[34m','\33[35m','\33[36m', '\33[31m','\33[32m','\33[33m','\33[34m','\33[35m','\33[36m', '\33[31m','\33[32m','\33[33m','\33[34m','\33[35m','\33[36m', '\33[31m','\33[32m','\33[33m','\33[34m','\33[35m','\33[36m', '\33[31m','\33[32m','\33[33m','\33[34m','\33[35m','\33[36m']

config = pd.read_csv('pegRNA_config.csv')

spacer_output = {"name": [], "seq": []}
peg_output = {"name": [], "seq": []}
nick_output = {"name": [], "seq": []}

with open('pegRNA_config.csv', mode='r') as infile:
    reader = csv.reader(infile)
    config_dict = {rows[0]:rows[2] for rows in reader}

# corrected:
seq = config_dict['sequence']
edits_conf = config_dict['edits'][:-1].split('.')

mindist = 0
spacer_size = 20

forward_conf = config_dict['forward']
reverse_conf = config_dict['reverse']

f_nick_conf = config_dict['forward_nick']
r_nick_conf = config_dict['reverse_nick']

premade_f_spacer = config_dict['premade_forward_spacer']
premade_r_spacer = config_dict['premade_reverse_spacer']

homology_maxsize = int(config_dict['homology_maxsize'])
homology_stepsize = int(config_dict['homology_stepsize'])
homology_variants = int(config_dict['homology_variants'])

annealing_maxsize = int(config_dict['annealing_maxsize'])
annealing_stepsize = int(config_dict['annealing_stepsize'])
annealing_variants = int(config_dict['annealing_variants'])

PAM_setting = config_dict['PAM']
nick_PAM_setting = config_dict['nick_PAM']
nick_dist_conf = int(config_dict['nick_max_dist'])


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '[': ']', ']': '['}

def split(word):
    return [char for char in word]

def reverse_complement(string):
    return "".join(complement.get(base, base) for base in reversed(string))

def analyze_5to3(seq):
    global spacer
    global filename
    global direcitonality
    global mindist
    global spacer_size
    global homology_maxsize
    global homology_stepsize
    global homology_variants
    global annealing_maxsize
    global annealing_stepsize
    global annealing_variants

    if(PAM_setting == 'NGG'):
        PAM_id_conf = 'GG'

    edit_start = seq.find('[')
    edit_end = seq.find(']')

    i = edit_start-1-mindist
    pams = []
    while i > 2:
        if seq[i-2:i] == PAM_id_conf:
            pams.append(i-3)
        i = i-1

    #spacer
    if((directionality == 'f1' and premade_f_spacer == 'n') or (directionality == 'r1' and premade_f_spacer == 'n')):
        visual = seq[0:edit_start+1]
        visual = split(visual)

        #create bucket of spacers
        all_spacers = []
        for i in range(len(pams)):
            tmpspacer_end = pams[i]
            pam_string = visual[tmpspacer_end - 20:pams[i]]
            tmpspacer_start = tmpspacer_end - 20
            all_spacers.append("".join(visual[tmpspacer_start:pams[i]]))

        for i in range(len(pams)):
            tmpspacer_end = pams[i]
            tmpspacer_start = tmpspacer_end - 20
            if(tmpspacer_start < 0):
                continue
            else:
                visual.insert(tmpspacer_end, str(i+1)+"\033[0m")
                if(i+1 < len(pams)):
                    if(tmpspacer_start < pams[i+1]):
                        visual.insert(pams[i+1], highlights[i])
                    else:
                        visual.insert(tmpspacer_start, highlights[i])
                else:
                    visual.insert(tmpspacer_start, highlights[i])

        string = "".join(visual)
        print(string)
        unsatisfactory = True

        while unsatisfactory:
            try:
                spacer_end = pams[i]
            except:
                print("Error. Check that edit tract is marked")
                sys.exit()
            while True:
                x = input("Input number to select corresponding spacer: ")
                try:
                    x = int(x)
                except ValueError:
                    print("Please input integer")
                    continue
                else:
                    if(len(all_spacers) < x or x < 1):
                        print("Number out of range")
                        continue
                    else:
                        break
            x = int(x)
            x = x - 1
            if (len(all_spacers) < x + 1 or x < 0):
                print("Error. Invalid number")
                continue
            spacer = all_spacers[x]
            if(spacer[0]=="G"):
                print("Added: ", spacer)
                Mod_spacer = spacer
                unsatisfactory = False
            else:
                print("Added: (+G)", spacer)
                Mod_spacer = "G" + spacer
                unsatisfactory = False

    if(directionality == 'f1' and premade_f_spacer != 'n'):
        spacer_start = seq.find(premade_f_spacer)
        spacer_end = spacer_start+len(premade_f_spacer)
        spacer = premade_f_spacer
        if(spacer[0] != "G"):
            print("(G will be added)")
            Mod_spacer = "G" + spacer
        else:
            Mod_spacer = spacer
    elif(directionality == 'r1' and premade_r_spacer != 'n'):
        spacer_start = seq.find(premade_r_spacer)
        spacer_end = spacer_start+len(premade_r_spacer)
        spacer = premade_r_spacer
        if(spacer[0] != "G"):
            print("(G will be added)")
            Mod_spacer = "G" + spacer
        else:
            Mod_spacer = spacer

    #annealing
    rc_spacer = reverse_complement(spacer)
    max_anneal = rc_spacer[3:]

    unsatisfactory = True
    annealing_strands = []
    annealing_lengths = []
    while unsatisfactory:
        i = 0
        annealing_tsize = annealing_maxsize
        while i < annealing_variants:
            while(max_anneal[annealing_tsize-1] == "C"):
                annealing_tsize -= 1
            current_anneal = max_anneal[0:annealing_tsize]
            annealing_strands.append(current_anneal)
            annealing_lengths.append(annealing_tsize)
            annealing_tsize -= 1
            i += 1
        i = 0
        while i < annealing_variants:
            print("length ", annealing_lengths[i], " ", annealing_strands[i])
            i += 1
        while True:
            x = input("Accept annealing strands? (y/n): ")
            if(x == "y"):
                unsatisfactory = False
                break
            elif(x == "n"):
                while True:
                    x = input("Set new max annealing strand length? (<17): ")
                    try:
                        x = int(x)
                    except ValueError:
                        print("Please input integer")
                        continue
                    else:
                        if(x > 17 or x < 1):
                            continue
                        else:
                            break
                annealing_maxsize = x
                annealing_strands = []
                annealing_lengths = []
                break

    #homology
    homology_strands = []
    homology_lengths = []
    i = 0
    homology_size = homology_maxsize
    while i < homology_variants:
        f_homology = seq[edit_end+1:edit_end+homology_size+1]
        homology_strand = reverse_complement(f_homology)
        homology_strands.append(homology_strand)
        homology_lengths.append(homology_size)
        print("homology_{}: {}".format(homology_size, homology_strand))
        homology_size -= homology_stepsize
        i += 1
    i = 0

    #synthesis
    edits_tmp = []
    if(directionality == "f1"):
        edits_tmp = edits_conf
    elif(directionality == "r1"):
        edits_tmp = r_edits_conf
    synthesis_strands = []
    i = 0
    for i in range(len(edits_tmp)):
        f_synthesis = seq[spacer_end-3:edit_start]+edits_tmp[i]
        synthesis_strands.append(reverse_complement(f_synthesis))
    print("synthesis strands: ", synthesis_strands)

    row_name = "{}_spacer".format(directionality)
    spacer_output["name"].append(row_name)
    spacer_output["seq"].append(Mod_spacer)

    k = 0
    while k < len(edits_tmp):
        i = 0
        while i < homology_variants:
            j = 0
            while j < annealing_variants:
                row_name = "{}_homology{}_{}edit_anneal{}".format(directionality, homology_lengths[i], edits_conf[k], annealing_lengths[j])
                row_value = homology_strands[i]+synthesis_strands[k]+annealing_strands[j]
                peg_output["name"].append(row_name)
                peg_output["seq"].append(row_value)
                j += 1
            i += 1
        k += 1

    file = "output/{}_{}_output.csv".format(filename, directionality)
    with open(file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        row_name = "{}_spacer".format(directionality)
        writer.writerow([row_name, Mod_spacer])

        k = 0
        while k < len(edits_tmp):
            i = 0
            while i < homology_variants:
                j = 0
                while j < annealing_variants:
                    row_name = "{}_homology{}_{}edit_anneal{}".format(directionality, homology_lengths[i], edits_conf[k], annealing_lengths[j])
                    row_value = homology_strands[i]+synthesis_strands[k]+annealing_strands[j]
                    writer.writerow([row_name, row_value])
                    j += 1
                i += 1
            k += 1

def nick_5to3(seq, spacer_SEQ):
    global filename
    global direcitonality
    global mindist
    global spacer_size
    global homology_maxsize
    global homology_stepsize
    global homology_variants
    global annealing_maxsize
    global annealing_stepsize
    global annealing_variants
    global nick_output

    if(nick_PAM_setting == 'NGG'):
        PAM_id_conf = 'GG'

    nicking_strands = []

    spacer_cut = seq.find(reverse_complement(spacer_SEQ))
    edit_start = seq.find('[')
    edit_end = seq.find(']')

    if((spacer_cut + 3 + (edit_start - spacer_cut) - nick_dist_conf - (edit_end - edit_start)) <= 0 or (spacer_cut + 3 + nick_dist_conf + 20) > len(seq)):
        print("ERROR: !!!Sequence is short, increase size or reduce nick_max_dist!!!")
        sys.exit()

    #nick
    visual = seq[spacer_cut + 3 + (edit_start - spacer_cut) - nick_dist_conf - (edit_end - edit_start) -20:spacer_cut + 3 + nick_dist_conf +20]

    vis_spacer_cut = visual.find(reverse_complement(spacer_SEQ))
    pams = []
    i = len(visual)
    while i > 3:
        if visual[i-2:i] == PAM_id_conf:
            pams.append(i-3)
        elif(visual[i-2] == "[" or visual[i-2] == "]"):
            if(visual[i-3] != "[" and visual[i-3] != "]"):
                if(visual[i-3]+visual[i-1] == PAM_id_conf):
                    pams.append(i-4)
            else:
                if(visual[i-4]+visual[i-1] == PAM_id_conf):
                    pams.append(i-5)

        i = i-1
    r_pams = pams

    visual = split(visual)
    #create bucket of nick guides
    all_nicks = []
    for i in range(len(r_pams)):
        if(r_pams[i] > 19):
            tmpspacer_end = r_pams[i]
            pam_string = visual[tmpspacer_end - 20:r_pams[i]]
            size_fix=20
            if("]" in pam_string):
                size_fix+=1
            pam_string = visual[tmpspacer_end - size_fix:r_pams[i]]
            if("[" in pam_string):
                size_fix+=1
            tmpspacer_start = tmpspacer_end - size_fix
            all_nicks.append("".join(visual[tmpspacer_start:r_pams[i]]).replace(']','').replace('[',''))

    before_spacer = True
    before_spacer_e = True
    for i in range(len(r_pams)):
        if(before_spacer and r_pams[i]<vis_spacer_cut + len(spacer_SEQ)):
            visual.insert(vis_spacer_cut + len(spacer_SEQ), "\33[44m<<\033[0m")
            before_spacer = False
        elif(before_spacer_e and r_pams[i]<vis_spacer_cut):
            visual.insert(vis_spacer_cut, "\33[44m<<\033[0m")
            before_spacer_e = False
        tmpspacer_end = r_pams[i]
        pam_string = visual[tmpspacer_end - 20:r_pams[i]]
        size_fix=20
        if("]" in pam_string):
            size_fix+=1
        pam_string = visual[tmpspacer_end - size_fix:r_pams[i]]
        if("[" in pam_string):
            size_fix+=1
        if(tmpspacer_end - size_fix < 0):
            break
        pam_string = visual[tmpspacer_end - size_fix:r_pams[i]]
        tmpspacer_start = tmpspacer_end - size_fix
        visual.insert(tmpspacer_end, highlights[i]+str(i+1)+"\033[0m")

    string = "".join(visual)
    print("Attention: this is the antisense strand")
    print(string)

    selected_nicks = []
    unsatisfactory = True
    while unsatisfactory:
        x='y'
        while True:
            x = input("Type number to add corresponding nick site, or 'n': ")
            if (x == "n"):
                break
            try:
                x = int(x)
            except ValueError:
                print("Please input integer or 'n'")
                continue
            else:
                if(len(all_nicks) < x or x < 1):
                    print("Number out of range")
                    continue
                else:
                    break
        if(x == "n"):
            unsatisfactory = False
            break
        x = int(x) - 1

        nick_seq = all_nicks[x]
        if(nick_seq[0]=="G"):
            print("Added: ", nick_seq)
            selected_nicks.append(nick_seq)
        else:
            print("Added: (+G)", nick_seq)
            selected_nicks.append("G"+nick_seq)
    for i in range(len(selected_nicks)):
        nick_output["name"].append("{}_{}_nick_guide_{}".format(filename, directionality, i))
        nick_output["seq"].append(selected_nicks[i])

filename = input("file name?:")

r_edits_conf = []
for i in range(len(edits_conf)):
    r_edits_conf.append(reverse_complement(edits_conf[i]))

if(forward_conf == "y"):
    print("sense spacer")
    directionality = "f1"
    analyze_5to3(seq)

if(f_nick_conf == "y"):
    print("sense nick guides")
    seq2 = reverse_complement(seq)
    directionality = "f1"
    nick_5to3(seq2, spacer)

if(reverse_conf == "y"):
    print("antisense spacer")
    seq2 = reverse_complement(seq)
    directionality = "r1"
    analyze_5to3(seq2)

if(r_nick_conf == "y"):
    print("antisense nick guides")
    seq2 = seq
    directionality = "r1"
    nick_5to3(seq2, spacer)

if(forward_conf == "y" or reverse_conf == "y"):
    tmp = pd.DataFrame.from_dict(spacer_output)
    spacer_buy_t = {"name": [], "seq": []}
    for index, row in tmp.iterrows():
        name = row["name"]
        seq = row["seq"]
        spacer_buy_t["name"].append(name)
        spacer_buy_t["seq"].append("cacc"+seq+"gtttt")

    spacer_buy_b = {"name": [], "seq": []}
    for index, row in tmp.iterrows():
        name = row["name"]
        seq = row["seq"]
        spacer_buy_b["name"].append(name+"_bottom")
        r_seq = reverse_complement(seq)
        spacer_buy_b["seq"].append("ctctaaaac"+r_seq)

    tmp = pd.DataFrame.from_dict(peg_output)
    peg_buy_tmp = {"name": [], "seq": []}
    for index, row in tmp.iterrows():
        name = row["name"]
        seq = row["seq"]

        peg_buy_tmp["name"].append(name)
        peg_buy_tmp["seq"].append("gtgc"+seq)
    for index, row in tmp.iterrows():
        name = row["name"]
        seq = row["seq"]

        peg_buy_tmp["name"].append(name+"_bottom")
        r_seq = reverse_complement(seq)
        peg_buy_tmp["seq"].append("aaaa"+r_seq)

    peg_buy_csv = pd.DataFrame.from_dict(peg_buy_tmp)
    peg_buy_csv.to_csv("output/peg_buy_output.csv")
    print("SAVED")

if(f_nick_conf == "y" or r_nick_conf == "y"):
    tmp = pd.DataFrame.from_dict(nick_output)
    nick_buy_tmp_t = {"name": [], "seq": []}
    for index, row in tmp.iterrows():
        name = row["name"]
        seq = row["seq"]

        nick_buy_tmp_t["name"].append(name)
        nick_buy_tmp_t["seq"].append("cacc"+seq)

    nick_buy_tmp_b = {"name": [], "seq": []}
    for index, row in tmp.iterrows():
        name = row["name"]
        seq = row["seq"]

        nick_buy_tmp_b["name"].append(name)
        r_seq = reverse_complement(seq)
        nick_buy_tmp_b["seq"].append("aaac"+r_seq)

if(spacer_buy_t and nick_buy_tmp_t):
    spacer_buy_t = pd.DataFrame.from_dict(spacer_buy_t)
    spacer_buy_b = pd.DataFrame.from_dict(spacer_buy_b)

    nick_buy_tmp_t = pd.DataFrame.from_dict(nick_buy_tmp_t)
    nick_buy_tmp_b = pd.DataFrame.from_dict(nick_buy_tmp_b)

    spacer_buy_t = spacer_buy_t.append(nick_buy_tmp_t, ignore_index=True)
    spacer_buy_t = spacer_buy_t.append(spacer_buy_b, ignore_index=True)
    spacer_buy_t = spacer_buy_t.append(nick_buy_tmp_b, ignore_index=True)
    spacer_buy_t.to_csv("output/spacers+guide_buy.csv")
    print("SAVED")
elif(spacer_buy_t):
    spacer_buy_t = pd.DataFrame.from_dict(spacer_buy_t)
    spacer_buy_b = pd.DataFrame.from_dict(spacer_buy_b)
    spacer_buy_t = spacer_buy_t.append(spacer_buy_b, ignore_index=True)
    spacer_buy_t.to_csv("output/spacers_buy.csv")
    print("SAVED")
elif(nick_buy_tmp_t):
    nick_buy_tmp_t = pd.DataFrame.from_dict(nick_buy_tmp_t)
    nick_buy_tmp_b = pd.DataFrame.from_dict(nick_buy_tmp_b)
    nick_buy_tmp_t = nick_buy_tmp_t.append(nick_buy_tmp_b, ignore_index=True)
    nick_buy_tmp_t.to_csv("output/nick_buy.csv")
    print("SAVED")
# print(max_anneal)

# i = annealing_maxsize
# while i > (annealing_stepsize*(annealing_variants+2)):
#     i = i - annealing_stepsize
#
#
# rc_spacer[]
