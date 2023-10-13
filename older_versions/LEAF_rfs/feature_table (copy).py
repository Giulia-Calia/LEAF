import numpy as np
import pandas as pd
import argparse
import os
import re
from Bio import SeqIO


def merge_dfs(df1, df2):
    return df1.merge(df2, on="seq_id", how="left").fillna(0)


def find_motif_start_pos(motif, sequence):
    # Find all indices of 'motif'
    indices_object = re.finditer(pattern=str(motif), string=str(sequence))
    indices = [index.start() for index in indices_object]
    return indices


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="",
                                     epilog="")
    parser.add_argument("-mi", "--mod_ids",
                        help=".txt file with a single column of name of effector "
                             "class associated with a protein id "
                             "e.g.: TENGU_A0A4P6MDK8")

    parser.add_argument("-i", "--input",
                        help="fasta file with protein sequences")
    parser.add_argument("-o", "--output_file",
                        help="path and name of the output file, a .tsv table of features ")
    parser.add_argument("-u", "--uniprot",
                        action="store_true")
    parser.add_argument("-p", "--prodigal",
                        action="store_true")
    parser.add_argument("-sp", "--signalP_input",
                        help="SignalP output file in .txt format")
    parser.add_argument("-tm", "--tm_input",
                        help="TMHMM output file in .txt format")
    parser.add_argument("-mb", "--mobi_input",
                        help="MobiDB-lite output file in .txt format")
    parser.add_argument("-pr", "--prosite_input",
                        help="Prosite output file in .fasta format")
    parser.add_argument("-prm", "--prosite_motifs",
                        help=".txt file used also as input for MOnSTER with a motif per line predicted in the "
                             "training positive_class")
    parser.add_argument("-ms", "--monster_scores",
                        help="MOnSTER summary output")
    parser.add_argument("-mmc", "--monster_mot_clump",
                        help="one of MOnSTER output file having as columns motif|CLUMP")
    parser.add_argument("--eff", help="if the fasta contains effector sequences",
                        action="store_true")
    parser.add_argument("--noneff", help="if the fasta contains non effector sequences",
                        action="store_true")

    args = parser.parse_args()

    # PARSING INPUT FASTA #
    input_proteins = {}  # a dictionary representation of a fasta file {protein_id: AA_sequence, ...}
    proteins_length = []
    proteins_name = []
    for record in SeqIO.parse(args.input, "fasta"):
        input_proteins[str(record.id)] = str(record.seq)
        proteins_length.append(len(str(record.seq)))
        proteins_name.append(str(record.id))

    name_list = []
    if args.eff:
        name_list = ["effector"] * len(input_proteins)
    elif args.noneff:
        name_list = ["non_effector"] * len(input_proteins)
    else:
        name_list = ["to_be_predicted"] * len(input_proteins)
    # 'returns' name_length_cols = a df with columns = protein ids, name of belonging class (effector or not)
    # and respective protein length
    name_length_cols = pd.DataFrame({"seq_id": proteins_name,
                                     "name": name_list,
                                     "sequence length": proteins_length})
    # SIGNALP #
    # 'returns' df_sp = a df with 3 columns = protein names, sequence length and signal peptide feature
    sp_predictions = {}
    with open(args.signalP_input, "r") as signalp_file:
        file = signalp_file.readlines()
        for line in file:
            # from line like Name=tr|A0A0A8JCG0|A0A0A8JCG0_9MOLU	SP='YES'    ...
            if line.startswith("Name="):
                tmp_id = line.split("\t")[0][line.split("\t")[0].find("=") + 1:]
                sp_in_line = line.split("\t")[1].split(" ")[0]
                sp_signal = str(sp_in_line[sp_in_line.find("=") + 2: -1])
                d_score_in_line = line.split("\t")[1][line.split("\t")[1].find("D=") + 2:]
                # if there is a signal peptide
                if sp_signal == "YES":
                    # e.g. Name=tr|A0A0A8JCG0|A0A0A8JCG0_9MOLU	SP='YES' Cleavage site pos. 34 and 35: VMG-MN D=0.696
                    d_score = float(d_score_in_line[:d_score_in_line.find(" ")])
                    sp_predictions[tmp_id] = d_score
                else:
                    pass
            else:
                pass

    df_sp_predictions = pd.DataFrame({"seq_id": list(sp_predictions.keys()),
                                      "signal peptide": list(sp_predictions.values())})
    # if not all the proteins are predicted to have a signal peptide, reconstruct the exact order of the proteins
    # merging the perv df with the new one having signalP predictions
    df_sp = merge_dfs(name_length_cols, df_sp_predictions)

    # TMHMM #
    # 'returns' df_sp_tm = a df with 8 columns = protein names, sequence length, signal peptide feature,
    # transmembrane domain,	aa in tr domain, first 60 aa, prob N-in, warning signal sequence

    tm_features = ["transmembrane domain", "aa in tr domain", "first 60 aa",
                   "prob N-in", "warning signal sequence"]
    tm_presence = {}
    aa_in_tm = {}
    first_sixty_aa = {}
    prob_n_in = {}
    warning_on_sp = {}

    # an entry in tmhmm output is like 
    # # tr|A0A0A8JCG0|A0A0A8JCG0_9MOLU Length: 125
    # # tr|A0A0A8JCG0|A0A0A8JCG0_9MOLU Number of predicted TMHs:  1
    # # tr|A0A0A8JCG0|A0A0A8JCG0_9MOLU Exp number of AAs in TMHs: 20.43649
    # ...
    # so the construction of the features is pretty straight forward
    with open(args.tm_input, "r") as tmhmm_file:
        file = tmhmm_file.readlines()
        for line in file:
            if line.startswith("#"):
                name = line.split(" ")[1]
                if "Number of predicted TMHs:" in line:
                    tm_presence[name] = int(line.split(" ")[-1].replace("\n", ""))
                if "Exp number of AAs in TMHs:" in line:
                    aa_in_tm[name] = float(line.split(" ")[-1].replace("\n", ""))
                if "Exp number, first 60 AAs:" in line:
                    first_sixty_aa[name] = float(line.split(" ")[-1].replace("\n", ""))
                if "Total prob of N-in:" in line:
                    prob_n_in[name] = float(line.split(" ")[-1].replace("\n", ""))
                # if there is a tm domain but there is also the possibility of a signal peptide, this
                # feature become 1 = True, 0 otherwise
                if "signal sequence" in line:
                    warning_on_sp[name] = 1
                else:
                    warning_on_sp[name] = 0

    df_tm_predictions = pd.DataFrame({"seq_id": list(tm_presence.keys()),
                                      tm_features[0]: list(tm_presence.values()),
                                      tm_features[1]: list(aa_in_tm.values()),
                                      tm_features[2]: list(first_sixty_aa.values()),
                                      tm_features[3]: list(prob_n_in.values()),
                                      tm_features[4]: list(warning_on_sp.values())
                                      })
    # the merging assure that even if there are some missing proteins in the tmhmm output,
    # they will be added to the dataframe in the same order of the input fasta file and their
    # values are filled with 0 value
    df_sp_tm = merge_dfs(df_sp, df_tm_predictions)

    # MOBIDB-LITE #
    #  # 'returns' df_sp_tm = a df with 9 columns = protein names, sequence length, signal peptide feature,
    #     # transmembrane domain,	aa in tr domain, first 60 aa, prob N-in, warning signal sequence, MobiDB-lite

    # if the file is empty, no disordered regions are found so return a 0 columns
    if os.stat(args.mobi_input).st_size == 0:
        mb_empty = pd.DataFrame({"seq_id": proteins_name, "MobiDB-lite": [0] * len(proteins_name)})
        df_sp_tm_mb = merge_dfs(df_sp_tm, mb_empty)

    # otherwise the entry would be like:
    # {
    # "acc": "tr|A0A1S9LYF9|A0A1S9LYF9_9MOLU",
    #     # "consensus": "SSSSSSSSSSSSSSSSSSSSSSSSDDDDDDDDDDDDDDDDDDDDD", ...
    else:
        # mb_predictions = {"seq_id": proteins_name, "MobiDB-lite": [0] * len(proteins_name)}
        k = []
        v = []
        with open(args.mobi_input, "r") as mobi_file:
            mobi_file = mobi_file.readlines()
            for line in mobi_file:
                line = line.strip()
                if line.startswith('"acc"'):
                    k.append(line.split(" ")[1].replace('"', '').replace(",", ""))
                if line.startswith('"consensus"'):
                    v.append(int(line.split(" ")[1].count("D")))

        # mb_predictions = dict(zip(k, v))
        df_mb_predictions = pd.DataFrame({"seq_id": k, "MobiDB-lite": v})
        df_sp_tm_mb = merge_dfs(df_sp_tm, df_mb_predictions)

    motif_list = list(pd.read_csv(args.prosite_motifs, header=None)[0])
    # PROSITE #
    # parse the PROSITE output file of effector class to build up a dictionary containing
    # {functional name of motifs: associated sequences}
    pros_mot = {}
    for record in SeqIO.parse(args.prosite_input, "fasta"):
        func_name = record.description.strip().split(" ")[-1]
        pros_mot[str(record.seq)] = func_name

    pros_mot_occ = {}
    for protein in list(input_proteins.keys()):
        pros_mot_occ[protein] = {}
        for m in list(pros_mot.keys()):
            mot_pres = find_motif_start_pos(m, input_proteins[protein])
            if mot_pres:
                if pros_mot[m] not in list(pros_mot_occ[protein].keys()):
                    pros_mot_occ[protein][pros_mot[m]] = len(mot_pres)
                else:
                    pros_mot_occ[protein][pros_mot[m]] += len(mot_pres)
        for fn in list(np.unique(list(pros_mot.values()))):
            if fn not in (pros_mot_occ[protein].keys()):
                pros_mot_occ[protein][fn] = 0

    df_pros_mot = pd.DataFrame(pros_mot_occ).transpose()
    df_pros_mot.reset_index(inplace=True)
    df_pros_mot = df_pros_mot.rename(columns={"index": "seq_id"})

    # CLUMPs #
    monster_sc = pd.read_csv(args.monster_scores, sep=",", header=0)[["CLUMP", "monster_score"]]
    selected_clumps = list(
        monster_sc["CLUMP"][monster_sc["monster_score"] > np.mean(list(monster_sc["monster_score"]))])
    selected_clumps.sort()
    monster_motsclumps = pd.read_csv(args.monster_mot_clump, sep=",", header=0)

    clump_features = {}
    for protein in list(input_proteins.keys()):
        clump_features[protein] = [0] * len(selected_clumps)
        for motif in motif_list:
            tmp_clump = list(monster_motsclumps["CLUMP"][monster_motsclumps["motif"] == motif])[0]
            # if motif in str(input_proteins[protein]) and tmp_clump in selected_clumps:
            if motif in str(input_proteins[protein]) and tmp_clump in selected_clumps:
                clump_features[protein][selected_clumps.index(tmp_clump)] += 1

    df_clump = pd.DataFrame(clump_features, index=[f"CLUMP{c}" for c in selected_clumps]).transpose()
    df_clump.reset_index(inplace=True)
    df_clump = df_clump.rename(columns={"index": "seq_id"})

    # CLUMPs positions #
    in_protein_bins = {}
    for protein in list(input_proteins.keys()):
        in_protein_bins[protein] = dict(zip([f"bin{i + 1}" for i in range(4)], [0 for i in range(4)]))

        # here comes the assign to the corresponding bin of each motifs belonging to each of the selected clumps
        for motif in motif_list:
            tmp_pos = find_motif_start_pos(motif, input_proteins[protein])

            for p in tmp_pos:  # all positions of motifs on seq
                tmp_norm_pos = int(p) / len(input_proteins[protein])  # normalized position with the seq len (0-100%)
                # find the motifs in the MOnSTER output to see at which CLUMP it belongs
                clump_belong = list(monster_motsclumps["CLUMP"][monster_motsclumps["motif"] == motif])
                # if the motifs finds a corresponding clump in the MOnSTER output and the clump is in the ones selected
                if clump_belong and clump_belong[0] in selected_clumps:
                    if tmp_norm_pos <= 0.25:  # 0-25% of the sequence = bin1
                        in_protein_bins[protein]["bin1"] += 1
                    elif 0.25 < tmp_norm_pos <= 0.50:  # 26-50% of the sequence = bin2
                        in_protein_bins[protein]["bin2"] += 1
                    elif 0.50 < tmp_norm_pos <= 0.75:  # 51-75% of the sequence = bin3
                        in_protein_bins[protein]["bin3"] += 1
                    else:  # 76-100% of the sequence = bin4
                        in_protein_bins[protein]["bin4"] += 1
                else:
                    pass

    df_clumps_pos = pd.DataFrame(in_protein_bins).transpose()
    df_clumps_pos.reset_index(inplace=True)
    df_clumps_pos = df_clumps_pos.rename(columns={"index": "seq_id"})
    # merge df_clumps into one clump feature and then all together
    df_clump_features = merge_dfs(df_clump, df_clumps_pos)

    # FINAL DF OF FEATURES #
    df_sp_tm_mb_pr = merge_dfs(df_sp_tm_mb, df_pros_mot)
    df_sp_tm_mb_pr_cl = merge_dfs(df_sp_tm_mb_pr, df_clump_features)
    # OUTPUT #
    df_sp_tm_mb_pr_cl.to_csv(f"{args.output_file}.tsv", sep="\t", index=False)
