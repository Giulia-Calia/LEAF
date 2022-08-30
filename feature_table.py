import numpy as np
import pandas as pd
import argparse
import os
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import binned_statistic
from scipy.stats import pearsonr


def parse_mod_ids(name_match):
    return pd.read_csv(name_match, header=None)[0].values.tolist()


def parse_input_seq(input_s, mod_names=None, uniprot=False, prodigal=False):
    # CREATE THE BASIC DF TO WHICH ADD ALL THE COLUMNS OF INTEREST
    # THE FILE MUST HAVE BEEN DOWNLOADED FROM UNIPROT OR PRODIGAL
    parsed_file = SeqIO.parse(input_s, "fasta")

    protein_names = []
    # proteins_associated_organism = []
    protein_length = []
    for record in parsed_file:
        if uniprot:
            name = record.id.split("|")[1]
            if mod_names is not None:
                for n in mod_names:
                    if len(name) == len(mod_names):
                        name = n
            protein_names.append(name)
            # record_desc = record.description.replace("[", "").replace("'", "").replace("]", "").replace(",", "")
            # if OS= or OX= are not present, the organism is the same as the record.id without the prefix
            # org = record_desc[record_desc.find("OS=") + 3:record_desc.rfind("OX=")]
            # proteins_associated_organism.append(org)
            protein_length.append(len(record.seq))

        elif prodigal:
            strain_name = record.id[:record.id.find("_")]
            protein_names.append(record.id)
            # proteins_associated_organism.append(f"Candidatus_phytoplasma_mali_({strain_name})")
            protein_length.append(len(record.seq))

        elif mod_names is None and uniprot is False and prodigal is False:
            protein_names.append(record.id)
            protein_length.append(len(record.seq))

    df_base = pd.DataFrame(list(zip(protein_names, protein_length)), columns=["name", "sequence length"])

    return df_base


def parse_signalp(input_f, uniprot=False, prodigal=False, mod_names=None):
    sp_predictions = {}
    with open(input_f, "r") as signalp_file:
        file = signalp_file.readlines()
        for el in file:
            if el.startswith("# Name="):
                if uniprot:
                    name = el.split("|")[1]
                elif prodigal:
                    name = el.split("=")[1].split("\t")[0]
                elif mod_names:
                    name = el.split("=")[1].split("\t")[0]
                else:
                    name = el.split("=")[1].split("\t")[0]

            if el.strip().startswith("D"):
                d_score = float(el.strip().split("  ")[-3])
                sp_predictions[name] = d_score
            else:
                pass
    sp_col = {"name": list(sp_predictions.keys()), "signal peptide": list(sp_predictions.values())}
    return sp_col


def parse_tmhmm(input_f, uniprot=False):
    tm_presence = {}
    aa_in_tm = {}
    first_sixty_aa = {}
    prob_n_in = {}
    warning_sp = {}
    with open(input_f, "r") as tmhmm_file:
        file = tmhmm_file.readlines()
        for el in file:
            if uniprot and "#" in el:
                name = el.split(" ")[1].split("|")[1]

            elif not uniprot and "#" in el:
                name = el.split(" ")[1]

            if "Number of predicted TMHs:" in el:
                tm_presence[name] = int(el.split(" ")[-1].replace("\n", ""))
            if "Exp number of AAs in TMHs:" in el:
                aa_in_tm[name] = float(el.split(" ")[-1].replace("\n", ""))
            if "Exp number, first 60 AAs:" in el:
                first_sixty_aa[name] = float(el.split(" ")[-1].replace("\n", ""))
            if "Total prob of N-in:" in el:
                prob_n_in[name] = float(el.split(" ")[-1].replace("\n", ""))
            if "signal sequence" in el:
                warning_sp[name] = "True"
    names = (tm_presence.keys())
    complete_warning_sp = {}
    for n in names:
        if n not in list(warning_sp.keys()):
            complete_warning_sp[n] = "False"
        else:
            complete_warning_sp[n] = "True"
    if len(names) == len(tm_presence) == len(aa_in_tm) == len(first_sixty_aa) == len(prob_n_in) == len(
            complete_warning_sp):
        tm_col = {"name": names, "transmembrane domain": list(tm_presence.values()),
                  "aa in tr domain": list(aa_in_tm.values()), "first 60 aa": list(first_sixty_aa.values()),
                  "prob N-in": list(prob_n_in.values()), "warning signal sequence": list(complete_warning_sp.values())}
        return tm_col


def parse_mobidb_lite(input_f, protein_names, uniprot=False):
    mb_empty = {"name": protein_names, "disordered regions": [0] * len(protein_names)}
    mb_predictions = dict(zip(protein_names, [0] * len(protein_names)))
    # se file e' vuoto, allora tutta la lista e' =FALSE()
    if os.stat(input_f).st_size == 0:
        return mb_empty
    # altrimenti:
    else:
        k = []
        v = []
        with open(input_f, "r") as mobi_file:
            mobi_file = mobi_file.readlines()
            for line in mobi_file:
                line = line.strip()
                if line.startswith('"acc"'):
                    if uniprot:
                        k.append(line.strip().split(" ")[1].split("|")[1])
                    else:
                        k.append(line.strip().split(" ")[1].replace('"', '').replace(",", ""))
                if line.startswith('"consensus"'):
                    v.append(line.strip().split(" ")[1].count("D"))

        mb_dict = dict(zip(k, v))
        for k in list(mb_dict.keys()):
            mb_predictions[k] = mb_dict[k]

        mb_col = {"name": list(mb_predictions.keys()), "MobiDB-lite": list(mb_predictions.values())}
        return mb_col


def parse_prosite(input_f, protein_names, protein_length, uniprot=False, prodigal=False, mod_names=None):
    # def of feature: any type of predicted motif, pattern or profile
    record_id = ""
    uniq_list_features = []
    # for each record=protein, all the feature assigned to it
    dict_record_feature = {}
    # --probable pattern or motifs in effectors--
    # ["CK2_PHOSPHO_SITE", "MYRISTYL", "PKC_PHOSPHO_SITE", "ASN_GLYCOSYLATION",
    # "CAMP_PHOSPHO_SITE", "AMIDATION", "TYR_PHOSPHO_SITE_1", "L=(-1)"]

    # list of any kind of feature predicted by prosite
    col_names_auto = []
    # list of lists with presence/absence of the feature for each protein
    pr_col = {}
    names = []
    pr_predictions = []
    prosite_file = SeqIO.parse(input_f, "fasta")
    tmp_record_id = ""
    for record in prosite_file:
        if uniprot:
            record_id = record.id[:record.id.rfind("/")].split("|")[1]
        elif prodigal:
            record_id = record.id[:record.id.rfind("/")]
        elif not uniprot and not prodigal:
            record_id = record.id[:record.id.rfind("/")]
        if record_id not in list(dict_record_feature.keys()):
            dict_record_feature[record_id] = []
        # record_id = record.id[record.id.find("|") + 1:record.id.rfind("|")]
        feature_name = record.description.split(" ")[-1]  # feature = motif name
        if feature_name.startswith("L"):
            feature_name = str(record.description.split(" ")[-2]) + "_" + str(record.description.split(" ")[-1])
        dict_record_feature[record_id].append(feature_name)
        if feature_name not in col_names_auto:
            col_names_auto.append(feature_name)
            pr_col[feature_name] = []
        else:
            pass
    print(dict_record_feature)
    # puo' essere che non tutte le sequenze contengano un pattern o un motivo quindi bisogna assicurarsi di non perdersi
    # quegli id per ricorstuire il file dopo
    complete_dict_record_feature = {}
    for i in range(len(protein_names)):
        # cioe' se le due liste non hanno stessa lunghezza aggiungi le posizioni(id) mancanti ad un nuovo dizionario che
        # le conterra' insieme a quelle gia' esistenti
        if protein_names[i] in list(dict_record_feature.keys()):
            complete_dict_record_feature[protein_names[i]] = dict_record_feature[protein_names[i]]
        else:
            complete_dict_record_feature[protein_names[i]] = []

    for el in complete_dict_record_feature:
        pr_predictions.append([])
        # row = [complete_dict_record_feature[el].count(cl)/protein_length[list(complete_dict_record_feature.keys()).index(el)] for cl in col_names_auto]
        # row = [complete_dict_record_feature[el].count(cl) for cl in col_names_auto]
        if len(complete_dict_record_feature[el]) == 0:
            row = [0 for cl in col_names_auto]
        else:
            row = [complete_dict_record_feature[el].count(cl) / len(complete_dict_record_feature[el]) for cl in col_names_auto]
        pr_predictions[list(complete_dict_record_feature.keys()).index(el)] += row


    pr_predictions_t = np.transpose(pr_predictions)
    pr_col["name"] = protein_names
    for i in range(len(list(pr_col.keys())[:-1])):
        pr_col[list(pr_col.keys())[i]] = list(pr_predictions_t[i])

    return pr_col


def parse_hydrophob_profile(input_f):
    h_bin_cols = pd.read_csv(input_f, sep="\t")
    return h_bin_cols


def parse_aac(input_f):
    aac_cols = pd.read_csv(input_f, sep="\t")
    return aac_cols


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="",
                                     epilog="")
    parser.add_argument("-mi", "--mod_ids",
                        help=".txt file with a single column of name of effector class associated with a protein id "
                             "e.g.: TENGU_A0A4P6MDK8")
    parser.add_argument("-i", "--input",
                        help="fasta file with protein sequences")
    parser.add_argument("-o", "--output_file",
                        help="path and name of the output file")
    parser.add_argument("-u", "--uniprot",
                        action="store_true")
    parser.add_argument("-p", "--prodigal",
                        action="store_true")
    parser.add_argument("-sp", "--signalP_input")
    parser.add_argument("-tm", "--tm_input")
    # parser.add_argument("-ph", "--phob_input")
    parser.add_argument("-mb", "--mobi_input")
    parser.add_argument("-pr", "--prosite_input")
    parser.add_argument("-hp", "--hydrophob_profile")
    parser.add_argument("-aac", "--aacomposition")

    args = parser.parse_args()

    if args.mod_ids:
        # if you want to modify the uniprot ids as you declared in the .txt passed with --mod_ids
        base_df = parse_input_seq(args.input,
                                  mod_names=parse_mod_ids(args.mod_ids),
                                  uniprot=args.uniprot,
                                  prodigal=args.prodigal)
        mb = parse_mobidb_lite(args.mobi_input, list(base_df["name"]),
                               uniprot=args.uniprot)
    else:
        # if you do not want to modify any ids
        base_df = parse_input_seq(args.input,
                                  uniprot=args.uniprot,
                                  prodigal=args.prodigal)
        mb = parse_mobidb_lite(args.mobi_input, list(base_df["name"]), uniprot=args.uniprot)

    sp = parse_signalp(args.signalP_input, uniprot=args.uniprot, prodigal=args.prodigal, mod_names=args.mod_ids)
    sp_df = base_df.merge(pd.DataFrame(sp), on="name")

    tm = parse_tmhmm(args.tm_input, uniprot=args.uniprot)
    sp_tm_df = sp_df.merge(pd.DataFrame(tm), on="name")

    sp_tm_mb_df = sp_tm_df.merge(pd.DataFrame(mb), on="name")

    pr = parse_prosite(args.prosite_input, list(base_df["name"]), list(base_df["sequence length"]), uniprot=args.uniprot, prodigal=args.prodigal,  mod_names=args.mod_ids)
    sp_tm_mb_pr_df = sp_tm_mb_df.merge(pd.DataFrame(pr), on="name")
    sp_tm_mb_pr_df.to_csv(args.output_file, sep="\t", index=False)

    # aac = parse_aac(args.aacomposition)
    # sp_tm_mb_pr_aac_df = sp_tm_mb_pr_df.merge(aac, on="name")
    # sp_tm_mb_pr_aac_df.to_csv(args.output_file, sep="\t", index=False)

    # hp = parse_hydrophob_profile(args.hydrophob_profile)
    # sp_tm_mb_pr_hp_df = pd.concat([sp_tm_mb_pr_df, hp], axis=1, join="inner")
    # sp_tm_mb_pr_hp_df.to_csv(args.output_file, sep="\t", index=False)

