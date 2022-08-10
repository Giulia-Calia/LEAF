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

        elif mod_names:
            protein_names.append(record.id)
            protein_length.append(len(record.seq))

        else:
            print("ERR: I don't know the input file format, thus I can't help you in parsing the feature prediction "
                  "results, please give me a fasta file coming from uniprot or prodigal")

    # df_base = pd.DataFrame(list(zip(protein_names, proteins_associated_organism, protein_length)), columns=["name", "organism", "sequence length"])
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

            if "SP='NO'" in el:
                sp_predictions[name] = "False"
            elif "SP='YES'" in el:
                sp_predictions[name] = "True"
            else:
                pass
    sp_col = {"name": list(sp_predictions.keys()), "signal peptide": list(sp_predictions.values())}
    return sp_col


def parse_tmhmm(input_f, uniprot=False):
    tm_predictions = {}
    with open(input_f, "r") as tmhmm_file:
        file = tmhmm_file.readlines()
        for el in file:
            if "Number of predicted TMHs:" in el and el.strip().split(" ")[-1] != "0":
                if uniprot:
                    name = el.split(" ")[1].split("|")[1]
                else:
                    name = el.split(" ")[1]
                tm_predictions[name] = "True"
            elif "Number of predicted TMHs:" in el and el.strip().split(" ")[-1] == "0":
                if uniprot:
                    name = el.split(" ")[1].split("|")[1]
                else:
                    name = el.split(" ")[1]
                tm_predictions[name] = "False"

    tm_col = {'name': list(tm_predictions.keys()), "transmembrane domain": list(tm_predictions.values())}
    return tm_col


def parse_phobius(input_f, uniprot=False):
    lines_list = [[]]
    ph_predictions_sp = {}
    ph_predictions_tm = {}
    i = 0

    with open(input_f, "r") as phob_file:
        file = phob_file.readlines()
        for el in file:
            el = el.strip().split()
            if el != ["//"] and el:
                lines_list[i] += el
            else:
                lines_list.append([])
                i += 1

    for chunk in lines_list[:-1]:
        if uniprot:
            name = str(chunk[1]).split("|")[1]
        else:
            name = chunk[1]
        if "SIGNAL" in chunk and "TRANSMEM" not in chunk:
            ph_predictions_sp[name] = "True"
            ph_predictions_tm[name] = "False"

        elif "SIGNAL" in chunk and "TRANSMEM" in chunk:
            ph_predictions_sp[name] = "True"
            ph_predictions_tm[name] = "True"

        elif "TRANSMEM" in chunk and "SIGNAL" not in chunk:
            ph_predictions_sp[name] = "False"
            ph_predictions_tm[name] = "True"

        elif "SIGNAL" not in chunk and "TRANSMEM" not in chunk:
            ph_predictions_sp[name] = "False"
            ph_predictions_tm[name] = "False"
        else:
            ph_predictions_sp[name] = "False"
            ph_predictions_tm[name] = "False"

    ph_cols = {"name": list(ph_predictions_sp.keys()),
               "Phob signal peptide": list(ph_predictions_sp.values()),
               "Phob transmembrane domain": list(ph_predictions_tm.values())}
    return ph_cols


def parse_mobidb_lite(input_f, protein_names, mod_names=None, uniprot=False):
    mb_empty = {"name": protein_names, "disordered regions": ["False"] * len(protein_names)}
    mb_predictions = dict(zip(protein_names, ["False"] * len(protein_names)))
    # se file e' vuoto, allora tutta la lista e' =FALSE()
    if os.stat(input_f).st_size == 0:
        return mb_empty
    # altrimenti:
    else:
        with open(input_f, "r") as mobi_file:
            lines = mobi_file.readlines()
            for el in lines:
                el = el.strip()
                if el.startswith('"acc"'):
                    for n in protein_names:
                        if n in el:
                            mb_predictions[n] = "True"

        mb_col = {"name": mb_predictions.keys(), "MobiDB-lite": mb_predictions.values()}
        return mb_col


def parse_prosite(input_f, protein_names, uniprot=False, prodigal=False, mod_names=None):
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

        if record_id not in list(dict_record_feature.keys()):
            dict_record_feature[record_id] = []
        # record_id = record.id[record.id.find("|") + 1:record.id.rfind("|")]
        feature_name = record.description.split(" ")[-1]  # feature = motif name
        if feature_name not in col_names_auto:
            if feature_name.startswith("L"):
                feature_name = str(record.description.split(" ")[-2]) + "_" + str(record.description.split(" ")[-1])
            col_names_auto.append(feature_name)
            pr_col[feature_name] = []
        else:
            pass

        if record_id != tmp_record_id:
            uniq_list_features = []
            if feature_name not in uniq_list_features:
                uniq_list_features.append(feature_name)
                dict_record_feature[record_id].append(feature_name)
            else:
                pass
        elif record_id == tmp_record_id:
            if feature_name not in uniq_list_features:
                uniq_list_features.append(feature_name)
                dict_record_feature[record_id].append(feature_name)
            else:
                pass
        else:
            pass
        tmp_record_id = record_id

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
        row = []
        pr_predictions.append([])
        for motif in col_names_auto:
            if motif in complete_dict_record_feature[el]:
                row.append("True")
            else:
                row.append("False")
        pr_predictions[list(complete_dict_record_feature.keys()).index(el)] += row

    pr_predictions_t = np.transpose(pr_predictions)
    pr_col["name"] = protein_names
    for i in range(len(list(pr_col.keys())[:-1])):
        pr_col[list(pr_col.keys())[i]] = list(pr_predictions_t[i])

    return pr_col


# TOBEADDED - hydrphobic profile
def parse_hydrophob_profile(input_f):
    h_bin_cols = pd.read_csv(input_f, sep="\t")
    return h_bin_cols


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
    parser.add_argument("-ph", "--phob_input")
    parser.add_argument("-mb", "--mobi_input")
    parser.add_argument("-pr", "--prosite_input")
    parser.add_argument("-hp", "--hydrophob_profile")

    args = parser.parse_args()

    # print(parse_tmhmm("/home/giulia/Workspace/PhytoPhD/effectors_analysis/classification/prediction_tools/pred_non_eff20220804/tmhmm_non_eff20220729.txt",
    #                   uniprot=True))
    # exit(1)
    if args.mod_ids:
        # if you want to modify the uniprot ids as you declared in the .txt passed with --mod_ids
        base_df = parse_input_seq(args.input,
                                  mod_names=parse_mod_ids(args.mod_ids),
                                  uniprot=args.uniprot,
                                  prodigal=args.prodigal)
        mb = parse_mobidb_lite(args.mobi_input, list(base_df["name"]), mod_names=parse_mod_ids(args.mod_ids),
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

    ph = parse_phobius(args.phob_input, uniprot=args.uniprot)
    sp_tm_ph_df = sp_tm_df.merge(pd.DataFrame(ph), on="name")

    sp_tm_ph_mb_df = sp_tm_ph_df.merge(pd.DataFrame(mb), on="name")

    pr = parse_prosite(args.prosite_input, list(base_df["name"]), uniprot=args.uniprot, prodigal=args.prodigal, mod_names=args.mod_ids)
    sp_tm_ph_mb_pr_df = sp_tm_ph_mb_df.merge(pd.DataFrame(pr), on="name")
    sp_tm_ph_mb_pr_df.to_csv(args.output_file, sep="\t", index=False)

    hp = parse_hydrophob_profile(args.hydrophob_profile)
    sp_tm_ph_mb_pr_hp_df = pd.concat([sp_tm_ph_mb_pr_df, hp], axis=1, join="inner")
    sp_tm_ph_mb_pr_hp_df.to_csv(args.output_file, sep="\t", index=False)
