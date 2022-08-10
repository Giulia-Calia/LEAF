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


def kyte_doolittle_score(protein_seq):
    biopy_analyzed_seq = ProteinAnalysis(str(protein_seq))
    kd_scores = (biopy_analyzed_seq.protein_scale(window=21, param_dict=ProtParamData.kd))
    kd_scores_70 = kd_scores[:70]
    return kd_scores_70

def bin_division(score_list, n_bins):
    freq = list(np.histogram(score_list, bins=n_bins)[0])
    tmp_start_pos = 0
    list_to_append = []
    for f in freq:
        tmp_end_pos = tmp_start_pos + f
        list_to_append.append(score_list[tmp_start_pos: tmp_end_pos])
        tmp_start_pos += f
    return list_to_append


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="",
                                     epilog="")

    parser.add_argument("-i", "--input_file",
                        help="input .fasta file with protein sequences on which calculate the hydrophobic profile")
    parser.add_argument("-o", "--output_file",
                        help="output .txt file containing dataframe of kyte_doolittle scores mean for each protein and "
                             "each bin")
    parser.add_argument("-b", "--n_bins",
                        type=int,
                        help="number of bins in which divide each hydrophobic profile")

    args = parser.parse_args()

    parse_input = SeqIO.parse(args.input_file, "fasta")
    columns = [f"h_bin{str(i+1)}" for i in list(range(args.n_bins))]
    bins = []

    # for record in parse_input:
    #     # for each sequence calculate the mean score fore each bin
    #     actual_bins_mean = []
    #     kds = kyte_doolittle_score(str(record.seq))
    #     actual_bin = bin_division(kds, args.n_bins)
    #
    #     for b in actual_bin:
    #         if not b:
    #             actual_bins_mean.append(0)
    #         else:
    #             actual_bins_mean.append(np.mean(b))
    #     bins.append(actual_bins_mean)
    # df_bins = pd.DataFrame(bins, columns=columns)
    uniform_bins = []
    for record in parse_input:
        kds = kyte_doolittle_score(str(record.seq))
        uniform_bins.append([kds[i:i+7] for i in range(0, len(kds), 7)])

    i = 0
    for el in uniform_bins:
        if len(el) < 10:
            el += [[0] for i in range(0, 10-len(el))]
        else:
            pass

    j = 0
    mean_uniform_bins = []
    for el in uniform_bins:
        mean_uniform_bins.append([])
        for l in el:
            if l == [0]:
                mean_uniform_bins[j].append(0)

            else:
                mean_uniform_bins[j].append(np.mean(l))

        j += 1

    print(mean_uniform_bins[0])
    df_bins = pd.DataFrame(mean_uniform_bins, columns=columns)
    df_bins.to_csv("/home/giulia/Workspace/PhytoPhD/effectors_analysis/classification/prediction_tools/pred_dataset_validazione20220805/uni_hydrophob_putative_eff_eff20220805.txt",
                   sep="\t",
                   index=False)


