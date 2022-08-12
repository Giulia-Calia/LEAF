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
import plotly.express as px
import plotly.graph_objects as go


def kyte_doolittle_score(protein_seq):
    if "X" in protein_seq:
        pass
    elif len(protein_seq) < 21:
        biopy_analyzed_seq = ProteinAnalysis(str(protein_seq))
        kd_scores = (biopy_analyzed_seq.protein_scale(window=9, param_dict=ProtParamData.kd))
        return kd_scores
    else:
        biopy_analyzed_seq = ProteinAnalysis(str(protein_seq))
        kd_scores = (biopy_analyzed_seq.protein_scale(window=21, param_dict=ProtParamData.kd))
        kd_scores_70 = kd_scores[:70]
        return kd_scores


def bin_division(score_list, n_bins):
    freq = list(np.histogram(score_list, bins=n_bins)[0])
    tmp_start_pos = 0
    list_to_append = []
    for f in freq:
        tmp_end_pos = tmp_start_pos + f
        list_to_append.append(score_list[tmp_start_pos: tmp_end_pos])
        tmp_start_pos += f
    return list_to_append


def asd(score_list, protein_sequence):
    np_score_list = np.array(score_list)
    # norm_score_list = ((np_score_list - np.min(np_score_list)) / (np.max(np_score_list) - np.min(np_score_list)) / len(score_list))
    norm_score_list = np_score_list / len(protein_sequence)
    return norm_score_list


def divide_asd(norm_score_list, n):
    p = len(norm_score_list) // n
    if len(norm_score_list) - p > 0:
        return [norm_score_list[:p]] + divide_asd(norm_score_list[p:], n - 1)
    else:
        return [norm_score_list]


def gravy_score(protein_seq):
    if "X" in protein_seq:
        pass
    else:
        biopy_analyzed_seq = ProteinAnalysis(str(protein_seq))
        return biopy_analyzed_seq.gravy()


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

    i = 0
    scores = []
    pre_scores = []
    gravy_list = []
    min_pre = []
    max_pre = []
    min_post = []
    max_post = []
    total_scores_pre = []
    total_scores_post = []

    for record in parse_input:
        if "X" in str(record.seq):
            pass
        else:
            kds = kyte_doolittle_score(str(record.seq))
            total_scores_pre += list(kds)
            if not min_pre and not max_pre:
                min_pre.append(np.min(kds))
                max_pre.append(np.max(kds))
            else:
                if np.min(kds) < min_pre[0]:
                    min_pre[0] = np.min(kds)
                if np.max(kds) > max_pre[0]:
                    max_pre[0] = np.max(kds)
            norm = asd(kds, str(record.seq))
            total_scores_post += list(norm)
            if not min_post and not max_post:
                min_post.append(np.min(norm))
                max_post.append(np.max(norm))
            else:
                if np.min(norm) < min_pre[0]:
                    min_pre[0] = np.min(norm)
                if np.max(norm) > max_pre[0]:
                    max_pre[0] = np.max(norm)
            div_norm = divide_asd(norm, 3)
            gravy_list.append(gravy_score(str(record.seq)))
            scores.append([])
            for el in div_norm:
                scores[i].append(np.mean(el))

            div_kds = divide_asd(kds, 3)
            pre_scores.append([])
            for el in div_kds:
                pre_scores[i].append(np.mean(el))
            i += 1

    norm_kds_df = pd.DataFrame(scores, columns=["begin", "center", "end"])
    # norm_kds_df["gravy_score"] = gravy_list
    print(norm_kds_df)
    print(f"global range pre: [{np.min(total_scores_pre)}, {np.max(total_scores_pre)}]\nglobal range post: [{np.min(total_scores_post)}, {np.max(total_scores_post)}]")

    kds_df = pd.DataFrame(pre_scores, columns=["begin", "center", "end"])
    print(kds_df)

    fig = px.histogram(total_scores_pre)
    fig.update_layout(title="GLOBAL SCORE PRE NORM")
    # fig.show()

    fig1 = px.histogram(total_scores_post)
    fig1.update_layout(title="GLOBAL SCORE POST NORM")
    # fig1.show()



    parse_input_noneff = SeqIO.parse("/home/giulia/Workspace/PhytoPhD/effectors_analysis/non_effectors_selection202207/non_effectors_selected20220729.fasta", "fasta")
    j = 0
    scores_noneff = []
    pre_scores_noneff = []
    gravy_list = []
    min_noneff_pre = []
    max_noneff_pre = []
    min_noneff_post = []
    max_noneff_post = []
    total_scores_noneff_pre = []
    total_scores_noneff_post = []

    for record in parse_input_noneff:
        if "X" in str(record.seq):
            pass
        else:
            kds_noneff = kyte_doolittle_score(str(record.seq))
            total_scores_noneff_pre += list(kds_noneff)
            if not min_noneff_pre and not max_noneff_pre:
                min_noneff_pre.append(np.min(kds_noneff))
                max_noneff_pre.append(np.max(kds_noneff))
            else:
                if np.min(kds_noneff) < min_noneff_pre[0]:
                    min_noneff_pre[0] = np.min(kds_noneff)
                if np.max(kds_noneff) > max_noneff_pre[0]:
                    max_noneff_pre[0] = np.max(kds_noneff)
            norm_noneff = asd(kds_noneff, str(record.seq))
            total_scores_noneff_post += list(norm_noneff)
            if not min_noneff_post and not max_noneff_post:
                min_noneff_post.append(np.min(norm_noneff))
                max_noneff_post.append(np.max(norm_noneff))
            else:
                if np.min(norm_noneff) < min_noneff_pre[0]:
                    min_noneff_pre[0] = np.min(norm_noneff)
                if np.max(norm_noneff) > max_noneff_pre[0]:
                    max_noneff_pre[0] = np.max(norm_noneff)
            div_norm_noneff = divide_asd(norm_noneff, 3)
            gravy_list.append(gravy_score(str(record.seq)))
            scores_noneff.append([])
            for el in div_norm_noneff:
                scores_noneff[j].append(np.mean(el))

            div_kds_noneff = divide_asd(kds_noneff, 3)
            pre_scores_noneff.append([])
            for el in div_kds_noneff:
                pre_scores_noneff[j].append(np.mean(el))
            j += 1

    kds_noneff_df = pd.DataFrame(pre_scores_noneff, columns=["begin", "center", "end"])
    norm_kds_noneff_df = pd.DataFrame(scores_noneff, columns=["begin", "center", "end"])

    fig5 = go.Figure()
    for s in kds_df:
        fig3 = px.histogram(kds_df[s])
        fig3.update_layout(title=f"PRE NORM {s} SCORES")
        # fig3.show()
        fig5.add_trace(go.Box(y=kds_df[s], boxpoints="all"))
        fig5.add_trace(go.Box(y=kds_noneff_df[s], boxpoints="all"))
        fig5.update_layout(title=f"PRE NORM {s} SCORES")
        # fig5.show()

    fig6 = go.Figure()

    for s in norm_kds_df:
        fig2 = px.histogram(norm_kds_df[s])
        fig2.update_layout(title=f"NORM {s} SCORES")
        # fig2.show()
        fig6.add_trace(go.Box(y=norm_kds_df[s], boxpoints="all"))
        fig6.add_trace(go.Box(y=norm_kds_noneff_df[s], boxpoints="all"))
        fig6.update_layout(title=f"NORM {s} SCORES")
        # fig6.show()

    fig7 = go.Figure()
    fig7.add_trace(go.Box(y=total_scores_pre, boxpoints="all"))
    fig7.add_trace(go.Box(y=total_scores_noneff_pre, boxpoints="all"))
    fig7.update_layout(title="PRE NORM GLOBAL SCORES DIST")

    fig8 = go.Figure()
    fig8.add_trace(go.Box(y=total_scores_post, boxpoints="all"))
    fig8.add_trace(go.Box(y=total_scores_noneff_post, boxpoints="all"))
    fig8.update_layout(title="NORM GLOBAL SCORES DIST")
    fig7.show()
    fig8.show()

    print(f"non_eff global range pre: [{np.min(total_scores_noneff_pre)}, {np.max(total_scores_noneff_pre)}]\nnon_eff global range post: [{np.min(total_scores_noneff_post)}, {np.max(total_scores_noneff_post)}]")
