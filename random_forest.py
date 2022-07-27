import pandas as pd
import numpy as np
import joblib
import argparse
from sklearn.ensemble import RandomForestClassifier  # random forest model
from sklearn.model_selection import GridSearchCV
from sklearn import metrics  # to calculate the accuracy of the model
from sklearn.metrics import (precision_recall_curve, PrecisionRecallDisplay, auc)
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedShuffleSplit
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from collections import Counter
from imblearn.under_sampling import RandomUnderSampler


def parse_table(input_file):
    # df of feature to be considered
    parsed_df = pd.read_csv(input_file, sep="\t", header=0)
    return parsed_df


def train_test_rfc(eff_pred_df, non_eff_pred_df, path_to_out_trained_model):
    original_common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    original_common_df.replace([True, False], [1, 0], inplace=True)

    # EFF PREDICTIONS DF for Random Forest
    eff_pred_df["name"] = ["effectors"] * len(eff_pred_df)

    # NON-EFF PREDICTIONS DF for Random Forest
    non_eff_pred_df["name"] = ["non_effectors"] * len(non_eff_pred_df)

    # COMBINE DFs (taking only effector motifs)
    for col in eff_pred_df.columns[3:]:
        if col not in non_eff_pred_df.columns[3:]:
            # if an effector motifs is not present in non effectors, a new column will be added to the non_effector_df
            non_eff_pred_df[col] = [0] * len(non_eff_pred_df)
        else:
            pass

    # "inner" will take only common cols
    common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    common_df.replace([True, False], [1, 0], inplace=True)

    # drop the column/s having the same value for all samples
    common_df = common_df.drop("EF_HAND_1_L=(-1)", 1)

    ### random under sampling
    # sampling_strategy = 0.8
    # rus = RandomUnderSampler(sampling_strategy=sampling_strategy)
    rus = RandomUnderSampler(random_state=42, replacement=True)
    # hot_encode_common_df = [0 if el == "non_effectors" else 1 for el in common_df["name"]]
    x_res, y_res = rus.fit_resample(common_df[list(common_df.columns)[1:]], common_df["name"])

    ## CREATE THE RANDOM FOREST CLASSIFIER
    ### K-fold CROSS-VALIDATION
    #### ITERATIVELY ESTIMATE THE BEST NUMBER OF DECISION TREES IN THE FOREST FOR EACH FOLD

    # shuffling is useful because positive and negative elements are ordered in the dataframe and this will help to
    # randomly pick the elements for each fold
    # kf5 = KFold(n_splits=5, shuffle=True)
    ssskf5 = StratifiedShuffleSplit(n_splits=5, test_size=0.3, random_state=42)
    i = 1
    acc = []
    all_auc = []
    all_prec = []
    all_rec = []

    trained_clf = RandomForestClassifier(random_state=42)  # RANDOM STATE IS IMPORTANT FOR FUTURE REPLICATION = SEED
    y_real = []
    y_pred_tot = []
    y_proba = []
    # SPLITTING
    for train_index, test_index in ssskf5.split(common_df, list(common_df["name"])):
    # for train_index, test_index in ssskf5.split(x_res, y_res):
        x_train = common_df.iloc[train_index].loc[:, list(common_df.columns)[1:]]
        x_test = common_df.iloc[test_index][list(common_df.columns)[1:]]
        y_train = common_df["name"].iloc[train_index]
        y_test = common_df["name"].iloc[test_index]

        # ith fold RANDOM FOREST
        clf = RandomForestClassifier(n_estimators=100, random_state=42)  # n_estimators = number of trees in the forest
        ## SELECT BEST n_estimators
        print("BEGIN best parameters selection for Random Forest")
        params_to_test = {"n_estimators": [50, 75, 100, 1000, 5000]}
        grid_search = GridSearchCV(clf, params_to_test, n_jobs=4)
        grid_search.fit(x_train, y_train)
        print("FINISH - best parameters selection")
        best_params = grid_search.best_params_
        clf = RandomForestClassifier(**best_params, random_state=42)

        ## TRAIN THE MODEL USING ith-FOLD TRAINING SETS
        clf.fit(x_train, y_train)

        ## ACCURACY
        y_pred = clf.predict(x_test)
        ### predict probability
        y_pred_prob = clf.predict_proba(x_test)
        tmp_accuracy = metrics.accuracy_score(y_test, y_pred)
        print(f"TEST set accuracy for the fold no. {i}: {tmp_accuracy}")
        acc.append(metrics.accuracy_score(y_test, y_pred))

        ## F-measure
        f_score = metrics.f1_score(y_test, y_pred, average="macro")
        print(f"\nF-Measure:\t{f_score}")
        y_test_hot_encode = [0 if el == "non_effectors" else 1 for el in y_test]
        y_pred_hot_encode = [0 if el == "non_effectors" else 1 for el in y_pred]

        precision, recall, _ = precision_recall_curve(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        print(f"Precision coordinates: {precision}\nRecall(Sensitivity) coordinates: {recall}")
        precision_perc = metrics.precision_score(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        recall_perc = metrics.recall_score(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        print(f"Precision (TP/TP+FP) on real data: {precision_perc}")
        print(f"Recall (TP/TP+FN) on real data: {recall_perc}")
        all_prec.append(precision_perc)
        all_rec.append(recall_perc)
        ## Area Under the Curve
        pr_auc = auc(recall, precision)
        all_auc.append(pr_auc)
        print(f"Area under the curve is: {pr_auc}")

        ## CONFUSION MATRIX
        conf_matrix = metrics.confusion_matrix(y_test, y_pred)
        print(f"\nconfusion matrix\n{conf_matrix}")

        ## FEATURE IMPORTANCE
        feat_imp = pd.Series(clf.feature_importances_, index=list(common_df.columns)[1:]).sort_values(ascending=False)
        print(f"Feature importance:\n{feat_imp}")

        trained_clf = clf
        # disp = PrecisionRecallDisplay(precision=precision, recall=recall)
        # disp.plot()
        y_real.append([0 if el == "non_effectors" else 1 for el in y_test])
        y_pred_tot.append([0 if el == "non_effectors" else 1 for el in y_pred])
        y_proba.append(y_pred_prob[:, 1])
        i += 1

    # AUC PLOT
    # y_real = np.concatenate(y_real)
    # y_pred = np.concatenate(y_pred_tot)
    # precision, recall, _ = precision_recall_curve(y_real, y_pred)
    # disp = PrecisionRecallDisplay(precision=precision, recall=recall)
    # disp.plot()
    # plt.show()

    print(f"Averaged accuracy for {i-1}-fold cross validation: {np.mean(acc)}")
    print(f"Averaged precision for {i-1}-fold cross validation: {np.mean(all_prec)}")
    print(f"Averaged recall for {i-1}-fold cross validation: {np.mean(all_rec)}")
    joblib.dump(trained_clf, f"{path_to_out_trained_model}/random_forest_model_acc{str(np.mean(acc))[:5]}_TMcorr_prof_undersampl_no_samp_strategy.pkl")
    return trained_clf


def train_test_rfc_no_signal_pep(eff_pred_df, non_eff_pred_df, path_to_out_trained_model):
    original_common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    original_common_df.replace([True, False], [1, 0], inplace=True)

    # EFF PREDICTIONS DF for Random Forest
    eff_pred_df["name"] = ["effectors"] * len(eff_pred_df)

    # NON-EFF PREDICTIONS DF for Random Forest
    non_eff_pred_df["name"] = ["non_effectors"] * len(non_eff_pred_df)

    # COMBINE DFs (taking only effector motifs)
    for col in eff_pred_df.columns[3:]:
        if col not in non_eff_pred_df.columns[3:]:
            # if an effector motifs is not present in non effectors, a new column will be added to the non_effector_df
            non_eff_pred_df[col] = [0] * len(non_eff_pred_df)
        else:
            pass

    # "inner" will take only common cols
    common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    common_df.replace([True, False], [1, 0], inplace=True)
    common_df = common_df.drop(["signal peptide", "Phob signal peptide"], 1)
    print(common_df.columns)
    ## CREATE THE RANDOM FOREST CLASSIFIER
    ### K-fold CROSS-VALIDATION
    #### ITERATIVELY ESTIMATE THE BEST NUMBER OF DECISION TREES IN THE FOREST FOR EACH FOLD

    # shuffling is useful because positive and negative elements are ordered in the dataframe and this will help to
    # randomly pick the elements for each fold
    kf5 = KFold(n_splits=5, shuffle=True)
    i = 1
    acc = []
    all_auc = []

    trained_clf = RandomForestClassifier(random_state=42)  # RANDOM STATE IS IMPORTANT FOR FUTURE REPLICATION = SEED
    y_real = []
    y_pred_tot = []
    y_proba = []
    # SPLITTING
    for train_index, test_index in kf5.split(common_df):
        x_train = common_df.iloc[train_index].loc[:, list(common_df.columns)[2:]]
        x_test = common_df.iloc[test_index][list(common_df.columns)[2:]]
        y_train = common_df.iloc[train_index].loc[:, "name"]
        y_test = common_df.iloc[test_index]["name"]

        # ith fold RANDOM FOREST
        clf = RandomForestClassifier(n_estimators=100, random_state=42)  # n_estimators = number of trees in the forest
        ## SELECT BEST n_estimators
        print("BEGIN best parameters selection for Random Forest")
        params_to_test = {"n_estimators": [50, 75, 100, 1000, 5000]}
        grid_search = GridSearchCV(clf, params_to_test, n_jobs=4)
        grid_search.fit(x_train, y_train)
        print("FINISH - best parameters selection")
        best_params = grid_search.best_params_
        clf = RandomForestClassifier(**best_params, random_state=42)

        ## TRAIN THE MODEL USING ith-FOLD TRAINING SETS
        clf.fit(x_train, y_train)

        ## ACCURACY
        y_pred = clf.predict(x_test)
        ### predict probability
        y_pred_prob = clf.predict_proba(x_test)
        tmp_accuracy = metrics.accuracy_score(y_test, y_pred)
        print(f"TEST set accuracy for the fold no. {i}: {tmp_accuracy}")
        acc.append(metrics.accuracy_score(y_test, y_pred))

        ## F-measure
        f_score = metrics.f1_score(y_test, y_pred, average="macro")
        print(f"\nF-Measure:\t{f_score}")
        y_test_hot_encode = [0 if el == "non_effectors" else 1 for el in y_test]
        y_pred_hot_encode = [0 if el == "non_effectors" else 1 for el in y_pred]

        precision, recall, _ = precision_recall_curve(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        print(f"Precision coordinates: {precision}\nRecall(Sensitivity) coordinates: {recall}")

        ## Area Under the Curve
        pr_auc = auc(recall, precision)
        all_auc.append(pr_auc)
        print(f"Area under the curve is: {pr_auc}")

        ## CONFUSION MATRIX
        conf_matrix = metrics.confusion_matrix(y_test, y_pred)
        print(f"\nconfusion matrix\n{conf_matrix}")

        ## FEATURE IMPORTANCE
        feat_imp = pd.Series(clf.feature_importances_, index=list(common_df.columns)[2:]).sort_values(ascending=False)
        print(f"Feature importance:\n{feat_imp}")

        trained_clf = clf
        # disp = PrecisionRecallDisplay(precision=precision, recall=recall)
        # disp.plot()
        y_real.append([0 if el == "non_effectors" else 1 for el in y_test])
        y_pred_tot.append([0 if el == "non_effectors" else 1 for el in y_pred])
        y_proba.append(y_pred_prob[:, 1])
        i += 1

    # AUC PLOT
    # y_real = np.concatenate(y_real)
    # y_pred = np.concatenate(y_pred_tot)
    # precision, recall, _ = precision_recall_curve(y_real, y_pred)
    # disp = PrecisionRecallDisplay(precision=precision, recall=recall)
    # disp.plot()
    # plt.show()

    print(f"Averaged accuracy for {i-1}-fold cross validation: {np.mean(acc)}")
    joblib.dump(trained_clf, f"{path_to_out_trained_model}/random_forest_model_acc{str(np.mean(acc))[:5]}.pkl")
    return trained_clf


def feature_extraction_new_data(data_pred_df, eff_pred_df):
    data_pred_df.replace([True, False], [1, 0], inplace=True)
    eff_pred_df.replace([True, False], [1, 0], inplace=True)
    cols_to_drop = []
    for col in list(data_pred_df.columns):
        if col not in list(eff_pred_df.columns)[3:]:
            cols_to_drop.append(col)
        else:
            pass
    # add columns from effectors that are not present in the input file and set the value to 0 since the only numerical
    # feature is the sequence length and the others are 0/1
    df_features = data_pred_df.drop(cols_to_drop, 1)

    for col in list(eff_pred_df.columns)[3:]:
        if col not in df_features:
            df_features[col] = [0] * len(df_features)
        else:
            pass
    df_features = df_features.drop("EF_HAND_1_L=(-1)", 1)
    # print(list(eff_pred_df.columns))
    # print(list(df_features.columns))
    # print(len(list(eff_pred_df.columns)[3:])==len(list(df_features.columns)))
    return df_features


def rf_positives_check(clf, annot_table):
    annot = pd.read_csv(annot_table, sep="\t")
    j = 0
    k = []
    v = []
    # add to the classification df only those annotation for the proteins RF classify as effectors
    for i in range(len(clf)):
        if clf.iloc[i]["RF_model_classification"] == "effectors":
            j += 1
            k.append(clf.iloc[i]["CaPm7_contig01_proteins"])
            v.append(str(list(annot[annot.columns[1:]].iloc[i])).replace("[", "").replace("]", ""))
        else:
            k.append(clf.iloc[i]["CaPm7_contig01_proteins"])
            v.append(" ")
    cols_to_add = pd.DataFrame(list(zip(k, v)), columns=[str(clf.columns[0]), "annotation"])
    check_df = clf.merge(cols_to_add, on=str(clf.columns[0]))
    return check_df


def parse_annot_table_to_eliminate_unknown(annot_table):
    at = pd.read_csv(annot_table, sep="\t")
    index_to_drop = []
    supervised_class = []
    for row_index in range(len(at)):
        if "Uncharacterized" in str(at.iloc[row_index]["Protein names"]) or "nan" in str(at.iloc[row_index]["Protein names"]):
            index_to_drop.append(row_index)
        elif "effector" in str(at.iloc[row_index]["Protein names"]) or "Effector" in str(at.iloc[row_index]["Protein names"]) or "SAP" in str(at.iloc[row_index]["Protein names"]) or "TENGU" in str(at.iloc[row_index]["Protein names"]):
            supervised_class.append("effectors")
        else:
            supervised_class.append("non_effectors")

    at = at.drop(index_to_drop)
    filt_annot = pd.DataFrame(list(zip(list(at["ID"]), supervised_class)), columns=["name", "annotation"])
    return filt_annot, index_to_drop


def parse_feature_table_to_eliminate_unknown(feature_table, index_to_drop):
    ft = feature_table.drop(index_to_drop)
    return ft


def evaluate_performances_on_real_data(parsed_annot_cl, real_data_clf):
    print(f"Accuracy on real data classification: {metrics.accuracy_score(parsed_annot_cl, real_data_clf)}")
    precision = metrics.precision_score(parsed_annot_cl, real_data_clf, labels=["effectors", "non-_effectors"], pos_label="effectors")
    recall = metrics.recall_score(parsed_annot_cl, real_data_clf, labels=["effectors", "non-_effectors"], pos_label="effectors")
    print(f"Precision (TP/TP+FP) on real data: {precision}")
    print(f"Recall (TP/TP+FN) on real data: {recall}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="",
                                     epilog="")

    parser.add_argument("-e", "--effector_predictions",
                        help="TSV file of feature predictions")
    parser.add_argument("-ne", "--non_effector_predictions",
                        help="TSV file of feature predictions")
    parser.add_argument("-orf", "--output_path_rf",
                        help="path to save the random forest model")
    parser.add_argument("-cl", "--classifier",
                        help="pickle file of the random forest model to classify proteins")
    parser.add_argument("--apply_model",
                        action="store_true",
                        help="when the training is already done, use --apply_model to classify new data")
    parser.add_argument("-d", "--data")
    parser.add_argument("-a", "--annotations_table")
    parser.add_argument("-o", "--output_classification")

    args = parser.parse_args()

    annot, index = parse_annot_table_to_eliminate_unknown(args.annotations_table)
    # data_features = parse_feature_table_to_eliminate_unknown(args.data, index)
    rf_classifier = joblib.load(args.classifier)
    data = feature_extraction_new_data(parse_table(args.data), parse_table(args.effector_predictions))
    classification = rf_classifier.predict(data[list(data.columns)])
    clf_df = pd.DataFrame(list(zip(list(parse_table(args.data)["name"]), classification)),
                          columns=["name", "RF_model_classification"])
    print(clf_df.value_counts())
    clf_df_filt = parse_feature_table_to_eliminate_unknown(clf_df, index)
    evaluate_performances_on_real_data(np.array(annot[annot.columns[1]]), np.array(clf_df_filt[clf_df_filt.columns[1]]))
    exit(1)
    if not args.apply_model:
        train_test_rfc(parse_table(args.effector_predictions),
                       parse_table(args.non_effector_predictions),
                       args.output_path_rf)
        # train_test_rfc_no_signal_pep(parse_table(args.effector_predictions),
        #                parse_table(args.non_effector_predictions),
        #                args.output_path_rf)
    else:
        rf_classifier = joblib.load(args.classifier)
        data = feature_extraction_new_data(parse_table(args.data), parse_table(args.effector_predictions))
        classification = rf_classifier.predict(data)
        clf_df = pd.DataFrame(list(zip(list(parse_table(args.data)["name"]), classification)),
                              columns=["name", "RF_model_classification"])

        print(clf_df["RF_model_classification"].value_counts())
        rf_positives_check(clf_df, args.annotations_table).to_csv(f"{args.output_classification}/check_on_classification.tsv", sep="\t")


