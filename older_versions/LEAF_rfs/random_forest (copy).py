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
from sklearn.model_selection import StratifiedKFold
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from collections import Counter
from imblearn.under_sampling import RandomUnderSampler
import random


def parse_table(input_file):
    # df of feature to be considered
    parsed_df = pd.read_csv(input_file, sep="\t", header=0)
    return parsed_df


def train_test_rfc(eff_pred_df, non_eff_pred_df, path_to_out_trained_model):
    original_common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    original_common_df.replace([True, False], [1, 0], inplace=True)

    # EFF PREDICTIONS DF for Random Forest
    eff_pred_df["name"] = ["effectors"] * len(eff_pred_df)
    eff_pred_df = eff_pred_df.drop(["ID", "organism"], 1)

    # NON-EFF PREDICTIONS DF for Random Forest
    non_eff_pred_df["name"] = ["non_effectors"] * len(non_eff_pred_df)

    # COMBINE DFs (taking only effector motifs)
    for col in eff_pred_df.columns[1:]:
        if col not in non_eff_pred_df.columns[1:]:
            # if an effector motifs is not present in non effectors, a new column will be added to the non_effector_df
            # having value False or 0 in this case
            non_eff_pred_df[col] = [0] * len(non_eff_pred_df)
        else:
            pass

    # CONCAT AND CLEAN UP THE DATAFRAME FOR THE TRAINING
    # "inner" will take only common cols
    common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    print(f"common_df: {common_df}")
    common_df.replace([True, False], [1, 0], inplace=True)

    # drop the column/s having the same value for all samples == ZERO VARIANCE / ALMOST ZERO VARIANCE
    common_df = common_df.drop("EF_HAND_1_L=(-1)", 1)

    ## CREATE THE RANDOM FOREST CLASSIFIER
    ### STRATIFIED SHUFFLE SPLIT K-fold CROSS-VALIDATION
    #### ITERATIVELY ESTIMATE THE BEST NUMBER OF DECISION TREES IN THE FOREST FOR EACH FOLD

    # shuffling is useful because positive and negative elements are ordered in the dataframe and this will help to
    # randomly pick the elements for each fold

    skf5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    i = 1

    # take trace of all evaluation parameters for each k-fold
    all_acc = []
    all_auc = []
    all_prec = []
    all_rec = []

    trained_clf = RandomForestClassifier(random_state=42)  # RANDOM STATE IS IMPORTANT FOR FUTURE REPLICATION = SEED
    y_real = []
    y_pred_tot = []
    y_proba = []

    all_best_params = []
    # SPLITTING
    for train_index, test_index in skf5.split(common_df, list(common_df["name"])):
        x_train = common_df.iloc[train_index].loc[:, list(common_df.columns)[1:]]  # all cols except the name_col
        x_test = common_df.iloc[test_index][list(common_df.columns)[1:]]
        y_train = common_df["name"].iloc[train_index]  # only the name_col
        y_test = common_df["name"].iloc[test_index]

        # kth fold RANDOM FOREST
        clf = RandomForestClassifier(n_estimators=100, random_state=42)  # n_estimators = number of trees in the forest
        ## SELECT BEST n_estimators
        print("BEGIN best parameters-n_trees selection for Random Forest")
        params_to_test = {"n_estimators": [50, 75, 100, 1000, 5000]}
        grid_search = GridSearchCV(clf, params_to_test, n_jobs=4)
        grid_search.fit(x_train, y_train)
        print("FINISH - best parameters-n_trees selection")
        best_params = grid_search.best_params_
        clf = RandomForestClassifier(**best_params, random_state=42)
        all_best_params.append(best_params["n_estimators"])
        ## TRAIN THE MODEL USING kth-FOLD TRAINING SETS
        clf.fit(x_train, y_train)
        ## APPLY THE MODEL
        y_pred = clf.predict(x_test)

        ## ACCURACY
        kth_accuracy = metrics.accuracy_score(y_test, y_pred)
        print(f"Training phase accuracy for the fold no. {i}: {kth_accuracy}")
        all_acc.append(kth_accuracy)

        ### predict probability
        y_pred_prob = clf.predict_proba(x_test)

        ## F-measure
        f_score = metrics.f1_score(y_test, y_pred, average="macro")
        print(f"\nF-Measure:\t{f_score}")

        ## PRECISION AND RECALL
        y_test_hot_encode = [0 if el == "non_effectors" else 1 for el in y_test]
        y_pred_hot_encode = [0 if el == "non_effectors" else 1 for el in y_pred]
        precision_perc = metrics.precision_score(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        recall_perc = metrics.recall_score(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        print(f"Precision (TP/TP+FP) on training phase: {precision_perc}")
        print(f"Recall (TP/TP+FN) on training phase: {recall_perc}")
        precision, recall, _ = precision_recall_curve(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        print(f"Precision coordinates: {precision}\nRecall(Sensitivity) coordinates: {recall}")
        all_prec.append(precision_perc)
        all_rec.append(recall_perc)
        ## Area Under the Curve
        pr_auc = auc(recall, precision)
        all_auc.append(pr_auc)
        # disp = PrecisionRecallDisplay(precision=precision, recall=recall)
        # disp.plot()
        print(f"Area under the curve is: {pr_auc}")

        ## CONFUSION MATRIX
        conf_matrix = metrics.confusion_matrix(y_test, y_pred)
        print(f"\nconfusion matrix on k-fold test set\n{conf_matrix}")

        ## FEATURE IMPORTANCE
        feat_imp = pd.Series(clf.feature_importances_, index=list(common_df.columns)[1:]).sort_values(ascending=False)
        print(f"Feature importance:\n{feat_imp}")

        trained_clf = clf
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

    print(f"Averaged accuracy for {i-1}-fold cross validation: {np.mean(all_acc)}")
    print(f"Averaged precision for {i-1}-fold cross validation: {np.mean(all_prec)}")
    print(f"Averaged recall for {i-1}-fold cross validation: {np.mean(all_rec)}")
    return np.mean(all_best_params), common_df


def best_model(n_trees, all_training_data, path_to_out_best_rf_model):
    best_rf = RandomForestClassifier(n_estimators=n_trees, random_state=42)
    best_rf.fit(all_training_data[list(all_training_data.columns)[1:]], all_training_data["name"])
    joblib.dump(best_rf, f"{path_to_out_best_rf_model}/random_forest_best_model.pkl")
    return best_rf


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


def test_best_model_random_labelling(clf, eff_pred_df, non_eff_pred_df):
    # EFF PREDICTIONS DF for Random Forest
    eff_pred_df["name"] = ["effectors"] * len(eff_pred_df)
    eff_pred_df = eff_pred_df.drop(["ID", "organism"], 1)

    # NON-EFF PREDICTIONS DF for Random Forest
    non_eff_pred_df["name"] = ["non_effectors"] * len(non_eff_pred_df)

    # COMBINE DFs (taking only effector motifs)
    for col in eff_pred_df.columns[1:]:
        if col not in non_eff_pred_df.columns[1:]:
            # if an effector motifs is not present in non effectors, a new column will be added to the non_effector_df
            # having value False or 0 in this case
            non_eff_pred_df[col] = [0] * len(non_eff_pred_df)
        else:
            pass

    # CONCAT AND CLEAN UP THE DATAFRAME FOR THE TRAINING
    # "inner" will take only common cols
    common_df = pd.concat([eff_pred_df, non_eff_pred_df], ignore_index=True, join="inner")
    common_df.replace([True, False], [1, 0], inplace=True)

    # drop the column/s having the same value for all samples == ZERO VARIANCE / ALMOST ZERO VARIANCE
    common_df = common_df.drop("EF_HAND_1_L=(-1)", 1)
    random.shuffle(x := list(common_df["name"]))
    pred = clf.predict(common_df[list(common_df)[1:]])
    acc = metrics.accuracy_score(x, pred)
    print(acc)


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

    if not args.apply_model:
        best_par = train_test_rfc(parse_table(args.effector_predictions),
                       parse_table(args.non_effector_predictions),
                       args.output_path_rf)

        best_model(best_par[0], best_par[1], args.output_path_rf)

    else:
        rf_classifier = joblib.load(args.classifier)
        data = feature_extraction_new_data(parse_table(args.data), parse_table(args.effector_predictions))
        classification = rf_classifier.predict(data)
        clf_df = pd.DataFrame(list(zip(list(parse_table(args.data)["name"]), classification)),
                              columns=["name", "RF_model_classification"])
        print(clf_df["RF_model_classification"].value_counts())
        test_best_model_random_labelling(rf_classifier, parse_table(args.effector_predictions),
                       parse_table(args.non_effector_predictions))



