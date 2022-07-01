import pandas as pd
import numpy as np
import joblib
import argparse
from sklearn.model_selection import train_test_split  # to split the dataset
from sklearn.ensemble import RandomForestClassifier  # random forest model
from sklearn.model_selection import GridSearchCV
from sklearn import metrics  # to calculate the accuracy of the model
from sklearn.metrics import (precision_recall_curve, PrecisionRecallDisplay)
from sklearn.model_selection import KFold


def parse_table(input_file, h=False):
    # df of feature to be considered
    if not h:
        parsed_df = pd.read_csv(input_file, sep="\t", header=1)
    else:
        parsed_df = pd.read_csv(input_file, sep="\t", header=0, index_col=0)
    return parsed_df


def parse_excell(input_file):
    parsed_df = pd.read_excel(input_file, sheet_name="protein_features_table", header=1)
    return parsed_df


def train_test_rfc(eff_df, non_eff_df):  #, hydrophob_df):
    original_common_df = pd.concat([eff_df, non_eff_df], ignore_index=True, join="inner")
    original_common_df.replace([True, False], [1, 0], inplace=True)

    # EFF DF for Random Forest
    eff_df["name"] = ["effectors"] * len(eff_df)
    # eff_df["sequence length"] = [True] * len(eff_df)
    ## length selection later
    # eff_min_len = np.min(eff_df["sequence length"])
    # eff_max_len = np.max(eff_df["sequence length"])

    # NON-EFF DF for Random Forest
    non_eff_df["name"] = ["non_effectors"] * len(non_eff_df)

    ## NO LENGTH SELECTION
    index_to_drop = []
    ### filter just to be sure to not have mistaken into negative set (non effectors)
    for i in range(len(non_eff_df)):
        desc = non_eff_df["description"].iloc[i]
        if "effector" in desc or "Effector" in desc:
            index_to_drop.append(i)

    non_eff_df = non_eff_df.drop(non_eff_df.index[index_to_drop])  # drop not considered rows
    non_eff_df = non_eff_df.drop("description", 1)  # to make the two dataframe compatible for the merging process

    # non_eff_df["sequence length"] = [True if eff_min_len > p_len > eff_max_len
    #                                  else False for p_len in non_eff_df["sequence length"]]
    # i = 0
    # for el in non_eff_df["sequence length"]:
    #     if eff_min_len < el < eff_max_len:
    #         print(f"TRUE\n{el}")
    #         i += 1
    #     else:
    #         print(".")

    # COMBINE DFs (taking only common features)
    common_df = pd.concat([eff_df, non_eff_df], ignore_index=True, join="inner")
    common_df.replace([True, False], [1, 0], inplace=True)

    # ADD HYDROPHOBICITY PREDICTION TO COMMON DF
    # for col in hydrophob_df:
    #     common_df[col] = list(hydrophob_df[col])

    # RANDOM FOREST
    x = common_df[list(common_df.columns)[1:]]  # features=variables
    y = common_df["name"]  # labels=output

    ## SPLIT DATASET (TRAIN/TEST - 70/30)
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3)
    # print(f"x_train\n{original_common_df.iloc[list(x_train.index)]['name']}")
    # print(f"x_test\n{x_test}")

    # ## CREATE A GAUSSIAN CLASSIFIER
    # clf = RandomForestClassifier(n_estimators=100)  # number of trees in the forest
    # ## SELECT BEST n_estimators
    # print("BEGIN best parameters selection for Random Forest")
    # params_to_test = {"n_estimators": [50, 75, 100, 1000, 5000]}
    # grid_search = GridSearchCV(clf, params_to_test, n_jobs=4)
    # grid_search.fit(x_train, y_train)
    # print("FINISH - best parameters selection")
    # best_params = grid_search.best_params_
    # print(best_params)
    # clf = RandomForestClassifier(**best_params)
    # clf.fit(x_train, y_train)
    ## K-fold CROSS-VALIDATION
    kf3 = KFold(n_splits=3, shuffle=True)
    i = 1
    acc = []
    # RANDOM STATE IS IMPORTANT FOR FUTURE REPLICATION
    trained_clf = RandomForestClassifier(random_state=42)
    for train_index, test_index in kf3.split(common_df):
        # print(f"train_set {common_df.iloc[train_index]}\ntest_set {common_df.iloc[test_index]}")
        x_train = common_df.iloc[train_index].loc[:, list(common_df.columns)[1:]]
        x_test = common_df.iloc[test_index][list(common_df.columns)[1:]]
        y_train = common_df.iloc[train_index].loc[:, "name"]
        y_test = common_df.iloc[test_index]["name"]

        # ## CREATE A GAUSSIAN CLASSIFIER
        clf = RandomForestClassifier(n_estimators=100, random_state=42)  # number of trees in the forest
        ## SELECT BEST n_estimators
        print("BEGIN best parameters selection for Random Forest")
        params_to_test = {"n_estimators": [50, 75, 100, 1000, 5000]}
        grid_search = GridSearchCV(clf, params_to_test, n_jobs=4)
        grid_search.fit(x_train, y_train)
        print("FINISH - best parameters selection")
        best_params = grid_search.best_params_
        print(best_params)
        clf = RandomForestClassifier(**best_params, random_state=42)

        ## TRAIN THE MODEL USING THE TRAINING SETS y_pred = clf.predict(x_test)
        # clf = RandomForestClassifier(n_estimators=100)
        clf.fit(x_train, y_train)
        print(f"Accuracy for the fold no. {i} on the test set: {metrics.accuracy_score(y_test, clf.predict(x_test))}")
        i += 1
        y_pred = clf.predict(x_test)

        ## ACCURACY
        print(f"Accuracy: {metrics.accuracy_score(y_test, y_pred)}")
        acc.append(metrics.accuracy_score(y_test, y_pred))
        # first iteration: Accuracy: 0.9894736842105263
        # second, third, forth, ..., : Accuracy: 1.0

        ## FEATURE IMPORTANCE
        feature_imp = pd.Series(clf.feature_importances_, index=list(common_df.columns)[1:]).sort_values(ascending=False)
        print(f"Feature importance:\n{feature_imp}")
        # sequence length 0.377685
        # signal peptide 0.297652
        # Phob signal peptide 0.231584
        # MYRISTYL 0.044473
        # disordered regions 0.015679
        # PKC_PHOSPHO_SITE 0.011141
        # ASN_GLYCOSYLATION 0.008138
        # CAMP_PHOSPHO_SITE 0.003997
        # CK2_PHOSPHO_SITE 0.003556
        # Phob transmembrane domain 0.003174
        # AMIDATION 0.002279
        # TYR_PHOSPHO_SITE_1 0.000642
        # transmembrane domain 0.000000

        ## CONFUSION MATRIX
        conf_matrix = metrics.confusion_matrix(y_test, y_pred)
        print(f"\nconfusion matrix\n{conf_matrix}")

        ## F-measure
        f_score = metrics.f1_score(y_test, y_pred, average="macro")
        print(f"\nF-Measure:\t{f_score}")
        y_test_hot_encode = [0 if el == "non_effectors" else 1 for el in y_test]
        y_pred_hot_encode = [0 if el == "non_effectors" else 1 for el in y_pred]

        precision, recall, _ = precision_recall_curve(y_test_hot_encode, y_pred_hot_encode, pos_label=1)
        print(f"Precision coordinates: {precision}\nRecall(Sensitivity) coordinates: {recall}")

        ## checks
        print(f"Y_test\n{original_common_df.iloc[list(y_test.index)]['name']}")
        print(f"Y_pred\n{y_pred}")
        trained_clf = clf

    print(f"final {trained_clf.predict(common_df.iloc[list(range(120, 300))].loc[:, list(common_df.columns)[1:]])}")
    print(np.mean(acc))
    joblib.dump(trained_clf, "random_forest_model.pkl")
    return trained_clf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="",
                                     epilog="")

    parser.add_argument("-e", "--effector_predictions",
                        help="TSV file of feature predictions")
    parser.add_argument("-ne", "--non_effector_predictions",
                        help="TSV file of feature predictions")
    parser.add_argument("-cl", "--classifier",
                        help="pickle file of the random forest model to classify proteins")
    parser.add_argument("--apply_model",
                        action="store_true",
                        help="when the training is already done, use --apply_model to classify new data")

    # eff = "...effectors/analysis/prediction_tools/effector_protein_feature_prediction.tsv"
    # non_eff = "...effectors/analysis/prediction_tools/062022non_effector_protein_feature_prediction.tsv"
    # hydrophob_pred = "/home/giulia/Workspace/PhytoPhD/effectors/analysis/KD/hydrophob_comparisons.tsv"

    args = parser.parse_args()
    if not args.apply_model:
        train_test_rfc()
    else:
        rfc = joblib.load(args.classifier)
        rfc.predict()


## RESULTS AFTER HYDROPHOBICITY PROFILE ADDITION
# {'n_estimators': 75}
# Accuracy: 1.0
# Feature importance:
# sequence length              0.308751
# Phob signal peptide          0.244183
# signal peptide               0.140814
# bin_2                        0.093630
# bin_1                        0.082008
# MYRISTYL                     0.042225
# bin_13                       0.011926
# disordered regions           0.011282
# bin_6                        0.010429
# bin_9                        0.008407
# PKC_PHOSPHO_SITE             0.005789
# bin_8                        0.005684
# bin_5                        0.005171
# bin_7                        0.005009
# bin_12                       0.004768
# bin_4                        0.004147
# bin_11                       0.003519
# bin_15                       0.002512
# ASN_GLYCOSYLATION            0.002429
# bin_3                        0.002303
# Phob transmembrane domain    0.001607
# bin_10                       0.001311
# CK2_PHOSPHO_SITE             0.001137
# TYR_PHOSPHO_SITE_1           0.000373
# AMIDATION                    0.000309
# bin_14                       0.000252
# CAMP_PHOSPHO_SITE            0.000024
# transmembrane domain         0.000000
# dtype: float64
#
# confusion matrix
# [[24  0]
#  [ 0 71]]
#
# F-Measure:      1.0
# Precision coordinates: [1. 1.]
# Recall(Sensitivity) coordinates: [1. 0.]
