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


def train_test_rfc(eff_df, non_eff_df, path_to_trained_model):  #, hydrophob_df):
    original_common_df = pd.concat([eff_df, non_eff_df], ignore_index=True, join="inner")
    original_common_df.replace([True, False], [1, 0], inplace=True)

    # EFF DF for Random Forest
    eff_df["name"] = ["effectors"] * len(eff_df)

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

    # COMBINE DFs (taking only common features)
    common_df = pd.concat([eff_df, non_eff_df], ignore_index=True, join="inner")
    common_df.replace([True, False], [1, 0], inplace=True)

    # ADD HYDROPHOBICITY PREDICTION TO COMMON DF (da rivedere)
    # for col in hydrophob_df:
    #     common_df[col] = list(hydrophob_df[col])

    # RANDOM FOREST
    x = common_df[list(common_df.columns)[1:]]  # features=variables
    y = common_df["name"]  # labels=output

    ## CREATE THE RANDOM FOREST CLASSIFIER
    ### K-fold CROSS-VALIDATION
    #### ITERATIVELY ESTIMATE THE BEST NUMBER OF DECISION TREES IN THE FOREST
    # shuffling is useful because positive and negative elements are ordered in the dataframe and this will help to
    # randomly pick the elements for each fold
    kf5 = KFold(n_splits=5, shuffle=True)
    i = 1
    acc = []

    trained_clf = RandomForestClassifier(random_state=42)  # RANDOM STATE IS IMPORTANT FOR FUTURE REPLICATION
    # SPLITTING
    for train_index, test_index in kf5.split(common_df):
        x_train = common_df.iloc[train_index].loc[:, list(common_df.columns)[1:]]
        x_test = common_df.iloc[test_index][list(common_df.columns)[1:]]
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

        ## CONFUSION MATRIX
        conf_matrix = metrics.confusion_matrix(y_test, y_pred)
        print(f"\nconfusion matrix\n{conf_matrix}")

        ## FEATURE IMPORTANCE
        feat_imp = pd.Series(clf.feature_importances_, index=list(common_df.columns)[1:]).sort_values(ascending=False)
        print(f"Feature importance:\n{feat_imp}")

        trained_clf = clf
        i += 1

    print(f"Averaged accuracy for {i-1}-fold cross validation: {np.mean(acc)}")
    joblib.dump(trained_clf, f"{path_to_trained_model}/random_forest_model.pkl")
    return trained_clf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="",
                                     epilog="")

    parser.add_argument("-e", "--effector_predictions",
                        help="TSV file of feature predictions")
    parser.add_argument("-ne", "--non_effector_predictions",
                        help="TSV file of feature predictions")
    parser.add_argument("-o", "--output_path",
                        help="path to save the random forest model")
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
        train_test_rfc(parse_table(args.effector_predictions),
                       parse_table(args.non_effector_predictions),
                       args.output_path)
    else:
        rf_classifier = joblib.load(args.classifier)
        rf_classifier.predict()


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
