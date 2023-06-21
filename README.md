# LEAF - phytopLasmas Effector prediction via rAndom Forest

Both positive (effector proteins) and negative (non-effector proteins) dataset are
retrieved from Uniprot database since **June 2022**.

## Protein Selection Criteria - Filtered Out

### POSITIVE DATASET
Proteins annotated in Uniprot as Effector, SAP05, SAP11, SAP54, TENGU, or PHYL1/phyllody, are selected. Comprising those experimentally validated and included in Uniprot.

### NGATIVE DATASET
The negative dataset at **July-August 2022** is composed by proteins belonging to ANY Candidatus phytoplasma whos annotation is MANUALLY CURATED, a publication is present, and the annotation seems to not be involved into the pathogenicity activity (no SecA-Y proteins, no SVM proteins, no trigger factors, no RuvA-B proteins, no FtsH), plus no putative, possible or recombinant proteins. 
Included ribosomal proteins

## Numbers Considered for the analysis
- Effector proteins: 184
- - Non-Effector proteins: 295 (some are fragments and some of them have and X in the AA sequence)

## Feature Considered
Protein features (presence/absence, except for sequence length) considered for
classification are:
- Protein length
- Signal Peptide (SignalP 4.1)
- Transmembrane Region (TMHMM)
 - AA in transmembrane domain
 - first 60 AA
 - Probability of N-in
 - Warning signal sequence 
-  Intrinsically disordered regions (MobiDB-lite)
- Motifs/Profiles in AA sequence (Prosite)
- CLUMPs (MOnSTER)
- bin of sequence and presence of CLUMPs


All features are considered as numerical variables, continuous or not.

## Classification
- Method: **Random Forest**
- Overall data-set partition and validation: **5-fold cross validation** 
- Goodness of model:
  - Accuracy
  - F-measure
  - Precision-Recall
  - Feature importance 
  - Random labeling
- number of variables (features): 13
Taking into consideration only features in common between Effectors and Non-
Effectors data-sets (to solve the problems due to motifs or profiles more abundant
in Non-Effectors compared to Effectors and also because of the classification objective to search for the positive class characteristics)

```
clf = RandomForestClassifier(**best_params, random_state=42, max_feature="sqrt")
```
****best_params** refers to number of trees in the forest and is calculated with the GridSearchCV with 4 iterations, **random_state=42** is used to set seed for future training of the model, to always start from the same point and **max_feature="sqrt"** is specified even thought is the default parameter, it allows a splitting of the features for each tree in the forest so that correlated features do not influence the decision, instead contribute more efficiently to the decision (some trees take decisions on some features, others on other features) this randomness and the combination of different decision trees will reduce the variance of the forest estimator --> reduce the overfitting typical tendency of single decision trees that presents a higher variance.  


The RF can reach an accuracy of 97% in effector prediction by the training process



