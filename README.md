# effectors_classifier/EFFECTOR_HUNTER

Both positive (effector proteins) and negative (non-effector proteins) dataset are
retrieved from Uniprot database.

## Protein Selection Criteria - Filtered Out
The negative dataset at **June 2022** is composed by proteins belonging to Candi-
datus Phytoplasma mali and other two closely related phytoplasmas, Candidatus
Phytoplasma pruni and Peanut Witches’ Broom phytoplasma.
A **first selection** was done by **sequence length** –> those proteins having a length
comprised between min(positive seq length) +- standard deviaton and max(positive
seq length) +- standard deviation were kept.
The **second selection** was done by **alignment** of the negative proteins to the
positive ones (BLAST). From this alignment a **threshold of protein similarity** was
set using the formula described here _Rost, B. Twilight zone of protein sequence
alignments. Protein Eng. 12, 85–94 (1999)_; correlating %identity and alignment
length for each **pairwise alignment** –> those proteins having %identity below the
similarity threshold were kept.
`%ide < pS (n = n + 420 ∗ L−0,335∗(1+e−L/2000))`

The negative dataset at **July-August 2022** is composed by proteins belonging to ANY Candidatus phytoplasma whos annotation is MANUALLY CURATED, a publication is present, and the annotation seems to not be involved into the pathogenicity activity (no SecA-Y proteins, no SVM proteins, no trigger factors, no RuvA-B proteins, no FtsH), plus no putative, possible or recombinant proteins. 
Included ribosomal proteins

## Numbers Considered for the analysis
- Effector proteins: 88
- Non-Effector proteins: 252

## Feature Considered
Protein features (presence/absence, except for sequence length) considered for
classification are:
 - Protein length
- Signal Peptide (SignalP 4.1 + Phobius)
- Transmembrane Region (TMHMM + Phobius)
-  Intrinsically disordered regions (MobiDB-lite)
- Motifs/Profiles in AA sequence (Prosite)
All features are considered as categorical variables, presence or absence, namely 1
or 0 values EXCEPT FOR SEQUENCE LENGTH.

## Classification
- Method: **Random Forest**
- Overall data-set partition and validation: **5-fold cross validation**
- Goodness of model:
  - Accuracy
  - F-measure
  - Area Under the Curve
  - Precision-Recall Curve
- number of variables (features): 13
Taking into consideration only features in common between Effectors and Non-
Effectors data-sets (to solve the problems due to motifs or profiles more abundant
in Non-Effectors compared to Effectors and also because of the classification objective to search for the positive class characteristics)

```
clf = RandomForestClassifier(**best_params, random_state=42, max_feature="sqrt")
```
****best_params** refers to number of trees in the forest and is calculated with the GridSearchCV with 4 iterations, **random_state=42** is used to set seed for future training of the model, to always start from the same point and **max_feature="sqrt"** is specified even thought is the default parameter, it allows a splitting of the features for each tree in the forest so that correlated features do not influence the decision, instead contribute more efficiently to the decision (some trees take decisions on some features, others on other features) this randomness and the combination of different decision trees will reduce the variance of the forest estimator --> reduce the overfitting typical tendency of single decision trees that presents a higher variance.  

## Results on CaPm7_contig01 proteins predictions
It classifies 14 proteins out of 462, as effectors.
As a first check the 462 sequences were aligned to the 88 effector dataset to see if some of them are effectively similar to some of the effectors used to train the model. NONE of the alignment was sufficiently long or have a sufficiently high % of identity with the effectors. This means that basically we cannot trast BLAST alignment to "validate" our findings. 
The second check was done intersecting the annotation of the proteins classified by the RF as effectors. These annotations are obtained aligning the predicted protein sequences of CaPm to a Uniprot database containing the complete set of known proteins in phytoplasmas. In this case we found that 3/14 classified effectors had a different annotation. At first side those annotations didn't have strong evidences against the RF classification of effectors but then searching with HMMER the domains present in the proteins into the 88 effector protein sequences, we couldn't find any matches. This means that none of the effectors considered for the training and testing of the model contains those domains. 

--> **The question is then why the RF classify those proteins as effectors? Do we have considered too few features in describing an effector? Or do we have to extend the negative proteins set?** 
### From some works on FtsH proteins, at least for CaPm7_contig01_0001 protein, we have some studies that evidence, for this protein, a role in the virulence of phytoplamsmas. Thus we are good in keep its effector classification. 

## Added the Hydrophobicity feature
Before the hydrophobicity scores were too complicated and applied in an unbalanced manner as features. Insteat we try (07/072022) to add as additional feature, the GRAVY score. GRAVY stands for GRand AVerage of hYdrpathy and is calculated on the basis of Kyte-Doolittle scores; it basically give a mean hydropathy score for all AA in the sequence. 
The gravy score in the non-effector set have a wider range of values with respect to the effector set, even when the non-effctors are picked randomly and in equal number of the effectors. 
```
effectors_g_mean:-0.53
effectors_g_median:-0.52
effectors_g_standard_deviation:0.20

NON-effectors_g_mean:-0.38
NON-effectors_g_median:-0.37
NON-effectors_g_standard_deviation:0.53

NON-effectors_g_mean(randomly taken 88):-0.37
NON-effectors_g_median(randomly taken 88):-0.38
NON-effectors_g_standard_deviation(randomly taken 88):0.55
```

## Results with Hydrophobicity feature on CaPm7_contig01 proteins predictions
We obtain 13 effector classification instead of 14. The missing classification is an uncharacterized protein from the Ca.P.mali AT strain. 
