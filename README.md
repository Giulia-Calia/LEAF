# effectors_classifier/EFFECTOR_HUNTER

Both positive (effector proteins) and negative (non-effector proteins) dataset are
retrieved from Uniprot database.

## Protein Selection Criteria - Filtered Out
The negative dataset at June 2022 is composed by proteins belonging to Candi-
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

